# Numerical accuracy of the nested donor solves in the ordinary-h beta root.
# These checks never replace the score/step, curvature or recession criteria.

hl_fo_inner_geometry <- function(state, K, beta_star) {
  S <- as.numeric(crossprod(K, state$S_beta))
  B <- crossprod(K, state$B_beta %*% K)
  eig <- tryCatch(eigen(-0.5 * (B + t(B)), symmetric = TRUE,
                        only.values = TRUE)$values, error = function(e) NA_real_)
  failure <- function(reason) list(success = FALSE, status = "weak_curvature",
                                  reason = reason, B = B, curvature = eig,
                                  beta = state$beta, eta = state$eta)
  if (any(!is.finite(eig)) || min(eig) <= 0 ||
      hl_fo_matrix_rank(B) < ncol(K)) {
    return(failure("inner_curvature_not_positive_definite"))
  }
  delta <- tryCatch(solve(B, -S), error = function(e) NULL)
  if (is.null(delta) || any(!is.finite(delta))) {
    return(failure("inner_linear_solve_failed"))
  }
  gamma <- as.numeric(crossprod(K, state$beta - beta_star))
  score_metric <- sqrt(max(0, sum(S * delta)))
  step_metric <- max(abs(delta)) / (1 + max(abs(gamma)))
  if (any(!is.finite(c(gamma, score_metric, step_metric)))) {
    return(failure("inner_nonfinite_newton_metric"))
  }
  list(success = TRUE, S = S, B = B, delta = delta, gamma = gamma,
       score_metric = score_metric, step_metric = step_metric,
       min_curvature = min(eig))
}

hl_fo_inner_mode_accuracy <- function(data, state, rho_full, K, geometry) {
  # A local first-order error estimate, not a substitute root certificate.
  # If A is donor curvature and D the donor-by-beta mixed derivative, polishing
  # changes the reduced score by approximately -K'D'A^{-1}h. Triangle bounds
  # retain individual donor contributions instead of relying on cancellation.
  diagnostics <- state$mode_diagnostics
  invalid <- list(success = FALSE, reason = "inner_mode_accuracy_invalid")
  if (length(diagnostics) != data$m || length(state$eta) != data$m) return(invalid)
  h <- tryCatch(vapply(diagnostics, function(x) x$h_eta, numeric(1)),
                error = function(e) NULL)
  a <- tryCatch(vapply(diagnostics, function(x) x$h_etaeta, numeric(1)),
                error = function(e) NULL)
  if (is.null(h) || is.null(a)) return(invalid)
  if (any(!is.finite(c(h, a, state$eta))) || any(a >= 0)) return(invalid)
  correction <- -h / a
  log_m <- data$log_offset + as.numeric(data$X %*% state$beta) +
    state$eta[data$subject]
  # Only the ordinary log-mean curvature is needed here. Use exactly its
  # primitive formula, without computing unrelated dispersion derivatives.
  xi <- rho_full[2L]
  phi <- exp(xi)
  curvature <- -(data$y + phi) * stats::plogis(log_m - xi) *
    stats::plogis(xi - log_m)
  if (!is.finite(phi) || phi <= 0 || any(!is.finite(curvature))) return(invalid)
  D <- rowsum(data$X * curvature, data$subject, reorder = FALSE) %*% K
  R <- tryCatch(chol(-geometry$B), error = function(e) NULL)
  if (is.null(R)) return(invalid)
  whitened <- forwardsolve(t(R), t(D))
  steps <- backsolve(R, whitened)
  score_error <- sum(sqrt(colSums(whitened^2)) * abs(correction))
  step_error <- max(as.numeric(abs(steps) %*% abs(correction))) /
    (1 + max(abs(geometry$gamma)))
  component_error <- as.numeric(crossprod(abs(D), abs(correction)))
  if (any(!is.finite(c(score_error, step_error, component_error)))) return(invalid)
  list(success = TRUE, score_error = score_error, step_error = step_error,
       component_error = component_error, correction = correction)
}

hl_fo_polish_inner_modes <- function(data, state, rho_full, K, beta_star,
                                     score_budget, step_budget,
                                     max_polish = 3L, geometry = NULL) {
  stopifnot(length(max_polish) == 1L, is.finite(max_polish),
            max_polish >= 0, max_polish == as.integer(max_polish),
            is.finite(score_budget), score_budget > 0,
            is.finite(step_budget), step_budget > 0)
  for (round in 0:as.integer(max_polish)) {
    if (!isTRUE(state$success)) return(state)
    if (round > 0L || is.null(geometry)) {
      geometry <- hl_fo_inner_geometry(state, K, beta_star)
    }
    if (!isTRUE(geometry$success)) return(geometry)
    accuracy <- hl_fo_inner_mode_accuracy(data, state, rho_full, K, geometry)
    failure <- function(reason) list(
      success = FALSE, status = "inner_no_convergence", reason = reason,
      beta = state$beta, eta = state$eta, mode_accuracy = accuracy,
      mode_polish_rounds = round, score_metric = geometry$score_metric,
      step_metric = geometry$step_metric
    )
    if (!isTRUE(accuracy$success)) return(failure(accuracy$reason))
    accuracy$rounds <- round
    if (accuracy$score_error <= score_budget && accuracy$step_error <= step_budget) {
      state$mode_accuracy <- accuracy
      return(list(success = TRUE, state = state, geometry = geometry))
    }
    if (round == max_polish) return(failure("inner_mode_accuracy_exhausted"))
    eta <- state$eta + accuracy$correction
    if (any(!is.finite(eta)) || identical(eta, state$eta)) {
      return(failure("inner_mode_accuracy_stalled"))
    }
    # Current donor roots already passed the original safeguarded solver.
    # Reevaluate every donor after this Newton polish; do not correct the beta
    # score algebraically or accept a linearized score as an actual solution.
    state <- hl_fo_inner_quantities(data, state$beta, rho_full, eta_start = eta)
  }
}

hl_fo_inner_armijo <- function(current, trial, K, displacement, midpoint = NULL) {
  reject <- list(accepted = FALSE, needs_integral = FALSE,
                 method = NA_character_, integrated_gain = NA_real_)
  if (!isTRUE(current$success) || !is.finite(current$P) ||
      !isTRUE(trial$success) || !is.finite(trial$P) ||
      any(!is.finite(displacement)) || !any(displacement != 0)) return(reject)
  direction <- function(s) sum(as.numeric(crossprod(K, s$S_beta)) * displacement)
  directional <- direction(current)
  if (!is.finite(directional) || directional <= 0) return(reject)
  gain <- trial$P - current$P
  cancellation <- 64 * .Machine$double.eps * (1 + abs(current$P) + abs(trial$P))
  if (abs(gain) > cancellation) {
    reject$accepted <- gain >= 1e-4 * directional
    reject$method <- "direct_profile_criterion"
    return(reject)
  }
  # Even a positive rounded difference is unreliable within this band.
  reject$needs_integral <- TRUE
  if (is.null(midpoint) || !isTRUE(midpoint$success)) return(reject)
  states <- list(current, midpoint, trial)
  if (!all(vapply(states, function(s) isTRUE(s$mode_accuracy$success), logical(1)))) {
    return(reject)
  }
  directions <- vapply(states, direction, numeric(1))
  errors <- vapply(states, function(s)
    sum(abs(displacement) * s$mode_accuracy$component_error), numeric(1))
  gain <- sum(c(1, 4, 1) * directions) / 6
  error <- sum(c(1, 4, 1) * errors) / 6
  reject$accepted <- all(is.finite(c(gain, error))) &&
    gain - error >= 1e-4 * (directional + errors[1L])
  reject$method <- "profile_score_integral"
  reject$integrated_gain <- gain
  reject$integrated_error <- error
  reject
}
