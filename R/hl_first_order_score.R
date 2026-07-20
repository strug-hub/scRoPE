# Separate full-chi null solver for the first-order nested HL Score test.
#
# The reduced constrained maximizer is deliberately kept separate from this
# algebraic nuisance root.  At every proposed active-dispersion vector, beta
# is restricted to beta_star + K gamma, all modes and the reduced ordinary
# h-score root are re-solved, and the full p-dimensional auxiliary equation is
# solved before g0 is evaluated.

.hl_fo_score_root_failure <- function(status, reason = NULL,
                                      contrast_rank = NA_integer_,
                                      m = NA_integer_, extras = list()) {
  c(list(
    method = "score_ee",
    success = FALSE,
    status = status,
    reason = if (is.null(reason)) status else reason,
    stat = NA_real_,
    statistic = NA_real_,
    df = as.integer(contrast_rank),
    p = NA_real_,
    p_value = NA_real_,
    m = as.integer(m),
    contrast_rank = as.integer(contrast_rank)
  ), extras)
}

.hl_fo_score_active <- function(active) {
  out <- if (is.character(active)) match(active, c("tau", "xi")) else
    as.integer(active)
  if (anyNA(out) || any(!out %in% 1:2) || anyDuplicated(out) ||
      is.unsorted(out)) {
    stop("active must be an ordered subset of tau and xi.")
  }
  as.integer(out)
}

.hl_fo_score_rho_template <- function(seed, fixed_natural) {
  rho <- as.numeric(seed$rho_full)
  if (length(rho) != 2L || any(!is.finite(rho))) {
    stop("The Score seed must contain two finite canonical dispersions.")
  }
  if (length(fixed_natural)) {
    fixed_names <- names(fixed_natural)
    fixed_idx <- match(fixed_names, c("sigma2", "phi"))
    numeric_names <- suppressWarnings(as.integer(fixed_names))
    fixed_idx[is.na(fixed_idx)] <- numeric_names[is.na(fixed_idx)]
    if (anyNA(fixed_idx) || any(!fixed_idx %in% 1:2) ||
        any(!is.finite(fixed_natural)) || any(fixed_natural <= 0)) {
      stop("fixed_natural must contain positive sigma2/phi values.")
    }
    rho[fixed_idx] <- log(as.numeric(fixed_natural))
  }
  rho
}

# Evaluate the restricted-null full-chi system at one active rho vector.
# This function is also the immediately lower-order analytic object used by
# the finite-difference validation of the root Jacobian.
hl_fo_score_state_at_rho <- function(data, rho_active, rho_template, active,
                                     beta_star, K, beta_start,
                                     eta_start = NULL, cold = FALSE,
                                     include_laplace_constant = TRUE,
                                     constraint_R = NULL,
                                     constraint_rhs = numeric(0),
                                     inner_acceptance_tier = "derivative") {
  active <- as.integer(active)
  rho_full <- as.numeric(rho_template)
  if (length(active)) rho_full[active] <- as.numeric(rho_active)
  p <- data$p
  q <- length(active)
  k <- ncol(K)

  inner <- hl_fo_solve_inner(
    data = data,
    rho_full = rho_full,
    beta_star = beta_star,
    K = K,
    beta_start = beta_start,
    eta_start = eta_start,
    cold = cold,
    acceptance_tier = inner_acceptance_tier
  )
  if (!isTRUE(inner$success)) return(inner)

  # Full beta and active-rho coordinates are required here.  Contracting with
  # K would incorrectly replace the full auxiliary equation by the reduced
  # constrained-maximizer equation.
  transform <- hl_fo_scientific_transform(p, diag(p), active)
  subjects <- vector("list", data$m)
  for (ii in seq_len(data$m)) {
    rows <- data$subject_rows[[ii]]
    raw <- hl_fo_subject_raw(
      y = data$y[rows],
      X = data$X[rows, , drop = FALSE],
      log_offset = data$log_offset[rows],
      beta = inner$beta,
      eta = inner$eta[ii],
      rho_full = rho_full,
      include_laplace_constant = include_laplace_constant
    )
    if (!isTRUE(raw$ok)) {
      return(list(success = FALSE, status = "mode_failure",
                  reason = raw$reason, failed_subject = ii))
    }
    subjects[[ii]] <- hl_fo_contract_subject(raw, transform)
  }
  augmented <- hl_fo_augmented_state(subjects, k_free = p, q = q)
  if (!isTRUE(augmented$success)) return(augmented)

  B <- augmented$B
  E <- augmented$E
  R_aux <- augmented$R
  U_aux <- augmented$U
  V_aux <- augmented$V
  beta_rho <- matrix(0, nrow = p, ncol = q)
  if (k && q) {
    B_reduced <- crossprod(K, B %*% K)
    E_reduced <- crossprod(K, E)
    gamma_rho <- tryCatch(
      solve(B_reduced, E_reduced),
      error = function(e) NULL
    )
    if (is.null(gamma_rho) || any(!is.finite(gamma_rho))) {
      return(list(success = FALSE, status = "sensitivity_singular",
                  reason = "score_reduced_beta_rho_solve_failed"))
    }
    beta_rho <- -K %*% gamma_rho
  }
  chi_rho <- matrix(0, nrow = p, ncol = q)
  if (q) {
    chi_rho <- tryCatch(
      -solve(t(B), R_aux %*% beta_rho + U_aux),
      error = function(e) NULL
    )
    if (is.null(chi_rho) || any(!is.finite(chi_rho))) {
      return(list(success = FALSE, status = "sensitivity_singular",
                  reason = "score_full_chi_rho_solve_failed"))
    }
  }
  G0 <- if (q) {
    V_aux - crossprod(U_aux, beta_rho) - crossprod(E, chi_rho)
  } else {
    matrix(0, 0L, 0L)
  }

  psi3 <- if (q) {
    augmented$psi_subject[2L * p + seq_len(q), , drop = FALSE]
  } else {
    matrix(0, nrow = 0L, ncol = data$m)
  }
  scale <- if (q) 1 + sqrt(rowSums(psi3^2)) else numeric(0)
  g0 <- as.numeric(augmented$g)
  scaled_root <- if (q) max(abs(g0) / scale) else 0
  constraint <- if (is.null(constraint_R)) NULL else
    hl_fo_constraint_residual(constraint_R, inner$beta, constraint_rhs)

  list(
    success = TRUE,
    status = "ok",
    reason = NULL,
    root_type = "full_chi_score_null",
    data = data,
    rho = as.numeric(rho_active),
    rho_full = rho_full,
    active = active,
    beta_star = as.numeric(beta_star),
    K = as.matrix(K),
    gamma = inner$gamma,
    beta = inner$beta,
    eta = inner$eta,
    inner = inner,
    subjects = subjects,
    augmented = augmented,
    B = B,
    E = E,
    chi = augmented$chi,
    g0 = g0,
    G0 = G0,
    beta_rho = beta_rho,
    chi_rho = chi_rho,
    root_scale = scale,
    root_scaled_residual = scaled_root,
    constraint = constraint,
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
}

hl_fo_score_make_evaluator <- function(data, rho_template, active, beta_star, K,
                                       seed_beta, include_laplace_constant = TRUE,
                                       constraint_R, constraint_rhs) {
  cache <- list()
  trace <- list()
  evaluate <- function(rho, cold = FALSE) {
    nearest <- NULL
    if (!isTRUE(cold) && length(cache)) {
      distances <- vapply(cache, function(old) {
        hl_fo_cache_distance(rho, old$rho)
      }, numeric(1))
      nearest <- cache[[which.min(distances)]]
    }
    state <- hl_fo_score_state_at_rho(
      data = data,
      rho_active = rho,
      rho_template = rho_template,
      active = active,
      beta_star = beta_star,
      K = K,
      beta_start = if (is.null(nearest)) seed_beta else nearest$beta,
      eta_start = if (is.null(nearest)) NULL else nearest$eta,
      cold = cold,
      include_laplace_constant = include_laplace_constant,
      constraint_R = constraint_R,
      constraint_rhs = constraint_rhs,
      inner_acceptance_tier = if (length(active)) "derivative" else "statistical"
    )
    trace[[length(trace) + 1L]] <<- list(
      rho = as.numeric(rho),
      cold = isTRUE(cold),
      success = isTRUE(state$success),
      status = state$status,
      reason = state$reason,
      root_scaled_residual = if (isTRUE(state$success))
        state$root_scaled_residual else NA_real_,
      constraint_residual = if (isTRUE(state$success))
        state$constraint$residual else NULL
    )
    if (isTRUE(state$success)) {
      # The root search only needs the nearest warm starts.  Retaining every
      # subject derivative/Jacobian state here scales poorly across ordered
      # multistarts and can exhaust memory for large subject counts.
      cache[[length(cache) + 1L]] <<- list(
        rho = state$rho,
        beta = state$beta,
        eta = state$eta
      )
    }
    state
  }
  list(
    evaluate = evaluate,
    count = function() length(trace),
    trace = function() trace,
    cache = function() cache
  )
}

hl_fo_score_optimize_start <- function(start, lower, upper, evaluator,
                                       root_tolerance = 1e-8) {
  radius <- 1
  invalid_consecutive <- 0L
  start_count <- evaluator$count()
  current <- evaluator$evaluate(start)
  if (!isTRUE(current$success)) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                reason = current$status, state = current, start = start,
                evaluations = evaluator$count() - start_count,
                trust_trace = list()))
  }
  trust_trace <- list()
  repeat {
    if (current$root_scaled_residual <= root_tolerance) {
      return(list(success = TRUE, status = "ok", state = current,
                  start = start,
                  evaluations = evaluator$count() - start_count,
                  trust_trace = trust_trace))
    }
    if (evaluator$count() - start_count >= 200L) {
      return(list(success = FALSE, status = "score_null_no_convergence",
                  reason = "score_root_evaluation_cap", state = current,
                  start = start,
                  evaluations = evaluator$count() - start_count,
                  trust_trace = trust_trace))
    }
    if (radius <= 1e-12 * max(1, sqrt(sum(current$rho^2)))) {
      return(list(success = FALSE, status = "score_null_no_convergence",
                  reason = "score_root_radius_collapse", state = current,
                  start = start,
                  evaluations = evaluator$count() - start_count,
                  trust_trace = trust_trace))
    }

    residual <- current$g0 / current$root_scale
    jacobian <- current$G0 / current$root_scale
    gradient <- as.numeric(crossprod(jacobian, residual))
    curvature <- crossprod(jacobian)
    step <- hl_fo_trust_step_max(-gradient, -curvature, radius)
    predicted <- if (all(is.finite(step))) {
      -sum(gradient * step) -
        0.5 * sum(step * (curvature %*% step))
    } else NA_real_
    entry <- list(
      radius = radius,
      rho = current$rho,
      residual = residual,
      step = step,
      predicted_reduction = predicted
    )
    if (!is.finite(predicted) || predicted <= 0) {
      radius <- radius / 4
      entry$accepted <- FALSE
      entry$reason <- "nonpositive_prediction"
      trust_trace[[length(trust_trace) + 1L]] <- entry
      next
    }

    proposal <- current$rho + step
    if (any(proposal <= lower) || any(proposal >= upper)) {
      radius <- radius / 4
      invalid_consecutive <- invalid_consecutive + 1L
      entry$accepted <- FALSE
      entry$reason <- "outside_guard"
      trust_trace[[length(trust_trace) + 1L]] <- entry
      if (invalid_consecutive >= 5L) {
        return(list(success = FALSE, status = "score_null_no_convergence",
                    reason = "five_invalid_trials", state = current,
                    start = start,
                    evaluations = evaluator$count() - start_count,
                    trust_trace = trust_trace))
      }
      next
    }

    trial <- evaluator$evaluate(proposal)
    if (!isTRUE(trial$success)) {
      radius <- radius / 4
      invalid_consecutive <- invalid_consecutive + 1L
      entry$accepted <- FALSE
      entry$reason <- trial$status
      trust_trace[[length(trust_trace) + 1L]] <- entry
      if (invalid_consecutive >= 5L) {
        return(list(success = FALSE, status = "score_null_no_convergence",
                    reason = "five_invalid_trials", state = current,
                    start = start,
                    evaluations = evaluator$count() - start_count,
                    trust_trace = trust_trace))
      }
      next
    }
    invalid_consecutive <- 0L
    # D0 is held fixed over this trust iteration, as required by the root
    # algorithm.  The trial's own scale is used only after acceptance.
    trial_residual_fixed_scale <- trial$g0 / current$root_scale
    actual <- 0.5 * sum(residual^2) -
      0.5 * sum(trial_residual_fixed_scale^2)
    ratio <- actual / predicted
    entry$actual_reduction <- actual
    entry$ratio <- ratio
    if (!is.finite(ratio) || ratio < 0.1) {
      radius <- radius / 4
      entry$accepted <- FALSE
      entry$reason <- "insufficient_reduction"
    } else {
      current <- trial
      entry$accepted <- TRUE
      entry$reason <- "accepted"
      if (ratio > 0.75 && sqrt(sum(step^2)) >= 0.9 * radius) {
        radius <- min(100, 2 * radius)
      }
    }
    trust_trace[[length(trust_trace) + 1L]] <- entry
  }
}

hl_fo_score_cluster_candidates <- function(candidates, tolerance = 1e-7) {
  kept <- list()
  for (candidate in candidates) {
    duplicate <- any(vapply(kept, function(old) {
      hl_fo_cache_distance(candidate$state$rho, old$state$rho) <= tolerance &&
        hl_fo_sdisc(candidate$state$beta, old$state$beta) <= tolerance
    }, logical(1)))
    if (!duplicate) kept[[length(kept) + 1L]] <- candidate
  }
  kept
}

hl_fo_score_select_candidates <- function(candidates) {
  distinct <- hl_fo_score_cluster_candidates(candidates, tolerance = 1e-7)
  if (!length(distinct)) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                candidates = distinct))
  }
  if (length(distinct) > 1L) {
    return(list(success = FALSE, status = "score_null_multiple_roots",
                candidates = distinct))
  }
  list(success = TRUE, status = "ok", selected = distinct[[1L]],
       candidates = distinct)
}

hl_fo_score_certify_candidate <- function(candidate, evaluator, R, rhs,
                                          C_R, K, root_tolerance = 1e-8) {
  warm <- candidate$state
  cold <- evaluator$evaluate(warm$rho, cold = TRUE)
  if (!isTRUE(cold$success)) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                reason = "score_cold_replay_failure", cold = cold))
  }
  cold_discrepancy <- list(
    beta = hl_fo_sdisc(cold$beta, warm$beta),
    eta = hl_fo_sdisc(cold$eta, warm$eta),
    g0 = hl_fo_sdisc(cold$g0, warm$g0)
  )
  if (cold_discrepancy$beta > 1e-8 || cold_discrepancy$eta > 1e-8 ||
      cold_discrepancy$g0 > 1e-8 ||
      cold$root_scaled_residual > root_tolerance) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                reason = "score_cold_replay_mismatch", state = cold,
                cold_discrepancy = cold_discrepancy))
  }

  constraint <- hl_fo_constraint_residual(R, cold$beta, rhs)
  if (!isTRUE(constraint$satisfied)) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                reason = "score_null_constraint_failure", state = cold,
                constraint = constraint))
  }
  mode_residual <- max(vapply(cold$inner$mode_diagnostics, function(x) {
    max(x$residual_scaled, x$step_scaled)
  }, numeric(1)))
  if (!is.finite(mode_residual) || mode_residual > 1e-10 ||
      cold$inner$score_metric > 1e-8 || cold$inner$step_metric > 1e-10) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                reason = "score_null_inner_not_certified", state = cold,
                mode_residual = mode_residual))
  }

  aux_residual <- as.numeric(
    t(cold$augmented$B) %*% cold$augmented$chi -
      cold$augmented$L_gamma
  )
  aux_backward_error <- max(abs(aux_residual)) /
    (1 + max(abs(t(cold$augmented$B) %*% cold$augmented$chi)) +
       max(abs(cold$augmented$L_gamma)))
  if (!is.finite(aux_backward_error) || aux_backward_error > 1e-10) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                reason = "score_null_auxiliary_backward_error", state = cold,
                auxiliary_backward_error = aux_backward_error))
  }

  p <- cold$data$p
  q <- length(cold$active)
  r <- nrow(R)
  k <- ncol(K)
  Q <- cbind(C_R, K)
  d <- 2L * p + q
  T <- matrix(0, d, d)
  T[seq_len(p), seq_len(p)] <- Q
  T[p + seq_len(p), p + seq_len(p)] <- Q
  if (q) {
    rho_idx <- 2L * p + seq_len(q)
    T[rho_idx, rho_idx] <- diag(q)
  }
  u_subject <- crossprod(T, cold$augmented$psi_subject)
  A <- crossprod(T, cold$augmented$A %*% T)
  beta_nuisance_idx <- if (k) r + seq_len(k) else integer(0)
  chi_idx <- p + seq_len(p)
  rho_idx <- if (q) 2L * p + seq_len(q) else integer(0)
  nuisance_idx <- c(beta_nuisance_idx, chi_idx, rho_idx)
  A_nn <- A[nuisance_idx, nuisance_idx, drop = FALSE]
  u_nuisance <- u_subject[nuisance_idx, , drop = FALSE]
  J_nn <- tcrossprod(u_nuisance)
  nuisance_rank <- hl_fo_num_rank(A_nn)$rank
  if (nuisance_rank < length(nuisance_idx)) {
    return(list(success = FALSE, status = "score_nuisance_singular",
                reason = "score_nuisance_singular", state = cold,
                nuisance_rank = nuisance_rank))
  }
  left <- tryCatch(solve(A_nn, J_nn), error = function(e) NULL)
  V_nuisance <- if (!is.null(left)) tryCatch(
    t(solve(A_nn, t(left))), error = function(e) NULL
  ) else NULL
  U_nuisance <- rowSums(u_nuisance)
  delta_nuisance <- tryCatch(solve(A_nn, U_nuisance), error = function(e) NULL)
  if (is.null(V_nuisance) || is.null(delta_nuisance) ||
      any(!is.finite(c(V_nuisance, delta_nuisance))) ||
      any(diag(V_nuisance) <= 0)) {
    return(list(success = FALSE, status = "score_nuisance_singular",
                reason = "score_nuisance_covariance_invalid", state = cold))
  }
  statistical_displacement <- abs(delta_nuisance) /
    sqrt(diag(V_nuisance))
  displacement_limit <- min(1e-3, cold$data$m^(-1 / 4))
  if (max(statistical_displacement) > displacement_limit) {
    return(list(success = FALSE, status = "score_null_no_convergence",
                reason = "score_nuisance_statistical_root_not_certified",
                state = cold,
                statistical_displacement = statistical_displacement,
                statistical_displacement_limit = displacement_limit))
  }

  cold$root_type <- "full_chi_score_null"
  cold$certification <- list(
    cold_discrepancy = cold_discrepancy,
    constraint = constraint,
    mode_residual = mode_residual,
    reduced_inner_score_metric = cold$inner$score_metric,
    reduced_inner_step_metric = cold$inner$step_metric,
    auxiliary_backward_error = aux_backward_error,
    root_scaled_residual = cold$root_scaled_residual,
    nuisance_rank = nuisance_rank,
    nuisance_indices = nuisance_idx,
    nuisance_sensitivity = A_nn,
    nuisance_meat = J_nn,
    nuisance_covariance = V_nuisance,
    nuisance_residual = U_nuisance,
    nuisance_displacement = delta_nuisance,
    statistical_displacement = statistical_displacement,
    statistical_displacement_limit = displacement_limit,
    transformed_subject = u_subject,
    transformed_sensitivity = A,
    T = T,
    Q = Q,
    C_R = C_R,
    K = K
  )
  list(success = TRUE, status = "ok", state = cold,
       start = candidate$start, trust_trace = candidate$trust_trace,
       evaluations = candidate$evaluations)
}

# Fit the full-chi Score nuisance root.  The returned constrained_fit is the
# reduced constrained maximizer and is retained only as a separate starting
# object; null_fit is the certified algebraic full-chi root.
hl_fo_score_null_fit <- function(data, R, rhs = NULL,
                                 unrestricted_fit = NULL,
                                 constrained_fit = NULL,
                                 active = c("tau", "xi"),
                                 fixed_natural = numeric(0),
                                 guard_min = c(1e-4, 1e-4),
                                 guard_max = c(10, 1000),
                                 include_laplace_constant = TRUE) {
  reduced <- .hl_fo_reduce_contrast(R, rhs)
  if (!isTRUE(reduced$ok)) {
    return(.hl_fo_score_root_failure(
      status = paste0("score_", reduced$status),
      reason = reduced$status,
      contrast_rank = reduced$rank,
      m = data$m,
      extras = list(contrast = reduced)
    ))
  }
  R_use <- reduced$R
  rhs_use <- reduced$rhs
  r <- reduced$rank
  active_idx <- .hl_fo_score_active(active)
  q <- length(active_idx)
  affine <- hl_fo_affine_constraint(R_use, rhs_use, p = data$p)
  if (!isTRUE(affine$success)) {
    return(.hl_fo_score_root_failure(
      status = "score_contrast_rhs_inconsistent",
      contrast_rank = r,
      m = data$m,
      extras = list(contrast = reduced, constraint = affine)
    ))
  }
  C_R <- t(R_use) %*% solve(tcrossprod(R_use))

  seed <- NULL
  if (is.list(unrestricted_fit) &&
      length(unrestricted_fit$beta) == data$p &&
      length(unrestricted_fit$rho_full) == 2L &&
      all(is.finite(c(unrestricted_fit$beta, unrestricted_fit$rho_full)))) {
    seed <- list(
      beta = as.numeric(unrestricted_fit$beta),
      rho_full = as.numeric(unrestricted_fit$rho_full),
      coordinate_tag = "canonical_conditional_beta_tau_xi",
      source = "unrestricted_fit_start_only"
    )
  }
  if (is.null(seed)) seed <- hl_fo_poisson_seed(data)
  rho_template <- .hl_fo_score_rho_template(seed, fixed_natural)
  seed$rho_full <- rho_template

  if (is.null(constrained_fit)) {
    constrained_fit <- hl_fo_fit(
      data = data,
      C = R_use,
      rhs = rhs_use,
      active = active_idx,
      fixed_natural = fixed_natural,
      guard_min = guard_min,
      guard_max = guard_max,
      seed = seed,
      include_laplace_constant = include_laplace_constant
    )
  }
  constrained_fit$object_type <- "reduced_constrained_maximizer"
  if (isTRUE(constrained_fit$point_success) &&
      length(constrained_fit$rho_full) == 2L &&
      all(is.finite(constrained_fit$rho_full))) {
    rho_template <- constrained_fit$rho_full
  }

  lower <- log(guard_min)[active_idx]
  upper <- log(guard_max)[active_idx]
  if (q && isTRUE(constrained_fit$point_success) &&
      length(constrained_fit$guard_bounds$lower) == q &&
      length(constrained_fit$guard_bounds$upper) == q) {
    lower <- constrained_fit$guard_bounds$lower
    upper <- constrained_fit$guard_bounds$upper
  }
  primary <- if (q && isTRUE(constrained_fit$point_success) &&
                 length(constrained_fit$rho) == q &&
                 all(is.finite(constrained_fit$rho))) {
    constrained_fit$rho
  } else {
    rho_template[active_idx]
  }
  if (q && (any(primary <= lower) || any(primary >= upper))) {
    primary <- (lower + upper) / 2
  }

  evaluator <- hl_fo_score_make_evaluator(
    data = data,
    rho_template = rho_template,
    active = active_idx,
    beta_star = affine$beta_star,
    K = affine$K,
    seed_beta = seed$beta,
    include_laplace_constant = include_laplace_constant,
    constraint_R = R_use,
    constraint_rhs = rhs_use
  )
  starts <- hl_fo_outer_starts(primary, lower, upper)
  runs <- lapply(starts, hl_fo_score_optimize_start,
                 lower = lower, upper = upper, evaluator = evaluator)
  algebraic <- Filter(function(x) isTRUE(x$success), runs)
  certification_attempts <- lapply(
    algebraic, hl_fo_score_certify_candidate,
    evaluator = evaluator, R = R_use, rhs = rhs_use,
    C_R = C_R, K = affine$K
  )
  certified <- Filter(function(x) isTRUE(x$success), certification_attempts)
  selection <- hl_fo_score_select_candidates(certified)
  if (!isTRUE(selection$success)) {
    certification_statuses <- vapply(
      certification_attempts,
      function(x) if (is.null(x$status)) "" else x$status,
      character(1)
    )
    failure_status <- if (!length(certified) &&
                          any(certification_statuses == "score_nuisance_singular")) {
      "score_nuisance_singular"
    } else {
      selection$status
    }
    return(.hl_fo_score_root_failure(
      status = failure_status,
      contrast_rank = r,
      m = data$m,
      extras = list(
        contrast = reduced,
        constraint = affine,
        constrained_fit = constrained_fit,
        constrained_fit_status = constrained_fit$fit_status,
        null_fit = NULL,
        null_runs = runs,
        certification_attempts = certification_attempts,
        starts = starts,
        root_trace = evaluator$trace(),
        root_candidates = selection$candidates
      )
    ))
  }

  null_fit <- selection$selected$state
  null_fit$object_type <- "full_chi_score_null"
  list(
    success = TRUE,
    status = "ok",
    contrast_rank = r,
    m = data$m,
    contrast = reduced,
    constraint = affine,
    constrained_fit = constrained_fit,
    constrained_fit_status = constrained_fit$fit_status,
    null_fit = null_fit,
    null_fit_status = "ok",
    null_runs = runs,
    starts = starts,
    root_trace = evaluator$trace(),
    root_candidates = selection$candidates
  )
}

hl_fo_score_test <- function(data, R, rhs = NULL,
                             unrestricted_fit = NULL,
                             guard_min = c(1e-4, 1e-4),
                             guard_max = c(10, 1000),
                             active = c("tau", "xi"),
                             fixed_natural = numeric(0),
                             constrained_fit = NULL,
                             include_laplace_constant = TRUE) {
  null <- hl_fo_score_null_fit(
    data = data,
    R = R,
    rhs = rhs,
    unrestricted_fit = unrestricted_fit,
    constrained_fit = constrained_fit,
    active = active,
    fixed_natural = fixed_natural,
    guard_min = guard_min,
    guard_max = guard_max,
    include_laplace_constant = include_laplace_constant
  )
  if (!isTRUE(null$success)) return(null)

  state <- null$null_fit
  post <- hl_fo_score_postroot(
    psi_subject = state$augmented$psi_subject,
    A_zeta = state$augmented$A,
    R = null$contrast$R,
    beta = state$beta,
    rhs = null$contrast$rhs,
    q = length(state$active)
  )
  if (identical(post$status, "score_adjusted_sensitivity_singular")) {
    post$algebra_status <- post$status
    post$status <- "score_nuisance_singular"
  }
  post$statistic <- post$stat
  post$p_value <- post$p
  post$constrained_fit <- null$constrained_fit
  post$constrained_fit_status <- null$constrained_fit_status
  post$null_fit <- state
  post$null_fit_status <- null$null_fit_status
  post$null_runs <- null$null_runs
  post$starts <- null$starts
  post$root_trace <- null$root_trace
  post$root_candidates <- null$root_candidates
  post$root_certification <- state$certification
  post
}
