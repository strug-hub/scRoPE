# Certified nested solvers for the smooth first-order HL estimator.

hl_fo_matrix_rank <- function(x) {
  x <- as.matrix(x)
  if (!length(x) || !nrow(x) || !ncol(x)) {
    return(0L)
  }
  ss <- svd(x, nu = 0, nv = 0)$d
  if (!length(ss) || ss[1] == 0) {
    return(0L)
  }
  as.integer(sum(ss > max(dim(x)) * .Machine$double.eps * ss[1]))
}

hl_fo_make_data <- function(y, X, log_offset, subject, intcol = 1L) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  log_offset <- as.numeric(log_offset)
  subject <- as.integer(factor(subject, levels = unique(subject)))
  n <- length(y)
  if (nrow(X) != n || length(log_offset) != n || length(subject) != n) {
    stop("y, X, log_offset, and subject must describe the same cells.")
  }
  if (any(!is.finite(y)) || any(y < 0) || any(abs(y - round(y)) > 0)) {
    stop("Counts must be finite nonnegative integers.")
  }
  if (any(!is.finite(X)) || any(!is.finite(log_offset))) {
    stop("Design and log offsets must be finite.")
  }
  if (intcol < 1L || intcol > ncol(X)) {
    stop("intcol is outside the design.")
  }
  rows <- split(seq_len(n), subject)
  if (any(lengths(rows) == 0L)) {
    stop("Every subject must contribute at least one cell.")
  }
  list(
    y = y,
    X = X,
    log_offset = log_offset,
    subject = subject,
    subject_rows = unname(rows),
    n = n,
    p = ncol(X),
    m = length(rows),
    intcol = as.integer(intcol),
    coordinate_order = c(paste0("beta", seq_len(ncol(X))), "tau", "xi")
  )
}

hl_fo_solve_mode <- function(y, log_m_base, tau, xi, start = NULL,
                             max_expand = 60L, max_iter = 100L) {
  v <- exp(tau)
  if (is.null(start) || length(start) != 1L || !is.finite(start)) {
    start <- v / 2
    start_type <- "cold_prior_mode"
  } else {
    start_type <- "warm"
  }
  eval_mode <- function(z) hl_fo_mode_components(y, log_m_base, z, tau, xi)
  first <- eval_mode(start)
  if (!isTRUE(first$ok)) {
    return(list(success = FALSE, status = "mode_failure", reason = first$reason,
                start_type = start_type))
  }

  width0 <- max(1, 0.25 * (1 + abs(start)))
  bracket <- NULL
  expansion <- NA_integer_
  for (ee in 0:max_expand) {
    width <- width0 * 2^ee
    lo <- start - width
    hi <- start + width
    flo <- eval_mode(lo)
    fhi <- eval_mode(hi)
    if (!isTRUE(flo$ok) || !isTRUE(fhi$ok)) {
      next
    }
    if (flo$h_eta > 0 && fhi$h_eta < 0) {
      bracket <- c(lo, hi)
      expansion <- ee
      break
    }
  }
  if (is.null(bracket)) {
    return(list(success = FALSE, status = "mode_failure",
                reason = "mode_bracket_failure", start_type = start_type,
                bracket_expansions = max_expand))
  }

  lo <- bracket[1]
  hi <- bracket[2]
  eta <- min(max(start, lo), hi)
  bisections <- 0L
  last <- NULL
  for (iter in seq_len(max_iter)) {
    cur <- eval_mode(eta)
    if (!isTRUE(cur$ok)) {
      return(list(success = FALSE, status = "mode_failure", reason = cur$reason,
                  iterations = iter, bracket_expansions = expansion))
    }
    a <- cur$h_etaeta
    if (!is.finite(a) || a >= 0) {
      return(list(success = FALSE, status = "mode_failure",
                  reason = "mode_nonnegative_curvature", eta = eta,
                  h_eta = cur$h_eta, h_etaeta = a, iterations = iter))
    }
    residual_scaled <- abs(cur$h_eta) / sqrt(-a)
    step_scaled <- abs(cur$h_eta / a) / (1 + abs(eta))
    if (is.finite(residual_scaled) && is.finite(step_scaled) &&
        residual_scaled <= 1e-10 && step_scaled <= 1e-10) {
      return(list(
        success = TRUE,
        status = "ok",
        reason = NULL,
        eta = eta,
        h_eta = cur$h_eta,
        h_etaeta = a,
        residual_scaled = residual_scaled,
        step_scaled = step_scaled,
        iterations = iter,
        bracket_expansions = expansion,
        bisections = bisections,
        bracket = c(lo, hi),
        start_type = start_type
      ))
    }

    if (cur$h_eta > 0) {
      lo <- eta
    } else {
      hi <- eta
    }
    proposal <- eta - cur$h_eta / a
    if (!is.finite(proposal) || proposal <= lo || proposal >= hi) {
      proposal <- (lo + hi) / 2
      bisections <- bisections + 1L
    }
    last <- cur
    eta <- proposal
  }
  list(
    success = FALSE,
    status = "mode_failure",
    reason = "mode_iteration_exhausted",
    eta = eta,
    h_eta = if (is.null(last)) NA_real_ else last$h_eta,
    h_etaeta = if (is.null(last)) NA_real_ else last$h_etaeta,
    iterations = max_iter,
    bracket_expansions = expansion,
    bisections = bisections,
    bracket = c(lo, hi),
    start_type = start_type
  )
}

hl_fo_mode_vector_eval <- function(data, log_m_base, eta, tau, xi) {
  phi <- exp(xi)
  v <- exp(tau)
  alpha <- 1 / expm1(v)
  d_prior <- eta - v / 2
  lambda_z <- alpha * exp(d_prior)
  log_m <- log_m_base + eta[data$subject]
  r <- stats::plogis(log_m - xi)
  omr <- stats::plogis(xi - log_m)
  s <- data$y * omr - phi * r
  tt <- -(data$y + phi) * r * omr
  if (any(!is.finite(c(s, tt, lambda_z)))) {
    return(list(ok = FALSE, reason = "nonfinite_mode_arithmetic"))
  }
  list(
    ok = TRUE,
    h_eta = as.numeric(rowsum(s, data$subject, reorder = FALSE)) -
      alpha * expm1(d_prior),
    h_etaeta = as.numeric(rowsum(tt, data$subject, reorder = FALSE)) - lambda_z
  )
}

hl_fo_mode_subject_eval <- function(data, log_m_base, subject_index, eta,
                                    tau, xi) {
  rows <- data$subject_rows[[subject_index]]
  phi <- exp(xi)
  v <- exp(tau)
  alpha <- 1 / expm1(v)
  d_prior <- eta - v / 2
  lambda_z <- alpha * exp(d_prior)
  log_m <- log_m_base[rows] + eta
  r <- stats::plogis(log_m - xi)
  omr <- stats::plogis(xi - log_m)
  s <- data$y[rows] * omr - phi * r
  tt <- -(data$y[rows] + phi) * r * omr
  if (any(!is.finite(c(s, tt, lambda_z)))) {
    return(list(ok = FALSE, reason = "nonfinite_mode_arithmetic"))
  }
  list(
    ok = TRUE,
    h_eta = sum(s) - alpha * expm1(d_prior),
    h_etaeta = sum(tt) - lambda_z
  )
}

hl_fo_solve_modes <- function(data, beta, rho_full, eta_start = NULL) {
  base <- as.numeric(data$log_offset + data$X %*% beta)
  v <- exp(rho_full[1])
  warm <- !is.null(eta_start) && length(eta_start) == data$m &&
    all(is.finite(eta_start))
  start <- if (warm) as.numeric(eta_start) else rep(v / 2, data$m)
  width0 <- pmax(1, 0.25 * (1 + abs(start)))
  lo <- hi <- rep(NA_real_, data$m)
  expansions <- rep(NA_integer_, data$m)
  unresolved <- rep(TRUE, data$m)
  for (ee in 0:60) {
    width <- width0 * 2^ee
    for (ii in which(unresolved)) {
      cand_lo <- start[ii] - width[ii]
      cand_hi <- start[ii] + width[ii]
      vlo <- hl_fo_mode_subject_eval(
        data, base, ii, cand_lo, rho_full[1], rho_full[2]
      )
      vhi <- hl_fo_mode_subject_eval(
        data, base, ii, cand_hi, rho_full[1], rho_full[2]
      )
    if (isTRUE(vlo$ok) && isTRUE(vhi$ok) &&
        vlo$h_eta >= 0 && vhi$h_eta <= 0) {
        lo[ii] <- cand_lo
        hi[ii] <- cand_hi
        expansions[ii] <- ee
        unresolved[ii] <- FALSE
      }
    }
    if (!any(unresolved)) break
  }
  if (any(unresolved)) {
    bad <- which(unresolved)[1]
    return(list(success = FALSE, status = "mode_failure",
                reason = "mode_bracket_failure", failed_subject = bad,
                diagnostics = list(unresolved_subjects = which(unresolved))))
  }

  eta <- pmin(pmax(start, lo), hi)
  bisections <- integer(data$m)
  iterations <- integer(data$m)
  converged <- rep(FALSE, data$m)
  last <- NULL
  for (iter in seq_len(100L)) {
    cur <- hl_fo_mode_vector_eval(data, base, eta, rho_full[1], rho_full[2])
    if (!isTRUE(cur$ok)) {
      return(list(success = FALSE, status = "mode_failure", reason = cur$reason))
    }
    if (any(!is.finite(cur$h_etaeta)) || any(cur$h_etaeta >= 0)) {
      bad <- which(!is.finite(cur$h_etaeta) | cur$h_etaeta >= 0)[1]
      return(list(success = FALSE, status = "mode_failure",
                  reason = "mode_nonnegative_curvature", failed_subject = bad))
    }
    residual <- abs(cur$h_eta) / sqrt(-cur$h_etaeta)
    step_resid <- abs(cur$h_eta / cur$h_etaeta) / (1 + abs(eta))
    newly <- !converged & residual <= 1e-10 & step_resid <= 1e-10
    iterations[newly] <- iter
    converged <- converged | newly
    if (all(converged)) {
      diagnostics <- lapply(seq_len(data$m), function(ii) list(
        success = TRUE, status = "ok", reason = NULL, eta = eta[ii],
        h_eta = cur$h_eta[ii], h_etaeta = cur$h_etaeta[ii],
        residual_scaled = residual[ii], step_scaled = step_resid[ii],
        iterations = iterations[ii], bracket_expansions = expansions[ii],
        bisections = bisections[ii], bracket = c(lo[ii], hi[ii]),
        start_type = if (warm) "warm" else "cold_prior_mode"
      ))
      return(list(success = TRUE, status = "ok", eta = eta,
                  diagnostics = diagnostics))
    }
    active <- !converged
    move_lo <- active & cur$h_eta > 0
    move_hi <- active & !move_lo
    lo[move_lo] <- eta[move_lo]
    hi[move_hi] <- eta[move_hi]
    proposal <- eta - cur$h_eta / cur$h_etaeta
    unsafe <- active & (!is.finite(proposal) | proposal <= lo | proposal >= hi)
    proposal[unsafe] <- (lo[unsafe] + hi[unsafe]) / 2
    bisections[unsafe] <- bisections[unsafe] + 1L
    eta[active] <- proposal[active]
    last <- cur
  }
  bad <- which(!converged)[1]
  list(success = FALSE, status = "mode_failure",
       reason = "mode_iteration_exhausted", failed_subject = bad,
       eta = eta, h_eta = if (is.null(last)) NA_real_ else last$h_eta[bad],
       h_etaeta = if (is.null(last)) NA_real_ else last$h_etaeta[bad])
}

hl_fo_inner_quantities_r <- function(data, beta, rho_full, eta_start = NULL) {
  modes <- hl_fo_solve_modes(data, beta, rho_full, eta_start = eta_start)
  if (!isTRUE(modes$success)) {
    return(modes)
  }
  p <- data$p
  log_m <- as.numeric(data$log_offset + data$X %*% beta + modes$eta[data$subject])
  nb <- hl_fo_nb_primitives(data$y, log_m, rho_full[2])
  if (!isTRUE(nb$ok)) {
    return(list(success = FALSE, status = "mode_failure", reason = nb$reason))
  }
  v <- exp(rho_full[1])
  alpha <- 1 / expm1(v)
  log_alpha <- -log(expm1(v))
  d_prior <- modes$eta - v / 2
  lambda_z <- alpha * exp(d_prior)
  normalizer <- hl_fo_gamma_normalizer_derivatives(alpha, log_alpha)$value
  prior_value <- vapply(
    d_prior,
    function(d) alpha * hl_fo_d_minus_expm1(d) + normalizer,
    numeric(1)
  )
  prior_a <- -lambda_z
  h_subject <- as.numeric(rowsum(nb$loglik, data$subject, reorder = FALSE)) + prior_value
  a_subject <- as.numeric(rowsum(nb$t, data$subject, reorder = FALSE)) + prior_a
  if (any(!is.finite(a_subject)) || any(a_subject >= 0)) {
    bad <- which(!is.finite(a_subject) | a_subject >= 0)[1]
    return(list(success = FALSE, status = "mode_failure",
                reason = "mode_nonnegative_curvature", failed_subject = bad))
  }
  score_beta <- as.numeric(crossprod(data$X, nb$s))
  h_bb <- crossprod(data$X, data$X * nb$t)
  d_subject <- rowsum(data$X * nb$t, data$subject, reorder = FALSE)
  B_beta <- h_bb - crossprod(d_subject, d_subject / a_subject)
  P <- sum(h_subject)
  subject <- lapply(seq_len(data$m), function(ii) list(
    p = h_subject[ii],
    h_beta = as.numeric(crossprod(
      data$X[data$subject_rows[[ii]], , drop = FALSE],
      nb$s[data$subject_rows[[ii]]]
    )),
    p_beta_beta = NA,
    a = a_subject[ii]
  ))
  list(
    success = TRUE,
    status = "ok",
    beta = beta,
    eta = modes$eta,
    mode_diagnostics = modes$diagnostics,
    P = P,
    S_beta = score_beta,
    B_beta = B_beta,
    subject = subject
  )
}

hl_fo_inner_quantities <- function(data, beta, rho_full, eta_start = NULL) {
  # Keep the literal R transaction as the executable reference and as a
  # source-only fallback.  Installed builds use the parity-tested compiled
  # transaction because every inner/outer trial must completely re-solve all
  # subject modes and the ordinary profiled h-score root.
  native_loaded <- isTRUE(tryCatch(
    is.loaded("_scRoPE_hl_fo_inner_eval_cpp", PACKAGE = "scRoPE"),
    error = function(e) FALSE
  ))
  if (!native_loaded) {
    return(hl_fo_inner_quantities_r(data, beta, rho_full, eta_start))
  }
  use_warm <- !is.null(eta_start) && length(eta_start) == data$m &&
    all(is.finite(eta_start))
  hl_fo_inner_eval_cpp(
    y = data$y,
    X = data$X,
    log_offset = data$log_offset,
    subject_zero = as.integer(data$subject - 1L),
    beta = beta,
    tau = rho_full[1L],
    xi = rho_full[2L],
    eta_start = if (use_warm) as.numeric(eta_start) else numeric(0),
    use_warm = use_warm
  )
}

hl_fo_inner_recession <- function(data, K) {
  # A beta recession exists exactly when a nonzero reduced direction leaves
  # every positive-count linear predictor unchanged and moves every zero-count
  # predictor weakly downward, with at least one strict move.  Subject-mode
  # directions cannot participate because either eta tail sends the complete
  # normalized Gamma contribution to -Inf.  The compact LP is used only when
  # the positive-count rows do not already span the reduced beta space.
  K <- as.matrix(K)
  k <- ncol(K)
  if (!k) return(list(recession = FALSE, direction = numeric(0)))
  Z <- data$X %*% K
  positive <- data$y > 0
  Z_positive <- Z[positive, , drop = FALSE]
  positive_rank <- hl_fo_num_rank(Z_positive)
  if (positive_rank$rank == k) {
    return(list(recession = FALSE, direction = numeric(k)))
  }
  reduced_positive <- hl_fo_reduce_linear_rows(
    Z_positive, numeric(nrow(Z_positive))
  )
  if (!isTRUE(reduced_positive$success)) {
    return(list(recession = FALSE, direction = NULL,
                status = "nonattainment_check_failure"))
  }
  if (reduced_positive$rank == k) {
    return(list(recession = FALSE, direction = numeric(k)))
  }
  N <- .hl_fo_null_basis(reduced_positive$matrix)
  Z_zero <- Z[!positive, , drop = FALSE]
  A <- Z_zero %*% N
  if (!nrow(A) || !ncol(A) || max(abs(A)) == 0) {
    return(list(recession = FALSE, direction = numeric(k)))
  }

  objective_gradient <- colSums(A)
  if (max(abs(objective_gradient)) == 0) {
    return(list(recession = FALSE, direction = numeric(k)))
  }
  lp <- tryCatch(
    nloptr::nloptr(
      x0 = numeric(ncol(A)),
      eval_f = function(x) list(
        objective = sum(A %*% x),
        gradient = objective_gradient
      ),
      eval_g_ineq = function(x) list(
        constraints = as.numeric(A %*% x),
        jacobian = A
      ),
      lb = rep(-1, ncol(A)),
      ub = rep(1, ncol(A)),
      opts = list(
        algorithm = "NLOPT_LD_SLSQP",
        xtol_rel = 1e-12,
        ftol_abs = 1e-12,
        maxeval = 500L,
        print_level = 0L
      )
    ),
    error = function(e) NULL
  )
  if (is.null(lp) || any(!is.finite(lp$solution))) {
    return(list(recession = FALSE, direction = NULL,
                status = "nonattainment_check_failure"))
  }
  movement <- as.numeric(A %*% lp$solution)
  tolerance <- 100 * max(nrow(A), ncol(A)) * .Machine$double.eps *
    max(1, max(abs(A)))
  recession <- max(movement) <= tolerance && sum(movement) < -tolerance
  direction <- if (recession) as.numeric(K %*% N %*% lp$solution) else numeric(k)
  list(
    recession = recession,
    direction = direction,
    reduced_direction = if (recession) as.numeric(N %*% lp$solution) else numeric(k),
    zero_cell_movement = movement,
    objective = lp$objective,
    optimizer_status = lp$status
  )
}

hl_fo_make_recession_check <- function(data, K) {
  ready <- FALSE
  result <- NULL
  function() {
    if (!ready) {
      result <<- hl_fo_inner_recession(data, K)
      ready <<- TRUE
    }
    result
  }
}

hl_fo_inner_acceptance_decision <- function(score_metric, step_metric,
                                            iteration, max_iter,
                                            acceptance_tier = c(
                                              "derivative", "statistical"
                                            )) {
  acceptance_tier <- match.arg(acceptance_tier)
  derivative_grade <- is.finite(score_metric) && is.finite(step_metric) &&
    score_metric <= 1e-10 && step_metric <= 1e-10
  statistical_grade <- is.finite(score_metric) && is.finite(step_metric) &&
    score_metric <= 1e-8 && step_metric <= 1e-10
  statistical_at_cap <- identical(acceptance_tier, "statistical") &&
    iteration >= max_iter && statistical_grade
  accepted <- derivative_grade || statistical_at_cap
  list(
    accepted = accepted,
    requested_tier = acceptance_tier,
    accepted_tier = if (derivative_grade) {
      "derivative"
    } else if (statistical_at_cap) {
      "statistical"
    } else {
      NA_character_
    },
    derivative_grade = derivative_grade,
    statistical_grade = statistical_grade,
    at_iteration_cap = iteration >= max_iter,
    score_metric = score_metric,
    step_metric = step_metric,
    score_tolerance = if (identical(acceptance_tier, "derivative")) 1e-10 else 1e-8,
    step_tolerance = 1e-10
  )
}

hl_fo_solve_inner <- function(data, rho_full, beta_star, K, beta_start,
                              eta_start = NULL, cold = FALSE,
                              max_iter = 100L,
                              acceptance_tier = c("derivative", "statistical"),
                              recession_check = NULL) {
  acceptance_tier <- match.arg(acceptance_tier)
  beta_star <- as.numeric(beta_star)
  K <- as.matrix(K)
  k_free <- ncol(K)
  if (length(beta_star) != data$p || nrow(K) != data$p) {
    stop("beta_star/K dimensions do not match the design.")
  }
  reduced_design <- data$X %*% K
  if (k_free > 0L && hl_fo_matrix_rank(reduced_design) < k_free) {
    return(list(success = FALSE, status = "inner_rank_deficient",
                reason = "reduced_design_rank_deficient"))
  }
  if (!is.null(recession_check) && !is.function(recession_check)) {
    stop("recession_check must be NULL or a zero-argument function.")
  }
  recession <- if (is.null(recession_check)) {
    hl_fo_inner_recession(data, K)
  } else {
    recession_check()
  }
  if (isTRUE(recession$recession)) {
    return(list(success = FALSE, status = "inner_no_finite_root",
                reason = if (sum(data$y) == 0)
                  "all_zero_intercept_recession" else
                  "separated_inner_recession",
                recession = recession))
  }
  if (identical(recession$status, "nonattainment_check_failure")) {
    return(list(success = FALSE, status = "inner_no_convergence",
                reason = "nonattainment_check_failure",
                recession = recession))
  }

  gamma <- if (k_free == 0L) numeric(0) else
    as.numeric(crossprod(K, as.numeric(beta_start) - beta_star))
  compose <- function(g) as.numeric(beta_star + if (k_free) K %*% g else 0)
  current <- hl_fo_inner_quantities(
    data, compose(gamma), rho_full,
    eta_start = if (isTRUE(cold)) NULL else eta_start
  )
  if (!isTRUE(current$success)) {
    return(current)
  }
  if (k_free == 0L) {
    current$gamma <- numeric(0)
    current$S <- numeric(0)
    current$B <- matrix(0, 0, 0)
    current$score_metric <- 0
    current$step_metric <- 0
    current$iterations <- 0L
    current$acceptance_tier_requested <- acceptance_tier
    current$accepted_tier <- "derivative"
    current$inner_acceptance <- hl_fo_inner_acceptance_decision(
      score_metric = 0, step_metric = 0,
      iteration = 0L, max_iter = max_iter,
      acceptance_tier = acceptance_tier
    )
    current$success <- TRUE
    return(current)
  }

  trace <- vector("list", max_iter)
  for (iter in seq_len(max_iter)) {
    S <- as.numeric(crossprod(K, current$S_beta))
    B <- crossprod(K, current$B_beta %*% K)
    B_sym <- 0.5 * (B + t(B))
    eig_neg_B <- tryCatch(eigen(-B_sym, symmetric = TRUE, only.values = TRUE)$values,
                          error = function(e) rep(NA_real_, k_free))
    if (any(!is.finite(eig_neg_B)) || min(eig_neg_B) <= 0 ||
        hl_fo_matrix_rank(B) < k_free) {
      return(list(success = FALSE, status = "weak_curvature",
                  reason = "inner_curvature_not_positive_definite",
                  gamma = gamma, beta = current$beta, eta = current$eta,
                  B = B, curvature = eig_neg_B, trace = Filter(Negate(is.null), trace)))
    }
    delta <- tryCatch(solve(B, -S), error = function(e) NULL)
    score_solve <- tryCatch(solve(-B, S), error = function(e) NULL)
    if (is.null(delta) || is.null(score_solve)) {
      return(list(success = FALSE, status = "weak_curvature",
                  reason = "inner_linear_solve_failed", B = B))
    }
    score_metric <- sqrt(max(0, sum(S * score_solve)))
    step_metric <- max(abs(delta)) / (1 + max(abs(gamma)))
    trace[[iter]] <- list(
      iteration = iter, P = current$P, score_metric = score_metric,
      step_metric = step_metric, min_curvature = min(eig_neg_B)
    )
    # Prefer derivative-grade polishing because total outer derivatives can
    # amplify an otherwise admissible beta-score residual.  The closed
    # statistical tier retains that target through the 100-iteration budget,
    # then accepts only the specification's 1e-8 score and 1e-10 step gate.
    # It is used only by fixed-dispersion profiles and separately cold-
    # certified full-chi Score states.
    acceptance <- hl_fo_inner_acceptance_decision(
      score_metric = score_metric,
      step_metric = step_metric,
      iteration = iter,
      max_iter = max_iter,
      acceptance_tier = acceptance_tier
    )
    trace[[iter]]$derivative_grade <- acceptance$derivative_grade
    trace[[iter]]$statistical_grade <- acceptance$statistical_grade
    trace[[iter]]$accepted_tier <- acceptance$accepted_tier
    if (isTRUE(acceptance$accepted)) {
      current$gamma <- gamma
      current$S <- S
      current$B <- B
      current$score_metric <- score_metric
      current$step_metric <- step_metric
      current$iterations <- iter
      current$trace <- Filter(Negate(is.null), trace)
      current$acceptance_tier_requested <- acceptance_tier
      current$accepted_tier <- acceptance$accepted_tier
      current$inner_acceptance <- acceptance
      current$success <- TRUE
      return(current)
    }

    directional <- sum(S * delta)
    if (!is.finite(directional) || directional <= 0) {
      return(list(success = FALSE, status = "weak_curvature",
                  reason = "inner_nonascent_newton_direction", B = B,
                  S = S, gamma = gamma))
    }
    accepted <- FALSE
    step <- 1
    armijo_method <- NA_character_
    integrated_gain <- NA_real_
    for (ls in 0:25) {
      gamma_trial <- gamma + step * delta
      trial <- hl_fo_inner_quantities(
        data, compose(gamma_trial), rho_full, eta_start = current$eta
      )
      if (isTRUE(trial$success) && is.finite(trial$P)) {
        direct_gain <- trial$P - current$P
        if (direct_gain >= 1e-4 * step * directional) {
          accepted <- TRUE
          armijo_method <- "direct_profile_criterion"
          break
        }
        # Near a root the true O(||delta||^2) ascent can be smaller than one
        # ulp of the O(m) criterion.  In that cancellation regime evaluate the
        # same profiled criterion difference as the line integral of its
        # ordinary analytic score, using a fully mode-reprofiled Simpson
        # midpoint.  This is an Armijo check on P, not a score-only acceptance
        # or an objective surrogate.
        cancellation_tolerance <- 64 * .Machine$double.eps *
          (1 + abs(current$P) + abs(trial$P))
        if (abs(direct_gain) <= cancellation_tolerance) {
          midpoint <- hl_fo_inner_quantities(
            data,
            compose(gamma + 0.5 * step * delta),
            rho_full,
            eta_start = current$eta
          )
          if (isTRUE(midpoint$success)) {
            directional_mid <- sum(as.numeric(crossprod(K, midpoint$S_beta)) * delta)
            directional_trial <- sum(as.numeric(crossprod(K, trial$S_beta)) * delta)
            integrated_gain <- step *
              (directional + 4 * directional_mid + directional_trial) / 6
            if (is.finite(integrated_gain) &&
                integrated_gain >= 1e-4 * step * directional) {
              accepted <- TRUE
              armijo_method <- "profile_score_integral"
              break
            }
          }
        }
      }
      step <- step / 2
    }
    trace[[iter]]$armijo_step <- if (accepted) step else NA_real_
    trace[[iter]]$armijo_reductions <- if (accepted) ls else 25L
    trace[[iter]]$armijo_method <- armijo_method
    trace[[iter]]$armijo_integrated_gain <- integrated_gain
    if (!accepted) {
      return(list(
        success = FALSE,
        status = "inner_no_convergence",
        reason = "inner_armijo_failure",
        gamma = gamma, beta = current$beta, eta = current$eta,
        trace = Filter(Negate(is.null), trace)
      ))
    }
    gamma <- gamma_trial
    current <- trial
  }

  list(
    success = FALSE,
    status = "inner_no_convergence",
    reason = "inner_iteration_exhausted",
    gamma = gamma,
    beta = current$beta,
    eta = current$eta,
    acceptance_tier_requested = acceptance_tier,
    trace = Filter(Negate(is.null), trace)
  )
}

hl_fo_poisson_seed <- function(data, v0 = 1, phi0 = 1) {
  beta_u <- numeric(data$p)
  converged <- FALSE
  for (iter in seq_len(100L)) {
    eta <- as.numeric(data$log_offset + data$X %*% beta_u)
    mu <- exp(eta)
    if (any(!is.finite(mu)) || any(mu <= 0)) {
      break
    }
    score <- as.numeric(crossprod(data$X, data$y - mu))
    info <- crossprod(data$X, data$X * mu)
    step <- tryCatch(solve(info, score), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      break
    }
    beta_new <- beta_u + step
    if (hl_fo_sdisc(beta_new, beta_u) <= 1e-10) {
      beta_u <- beta_new
      converged <- TRUE
      break
    }
    beta_u <- beta_new
  }
  method <- "poisson_irls"
  if (!converged || any(!is.finite(beta_u))) {
    beta_u <- numeric(data$p)
    beta_u[data$intcol] <- log((sum(data$y) + 0.5) / sum(exp(data$log_offset)))
    method <- "pseudocount_intercept"
  }
  beta <- beta_u
  beta[data$intcol] <- beta[data$intcol] - v0 / 2
  list(
    beta = beta,
    rho_full = log(c(v0, phi0)),
    coordinate_tag = "canonical_conditional_beta_tau_xi",
    source_coordinate = "mean_one_beta_natural_dispersion",
    method = method,
    poisson_converged = converged
  )
}

hl_fo_scientific_transform <- function(p, K, active) {
  K <- as.matrix(K)
  active <- as.integer(active)
  q <- length(active)
  out <- matrix(0, p + 2L, ncol(K) + q)
  if (ncol(K)) {
    out[seq_len(p), seq_len(ncol(K))] <- K
  }
  if (q) {
    for (jj in seq_len(q)) {
      out[p + active[jj], ncol(K) + jj] <- 1
    }
  }
  out
}

hl_fo_augmented_state <- function(subjects, k_free, q) {
  m <- length(subjects)
  gamma_idx <- if (k_free) seq_len(k_free) else integer(0)
  rho_idx <- if (q) k_free + seq_len(q) else integer(0)

  B <- matrix(0, k_free, k_free)
  E <- matrix(0, k_free, q)
  L_gamma <- numeric(k_free)
  L_rho <- numeric(q)
  for (ss in subjects) {
    if (k_free) {
      B <- B + ss$p_hess[gamma_idx, gamma_idx, drop = FALSE]
      L_gamma <- L_gamma + ss$ell_grad[gamma_idx]
    }
    if (q) {
      L_rho <- L_rho + ss$ell_grad[rho_idx]
      if (k_free) {
        E <- E + ss$p_hess[gamma_idx, rho_idx, drop = FALSE]
      }
    }
  }
  chi <- if (k_free) tryCatch(solve(t(B), L_gamma), error = function(e) NULL) else numeric(0)
  if (k_free && is.null(chi)) {
    return(list(success = FALSE, status = "sensitivity_singular",
                reason = "auxiliary_beta_solve_failed", B = B))
  }
  g <- L_rho - if (k_free && q) as.numeric(crossprod(E, chi)) else numeric(q)

  d_aug <- 2L * k_free + q
  psi_subject <- matrix(0, d_aug, m)
  jac_subject <- vector("list", m)
  R_sum <- matrix(0, k_free, k_free)
  U_sum <- matrix(0, k_free, q)
  V_sum <- matrix(0, q, q)

  for (ii in seq_along(subjects)) {
    ss <- subjects[[ii]]
    Bi <- if (k_free) ss$p_hess[gamma_idx, gamma_idx, drop = FALSE] else
      matrix(0, 0, 0)
    Ei <- if (k_free && q) ss$p_hess[gamma_idx, rho_idx, drop = FALSE] else
      matrix(0, k_free, q)
    Ri <- matrix(0, k_free, k_free)
    Ui <- matrix(0, k_free, q)
    Vi <- if (q) ss$ell_hess[rho_idx, rho_idx, drop = FALSE] else matrix(0, 0, 0)
    if (k_free) {
      Ri <- -ss$ell_hess[gamma_idx, gamma_idx, drop = FALSE]
      for (aa in seq_len(k_free)) {
        for (bb in seq_len(k_free)) {
          Ri[aa, bb] <- Ri[aa, bb] +
            sum(chi * ss$p_third[gamma_idx, aa, bb])
        }
      }
      if (q) {
        Ui <- -ss$ell_hess[gamma_idx, rho_idx, drop = FALSE]
        for (aa in seq_len(k_free)) {
          for (rr in seq_len(q)) {
            Ui[aa, rr] <- Ui[aa, rr] +
              sum(chi * ss$p_third[gamma_idx, aa, rho_idx[rr]])
          }
        }
        for (rr in seq_len(q)) {
          for (ss2 in seq_len(q)) {
            Vi[rr, ss2] <- Vi[rr, ss2] -
              sum(chi * ss$p_third[gamma_idx, rho_idx[rr], rho_idx[ss2]])
          }
        }
      }
    }

    psi <- numeric(d_aug)
    if (k_free) {
      psi[seq_len(k_free)] <- ss$p_grad[gamma_idx]
      psi[k_free + seq_len(k_free)] <- as.numeric(t(Bi) %*% chi) -
        ss$ell_grad[gamma_idx]
    }
    if (q) {
      psi[2L * k_free + seq_len(q)] <- ss$ell_grad[rho_idx] -
        if (k_free) as.numeric(crossprod(Ei, chi)) else 0
    }
    psi_subject[, ii] <- psi

    Ji <- matrix(0, d_aug, d_aug)
    if (k_free) {
      i1 <- seq_len(k_free)
      i2 <- k_free + seq_len(k_free)
      Ji[i1, i1] <- Bi
      Ji[i2, i1] <- Ri
      Ji[i2, i2] <- t(Bi)
      if (q) {
        i3 <- 2L * k_free + seq_len(q)
        Ji[i1, i3] <- Ei
        Ji[i2, i3] <- Ui
        Ji[i3, i1] <- -t(Ui)
        Ji[i3, i2] <- -t(Ei)
        Ji[i3, i3] <- Vi
      }
    } else if (q) {
      Ji[,] <- Vi
    }
    jac_subject[[ii]] <- Ji
    R_sum <- R_sum + Ri
    U_sum <- U_sum + Ui
    V_sum <- V_sum + Vi
  }

  H_F <- if (!q) {
    matrix(0, 0, 0)
  } else if (!k_free) {
    V_sum
  } else {
    BinvE <- tryCatch(solve(B, E), error = function(e) NULL)
    BtinvU <- tryCatch(solve(t(B), U_sum), error = function(e) NULL)
    BtinvRBinvE <- if (!is.null(BinvE))
      tryCatch(solve(t(B), R_sum %*% BinvE), error = function(e) NULL) else NULL
    if (is.null(BinvE) || is.null(BtinvU) || is.null(BtinvRBinvE)) {
      return(list(success = FALSE, status = "sensitivity_singular",
                  reason = "outer_hessian_solve_failed", B = B, E = E))
    }
    V_sum + t(U_sum) %*% BinvE + t(E) %*% BtinvU -
      t(E) %*% BtinvRBinvE
  }

  A <- if (d_aug) -Reduce(`+`, jac_subject) else matrix(0, 0, 0)
  J <- if (d_aug) tcrossprod(psi_subject) else matrix(0, 0, 0)
  list(
    success = TRUE,
    status = "ok",
    B = B,
    E = E,
    L_gamma = L_gamma,
    L_rho = L_rho,
    chi = as.numeric(chi),
    g = as.numeric(g),
    R = R_sum,
    U = U_sum,
    V = V_sum,
    H_F = H_F,
    psi_subject = psi_subject,
    jacobian_subject = jac_subject,
    A = A,
    J = J,
    augmented_residual = if (d_aug) rowSums(psi_subject) else numeric(0)
  )
}

hl_fo_profile_evaluate <- function(data, rho_active, rho_template, active,
                                   beta_star, K, beta_start, eta_start = NULL,
                                   cold = FALSE,
                                   include_laplace_constant = TRUE,
                                   inner_acceptance_tier = "derivative",
                                   recession_check = NULL) {
  rho_full <- as.numeric(rho_template)
  active <- as.integer(active)
  if (length(active)) {
    rho_full[active] <- rho_active
  }
  inner <- hl_fo_solve_inner(
    data = data,
    rho_full = rho_full,
    beta_star = beta_star,
    K = K,
    beta_start = beta_start,
    eta_start = eta_start,
    cold = cold,
    acceptance_tier = inner_acceptance_tier,
    recession_check = recession_check
  )
  if (!isTRUE(inner$success)) {
    return(inner)
  }
  transform <- hl_fo_scientific_transform(data$p, K, active)
  subjects <- vector("list", data$m)
  for (ii in seq_len(data$m)) {
    rows <- data$subject_rows[[ii]]
    ss <- hl_fo_subject_raw(
      y = data$y[rows],
      X = data$X[rows, , drop = FALSE],
      log_offset = data$log_offset[rows],
      beta = inner$beta,
      eta = inner$eta[ii],
      rho_full = rho_full,
      include_laplace_constant = include_laplace_constant
    )
    if (!isTRUE(ss$ok)) {
      return(list(success = FALSE, status = "mode_failure",
                  reason = ss$reason, failed_subject = ii))
    }
    subjects[[ii]] <- hl_fo_contract_subject(ss, transform)
  }
  augmented <- hl_fo_augmented_state(subjects, ncol(K), length(active))
  if (!isTRUE(augmented$success)) {
    return(augmented)
  }
  list(
    success = TRUE,
    status = "ok",
    data = data,
    rho = as.numeric(rho_active),
    rho_full = rho_full,
    active = active,
    beta_star = as.numeric(beta_star),
    K = K,
    gamma = inner$gamma,
    beta = inner$beta,
    eta = inner$eta,
    inner = inner,
    subjects = subjects,
    F = sum(vapply(subjects, `[[`, numeric(1), "ell")),
    P = sum(vapply(subjects, `[[`, numeric(1), "p")),
    g = augmented$g,
    H_F = augmented$H_F,
    augmented = augmented,
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
}

hl_fo_cache_distance <- function(a, b) {
  if (!length(a) && !length(b)) return(0)
  max(abs(a - b) / (1 + pmax(abs(a), abs(b))))
}

hl_fo_outer_starts <- function(primary, lower, upper) {
  q <- length(lower)
  if (!q) return(list(numeric(0)))
  starts <- list(as.numeric(primary), (lower + upper) / 2)
  grid <- expand.grid(rep(list(c(0.2, 0.5, 0.8)), q), KEEP.OUT.ATTRS = FALSE)
  for (ii in seq_len(nrow(grid))) {
    frac <- as.numeric(grid[ii, ])
    starts[[length(starts) + 1L]] <- lower + frac * (upper - lower)
  }
  kept <- list()
  for (candidate in starts) {
    if (length(candidate) != q || any(!is.finite(candidate)) ||
        any(candidate <= lower) || any(candidate >= upper)) next
    duplicate <- any(vapply(kept, function(old) {
      hl_fo_cache_distance(candidate, old) <= 1e-12
    }, logical(1)))
    if (!duplicate) kept[[length(kept) + 1L]] <- candidate
  }
  kept
}

hl_fo_trust_step_max <- function(g, H, radius) {
  q <- length(g)
  if (!q) return(numeric(0))
  Hs <- 0.5 * (H + t(H))
  ee <- eigen(Hs, symmetric = TRUE)
  lam <- ee$values
  Q <- ee$vectors
  d <- as.numeric(crossprod(Q, g))
  tol <- 100 * .Machine$double.eps * max(1, max(abs(lam)))

  if (max(lam) < -tol) {
    unconstrained <- tryCatch(solve(Hs, -g), error = function(e) NULL)
    if (!is.null(unconstrained) && all(is.finite(unconstrained)) &&
        sqrt(sum(unconstrained^2)) <= radius) {
      return(as.numeric(unconstrained))
    }
  }

  lower <- max(0, max(lam))
  top <- which(abs(lam - max(lam)) <= tol)
  non_top <- setdiff(seq_len(q), top)
  top_signal <- sqrt(sum(d[top]^2))
  if (top_signal <= tol) {
    base_eig <- numeric(q)
    if (length(non_top)) {
      base_eig[non_top] <- d[non_top] / (lower - lam[non_top])
    }
    base <- as.numeric(Q %*% base_eig)
    remaining <- radius^2 - sum(base^2)
    if (remaining >= -tol) {
      remaining <- max(0, remaining)
      if (length(top) == q) {
        direction <- numeric(q)
        direction[1] <- -1
      } else if (length(top) == 1L) {
        direction <- Q[, top]
        plus <- base + sqrt(remaining) * direction
        minus <- base - sqrt(remaining) * direction
        first_difference <- which(plus != minus)[1L]
        direction <- if (is.na(first_difference) ||
                         plus[first_difference] <= minus[first_difference]) {
          direction
        } else {
          -direction
        }
      } else {
        # q <= 2: a repeated top eigenspace is the whole space.
        direction <- numeric(q)
        direction[1] <- -1
      }
      return(as.numeric(base + sqrt(remaining) * direction))
    }
  }

  # Solve the secular equation in the distance from the lower multiplier
  # bound.  Adding a tiny root to a positive leading eigenvalue can round back
  # to that eigenvalue; an ordinary uniroot bracket in the multiplier then has
  # equal-sign endpoints even though the boundary solution exists.  Working
  # in `delta` keeps the top-eigenspace denominator representable and gives a
  # deterministic midpoint fallback all the way to the trust boundary.
  gaps <- lower - lam
  norm_at_delta <- function(delta) {
    denominators <- gaps + delta
    terms <- d / denominators
    sqrt(sum(terms^2))
  }
  hi <- max(1, top_signal / radius)
  while (norm_at_delta(hi) > radius && hi < 1e16) hi <- hi * 2
  hi_norm <- norm_at_delta(hi)
  if (!is.finite(hi_norm) || hi_norm > radius) {
    return(rep(NA_real_, q))
  }

  lo <- 0
  for (iteration in seq_len(200L)) {
    mid <- lo + (hi - lo) / 2
    if (mid <= lo || mid >= hi) break
    mid_norm <- norm_at_delta(mid)
    if (!is.finite(mid_norm) || mid_norm > radius) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  denominators <- gaps + hi
  if (any(denominators <= 0) || any(!is.finite(denominators))) {
    return(rep(NA_real_, q))
  }
  as.numeric(Q %*% (d / denominators))
}

hl_fo_outer_stationarity <- function(state, tolerance = 1e-7) {
  q <- length(state$g)
  if (!q) {
    return(list(success = TRUE, scaled_score = 0, quadratic = 0,
                displacement = 0, maximum = TRUE))
  }
  H <- state$H_F
  Hs <- 0.5 * (H + t(H))
  eig <- tryCatch(eigen(Hs, symmetric = TRUE, only.values = TRUE)$values,
                  error = function(e) rep(NA_real_, q))
  maximum <- all(is.finite(eig)) && max(eig) < 0
  psi3 <- state$augmented$psi_subject[
    2L * ncol(state$K) + seq_len(q), , drop = FALSE
  ]
  scale <- 1 + sqrt(rowSums(psi3^2))
  scaled_score <- max(abs(state$g) / scale)
  quadratic <- Inf
  displacement <- Inf
  if (maximum) {
    qsolve <- tryCatch(solve(-Hs, state$g), error = function(e) NULL)
    dsolve <- tryCatch(solve(Hs, state$g), error = function(e) NULL)
    if (!is.null(qsolve)) quadratic <- sqrt(max(0, sum(state$g * qsolve)))
    if (!is.null(dsolve)) {
      displacement <- max(abs(dsolve)) / (1 + max(abs(state$rho)))
    }
  }
  list(
    success = maximum && scaled_score <= tolerance && quadratic <= tolerance &&
      displacement <= 1e-8,
    scaled_score = scaled_score,
    quadratic = quadratic,
    displacement = displacement,
    maximum = maximum,
    eigenvalues = eig
  )
}

hl_fo_make_outer_evaluator <- function(data, rho_template, active, beta_star, K,
                                       seed_beta, include_laplace_constant = TRUE,
                                       constraint_C = NULL,
                                       constraint_rhs = numeric(0),
                                       inner_acceptance_tier = "derivative",
                                       recession_check = NULL) {
  if (is.null(recession_check)) {
    recession_check <- hl_fo_make_recession_check(data, K)
  }
  cache <- list()
  trace <- list()
  evaluate <- function(rho, cold = FALSE) {
    nearest <- NULL
    if (!cold && length(cache)) {
      distances <- vapply(cache, function(x) hl_fo_cache_distance(rho, x$rho), numeric(1))
      nearest <- cache[[which.min(distances)]]
    }
    beta_start <- if (is.null(nearest)) seed_beta else nearest$beta
    eta_start <- if (is.null(nearest)) NULL else nearest$eta
    ans <- hl_fo_profile_evaluate(
      data = data,
      rho_active = rho,
      rho_template = rho_template,
      active = active,
      beta_star = beta_star,
      K = K,
      beta_start = beta_start,
      eta_start = eta_start,
      cold = cold,
      include_laplace_constant = include_laplace_constant,
      inner_acceptance_tier = inner_acceptance_tier,
      recession_check = recession_check
    )
    reconstruction <- if (isTRUE(ans$success))
      as.numeric(ans$beta - beta_star - K %*% ans$gamma) else NULL
    true_residual <- if (isTRUE(ans$success) && !is.null(constraint_C))
      as.numeric(constraint_C %*% ans$beta - constraint_rhs) else reconstruction
    true_scaled <- if (isTRUE(ans$success) && !is.null(constraint_C)) {
      max(abs(true_residual)) /
        (1 + norm(constraint_C, "I") * max(abs(ans$beta)) +
           max(abs(constraint_rhs)))
    } else if (isTRUE(ans$success) && length(true_residual)) {
      max(abs(true_residual)) / (1 + max(abs(ans$beta)))
    } else if (isTRUE(ans$success)) 0 else NA_real_
    trace[[length(trace) + 1L]] <<- list(
      rho = as.numeric(rho), success = isTRUE(ans$success),
      status = ans$status, reason = ans$reason,
      constraint_residual = true_residual,
      constraint_residual_scaled = true_scaled,
      affine_reconstruction_residual = reconstruction
    )
    if (isTRUE(ans$success)) {
      # The cache is only a deterministic warm-start index.  Retaining every
      # subject derivative tensor from every outer trial can otherwise grow to
      # gigabytes on multistart constrained fits without changing a single
      # numerical decision.
      cache[[length(cache) + 1L]] <<- list(
        rho = ans$rho,
        beta = ans$beta,
        eta = ans$eta
      )
    }
    ans
  }
  list(
    evaluate = evaluate,
    cache = function() cache,
    trace = function() trace,
    count = function() length(trace)
  )
}

hl_fo_optimize_outer_start <- function(start, lower, upper, evaluator) {
  radius <- 1
  invalid_consecutive <- 0L
  invalid_guard_faces <- character(0)
  start_count <- evaluator$count()
  current <- evaluator$evaluate(start)
  if (!isTRUE(current$success)) {
    terminal_inner <- current$status %in% c(
      "mode_failure", "inner_rank_deficient", "inner_no_finite_root",
      "inner_no_convergence", "weak_curvature"
    )
    return(list(success = FALSE,
                status = if (terminal_inner) current$status else
                  "outer_no_stationary_root",
                reason = if (terminal_inner) current$reason else "invalid_outer_start",
                inner_failure = if (terminal_inner) current else NULL,
                start = start,
                evaluations = evaluator$count() - start_count))
  }
  trust_trace <- list()
  polish_attempts <- 0L
  repeat {
    station <- hl_fo_outer_stationarity(current)
    if (isTRUE(station$success)) {
      # A warm inner cache may satisfy every root tolerance while leaving a
      # delicate total derivative different at the required replay scale.
      # Reprofile the stationary candidate from the declared cold starts
      # before it is admitted to the candidate set.  Selection therefore uses
      # the same fully cold-replayed derivative state that final certification
      # must reproduce, without relaxing the 1e-8 replay discrepancy gate.
      cold_candidate <- evaluator$evaluate(current$rho, cold = TRUE)
      if (!isTRUE(cold_candidate$success)) {
        return(list(
          success = FALSE,
          status = cold_candidate$status,
          reason = "stationary_candidate_cold_replay_failure",
          inner_failure = cold_candidate,
          state = current,
          start = start,
          evaluations = evaluator$count() - start_count,
          trust_trace = trust_trace
        ))
      }
      cold_station <- hl_fo_outer_stationarity(cold_candidate)
      trust_trace[[length(trust_trace) + 1L]] <- list(
        radius = radius,
        rho = current$rho,
        step = rep(0, length(current$rho)),
        pred = 0,
        ared = cold_candidate$F - current$F,
        accepted = isTRUE(cold_station$success),
        reason = "stationary_candidate_cold_replay",
        discrepancy = list(
          beta = hl_fo_sdisc(cold_candidate$beta, current$beta),
          eta = hl_fo_sdisc(cold_candidate$eta, current$eta),
          F = hl_fo_sdisc(cold_candidate$F, current$F),
          g = hl_fo_sdisc(cold_candidate$g, current$g),
          H_F = hl_fo_sdisc(cold_candidate$H_F, current$H_F)
        )
      )
      if (isTRUE(cold_station$success)) {
        return(list(success = TRUE, status = "ok", state = cold_candidate,
                    stationarity = cold_station, start = start,
                    evaluations = evaluator$count() - start_count,
                    trust_trace = trust_trace))
      }
      current <- cold_candidate
      station <- cold_station
    }
    # Once the trust iterations have identified the local maximum basin, use
    # the specification's total-gradient polishing step.  Every proposal is a
    # complete reprofiled evaluation; acceptance is by a strict reduction of
    # the analytic total-gradient residual, not by an objective-only code.
    if (station$maximum && station$scaled_score <= 1e-4 && polish_attempts < 10L) {
      polish_step <- tryCatch(-solve(current$H_F, current$g), error = function(e) NULL)
      polish_attempts <- polish_attempts + 1L
      if (!is.null(polish_step) && all(is.finite(polish_step))) {
        polish_rho <- current$rho + polish_step
        if (all(polish_rho > lower) && all(polish_rho < upper)) {
          polished <- evaluator$evaluate(polish_rho)
          if (isTRUE(polished$success)) {
            polished_station <- hl_fo_outer_stationarity(polished)
            trust_trace[[length(trust_trace) + 1L]] <- list(
              radius = radius, rho = current$rho, step = polish_step,
              pred = NA_real_, ared = polished$F - current$F,
              ratio = NA_real_, accepted =
                polished_station$scaled_score < station$scaled_score,
              reason = "total_gradient_polish"
            )
            if (polished_station$scaled_score < station$scaled_score) {
              current <- polished
              next
            }
          }
        }
      }
    }
    if (evaluator$count() - start_count >= 200L) {
      return(list(success = FALSE, status = "outer_no_stationary_root",
                  reason = "outer_evaluation_cap", state = current,
                  start = start, evaluations = evaluator$count() - start_count,
                  trust_trace = trust_trace,
                  invalid_guard_faces = unique(invalid_guard_faces)))
    }
    if (radius <= 1e-12 * max(1, sqrt(sum(current$rho^2)))) {
      return(list(success = FALSE,
                  status = if (station$scaled_score <= 1e-7 && !station$maximum)
                    "outer_not_maximum" else "outer_no_stationary_root",
                  reason = "outer_radius_collapse", state = current,
                  start = start, evaluations = evaluator$count() - start_count,
                  trust_trace = trust_trace,
                  invalid_guard_faces = unique(invalid_guard_faces)))
    }

    step <- hl_fo_trust_step_max(current$g, current$H_F, radius)
    pred <- if (all(is.finite(step)))
      sum(current$g * step) + 0.5 * sum(step * (current$H_F %*% step)) else NA_real_
    entry <- list(radius = radius, rho = current$rho, step = step, pred = pred)
    if (!is.finite(pred) || pred <= 0) {
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
      lower_faces <- which(proposal <= lower)
      upper_faces <- which(proposal >= upper)
      faces <- c(
        if (length(lower_faces)) paste0(lower_faces, ":lower") else character(0),
        if (length(upper_faces)) paste0(upper_faces, ":upper") else character(0)
      )
      invalid_guard_faces <- c(invalid_guard_faces, faces)
      entry$accepted <- FALSE
      entry$reason <- "outside_guard"
      entry$guard_faces <- faces
      trust_trace[[length(trust_trace) + 1L]] <- entry
      if (invalid_consecutive >= 5L) {
        return(list(success = FALSE, status = "outer_no_stationary_root",
                    reason = "five_invalid_trials", state = current,
                    start = start, evaluations = evaluator$count() - start_count,
                    trust_trace = trust_trace,
                    invalid_guard_faces = unique(invalid_guard_faces)))
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
        return(list(success = FALSE, status = "outer_no_stationary_root",
                    reason = "five_invalid_trials", state = current,
                    start = start, evaluations = evaluator$count() - start_count,
                    trust_trace = trust_trace,
                    invalid_guard_faces = unique(invalid_guard_faces)))
      }
      next
    }
    invalid_consecutive <- 0L
    ared <- trial$F - current$F
    ratio <- ared / pred
    entry$ared <- ared
    entry$ratio <- ratio
    if (!is.finite(ratio) || ratio < 0.1) {
      radius <- radius / 4
      entry$accepted <- FALSE
    } else {
      current <- trial
      entry$accepted <- TRUE
      if (ratio > 0.75 && sqrt(sum(step^2)) >= 0.9 * radius) {
        radius <- min(100, radius * 2)
      }
    }
    trust_trace[[length(trust_trace) + 1L]] <- entry
  }
}

hl_fo_cluster_candidates <- function(candidates) {
  kept <- list()
  for (candidate in candidates) {
    duplicate <- any(vapply(kept, function(old) {
      hl_fo_cache_distance(candidate$state$rho, old$state$rho) <= 1e-8 &&
        hl_fo_sdisc(candidate$state$beta, old$state$beta) <= 1e-8
    }, logical(1)))
    if (!duplicate) kept[[length(kept) + 1L]] <- candidate
  }
  kept
}

hl_fo_terminal_boundary_faces <- function(run, lower, upper) {
  recorded <- if (is.null(run$invalid_guard_faces)) character(0) else
    unique(run$invalid_guard_faces)
  if (identical(run$reason, "five_invalid_trials") && length(recorded)) {
    return(recorded)
  }
  state <- run$state
  if (!is.list(state) || !isTRUE(state$success) ||
      length(state$rho) != length(lower) || length(state$g) != length(lower) ||
      any(!is.finite(c(state$rho, state$g)))) {
    return(character(0))
  }
  # The 1e-6 face-distance rule applies to a selected, valid candidate.  This
  # helper implements the specification's separate alternative in which all
  # starts require invalid evaluations beyond a guard face.  Once such a face
  # has actually been crossed by a rejected proposal, retain it whenever the
  # terminal total gradient still points outward; requiring the terminal state
  # itself to remain within 1e-6 of the face can otherwise erase direct guard-
  # crossing evidence after the radius has been quartered.
  outward <- vapply(recorded, function(face) {
    pieces <- strsplit(face, ":", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !pieces[[2L]] %in% c("lower", "upper")) {
      return(FALSE)
    }
    index <- suppressWarnings(as.integer(pieces[[1L]]))
    if (is.na(index) || index < 1L || index > length(state$g)) {
      return(FALSE)
    }
    if (identical(pieces[[2L]], "lower")) state$g[index] < 0 else
      state$g[index] > 0
  }, logical(1))
  recorded[outward]
}

hl_fo_polish_cold_outer <- function(state, evaluator, lower, upper,
                                    max_iter = 10L) {
  trace <- list()
  current <- state
  if (!length(current$g)) {
    return(list(success = TRUE, state = current,
                stationarity = hl_fo_outer_stationarity(current), trace = trace))
  }
  for (iteration in 0:max_iter) {
    stationarity <- hl_fo_outer_stationarity(current)
    if (isTRUE(stationarity$success)) {
      return(list(success = TRUE, state = current,
                  stationarity = stationarity, trace = trace))
    }
    if (iteration == max_iter) break
    newton <- tryCatch(-solve(current$H_F, current$g), error = function(e) NULL)
    if (is.null(newton) || any(!is.finite(newton))) break
    accepted <- FALSE
    factor <- 1
    for (reduction in 0:12) {
      proposal <- current$rho + factor * newton
      if (all(proposal > lower) && all(proposal < upper)) {
        trial <- evaluator$evaluate(proposal, cold = TRUE)
        if (isTRUE(trial$success)) {
          trial_stationarity <- hl_fo_outer_stationarity(trial)
          improved <- trial_stationarity$scaled_score < stationarity$scaled_score ||
            trial_stationarity$quadratic < stationarity$quadratic
          if (isTRUE(improved)) {
            trace[[length(trace) + 1L]] <- list(
              iteration = iteration + 1L,
              factor = factor,
              step = factor * newton,
              before = stationarity,
              after = trial_stationarity
            )
            current <- trial
            accepted <- TRUE
            break
          }
        }
      }
      factor <- factor / 2
    }
    if (!accepted) break
  }
  list(success = FALSE, state = current,
       stationarity = hl_fo_outer_stationarity(current), trace = trace)
}

hl_fo_sandwich <- function(A, J) {
  if (!nrow(A)) {
    return(list(success = TRUE, V = matrix(0, 0, 0), rcond = Inf,
                asymmetry = 0, eigenvalues = numeric(0)))
  }
  rc <- tryCatch(rcond(A), error = function(e) 0)
  if (!is.finite(rc) || rc < sqrt(.Machine$double.eps)) {
    return(list(success = FALSE, status = "sensitivity_singular",
                reason = "augmented_sensitivity_condition", rcond = rc))
  }
  left <- tryCatch(solve(A, J), error = function(e) NULL)
  V <- if (!is.null(left)) tryCatch(t(solve(A, t(left))), error = function(e) NULL) else NULL
  if (is.null(V) || any(!is.finite(V))) {
    return(list(success = FALSE, status = "covariance_invalid",
                reason = "sandwich_linear_solve_failed", rcond = rc))
  }
  asym <- hl_fo_asymmetry(V)
  eig <- tryCatch(eigen(0.5 * (V + t(V)), symmetric = TRUE, only.values = TRUE)$values,
                  error = function(e) NA_real_)
  psd_ok <- all(is.finite(eig)) && min(eig) >= -1e-8 * max(1, max(eig))
  if (asym > 1e-10 || !psd_ok) {
    return(list(success = FALSE, status = "covariance_invalid",
                reason = if (asym > 1e-10) "sandwich_asymmetry" else "sandwich_not_psd",
                V = V, rcond = rc, asymmetry = asym, eigenvalues = eig))
  }
  list(success = TRUE, status = "ok", V = V, rcond = rc,
       asymmetry = asym, eigenvalues = eig)
}

hl_fo_rank_limit_diagnostics <- function(data, k_free, q,
                                         J = NULL) {
  augmented_dimension <- 2L * as.integer(k_free) + as.integer(q)
  list(
    subject_count = data$m,
    augmented_dimension = augmented_dimension,
    J_rank = if (is.null(J)) NA_integer_ else hl_fo_matrix_rank(J),
    J_rank_bound = min(augmented_dimension, data$m - 1L),
    asymptotic_calibration = "chi_square_subject_clusters",
    small_subject_rank_limited = data$m - 1L < augmented_dimension
  )
}

hl_fo_finalize_state <- function(state, constraint_C = NULL, constraint_rhs = NULL) {
  k_free <- ncol(state$K)
  q <- length(state$active)
  B <- state$augmented$B
  H <- state$H_F
  diagnostics <- hl_fo_rank_limit_diagnostics(
    state$data, k_free, q, J = state$augmented$J
  )

  if (k_free) {
    diagnostics$B_asymmetry <- hl_fo_asymmetry(B)
    diagnostics$B_rcond <- rcond(-B)
    diagnostics$B_eigenvalues <- eigen(-0.5 * (B + t(B)), symmetric = TRUE,
                                        only.values = TRUE)$values
    b_gate <- diagnostics$B_asymmetry <= 1e-8 &&
      diagnostics$B_rcond >= sqrt(.Machine$double.eps) &&
      min(diagnostics$B_eigenvalues) > sqrt(.Machine$double.eps) *
        max(1, norm(B, "2"))
    if (!b_gate) {
      state$success <- FALSE
      state$point_success <- FALSE
      state$fit_status <- "weak_curvature"
      state$fit_failure_reason <- "inner_curvature_certification"
      state$certification <- diagnostics
      return(state)
    }
  }
  if (q) {
    diagnostics$H_asymmetry <- hl_fo_asymmetry(H)
    diagnostics$H_rcond <- rcond(-H)
    diagnostics$H_eigenvalues <- eigen(-0.5 * (H + t(H)), symmetric = TRUE,
                                        only.values = TRUE)$values
    h_gate <- diagnostics$H_asymmetry <= 1e-8 &&
      diagnostics$H_rcond >= sqrt(.Machine$double.eps) &&
      min(diagnostics$H_eigenvalues) > sqrt(.Machine$double.eps) *
        max(1, norm(H, "2"))
    if (!h_gate) {
      state$success <- FALSE
      state$point_success <- FALSE
      state$fit_status <- "weak_curvature"
      state$fit_failure_reason <- "outer_curvature_certification"
      state$certification <- diagnostics
      return(state)
    }
  }

  if (!is.null(constraint_C)) {
    resid <- as.numeric(constraint_C %*% state$beta - constraint_rhs)
    denom <- 1 + norm(constraint_C, "I") * max(abs(state$beta)) +
      max(abs(constraint_rhs))
    diagnostics$constraint_residual <- resid
    diagnostics$constraint_residual_scaled <- max(abs(resid)) / denom
    if (diagnostics$constraint_residual_scaled > 100 * .Machine$double.eps) {
      state$success <- FALSE
      state$point_success <- FALSE
      state$fit_status <- "constraint_infeasible"
      state$fit_failure_reason <- "constraint_residual"
      state$certification <- diagnostics
      return(state)
    }
  }

  if (k_free == 0L && q == 0L) {
    state$point_success <- TRUE
    state$fit_status <- "ok"
    state$fit_failure_reason <- NULL
    state$covariance_success <- FALSE
    state$covariance_status <- "not_applicable_all_parameters_fixed"
    state$covariance_failure_reason <- NULL
    state$covariance <- NULL
    state$certification <- diagnostics
    return(state)
  }

  sandwich <- hl_fo_sandwich(state$augmented$A, state$augmented$J)
  state$point_success <- TRUE
  state$fit_status <- if (identical(state$status, "multiple_near_tied_maxima"))
    "multiple_near_tied_maxima" else "ok"
  state$fit_failure_reason <- NULL
  state$covariance_success <- isTRUE(sandwich$success)
  state$covariance_status <- if (isTRUE(sandwich$success)) "ok" else sandwich$status
  state$covariance_failure_reason <- if (isTRUE(sandwich$success)) NULL else sandwich$reason
  diagnostics$sensitivity_rcond <- sandwich$rcond

  if (!isTRUE(sandwich$success)) {
    state$certification <- diagnostics
    return(state)
  }
  Vz <- sandwich$V
  if (nrow(Vz)) {
    if (any(!is.finite(diag(Vz))) || any(diag(Vz) <= 0)) {
      state$covariance_success <- FALSE
      state$covariance_status <- "covariance_invalid"
      state$covariance_failure_reason <- "nonpositive_raw_sandwich_diagonal"
      state$certification <- diagnostics
      return(state)
    }
    delta <- tryCatch(solve(state$augmented$A, state$augmented$augmented_residual),
                      error = function(e) NULL)
    if (is.null(delta)) {
      state$covariance_success <- FALSE
      state$covariance_status <- "sensitivity_singular"
      state$covariance_failure_reason <- "statistical_displacement_solve"
      state$certification <- diagnostics
      return(state)
    }
    displacement <- abs(delta) / sqrt(diag(Vz))
    diagnostics$statistical_displacement <- displacement
    diagnostics$statistical_displacement_limit <- min(1e-3, state$data$m^(-1 / 4))
    if (max(displacement) > diagnostics$statistical_displacement_limit) {
      state$covariance_success <- FALSE
      state$covariance_status <- "covariance_invalid"
      state$covariance_failure_reason <- "statistical_root_not_certified"
      state$certification <- diagnostics
      return(state)
    }
  }
  if (k_free) {
    aux_resid <- as.numeric(t(B) %*% state$augmented$chi - state$augmented$L_gamma)
    diagnostics$auxiliary_backward_error <- max(abs(aux_resid)) /
      (1 + max(abs(t(B) %*% state$augmented$chi)) +
         max(abs(state$augmented$L_gamma)))
    if (diagnostics$auxiliary_backward_error > 1e-10) {
      state$covariance_success <- FALSE
      state$covariance_status <- "covariance_invalid"
      state$covariance_failure_reason <- "auxiliary_backward_error"
      state$certification <- diagnostics
      return(state)
    }
  }

  scientific_idx <- c(if (k_free) seq_len(k_free) else integer(0),
                      if (q) 2L * k_free + seq_len(q) else integer(0))
  Vgr <- Vz[scientific_idx, scientific_idx, drop = FALSE]
  map <- matrix(0, state$data$p + q, k_free + q)
  if (k_free) map[seq_len(state$data$p), seq_len(k_free)] <- state$K
  if (q) {
    map[state$data$p + seq_len(q), k_free + seq_len(q)] <- diag(q)
  }
  Vbr <- map %*% Vgr %*% t(map)
  map_full <- matrix(0, state$data$p + 2L, k_free + q)
  if (k_free) map_full[seq_len(state$data$p), seq_len(k_free)] <- state$K
  if (q) {
    for (jj in seq_len(q)) {
      map_full[state$data$p + state$active[jj], k_free + jj] <- 1
    }
  }
  state$covariance <- list(
    V_zeta_reduced = Vz,
    V_gamma_rho = Vgr,
    V_beta_rho = Vbr,
    V_canonical_full = map_full %*% Vgr %*% t(map_full),
    V_beta = Vbr[seq_len(state$data$p), seq_len(state$data$p), drop = FALSE],
    sandwich = sandwich,
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  state$certification <- diagnostics
  state
}

hl_fo_fit_affine <- function(data, beta_star = rep(0, data$p), K = diag(data$p),
                             active = c(1L, 2L), fixed_natural = numeric(0),
                             guard_min = c(1e-4, 1e-4),
                             guard_max = c(10, 1000), seed = NULL,
                             include_laplace_constant = TRUE,
                             constraint_C = NULL, constraint_rhs = NULL) {
  active <- as.integer(active)
  if (any(!active %in% 1:2) || anyDuplicated(active) || is.unsorted(active)) {
    stop("active dispersion indices must be an ordered subset of c(1, 2).")
  }
  if (length(guard_min) != 2L || length(guard_max) != 2L ||
      any(!is.finite(c(guard_min, guard_max))) || any(guard_min <= 0) ||
      any(guard_max <= guard_min)) {
    stop("Dispersion guards must be finite positive ordered pairs.")
  }
  beta_star <- as.numeric(beta_star)
  K <- as.matrix(K)
  if (length(beta_star) != data$p || any(!is.finite(beta_star)) ||
      nrow(K) != data$p || any(!is.finite(K))) {
    stop("beta_star and K must be finite and match the analysis design.")
  }
  k_free <- ncol(K)
  orthogonality <- hl_fo_sdisc(crossprod(K), diag(k_free))
  if (!is.finite(orthogonality) || orthogonality > 1e-10) {
    stop("K must have deterministic orthonormal columns.")
  }
  if (!is.null(constraint_C)) {
    constraint_C <- if (is.null(dim(constraint_C)))
      matrix(constraint_C, nrow = 1L) else as.matrix(constraint_C)
    if (ncol(constraint_C) != data$p || any(!is.finite(constraint_C))) {
      stop("constraint_C must be finite with one column per beta coordinate.")
    }
    if (is.null(constraint_rhs) || !length(constraint_rhs)) {
      constraint_rhs <- numeric(nrow(constraint_C))
    }
    if (length(constraint_rhs) != nrow(constraint_C) ||
        any(!is.finite(constraint_rhs))) {
      stop("constraint_rhs must be finite with one entry per constraint row.")
    }
    anchor_check <- hl_fo_constraint_residual(
      constraint_C, beta_star, constraint_rhs
    )
    null_defect <- if (k_free) max(abs(constraint_C %*% K)) /
      max(1, norm(constraint_C, "I")) else 0
    if (!isTRUE(anchor_check$satisfied) || null_defect > 1e-10) {
      stop("beta_star and K do not satisfy the declared affine constraint.")
    }
  }
  if (hl_fo_matrix_rank(data$X) < data$p &&
      ncol(K) == data$p && max(abs(K - diag(data$p))) < 1e-14) {
    return(list(success = FALSE, point_success = FALSE,
                status = "inner_rank_deficient", fit_status = "inner_rank_deficient",
                reason = "analysis_design_rank_deficient", data = data))
  }
  if (is.null(seed)) seed <- hl_fo_poisson_seed(data)
  if (!identical(seed$coordinate_tag, "canonical_conditional_beta_tau_xi")) {
    canonical_seed <- tryCatch(
      hl_fo_canonicalize_beta_seed(seed, intercept_idx = data$intcol,
                                   expected_p = data$p),
      error = function(e) stop(conditionMessage(e), call. = FALSE)
    )
    seed <- list(
      beta = canonical_seed$beta,
      rho_full = c(
        if (is.null(canonical_seed$tau)) 0 else canonical_seed$tau,
        if (is.null(canonical_seed$xi)) 0 else canonical_seed$xi
      ),
      coordinate_tag = "canonical_conditional_beta_tau_xi",
      source = canonical_seed$source
    )
  }
  if (!is.list(seed) || length(seed$beta) != data$p ||
      any(!is.finite(seed$beta)) || length(seed$rho_full) != 2L ||
      any(!is.finite(seed$rho_full)) ||
      !identical(seed$coordinate_tag, "canonical_conditional_beta_tau_xi")) {
    stop("The canonical seed must contain finite beta and (tau, xi) values.")
  }
  rho_template <- as.numeric(seed$rho_full)
  if (length(fixed_natural)) {
    fixed_names <- names(fixed_natural)
    fixed_idx <- match(fixed_names, c("sigma2", "phi"))
    numeric_names <- suppressWarnings(as.integer(fixed_names))
    fixed_idx[is.na(fixed_idx)] <- numeric_names[is.na(fixed_idx)]
    if (anyNA(fixed_idx) || any(!fixed_idx %in% 1:2) ||
        anyDuplicated(fixed_idx) || any(!is.finite(fixed_natural)) ||
        any(fixed_natural <= 0)) {
      stop("fixed_natural must be named sigma2/phi or by indices 1/2.")
    }
    if (any(fixed_idx %in% active)) {
      stop("A dispersion coordinate cannot be both active and fixed.")
    }
    rho_template[fixed_idx] <- log(as.numeric(fixed_natural))
  }
  lower_full <- log(guard_min)
  upper_full <- log(guard_max)
  expansion_trace <- list()
  boundary_rounds <- list()
  recession_check <- hl_fo_make_recession_check(data, K)

  if (!length(active)) {
    evaluator <- hl_fo_make_outer_evaluator(
      data, rho_template, active, beta_star, K, seed$beta,
      include_laplace_constant, constraint_C, constraint_rhs,
      inner_acceptance_tier = "statistical",
      recession_check = recession_check
    )
    state <- evaluator$evaluate(numeric(0), cold = TRUE)
    if (!isTRUE(state$success)) {
      state$point_success <- FALSE
      state$fit_status <- state$status
      return(state)
    }
    state$outer_trace <- evaluator$trace()
    state$guard_expansions <- expansion_trace
    state$cold_replay <- list(beta_discrepancy = 0, F_discrepancy = 0)
    return(hl_fo_finalize_state(state, constraint_C, constraint_rhs))
  }

  selected <- NULL
  last_runs <- NULL
  for (round in 0:8) {
    lower <- lower_full[active]
    upper <- upper_full[active]
    primary <- rho_template[active]
    if (any(primary <= lower) || any(primary >= upper)) primary <- (lower + upper) / 2
    evaluator <- hl_fo_make_outer_evaluator(
      data, rho_template, active, beta_star, K, seed$beta,
      include_laplace_constant, constraint_C, constraint_rhs,
      recession_check = recession_check
    )
    starts <- hl_fo_outer_starts(primary, lower, upper)
    runs <- lapply(starts, hl_fo_optimize_outer_start,
                   lower = lower, upper = upper, evaluator = evaluator)
    valid <- Filter(function(x) isTRUE(x$success), runs)
    if (!length(valid)) {
      run_status <- vapply(runs, `[[`, character(1), "status")
      terminal <- c("mode_failure", "inner_rank_deficient", "inner_no_finite_root",
                    "inner_no_convergence", "weak_curvature")
      run_boundary_faces <- lapply(
        runs, hl_fo_terminal_boundary_faces, lower = lower, upper = upper
      )
      boundary_only <- all(lengths(run_boundary_faces) > 0L)
      boundary_faces <- unique(unlist(run_boundary_faces, use.names = FALSE))
      if (boundary_only) {
        boundary_rounds[[length(boundary_rounds) + 1L]] <- list(
          round = round,
          bounds = list(lower = lower, upper = upper),
          faces = boundary_faces,
          runs = runs,
          trace = evaluator$trace()
        )
      }
      if (boundary_only && round < 8L) {
        parsed_faces <- do.call(rbind, lapply(boundary_faces, function(face) {
          pieces <- strsplit(face, ":", fixed = TRUE)[[1L]]
          c(index = as.integer(pieces[[1L]]),
            side = match(pieces[[2L]], c("lower", "upper")))
        }))
        parsed_faces <- unique(parsed_faces[order(
          parsed_faces[, "index"], parsed_faces[, "side"]
        ), , drop = FALSE])
        for (ff in seq_len(nrow(parsed_faces))) {
          jj <- parsed_faces[ff, "index"]
          side <- c("lower", "upper")[[parsed_faces[ff, "side"]]]
          if (identical(side, "lower")) {
            lower_full[active[jj]] <- lower_full[active[jj]] - 2
          } else {
            upper_full[active[jj]] <- upper_full[active[jj]] + 2
          }
          expansion_trace[[length(expansion_trace) + 1L]] <- list(
            round = round + 1L,
            coordinate = active[jj],
            face = side,
            range = c(lower_full[active[jj]], upper_full[active[jj]])
          )
        }
        terminal_states <- Filter(
          function(run) is.list(run$state) && isTRUE(run$state$success),
          runs
        )
        if (length(terminal_states)) {
          terminal_values <- vapply(
            terminal_states, function(run) run$state$F, numeric(1)
          )
          rho_template[active] <- terminal_states[[which.max(terminal_values)]]$state$rho
        }
        next
      }
      chosen_status <- if (boundary_only) {
        "outer_boundary"
      } else if (length(unique(run_status)) == 1L && run_status[1] %in% terminal) {
        run_status[1]
      } else {
        "outer_no_stationary_root"
      }
      return(list(success = FALSE, point_success = FALSE,
                  status = chosen_status,
                  fit_status = chosen_status,
                  reason = if (chosen_status == "outer_boundary")
                    "all_starts_exhausted_beyond_guard" else if (
                      chosen_status == "outer_no_stationary_root"
                    ) "no_certified_multistart_candidate" else runs[[1]]$reason,
                  boundary_faces = boundary_faces,
                  guard_bounds = list(lower = lower, upper = upper),
                  guard_expansions = expansion_trace,
                  boundary_rounds = boundary_rounds,
                  outer_runs = runs, outer_trace = evaluator$trace(), data = data))
    }
    distinct <- hl_fo_cluster_candidates(valid)
    values <- vapply(distinct, function(x) x$state$F, numeric(1))
    best_value <- max(values)
    tied <- which(abs(values - best_value) <=
                    1e-10 * (1 + pmax(abs(values), abs(best_value))))
    best <- tied[1]
    selected <- distinct[[best]]$state
    selected$status <- if (length(tied) > 1L) "multiple_near_tied_maxima" else "ok"
    selected$selection_status <- selected$status
    selected$outer_runs <- runs
    selected$outer_trace <- evaluator$trace()
    selected$starts <- starts
    selected$candidates <- lapply(distinct, `[[`, "state")
    selected$guard_bounds <- list(lower = lower, upper = upper)
    last_runs <- runs

    rel_lower <- (selected$rho - lower) / (upper - lower)
    rel_upper <- (upper - selected$rho) / (upper - lower)
    faces_lower <- which(rel_lower <= 1e-6)
    faces_upper <- which(rel_upper <= 1e-6)
    if (!length(faces_lower) && !length(faces_upper)) break
    if (round == 8L) {
      selected$success <- FALSE
      selected$point_success <- FALSE
      selected$fit_status <- "outer_boundary"
      selected$fit_failure_reason <- "expanded_guard_still_active"
      selected$guard_expansions <- expansion_trace
      return(selected)
    }
    for (jj in faces_lower) {
      lower_full[active[jj]] <- lower_full[active[jj]] - 2
      expansion_trace[[length(expansion_trace) + 1L]] <-
        list(round = round + 1L, coordinate = active[jj], face = "lower",
             range = c(lower_full[active[jj]], upper_full[active[jj]]))
    }
    for (jj in faces_upper) {
      upper_full[active[jj]] <- upper_full[active[jj]] + 2
      expansion_trace[[length(expansion_trace) + 1L]] <-
        list(round = round + 1L, coordinate = active[jj], face = "upper",
             range = c(lower_full[active[jj]], upper_full[active[jj]]))
    }
    rho_template[active] <- selected$rho
  }

  cold_evaluator <- hl_fo_make_outer_evaluator(
    data, selected$rho_full, active, beta_star, K, seed$beta,
    include_laplace_constant, constraint_C, constraint_rhs,
    recession_check = recession_check
  )
  cold <- cold_evaluator$evaluate(selected$rho, cold = TRUE)
  if (isTRUE(cold$success)) {
    polished_cold <- hl_fo_polish_cold_outer(
      cold, cold_evaluator,
      lower = selected$guard_bounds$lower,
      upper = selected$guard_bounds$upper
    )
    cold <- polished_cold$state
  } else {
    polished_cold <- list(success = FALSE, state = cold,
                          stationarity = NULL, trace = list())
  }
  cold_discrepancy <- if (isTRUE(cold$success)) list(
    beta = hl_fo_sdisc(cold$beta, selected$beta),
    eta = hl_fo_sdisc(cold$eta, selected$eta),
    F = hl_fo_sdisc(cold$F, selected$F),
    g = hl_fo_sdisc(cold$g, selected$g),
    H_F = hl_fo_sdisc(cold$H_F, selected$H_F),
    A = hl_fo_sdisc(cold$augmented$A, selected$augmented$A),
    J = hl_fo_sdisc(cold$augmented$J, selected$augmented$J)
  ) else NULL
  cold_mismatch <- !isTRUE(cold$success) || !isTRUE(polished_cold$success) ||
    any(unlist(cold_discrepancy, use.names = FALSE) > 1e-8)
  if (cold_mismatch) {
    selected$success <- FALSE
    selected$point_success <- FALSE
    selected$fit_status <- if (isTRUE(cold$success) &&
                               !isTRUE(polished_cold$success))
      "outer_no_stationary_root" else "inner_no_convergence"
    selected$fit_failure_reason <- if (identical(selected$fit_status,
                                                  "outer_no_stationary_root"))
      "cold_outer_stationarity_not_certified" else "cold_replay_mismatch"
    selected$cold_replay <- list(
      state = cold,
      discrepancy = cold_discrepancy,
      stationarity = polished_cold$stationarity,
      polish_trace = polished_cold$trace
    )
    return(selected)
  }
  # Use the cold replay as the authoritative final derivative state.
  cold$status <- selected$status
  cold$selection_status <- selected$selection_status
  cold$outer_runs <- selected$outer_runs
  cold$outer_trace <- selected$outer_trace
  cold$starts <- selected$starts
  cold$candidates <- selected$candidates
  cold$guard_bounds <- selected$guard_bounds
  cold$guard_expansions <- expansion_trace
  cold$cold_replay <- c(list(
    beta_discrepancy = cold_discrepancy$beta,
    F_discrepancy = cold_discrepancy$F
  ), cold_discrepancy, list(
    stationarity = polished_cold$stationarity,
    polish_trace = polished_cold$trace,
    evaluator_trace = cold_evaluator$trace()
  ))
  hl_fo_finalize_state(cold, constraint_C, constraint_rhs)
}

hl_fo_fit <- function(data, C = NULL, rhs = numeric(0),
                      active = c("tau", "xi"), fixed_natural = numeric(0),
                      guard_min = c(1e-4, 1e-4),
                      guard_max = c(10, 1000), seed = NULL,
                      include_laplace_constant = TRUE) {
  active_idx <- if (is.character(active)) match(active, c("tau", "xi")) else
    as.integer(active)
  if (anyNA(active_idx)) stop("active must select tau and/or xi.")
  affine <- hl_fo_affine_constraint(C = C, rhs = rhs, p = data$p)
  if (!isTRUE(affine$success)) {
    return(list(
      success = FALSE,
      point_success = FALSE,
      status = "constraint_infeasible",
      fit_status = "constraint_infeasible",
      reason = affine$status,
      constraint = affine,
      data = data
    ))
  }
  fit <- hl_fo_fit_affine(
    data = data,
    beta_star = affine$beta_star,
    K = affine$K,
    active = active_idx,
    fixed_natural = fixed_natural,
    guard_min = guard_min,
    guard_max = guard_max,
    seed = seed,
    include_laplace_constant = include_laplace_constant,
    constraint_C = if (is.null(C)) NULL else if (is.null(dim(C)))
      matrix(C, nrow = 1L) else as.matrix(C),
    constraint_rhs = rhs
  )
  if (is.null(fit$certification)) {
    fit$certification <- hl_fo_rank_limit_diagnostics(
      data, ncol(affine$K), length(active_idx)
    )
  }
  fit$constraint <- affine
  fit
}
