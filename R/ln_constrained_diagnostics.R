ln_reduced_mapping <- function(nb, K) {
  K <- as.matrix(K)
  if (nrow(K) != nb) {
    stop("K must have nb rows for the LN reduced-coordinate map.")
  }

  free_dim <- ncol(K)
  map <- matrix(0, nrow = nb + 2L, ncol = free_dim + 2L)
  if (free_dim > 0L) {
    map[seq_len(nb), seq_len(free_dim)] <- K
  }
  map[nb + 1L, free_dim + 1L] <- 1
  map[nb + 2L, free_dim + 2L] <- 1
  map
}

ln_reduced_derivative_state <- function(full_gradient,
                                        full_hessian,
                                        K,
                                        theta_red,
                                        per_subject_gradients = NULL) {
  full_gradient <- as.numeric(full_gradient)
  full_hessian <- as.matrix(full_hessian)
  K <- as.matrix(K)
  theta_red <- as.numeric(theta_red)
  nb <- nrow(K)
  map <- ln_reduced_mapping(nb, K)

  if (length(full_gradient) != nb + 2L) {
    stop("LN full gradient must have length nb + 2.")
  }
  if (!all(dim(full_hessian) == c(nb + 2L, nb + 2L))) {
    stop("LN full Hessian must be a square (nb + 2) matrix.")
  }
  if (length(theta_red) != ncol(map)) {
    stop("theta_red must contain gamma, log(sigma^2), and log(phi).")
  }
  if (!is.null(per_subject_gradients) &&
    (!is.matrix(per_subject_gradients) || nrow(per_subject_gradients) != nb + 2L)) {
    stop("LN per-subject gradients must have nb + 2 rows.")
  }

  g_red <- as.numeric(crossprod(map, full_gradient))
  H_red <- crossprod(map, full_hessian %*% map)
  H_red <- 0.5 * (H_red + t(H_red))
  per_subject_red <- if (is.null(per_subject_gradients)) {
    NULL
  } else {
    crossprod(map, per_subject_gradients)
  }
  newton <- reduced_newton_diagnostics(H_red, g_red, theta_red)

  finite_gradient <- is.finite(g_red)
  list(
    parameterization = "gamma, log(sigma^2), log(phi)",
    map = map,
    theta_red = theta_red,
    g_red = g_red,
    H_red = H_red,
    per_subject_gradients_red = per_subject_red,
    diagnostics = list(
      reduced_gradient_max_abs = if (any(finite_gradient)) {
        max(abs(g_red[finite_gradient]))
      } else {
        NA_real_
      },
      reduced_gradient_l2 = if (all(finite_gradient)) sqrt(sum(g_red^2)) else NA_real_,
      H_red_positive_definite = newton$H_red_positive_definite,
      H_red_rcond = newton$H_red_rcond,
      newton_step_red = newton$newton_step_red,
      newton_step_red_max_abs = newton$newton_step_red_max_abs,
      max_abs_newton_step_scaled = newton$max_abs_newton_step_scaled,
      gradient_quadratic_lr_bound = newton$gradient_quadratic_lr_bound
    )
  )
}

ln_assess_lrt_convergence <- function(fit,
                                      use_information_scaled_convergence = TRUE,
                                      newton_step_scaled_tol = 1e-5,
                                      gradient_quadratic_lr_bound_tol = 1e-6) {
  diagnostics <- if (is.null(fit$diagnostics)) list() else fit$diagnostics
  legacy_converged <- isTRUE(fit$lrt_gate$converged)
  conv_code <- if (!is.null(fit$debug$conv_code) && length(fit$debug$conv_code)) {
    as.integer(fit$debug$conv_code[[1L]])
  } else {
    NA_integer_
  }

  H_red_positive_definite <- isTRUE(diagnostics$H_red_positive_definite)
  newton_step_scaled <- if (!is.null(diagnostics$max_abs_newton_step_scaled)) {
    as.numeric(diagnostics$max_abs_newton_step_scaled[[1L]])
  } else {
    NA_real_
  }
  quadratic_lr_bound <- if (!is.null(diagnostics$gradient_quadratic_lr_bound)) {
    as.numeric(diagnostics$gradient_quadratic_lr_bound[[1L]])
  } else {
    NA_real_
  }
  reduced_gradient <- if (!is.null(diagnostics$reduced_gradient_max_abs)) {
    as.numeric(diagnostics$reduced_gradient_max_abs[[1L]])
  } else {
    NA_real_
  }
  H_red_rcond <- if (!is.null(diagnostics$H_red_rcond)) {
    as.numeric(diagnostics$H_red_rcond[[1L]])
  } else {
    NA_real_
  }
  diagnostics_finite <- all(is.finite(c(
    reduced_gradient,
    H_red_rcond,
    newton_step_scaled,
    quadratic_lr_bound
  )))

  newton_step_scaled_ok <- isTRUE(use_information_scaled_convergence) &&
    H_red_positive_definite &&
    diagnostics_finite &&
    is.finite(newton_step_scaled) &&
    newton_step_scaled >= 0 &&
    newton_step_scaled <= newton_step_scaled_tol
  gradient_quadratic_lr_bound_ok <- isTRUE(use_information_scaled_convergence) &&
    H_red_positive_definite &&
    diagnostics_finite &&
    is.finite(quadratic_lr_bound) &&
    quadratic_lr_bound >= 0 &&
    quadratic_lr_bound <= gradient_quadratic_lr_bound_tol
  information_scaled_convergence_ok <- newton_step_scaled_ok ||
    gradient_quadratic_lr_bound_ok

  finite_likelihood_state <- is.finite(fit$loglik_ln) &&
    !is.null(fit$theta_ln$beta) &&
    all(is.finite(c(fit$theta_ln$beta, fit$theta_ln$subVar, fit$theta_ln$cellVar)))
  # Only legacy PML iteration/damping failures are eligible. Boundary,
  # curvature, nonfinite-likelihood, and robust-adjustment failures remain hard.
  information_scaled_eligible <- isTRUE(finite_likelihood_state) &&
    !is.na(conv_code) &&
    conv_code %in% c(-10L, -20L, -40L)
  information_scaled_accepted <- information_scaled_eligible &&
    information_scaled_convergence_ok
  accepted <- legacy_converged || information_scaled_accepted

  list(
    accepted = accepted,
    convergence_basis = if (legacy_converged) {
      "legacy"
    } else if (information_scaled_accepted) {
      "information_scaled"
    } else {
      NA_character_
    },
    legacy_converged = legacy_converged,
    legacy_conv_code = conv_code,
    information_scaled_eligible = information_scaled_eligible,
    information_scaled_convergence_ok = information_scaled_convergence_ok,
    information_scaled_accepted = information_scaled_accepted,
    H_red_positive_definite = H_red_positive_definite,
    diagnostics_finite = diagnostics_finite,
    H_red_rcond = H_red_rcond,
    reduced_gradient_max_abs = reduced_gradient,
    max_abs_newton_step_scaled = newton_step_scaled,
    gradient_quadratic_lr_bound = quadratic_lr_bound,
    newton_step_scaled_ok = newton_step_scaled_ok,
    gradient_quadratic_lr_bound_ok = gradient_quadratic_lr_bound_ok,
    use_information_scaled_convergence = isTRUE(use_information_scaled_convergence),
    newton_step_scaled_tol = newton_step_scaled_tol,
    gradient_quadratic_lr_bound_tol = gradient_quadratic_lr_bound_tol
  )
}

ln_lrt_gate_override <- function(fit,
                                 L,
                                 b,
                                 use_information_scaled_convergence = TRUE,
                                 newton_step_scaled_tol = 1e-5,
                                 gradient_quadratic_lr_bound_tol = 1e-6) {
  convergence <- ln_assess_lrt_convergence(
    fit,
    use_information_scaled_convergence = use_information_scaled_convergence,
    newton_step_scaled_tol = newton_step_scaled_tol,
    gradient_quadratic_lr_bound_tol = gradient_quadratic_lr_bound_tol
  )
  beta_ln <- fit$theta_ln$beta
  constraint_residual <- if (is.null(beta_ln)) {
    NA_real_
  } else {
    as.numeric(L %*% beta_ln - b)
  }
  constraint_residual_max_abs <- if (length(constraint_residual) &&
    all(is.finite(constraint_residual))) {
    max(abs(constraint_residual))
  } else {
    NA_real_
  }
  constraint_satisfied <- is.finite(constraint_residual_max_abs) &&
    constraint_residual_max_abs <= sqrt(.Machine$double.eps)

  list(
    converged = convergence$accepted,
    convergence_basis = convergence$convergence_basis,
    optimizer_converged = convergence$legacy_converged,
    legacy_converged = convergence$legacy_converged,
    legacy_conv_code = convergence$legacy_conv_code,
    information_scaled_eligible = convergence$information_scaled_eligible,
    information_scaled_convergence_ok = convergence$information_scaled_convergence_ok,
    information_scaled_accepted = convergence$information_scaled_accepted,
    H_red_positive_definite = convergence$H_red_positive_definite,
    diagnostics_finite = convergence$diagnostics_finite,
    reduced_gradient_max_abs = convergence$reduced_gradient_max_abs,
    max_abs_newton_step_scaled = convergence$max_abs_newton_step_scaled,
    gradient_quadratic_lr_bound = convergence$gradient_quadratic_lr_bound,
    newton_step_scaled_ok = convergence$newton_step_scaled_ok,
    gradient_quadratic_lr_bound_ok = convergence$gradient_quadratic_lr_bound_ok,
    use_information_scaled_convergence = convergence$use_information_scaled_convergence,
    newton_step_scaled_tol = convergence$newton_step_scaled_tol,
    gradient_quadratic_lr_bound_tol = convergence$gradient_quadratic_lr_bound_tol,
    H_red_rcond = convergence$H_red_rcond,
    reduced_parameterization = if (!is.null(fit$reduced$parameterization)) {
      fit$reduced$parameterization
    } else {
      NA_character_
    },
    constraint_residual_max_abs = constraint_residual_max_abs,
    constraint_satisfied = constraint_satisfied
  )
}
