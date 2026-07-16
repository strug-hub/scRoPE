reduced_newton_diagnostics <- function(H_red, g_red, theta_red) {
  H_red <- as.matrix(H_red)
  g_red <- as.numeric(g_red)
  theta_red <- as.numeric(theta_red)

  empty <- list(
    H_red_positive_definite = FALSE,
    H_red_rcond = NA_real_,
    newton_step_red = rep(NA_real_, length(g_red)),
    newton_step_red_max_abs = NA_real_,
    max_abs_newton_step_scaled = NA_real_,
    gradient_quadratic_lr_bound = NA_real_
  )

  if (!length(g_red) ||
    !all(dim(H_red) == c(length(g_red), length(g_red))) ||
    length(theta_red) != length(g_red) ||
    !all(is.finite(H_red)) ||
    !all(is.finite(g_red)) ||
    !all(is.finite(theta_red))) {
    return(empty)
  }

  H_red <- 0.5 * (H_red + t(H_red))
  empty$H_red_rcond <- tryCatch(
    rcond(H_red),
    error = function(e) NA_real_
  )

  chol_H <- tryCatch(
    chol(H_red),
    error = function(e) NULL
  )
  if (is.null(chol_H)) {
    return(empty)
  }

  step <- tryCatch(
    backsolve(chol_H, forwardsolve(t(chol_H), g_red)),
    error = function(e) NULL
  )
  if (is.null(step) || !all(is.finite(step))) {
    return(empty)
  }

  quad <- as.numeric(crossprod(g_red, step))
  if (!is.finite(quad)) {
    return(empty)
  }

  list(
    H_red_positive_definite = TRUE,
    H_red_rcond = empty$H_red_rcond,
    newton_step_red = as.numeric(step),
    newton_step_red_max_abs = max(abs(step)),
    max_abs_newton_step_scaled = max(abs(step) / pmax(1, abs(theta_red))),
    gradient_quadratic_lr_bound = quad
  )
}

assess_constrained_convergence <- function(diagnostics,
                                           optimizer_converged,
                                           n_subjects = NA_real_,
                                           accept_nonzero_optimizer_code = TRUE,
                                           gradient_tol = 1e-2,
                                           average_gradient_tol = 1e-4,
                                           use_information_scaled_convergence = TRUE,
                                           newton_step_scaled_tol = 1e-5,
                                           gradient_quadratic_lr_bound_tol = 1e-6) {
  if (is.null(diagnostics)) {
    diagnostics <- list()
  }

  reduced_gradient <- if (!is.null(diagnostics$reduced_gradient_max_abs)) {
    as.numeric(diagnostics$reduced_gradient_max_abs)
  } else {
    NA_real_
  }
  n_subjects <- as.numeric(n_subjects)
  average_reduced_gradient <- if (length(reduced_gradient) == 1L &&
    is.finite(reduced_gradient) &&
    length(n_subjects) == 1L &&
    is.finite(n_subjects) &&
    n_subjects > 0) {
    reduced_gradient / n_subjects
  } else {
    NA_real_
  }

  raw_gradient_ok <- length(reduced_gradient) == 1L &&
    is.finite(reduced_gradient) &&
    reduced_gradient <= gradient_tol
  average_gradient_ok <- length(average_reduced_gradient) == 1L &&
    is.finite(average_reduced_gradient) &&
    average_reduced_gradient <= average_gradient_tol
  gradient_ok <- raw_gradient_ok || average_gradient_ok

  H_red_positive_definite <- isTRUE(diagnostics$H_red_positive_definite)
  newton_step_scaled <- if (!is.null(diagnostics$max_abs_newton_step_scaled)) {
    as.numeric(diagnostics$max_abs_newton_step_scaled)
  } else {
    NA_real_
  }
  quadratic_lr_bound <- if (!is.null(diagnostics$gradient_quadratic_lr_bound)) {
    as.numeric(diagnostics$gradient_quadratic_lr_bound)
  } else {
    NA_real_
  }

  newton_step_scaled_ok <- isTRUE(use_information_scaled_convergence) &&
    H_red_positive_definite &&
    length(newton_step_scaled) == 1L &&
    is.finite(newton_step_scaled) &&
    newton_step_scaled >= 0 &&
    newton_step_scaled <= newton_step_scaled_tol
  gradient_quadratic_lr_bound_ok <- isTRUE(use_information_scaled_convergence) &&
    H_red_positive_definite &&
    length(quadratic_lr_bound) == 1L &&
    is.finite(quadratic_lr_bound) &&
    quadratic_lr_bound >= 0 &&
    quadratic_lr_bound <= gradient_quadratic_lr_bound_tol
  information_scaled_convergence_ok <- newton_step_scaled_ok ||
    gradient_quadratic_lr_bound_ok
  numerical_convergence_ok <- gradient_ok || information_scaled_convergence_ok

  optimizer_converged <- isTRUE(optimizer_converged)
  accepted <- optimizer_converged ||
    (isTRUE(accept_nonzero_optimizer_code) && numerical_convergence_ok)
  convergence_basis <- if (optimizer_converged) {
    "optimizer"
  } else if (isTRUE(accept_nonzero_optimizer_code) && gradient_ok) {
    "raw_gradient"
  } else if (isTRUE(accept_nonzero_optimizer_code) && information_scaled_convergence_ok) {
    "information_scaled"
  } else {
    NA_character_
  }

  list(
    accepted = accepted,
    convergence_basis = convergence_basis,
    optimizer_converged = optimizer_converged,
    optimizer_warning = accepted && !optimizer_converged,
    raw_gradient_ok = raw_gradient_ok,
    average_gradient_ok = average_gradient_ok,
    gradient_ok = gradient_ok,
    information_scaled_convergence_ok = information_scaled_convergence_ok,
    numerical_convergence_ok = numerical_convergence_ok,
    H_red_positive_definite = H_red_positive_definite,
    reduced_gradient_max_abs = reduced_gradient,
    average_reduced_gradient_max_abs = average_reduced_gradient,
    max_abs_newton_step_scaled = newton_step_scaled,
    gradient_quadratic_lr_bound = quadratic_lr_bound,
    newton_step_scaled_ok = newton_step_scaled_ok,
    gradient_quadratic_lr_bound_ok = gradient_quadratic_lr_bound_ok,
    use_information_scaled_convergence = isTRUE(use_information_scaled_convergence),
    gradient_tol = gradient_tol,
    average_gradient_tol = average_gradient_tol,
    newton_step_scaled_tol = newton_step_scaled_tol,
    gradient_quadratic_lr_bound_tol = gradient_quadratic_lr_bound_tol,
    failure_reason = if (accepted) NA_character_ else "nonconverged_constrained"
  )
}
