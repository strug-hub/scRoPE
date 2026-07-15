pmg_build_contrast_def <- function(def, idx) {
  if (is.null(def$L)) {
    stop("Each contrast definition must include an L matrix.")
  }
  L <- as.matrix(def$L)
  label <- if (!is.null(def$label)) {
    as.character(def$label)
  } else {
    paste0("contrast_", idx)
  }
  b_vec <- if (!is.null(def$b)) {
    as.numeric(def$b)
  } else {
    rep(0, nrow(L))
  }
  list(L = L, label = label, b = b_vec)
}

pmg_robust_wald <- function(unconstrained,
                            L,
                            b,
                            ridge_factor = 1e-8,
                            cond_limit = 1e8) {
  beta_hat <- unconstrained$theta$beta
  cov_beta <- unconstrained$robust_cov

  if (is.null(beta_hat) || is.null(cov_beta)) {
    return(list(
      stat = NA_real_,
      delta = matrix(NA_real_, nrow(L), 1),
      variance = matrix(NA_real_, nrow(L), nrow(L)),
      ridge = NA_real_,
      method = NA_character_,
      success = FALSE,
      failure_reason = normalise_failure("missing_godambe", "wald")
    ))
  }

  delta <- matrix(as.numeric(L %*% beta_hat - b), ncol = 1)
  cov_psi <- L %*% cov_beta %*% t(L)
  cov_psi <- 0.5 * (cov_psi + t(cov_psi))

  cov_inv_info <- invert_with_ridge(
    cov_psi,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )

  stat <- NA_real_
  method <- NA_character_
  failure_reason <- NA_character_
  success <- FALSE

  if (cov_inv_info$success) {
    stat <- as.numeric(t(delta) %*% cov_inv_info$inverse %*% delta)
    if (!is.na(stat) && stat < 0 && abs(stat) <= 1e-8) {
      stat <- 0
    }
    method <- if (cov_inv_info$ridge > 0) "wald_ridge" else "wald"
    success <- is.finite(stat)
    if (!success) {
      failure_reason <- normalise_failure("undefined_stat", "wald")
    }
  } else {
    failure_reason <- normalise_failure("singular_cov", "wald")
  }

  list(
    stat = stat,
    delta = delta,
    variance = cov_psi,
    ridge = cov_inv_info$ridge,
    method = method,
    success = success,
    failure_reason = failure_reason
  )
}

pmg_adjusted_contrasts <- function(ctx,
                                   gene_index,
                                   posv,
                                   contrasts,
                                   unconstrained = NULL,
                                   ridge_factor = 1e-8,
                                   cond_limit = 1e8,
                                   retry_constrained = TRUE,
                                   use_average_objective = FALSE,
                                   accept_nonzero_optimizer_code = TRUE,
                                   gradient_tol = 1e-2,
                                   average_gradient_tol = 1e-4,
                                   use_information_scaled_convergence = TRUE,
                                   newton_step_scaled_tol = 1e-5,
                                   gradient_quadratic_lr_bound_tol = 1e-6) {
  if (is.null(contrasts) || length(contrasts) == 0) {
    return(vector("list", 0))
  }

  if (is.null(unconstrained)) {
    unconstrained <- fit_gene_pmg_unconstrained(
      gene_index = gene_index,
      posv = posv,
      ctx = ctx
    )
  }

  defs <- vector("list", length(contrasts))
  for (idx in seq_along(contrasts)) {
    defs[[idx]] <- pmg_build_contrast_def(contrasts[[idx]], idx)
  }

  results <- vector("list", length(defs))
  names_vec <- character(length(defs))

  for (idx in seq_along(defs)) {
    def <- defs[[idx]]
    assess_constrained_fit <- function(fit) {
      pmg_assess_constrained_convergence(
        diagnostics = fit$diagnostics,
        optimizer_converged = isTRUE(fit$convergence > 0L),
        n_subjects = ctx$k,
        accept_nonzero_optimizer_code = accept_nonzero_optimizer_code,
        gradient_tol = gradient_tol,
        average_gradient_tol = average_gradient_tol,
        use_information_scaled_convergence = use_information_scaled_convergence,
        newton_step_scaled_tol = newton_step_scaled_tol,
        gradient_quadratic_lr_bound_tol = gradient_quadratic_lr_bound_tol
      )
    }

    fit_pair <- pmg_fit_constrained_pair(
      gene_index = gene_index,
      posv = posv,
      ctx = ctx,
      L = def$L,
      b = def$b,
      unconstrained_fit = unconstrained,
      warm_start = unconstrained,
      retry = FALSE,
      use_average_objective = use_average_objective
    )
    fit_unconstrained <- fit_pair$unconstrained
    fit_con <- fit_pair$constrained
    convergence <- assess_constrained_fit(fit_con)

    negative_raw_lr_initial <- is.finite(fit_pair$diagnostics$raw_lr_final) &&
      fit_pair$diagnostics$raw_lr_final < -1e-6
    retry_triggered <- isTRUE(retry_constrained) &&
      (!isTRUE(convergence$accepted) || negative_raw_lr_initial)

    if (retry_triggered) {
      fit_pair <- pmg_fit_constrained_pair(
        gene_index = gene_index,
        posv = posv,
        ctx = ctx,
        L = def$L,
        b = def$b,
        unconstrained_fit = unconstrained,
        warm_start = fit_con,
        retry = TRUE,
        use_average_objective = use_average_objective
      )
      fit_unconstrained <- fit_pair$unconstrained
      fit_con <- fit_pair$constrained
      convergence <- assess_constrained_fit(fit_con)
    }

    constraint_residual <- as.numeric(def$L %*% fit_con$theta$beta - def$b)
    constraint_residual_max_abs <- if (length(constraint_residual) &&
      all(is.finite(constraint_residual))) {
      max(abs(constraint_residual))
    } else {
      NA_real_
    }
    constraint_satisfied <- is.finite(constraint_residual_max_abs) &&
      constraint_residual_max_abs <= sqrt(.Machine$double.eps)

    lrt_gate_override <- list(
      converged = convergence$accepted,
      convergence_basis = convergence$convergence_basis,
      optimizer_converged = convergence$optimizer_converged,
      raw_gradient_ok = convergence$raw_gradient_ok,
      average_gradient_ok = convergence$average_gradient_ok,
      information_scaled_convergence_ok = convergence$information_scaled_convergence_ok,
      H_red_positive_definite = convergence$H_red_positive_definite,
      reduced_gradient_max_abs = convergence$reduced_gradient_max_abs,
      average_reduced_gradient_max_abs = convergence$average_reduced_gradient_max_abs,
      max_abs_newton_step_scaled = convergence$max_abs_newton_step_scaled,
      gradient_quadratic_lr_bound = convergence$gradient_quadratic_lr_bound,
      retry_triggered = retry_triggered,
      retry_used = isTRUE(fit_con$diagnostics$retry_used),
      optimizer_attempts = fit_con$diagnostics$optimizer_attempts,
      use_average_objective = isTRUE(fit_con$diagnostics$use_average_objective),
      constraint_residual_max_abs = constraint_residual_max_abs,
      constraint_satisfied = constraint_satisfied,
      raw_lr_initial = fit_pair$diagnostics$raw_lr_initial,
      raw_lr_final = fit_pair$diagnostics$raw_lr_final,
      unconstrained_repair_attempted = fit_pair$diagnostics$unconstrained_repair_attempted,
      unconstrained_repair_applied = fit_pair$diagnostics$unconstrained_repair_applied
    )

    score_res <- compute_adj_score(
      constrained = fit_con,
      L = def$L,
      ridge_factor = ridge_factor
    )
    lrt_res <- compute_profile_lrt(
      unconstrained = fit_unconstrained,
      constrained = fit_con,
      L = def$L,
      b = def$b,
      ridge_factor = ridge_factor,
      cond_limit = cond_limit,
      lrt_gate_override = lrt_gate_override
    )
    wald_res <- pmg_robust_wald(
      unconstrained = fit_unconstrained,
      L = def$L,
      b = def$b,
      ridge_factor = ridge_factor,
      cond_limit = cond_limit
    )
    df <- nrow(def$L)

    results[[idx]] <- list(
      label = def$label,
      df = df,
      contrast = def$L,
      rhs = def$b,
      unconstrained = fit_unconstrained,
      constrained = fit_con,
      fit_pair_diagnostics = fit_pair$diagnostics,
      adj = list(
        wald = wald_res,
        score = score_res,
        lr = lrt_res
      ),
      p_values = list(
        lr = if (is.finite(lrt_res$wP_adjusted)) pchisq(max(lrt_res$wP_adjusted, 0), df = df, lower.tail = FALSE) else NA_real_,
        wald = if (is.finite(wald_res$stat)) pchisq(max(wald_res$stat, 0), df = df, lower.tail = FALSE) else NA_real_,
        score = if (is.finite(score_res$stat)) pchisq(max(score_res$stat, 0), df = df, lower.tail = FALSE) else NA_real_
      )
    )
    names_vec[[idx]] <- def$label
  }

  names(results) <- names_vec
  results
}
