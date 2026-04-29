# Internal developer note:
# These helpers add a profile-likelihood computation layer for the exact PGMM
# path without extending the current Wald / score / LRT wrappers directly.
# The layer is intentionally narrow in v1:
# - exact PGMM only;
# - scalar fixed-effect targets;
# - pointwise or grid evaluation at user-supplied profile values.
#
# Workflow:
# 1. `pgmm_fit_unconstrained()` caches the exact unconstrained PGMM fit and the
#    scalar target definition.
# 2. `pgmm_relative_profile_coef_grid()` is a thin internal convenience wrapper
#    for the common "one coefficient index plus supplied grid" workflow.
# 3. `pgmm_profile_point()` evaluates one constrained PGMM fit at a fixed
#    profile value and returns score / observed-information blocks in the
#    scalar `(psi, eta)` parametrization.
# 4. `pgmm_relative_profile_point()` turns one constrained point into ordinary
#    and robust relative-profile quantities.
# 5. `pgmm_relative_profile_grid()` evaluates a supplied grid and standardizes
#    the ordinary and robust curves so their empirical maxima over the grid are
#    equal to 1.
#
# Relation to the existing robust PGMM LRT code:
# `compute_profile_lrt()` is test-oriented and works on variance-scale objects
# built from the same constrained exact PGMM derivatives. The helpers below
# instead return the scalar profile-sensitivity / variability quantities from
# the implementation note. For scalar targets, the resulting robust adjusted
# LR quantity `wP_star` agrees with `compute_profile_lrt(... )$wP_adjusted`
# whenever the same constrained point is evaluated and the adjustment is
# well-defined. No exported user-facing wrapper is added here yet.

pgmm_profile_target <- function(ctx, psi_index = NULL, c = NULL) {
  if (!is.null(psi_index) && !is.null(c)) {
    stop("Supply only one of psi_index or c.")
  }

  nb <- ctx$nb
  if (is.null(c)) {
    if (is.null(psi_index) || length(psi_index) != 1L || !is.finite(psi_index)) {
      stop("psi_index must be a single finite index.")
    }
    psi_index <- as.integer(psi_index)
    if (psi_index < 1L || psi_index > nb) {
      stop("psi_index is out of range for the PGMM fixed-effect vector.")
    }
    c_vec <- rep(0, nb)
    c_vec[psi_index] <- 1
  } else {
    c_vec <- as.numeric(c)
    if (length(c_vec) != nb) {
      stop("Contrast c must have length ctx$nb.")
    }
    if (!all(is.finite(c_vec))) {
      stop("Contrast c must be finite.")
    }
    if (sum(abs(c_vec)) <= 0) {
      stop("Contrast c must not be the zero vector.")
    }
    psi_index <- NA_integer_
    nz <- which(abs(c_vec) > 0)
    if (length(nz) == 1L && abs(c_vec[nz] - 1) <= 1e-12) {
      psi_index <- nz
    }
  }

  target <- scalar_profile_reparameterization(matrix(c_vec, nrow = 1L))

  beta_names <- colnames(ctx$pred)
  label <- "contrast"
  if (!is.na(psi_index) && !is.null(beta_names) && length(beta_names) >= psi_index) {
    label <- beta_names[[psi_index]]
  } else if (!is.na(psi_index)) {
    label <- paste0("beta", psi_index)
  }

  list(
    c = target$c,
    L = target$L,
    a = target$a,
    K = target$K,
    free_dim = target$free_dim,
    psi_index = psi_index,
    label = label
  )
}

pgmm_profile_wrap_unconstrained <- function(fit, target) {
  if (is.null(fit$theta) || is.null(fit$theta$beta)) {
    stop("Unconstrained PGMM fit must provide theta$beta.")
  }

  psi_hat <- as.numeric(crossprod(target$c, fit$theta$beta))

  list(
    fit = fit,
    theta = fit$theta,
    psi_hat = psi_hat,
    loglik = fit$loglik,
    target = target,
    convergence = fit$convergence,
    optimizer = fit$optimizer
  )
}

pgmm_profile_normalize_unconstrained <- function(unconstrained,
                                                 ctx,
                                                 gene_index = NULL,
                                                 posv = NULL,
                                                 psi_index = NULL,
                                                 c = NULL,
                                                 start = NULL,
                                                 control = NULL,
                                                 retry_failed_profiles = FALSE,
                                                 use_average_objective = FALSE) {
  if (!is.null(unconstrained) &&
    !is.null(unconstrained$fit) &&
    !is.null(unconstrained$target) &&
    !is.null(unconstrained$psi_hat)) {
    if (!is.null(psi_index) || !is.null(c)) {
      target_check <- pgmm_profile_target(ctx, psi_index = psi_index, c = c)
      if (!isTRUE(all.equal(as.numeric(unconstrained$target$c), as.numeric(target_check$c), tolerance = 0))) {
        stop("Cached unconstrained profile fit target does not match psi_index / c.")
      }
    }
    return(unconstrained)
  }

  target <- pgmm_profile_target(ctx, psi_index = psi_index, c = c)
  if (is.null(unconstrained)) {
    if (is.null(gene_index) || is.null(posv)) {
      stop("gene_index and posv are required when unconstrained is not supplied.")
    }
    unconstrained <- fit_gene_pmg_unconstrained(
      gene_index = gene_index,
      posv = posv,
      ctx = ctx,
      start = start,
      control = control,
      retry = retry_failed_profiles,
      use_average_objective = use_average_objective
    )
  }

  pgmm_profile_wrap_unconstrained(unconstrained, target)
}

pgmm_profile_standardize_curve <- function(stat, tol = 1e-8) {
  out <- rep(NA_real_, length(stat))
  valid <- is.finite(stat) & (stat >= -tol)
  if (!any(valid)) {
    return(out)
  }

  stat_adj <- stat
  stat_adj[valid & stat_adj < 0] <- 0
  min_stat <- min(stat_adj[valid])
  out[valid] <- exp(-0.5 * (stat_adj[valid] - min_stat))
  out
}

pgmm_clip_raw_lr <- function(raw_lr, tol = 1e-8) {
  lr_tol <- max(tol, 1e-6)
  clipped <- is.finite(raw_lr) && raw_lr < 0 && abs(raw_lr) <= lr_tol
  list(
    wP = if (isTRUE(clipped)) 0 else raw_lr,
    raw_lr_clipped = clipped,
    tol = lr_tol
  )
}

pgmm_profile_accept_constrained <- function(point,
                                            wP_nonnegative,
                                            accept_nonzero_optimizer_code = TRUE,
                                            gradient_tol = 1e-2,
                                            average_gradient_tol = 1e-4) {
  diag <- point$diagnostics
  fit_c <- point$constrained
  theta <- fit_c$theta
  full_par <- c(theta$beta, theta$subVar)
  reduced_grad <- if (!is.null(diag$reduced_gradient_max_abs)) diag$reduced_gradient_max_abs else NA_real_
  n_subjects <- if (!is.null(diag$n_subjects)) diag$n_subjects else NA_real_
  objective <- if (!is.null(diag$objective)) diag$objective else NA_real_
  avg_reduced_grad <- if (is.finite(reduced_grad) && is.finite(n_subjects) && n_subjects > 0) {
    reduced_grad / n_subjects
  } else {
    NA_real_
  }

  finite_fit <- is.finite(point$loglik_profile) &&
    all(is.finite(full_par)) &&
    is.finite(objective)
  constraint_ok <- isTRUE(diag$constraint_satisfied)
  h_eta_ok <- isTRUE(diag$H_etaeta_invertible)
  boundary_ok <- !isTRUE(diag$boundary_reduced_any)
  lr_ok <- isTRUE(wP_nonnegative)
  gradient_ok <- is.finite(reduced_grad) &&
    (reduced_grad <= gradient_tol ||
      (is.finite(avg_reduced_grad) && avg_reduced_grad <= average_gradient_tol))
  optimizer_ok <- isTRUE(diag$converged) ||
    (isTRUE(accept_nonzero_optimizer_code) && gradient_ok)

  accepted <- finite_fit &&
    constraint_ok &&
    h_eta_ok &&
    boundary_ok &&
    lr_ok &&
    optimizer_ok

  failure_reason <- NA_character_
  if (!accepted) {
    failure_reason <- if (!finite_fit) {
      "nonfinite_profile_fit"
    } else if (!constraint_ok) {
      "constraint_not_satisfied"
    } else if (!lr_ok) {
      "negative_raw_lr"
    } else if (!h_eta_ok) {
      "singular_H_etaeta"
    } else if (!boundary_ok) {
      "boundary_constrained"
    } else if (!optimizer_ok) {
      "nonconverged_constrained"
    } else {
      "profile_point_failed"
    }
  }

  list(
    accepted = accepted,
    optimizer_warning = accepted && !isTRUE(diag$converged),
    finite_fit = finite_fit,
    constraint_ok = constraint_ok,
    h_eta_ok = h_eta_ok,
    boundary_ok = boundary_ok,
    lr_ok = lr_ok,
    gradient_ok = gradient_ok,
    avg_reduced_gradient_max_abs = avg_reduced_grad,
    failure_reason = failure_reason
  )
}

pgmm_fit_unconstrained <- function(gene_index,
                                   posv,
                                   ctx,
                                   psi_index = NULL,
                                   c = NULL,
                                   start = NULL,
                                   control = NULL,
                                   retry_failed_profiles = FALSE,
                                   use_average_objective = FALSE) {
  target <- pgmm_profile_target(ctx, psi_index = psi_index, c = c)
  fit <- fit_gene_pmg_unconstrained(
    gene_index = gene_index,
    posv = posv,
    ctx = ctx,
    start = start,
    control = control,
    retry = retry_failed_profiles,
    use_average_objective = use_average_objective
  )
  pgmm_profile_wrap_unconstrained(fit, target)
}

pgmm_relative_profile_coef_grid <- function(ctx,
                                            gene_index,
                                            posv,
                                            psi_index,
                                            psi_values,
                                            compute_robust = TRUE,
                                            ridge_factor = 1e-8,
                                            cond_limit = 1e8,
                                            tol = 1e-8,
                                            accept_nonzero_optimizer_code = TRUE,
                                            retry_failed_profiles = TRUE,
                                            use_average_objective = FALSE,
                                            profile_gradient_tol = 1e-2,
                                            profile_average_gradient_tol = 1e-4) {
  unc <- pgmm_fit_unconstrained(
    gene_index = gene_index,
    posv = posv,
    ctx = ctx,
    psi_index = psi_index,
    retry_failed_profiles = retry_failed_profiles,
    use_average_objective = use_average_objective
  )

  pgmm_relative_profile_grid(
    ctx = ctx,
    gene_index = gene_index,
    posv = posv,
    psi_values = psi_values,
    unconstrained = unc,
    compute_robust = compute_robust,
    ridge_factor = ridge_factor,
    cond_limit = cond_limit,
    tol = tol,
    accept_nonzero_optimizer_code = accept_nonzero_optimizer_code,
    retry_failed_profiles = retry_failed_profiles,
    use_average_objective = use_average_objective,
    profile_gradient_tol = profile_gradient_tol,
    profile_average_gradient_tol = profile_average_gradient_tol
  )
}

pgmm_profile_point <- function(ctx,
                               gene_index,
                               posv,
                               psi_value,
                               psi_index = NULL,
                               c = NULL,
                               unconstrained = NULL,
                               include_robust = TRUE,
                               ridge_factor = 1e-8,
                               cond_limit = 1e8,
                               tol = 1e-8,
                               retry_failed_profiles = TRUE,
                               use_average_objective = FALSE,
                               warm_start = NULL) {
  psi_value <- as.numeric(psi_value)
  if (length(psi_value) != 1L || !is.finite(psi_value)) {
    stop("psi_value must be a single finite number.")
  }

  unc <- pgmm_profile_normalize_unconstrained(
    unconstrained = unconstrained,
    ctx = ctx,
    gene_index = gene_index,
    posv = posv,
    psi_index = psi_index,
    c = c,
    retry_failed_profiles = retry_failed_profiles,
    use_average_objective = use_average_objective
  )
  target <- unc$target

  fit_c <- fit_gene_pmg_constrained(
    gene_index = gene_index,
    posv = posv,
    ctx = ctx,
    L = target$L,
    b = psi_value,
    unconstrained_fit = unc$fit,
    warm_start = warm_start,
    retry = retry_failed_profiles,
    use_average_objective = use_average_objective
  )

  lr_tol <- max(tol, 1e-6)
  raw_lr_initial <- 2 * (unc$loglik - fit_c$loglik)
  if (isTRUE(retry_failed_profiles) &&
    is.finite(raw_lr_initial) &&
    raw_lr_initial < -lr_tol) {
    retry_start <- c(fit_c$theta$beta, fit_c$theta$subVar)
    fit_u_retry <- fit_gene_pmg_unconstrained(
      gene_index = gene_index,
      posv = posv,
      ctx = ctx,
      start = retry_start,
      retry = TRUE,
      use_average_objective = use_average_objective
    )
    if (is.finite(fit_u_retry$loglik) && fit_u_retry$loglik > unc$loglik + lr_tol / 2) {
      unc <- pgmm_profile_wrap_unconstrained(fit_u_retry, target)
      fit_c <- fit_gene_pmg_constrained(
        gene_index = gene_index,
        posv = posv,
        ctx = ctx,
        L = target$L,
        b = psi_value,
        unconstrained_fit = unc$fit,
        warm_start = fit_c,
        retry = TRUE,
        use_average_objective = use_average_objective
      )
    }
  }

  profile_state <- compute_scalar_profile_quantities(
    full_score = -fit_c$score,
    full_per_subject_scores = -fit_c$per_subject_gradients,
    observed_info_full = fit_c$observed_info,
    L = target$L,
    ridge_factor = ridge_factor,
    cond_limit = cond_limit,
    tol = tol
  )

  constraint_residual <- as.numeric(target$L %*% fit_c$theta$beta - psi_value)
  constraint_satisfied <- is.finite(constraint_residual) && abs(constraint_residual) <= sqrt(.Machine$double.eps)
  constrained_diagnostics <- fit_c$diagnostics
  if (is.null(constrained_diagnostics)) {
    constrained_diagnostics <- list()
  }

  list(
    psi = psi_value,
    psi_hat = unc$psi_hat,
    target = target,
    unconstrained = unc,
    constrained = fit_c,
    loglik_unconstrained = unc$loglik,
    loglik_profile = fit_c$loglik,
    score_blocks = profile_state$score_blocks,
    per_subject_scores = profile_state$per_subject_scores,
    observed_info_blocks = profile_state$observed_info_blocks,
    profile = if (isTRUE(include_robust)) {
      list(
        subject_scores = profile_state$profile$subject_scores,
        score = profile_state$profile$score,
        variability = profile_state$profile$variability,
        sensitivity = profile_state$profile$sensitivity,
        godambe = profile_state$profile$godambe,
        scale = profile_state$profile$scale
      )
    } else {
      NULL
    },
    diagnostics = list(
      converged = isTRUE(fit_c$convergence > 0L),
      convergence = fit_c$convergence,
      optimizer = fit_c$optimizer,
      raw_convergence = constrained_diagnostics$raw_convergence,
      raw_message = constrained_diagnostics$raw_message,
      objective = constrained_diagnostics$objective,
      full_gradient_max_abs = constrained_diagnostics$full_gradient_max_abs,
      full_gradient_l2 = constrained_diagnostics$full_gradient_l2,
      reduced_gradient_max_abs = constrained_diagnostics$reduced_gradient_max_abs,
      reduced_gradient_l2 = constrained_diagnostics$reduced_gradient_l2,
      active_reduced_all = constrained_diagnostics$active_reduced_all,
      boundary_reduced_any = constrained_diagnostics$boundary_reduced_any,
      retry_used = constrained_diagnostics$retry_used,
      use_average_objective = constrained_diagnostics$use_average_objective,
      n_subjects = ctx$k,
      active_reduced = constrained_diagnostics$active_reduced,
      at_lower_reduced = constrained_diagnostics$at_lower_reduced,
      at_upper_reduced = constrained_diagnostics$at_upper_reduced,
      constraint_residual = constraint_residual,
      constraint_satisfied = constraint_satisfied,
      H_etaeta_invertible = isTRUE(profile_state$diagnostics$H_etaeta_invertible),
      H_etaeta_ridge = profile_state$diagnostics$H_etaeta_ridge,
      Hpsi_positive = isTRUE(profile_state$diagnostics$Hpsi_positive),
      Jpsi_positive = isTRUE(profile_state$diagnostics$Jpsi_positive)
    )
  )
}

pgmm_relative_profile_point <- function(ctx,
                                        gene_index,
                                        posv,
                                        psi_value,
                                        psi_index = NULL,
                                        c = NULL,
                                        unconstrained = NULL,
                                        compute_robust = TRUE,
                                        ridge_factor = 1e-8,
                                        cond_limit = 1e8,
                                        tol = 1e-8,
                                        accept_nonzero_optimizer_code = TRUE,
                                        retry_failed_profiles = TRUE,
                                        use_average_objective = FALSE,
                                        warm_start = NULL,
                                        profile_gradient_tol = 1e-2,
                                        profile_average_gradient_tol = 1e-4) {
  point <- pgmm_profile_point(
    ctx = ctx,
    gene_index = gene_index,
    posv = posv,
    psi_value = psi_value,
    psi_index = psi_index,
    c = c,
    unconstrained = unconstrained,
    include_robust = compute_robust,
    ridge_factor = ridge_factor,
    cond_limit = cond_limit,
    tol = tol,
    retry_failed_profiles = retry_failed_profiles,
    use_average_objective = use_average_objective,
    warm_start = warm_start
  )

  raw_lr <- 2 * (point$loglik_unconstrained - point$loglik_profile)
  raw_lr_info <- pgmm_clip_raw_lr(raw_lr, tol = tol)
  wP <- raw_lr_info$wP
  lr_tol <- raw_lr_info$tol
  raw_lr_clipped <- raw_lr_info$raw_lr_clipped
  wP_nonnegative <- is.finite(wP) && wP >= -lr_tol
  relLik <- if (wP_nonnegative) exp(-0.5 * max(wP, 0)) else NA_real_
  acceptance <- pgmm_profile_accept_constrained(
    point = point,
    wP_nonnegative = wP_nonnegative,
    accept_nonzero_optimizer_code = accept_nonzero_optimizer_code,
    gradient_tol = profile_gradient_tol,
    average_gradient_tol = profile_average_gradient_tol
  )

  UP <- NA_real_
  Jpsi <- NA_real_
  Hpsi <- NA_real_
  Gpsi <- NA_real_
  wU <- NA_real_
  qP <- NA_real_
  scale <- NA_real_
  wP_star <- NA_real_
  relLik_star <- NA_real_
  wP_star_nonnegative <- FALSE
  robust_available <- FALSE
  failure_reason <- NA_character_

  if (isTRUE(compute_robust) && !is.null(point$profile)) {
    UP <- point$profile$score
    Jpsi <- point$profile$variability
    Hpsi <- point$profile$sensitivity
    Gpsi <- point$profile$godambe
    scale <- point$profile$scale

    if (is.finite(UP) && is.finite(Jpsi) && Jpsi > tol) {
      wU <- (UP^2) / Jpsi
      if (wU < 0 && abs(wU) <= tol) {
        wU <- 0
      }
    }
    if (is.finite(UP) && is.finite(Hpsi) && Hpsi > tol) {
      qP <- (UP^2) / Hpsi
      if (qP < 0 && abs(qP) <= tol) {
        qP <- 0
      }
    }

    if ((!is.finite(scale) || scale <= 0) &&
      is.finite(wU) &&
      is.finite(qP) &&
      abs(qP) > tol) {
      scale <- wU / qP
    }

    if (is.finite(scale) && scale > 0 && wP_nonnegative) {
      wP_star <- scale * max(wP, 0)
      if (is.finite(wP_star) && wP_star < 0 && abs(wP_star) <= lr_tol) {
        wP_star <- 0
      }
    }

    wP_star_nonnegative <- is.finite(wP_star) && wP_star >= -lr_tol
    relLik_star <- if (wP_star_nonnegative) exp(-0.5 * max(wP_star, 0)) else NA_real_

    robust_available <- isTRUE(acceptance$accepted) &&
      isTRUE(point$diagnostics$H_etaeta_invertible) &&
      isTRUE(point$diagnostics$Hpsi_positive) &&
      isTRUE(point$diagnostics$Jpsi_positive) &&
      wP_nonnegative &&
      wP_star_nonnegative

    if (!robust_available) {
      if (!isTRUE(acceptance$accepted)) {
        failure_reason <- acceptance$failure_reason
      } else if (!wP_nonnegative) {
        failure_reason <- "negative_raw_lr"
      } else if (!isTRUE(point$diagnostics$H_etaeta_invertible)) {
        failure_reason <- "singular_H_etaeta"
      } else if (!isTRUE(point$diagnostics$Hpsi_positive)) {
        failure_reason <- "nonpositive_Hpsi"
      } else if (!isTRUE(point$diagnostics$Jpsi_positive)) {
        failure_reason <- "nonpositive_Jpsi"
      } else if (!is.finite(scale) || scale <= 0) {
        failure_reason <- "nonpositive_scale"
      } else if (!is.finite(wP_star)) {
        failure_reason <- "undefined_adjustment"
      } else {
        failure_reason <- "profile_point_failed"
      }
    }
  }

  point$relative_profile <- list(
    raw_lr = raw_lr,
    wP = wP,
    raw_lr_clipped = raw_lr_clipped,
    relLik = relLik,
    UP = UP,
    Jpsi = Jpsi,
    Hpsi = Hpsi,
    Gpsi = Gpsi,
    wU = wU,
    qP = qP,
    scale = scale,
    wP_star = wP_star,
    relLik_star = relLik_star,
    wP_nonnegative = wP_nonnegative,
    wP_star_nonnegative = wP_star_nonnegative,
    robust_available = robust_available,
    profile_accepted = acceptance$accepted,
    optimizer_warning = acceptance$optimizer_warning,
    acceptance_failure_reason = acceptance$failure_reason,
    avg_reduced_gradient_max_abs = acceptance$avg_reduced_gradient_max_abs,
    failure_reason = failure_reason
  )
  point
}

pgmm_relative_profile_grid <- function(ctx,
                                       gene_index,
                                       posv,
                                       psi_values,
                                       psi_index = NULL,
                                       c = NULL,
                                       unconstrained = NULL,
                                       compute_robust = TRUE,
                                       ridge_factor = 1e-8,
                                       cond_limit = 1e8,
                                       tol = 1e-8,
                                       accept_nonzero_optimizer_code = TRUE,
                                       retry_failed_profiles = TRUE,
                                       use_average_objective = FALSE,
                                       warm_start = TRUE,
                                       profile_gradient_tol = 1e-2,
                                       profile_average_gradient_tol = 1e-4) {
  psi_values <- as.numeric(psi_values)
  if (length(psi_values) == 0L) {
    stop("psi_values must contain at least one profile point.")
  }

  unc <- pgmm_profile_normalize_unconstrained(
    unconstrained = unconstrained,
    ctx = ctx,
    gene_index = gene_index,
    posv = posv,
    psi_index = psi_index,
    c = c,
    retry_failed_profiles = retry_failed_profiles,
    use_average_objective = use_average_objective
  )

  rows <- vector("list", length(psi_values))
  eval_order <- if (isTRUE(warm_start)) {
    order(abs(psi_values - unc$psi_hat), psi_values)
  } else {
    seq_along(psi_values)
  }
  current_warm_start <- NULL
  for (idx in eval_order) {
    point <- pgmm_relative_profile_point(
      ctx = ctx,
      gene_index = gene_index,
      posv = posv,
      psi_value = psi_values[[idx]],
      unconstrained = unc,
      compute_robust = compute_robust,
      ridge_factor = ridge_factor,
      cond_limit = cond_limit,
      tol = tol,
      accept_nonzero_optimizer_code = accept_nonzero_optimizer_code,
      retry_failed_profiles = retry_failed_profiles,
      use_average_objective = use_average_objective,
      warm_start = current_warm_start,
      profile_gradient_tol = profile_gradient_tol,
      profile_average_gradient_tol = profile_average_gradient_tol
    )
    if (isTRUE(point$relative_profile$profile_accepted) && !is.null(point$constrained)) {
      current_warm_start <- point$constrained
    }

    rel <- point$relative_profile
    rows[[idx]] <- data.frame(
      psi = point$psi,
      psi_hat = point$psi_hat,
      target_label = point$target$label,
      psi_index = point$target$psi_index,
      loglik_unconstrained = point$loglik_unconstrained,
      loglik_profile = point$loglik_profile,
      raw_lr = rel$raw_lr,
      wP = rel$wP,
      raw_lr_clipped = rel$raw_lr_clipped,
      relLik = rel$relLik,
      UP = rel$UP,
      Jpsi = rel$Jpsi,
      Hpsi = rel$Hpsi,
      Gpsi = rel$Gpsi,
      wU = rel$wU,
      qP = rel$qP,
      scale = rel$scale,
      wP_star = rel$wP_star,
      relLik_star = rel$relLik_star,
      converged = point$diagnostics$converged,
      convergence = point$diagnostics$convergence,
      optimizer = point$diagnostics$optimizer,
      raw_convergence = point$diagnostics$raw_convergence,
      raw_message = point$diagnostics$raw_message,
      objective = point$diagnostics$objective,
      full_gradient_max_abs = point$diagnostics$full_gradient_max_abs,
      full_gradient_l2 = point$diagnostics$full_gradient_l2,
      reduced_gradient_max_abs = point$diagnostics$reduced_gradient_max_abs,
      reduced_gradient_l2 = point$diagnostics$reduced_gradient_l2,
      avg_reduced_gradient_max_abs = rel$avg_reduced_gradient_max_abs,
      active_reduced_all = point$diagnostics$active_reduced_all,
      boundary_reduced_any = point$diagnostics$boundary_reduced_any,
      retry_used = point$diagnostics$retry_used,
      use_average_objective = point$diagnostics$use_average_objective,
      constraint_residual = point$diagnostics$constraint_residual,
      constraint_satisfied = point$diagnostics$constraint_satisfied,
      H_etaeta_invertible = point$diagnostics$H_etaeta_invertible,
      H_etaeta_ridge = point$diagnostics$H_etaeta_ridge,
      Hpsi_positive = point$diagnostics$Hpsi_positive,
      Jpsi_positive = point$diagnostics$Jpsi_positive,
      wP_nonnegative = rel$wP_nonnegative,
      wP_star_nonnegative = rel$wP_star_nonnegative,
      robust_available = rel$robust_available,
      profile_accepted = rel$profile_accepted,
      optimizer_warning = rel$optimizer_warning,
      acceptance_failure_reason = rel$acceptance_failure_reason,
      failure_reason = rel$failure_reason,
      stringsAsFactors = FALSE
    )
  }

  grid <- do.call(rbind, rows)
  grid$relLik_std <- pgmm_profile_standardize_curve(grid$wP, tol = tol)
  grid$relLik_star_std <- pgmm_profile_standardize_curve(grid$wP_star, tol = tol)

  attr(grid, "target") <- unc$target
  attr(grid, "unconstrained") <- unc
  grid
}
