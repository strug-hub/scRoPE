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
                                                 c = NULL) {
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
      ctx = ctx
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

pgmm_fit_unconstrained <- function(gene_index,
                                   posv,
                                   ctx,
                                   psi_index = NULL,
                                   c = NULL) {
  target <- pgmm_profile_target(ctx, psi_index = psi_index, c = c)
  fit <- fit_gene_pmg_unconstrained(
    gene_index = gene_index,
    posv = posv,
    ctx = ctx
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
                                            tol = 1e-8) {
  unc <- pgmm_fit_unconstrained(
    gene_index = gene_index,
    posv = posv,
    ctx = ctx,
    psi_index = psi_index
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
    tol = tol
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
                               tol = 1e-8) {
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
    c = c
  )
  target <- unc$target

  fit_c <- fit_gene_pmg_constrained(
    gene_index = gene_index,
    posv = posv,
    ctx = ctx,
    L = target$L,
    b = psi_value,
    unconstrained_fit = unc$fit
  )

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
                                        tol = 1e-8) {
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
    tol = tol
  )

  wP <- 2 * (point$loglik_unconstrained - point$loglik_profile)
  if (is.finite(wP) && wP < 0 && abs(wP) <= tol) {
    wP <- 0
  }
  wP_nonnegative <- is.finite(wP) && wP >= -tol
  relLik <- if (wP_nonnegative) exp(-0.5 * max(wP, 0)) else NA_real_

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
      if (is.finite(wP_star) && wP_star < 0 && abs(wP_star) <= tol) {
        wP_star <- 0
      }
    }

    wP_star_nonnegative <- is.finite(wP_star) && wP_star >= -tol
    relLik_star <- if (wP_star_nonnegative) exp(-0.5 * max(wP_star, 0)) else NA_real_

    robust_available <- isTRUE(point$diagnostics$converged) &&
      isTRUE(point$diagnostics$H_etaeta_invertible) &&
      isTRUE(point$diagnostics$Hpsi_positive) &&
      isTRUE(point$diagnostics$Jpsi_positive) &&
      wP_nonnegative &&
      wP_star_nonnegative

    if (!robust_available) {
      if (!isTRUE(point$diagnostics$converged)) {
        failure_reason <- "nonconverged_constrained"
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
    wP = wP,
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
                                       tol = 1e-8) {
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
    c = c
  )

  rows <- vector("list", length(psi_values))
  for (idx in seq_along(psi_values)) {
    point <- pgmm_relative_profile_point(
      ctx = ctx,
      gene_index = gene_index,
      posv = posv,
      psi_value = psi_values[[idx]],
      unconstrained = unc,
      compute_robust = compute_robust,
      ridge_factor = ridge_factor,
      cond_limit = cond_limit,
      tol = tol
    )

    rel <- point$relative_profile
    rows[[idx]] <- data.frame(
      psi = point$psi,
      psi_hat = point$psi_hat,
      target_label = point$target$label,
      psi_index = point$target$psi_index,
      loglik_unconstrained = point$loglik_unconstrained,
      loglik_profile = point$loglik_profile,
      wP = rel$wP,
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
      constraint_residual = point$diagnostics$constraint_residual,
      constraint_satisfied = point$diagnostics$constraint_satisfied,
      H_etaeta_invertible = point$diagnostics$H_etaeta_invertible,
      H_etaeta_ridge = point$diagnostics$H_etaeta_ridge,
      Hpsi_positive = point$diagnostics$Hpsi_positive,
      Jpsi_positive = point$diagnostics$Jpsi_positive,
      wP_nonnegative = rel$wP_nonnegative,
      wP_star_nonnegative = rel$wP_star_nonnegative,
      robust_available = rel$robust_available,
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
