hl_fit_unconstrained <- function(ctx, gene_index, count, posv = NULL) {
  if (is.null(posv)) {
    posv <- call_posindy(count, gene_index - 1L, ctx$nind)
  }
  core <- hl_fit_gene_core(gene_index, posv, ctx)
  ps <- core$hl_state$per_subject_stats
  if (!is.null(ps)) {
    score <- hl_score_from_stats(ps, ctx$intcol, ctx$nb, ctx$k)
    core$hl_state$score <- score$U
    core$hl_state$score_subject <- score$per_subject
    core$hl_state$meat <- hl_meat_from_scores(score$per_subject)
    hess_blocks <- hl_hessian_blocks_from_stats(ps)
    eff <- hl_efficient_components(hess_blocks$H_bb, hess_blocks$H_bl, hess_blocks$H_ll, score$U, ctx$nb)
    core$hl_state$H_blocks <- hess_blocks
    core$hl_state$efficient_score <- eff$score_beta_eff
    core$hl_state$H_schur <- eff$H_beta_beta_schur
    core$hl_state$H_obs <- hl_hessian_from_stats(ps)
  }
  core$beta <- core$repml$beta
  core$dispersion <- core$vare
  core$nind <- ctx$nind
  core$k <- ctx$k
  core$posv <- posv
  core
}

hl_fit_constrained <- function(ctx,
                               gene_index,
                               count,
                               C,
                               posv = NULL,
                               unconstrained_fit = NULL,
                               compute_adj = FALSE,
                               ridge_factor = 1e-8,
                               cond_limit = 1e8,
                               rhs = NULL) {
  if (is.null(posv)) {
    posv <- call_posindy(count, gene_index - 1L, ctx$nind)
  }
  if (is.null(dim(C))) {
    C <- matrix(C, nrow = 1)
  }
  if (!is.matrix(C)) {
    stop("C must be a matrix.")
  }
  if (ncol(C) != ctx$nb) {
    stop("Constraint matrix column count must match the number of predictors.")
  }
  if (nrow(C) == 0) {
    return(hl_fit_unconstrained(ctx, gene_index, count, posv))
  }

  if (is.null(rhs)) {
    rhs <- rep(0, nrow(C))
  }
  if (length(rhs) != nrow(C)) {
    stop("Length of rhs must match number of rows of C.")
  }

  constr <- constraint_reparameterization(C, rhs)
  beta_star <- constr$beta_star
  K <- constr$null_basis
  if (is.null(K)) {
    K <- matrix(0, nrow = ctx$nb, ncol = 0)
  }
  free_dim <- ncol(K)
  intercept_idx <- ctx$intcol

  compose_beta <- function(gamma) {
    if (free_dim == 0) {
      beta_star
    } else {
      as.vector(beta_star + K %*% gamma)
    }
  }

  build_theta_full <- function(theta) {
    if (free_dim == 0) {
      gamma <- numeric(0)
      disp <- theta
    } else {
      gamma <- theta[seq_len(free_dim)]
      disp <- theta[-seq_len(free_dim)]
    }
    c(compose_beta(gamma), disp)
  }

  lbfgs_obj <- function(theta) {
    full_para <- build_theta_full(theta)
    res <- ptmg_ll_der(
      full_para,
      X = ctx$pred,
      offset = ctx$offset,
      Y = posv$Y,
      n_one = posv$n_onetwo,
      ytwo = posv$ytwo,
      fid = ctx$fid,
      cumsumy = ctx$cumsumy[gene_index, ],
      posind = ctx$posind[[gene_index]],
      posindy = posv$posindy,
      nb = ctx$nb,
      nind = ctx$nind,
      k = ctx$k
    )
    grad_full <- res$gradient
    grad_beta <- grad_full[seq_len(ctx$nb)]
    grad_gamma <- if (free_dim == 0) {
      numeric(0)
    } else {
      as.vector(t(K) %*% grad_beta)
    }
    sigma_idx <- ctx$nb + 1L
    grad_full[sigma_idx] <- grad_full[sigma_idx] + 0.5 * grad_beta[intercept_idx]
    grad_disp <- grad_full[(ctx$nb + 1):(ctx$nb + 2)]
    list(
      objective = res$objective,
      gradient = c(grad_gamma, grad_disp)
    )
  }

  lower_theta <- c(if (free_dim == 0) numeric(0) else rep(-100, free_dim), ctx$min[1], ctx$min[2])
  upper_theta <- c(if (free_dim == 0) numeric(0) else rep(100, free_dim), ctx$max[1], ctx$max[2])
  gamma_init <- if (free_dim == 0) numeric(0) else rep(0, free_dim)
  theta_init <- c(gamma_init, 1, ctx$cell_init)

  opt_res <- tryCatch(
    {
      ref <- lbfgs(
        theta_init,
        lbfgs_obj,
        lower = lower_theta,
        upper = upper_theta,
        control = list(ftol_abs = ctx$eps)
      )
      list(theta = ref$par, conv_flag = 1)
    },
    error = function(e) {
      obj_fn <- function(theta) {
        full_para <- build_theta_full(theta)
        ptmg_ll(
          full_para,
          X = ctx$pred,
          offset = ctx$offset,
          Y = posv$Y,
          n_one = posv$n_onetwo,
          ytwo = posv$ytwo,
          fam = ctx$id,
          fid = ctx$fid,
          cumsumy = ctx$cumsumy[gene_index, ],
          posind = ctx$posind[[gene_index]],
          posindy = posv$posindy,
          nb = ctx$nb,
          nind = ctx$nind,
          k = ctx$k
        )
      }
      grad_fn <- function(theta) {
        full_para <- build_theta_full(theta)
        gr_full <- ptmg_der(
          full_para,
          X = ctx$pred,
          offset = ctx$offset,
          Y = posv$Y,
          n_one = posv$n_onetwo,
          ytwo = posv$ytwo,
          fam = ctx$id,
          fid = ctx$fid,
          cumsumy = ctx$cumsumy[gene_index, ],
          posind = ctx$posind[[gene_index]],
          posindy = posv$posindy,
          nb = ctx$nb,
          nind = ctx$nind,
          k = ctx$k
        )
        grad_beta <- gr_full[seq_len(ctx$nb)]
        grad_gamma <- if (free_dim == 0) {
          numeric(0)
        } else {
          as.vector(t(K) %*% grad_beta)
        }
        sigma_idx <- ctx$nb + 1L
        gr_full[sigma_idx] <- gr_full[sigma_idx] + 0.5 * grad_beta[intercept_idx]
        grad_disp <- gr_full[(ctx$nb + 1):(ctx$nb + 2)]
        c(grad_gamma, grad_disp)
      }
      ref2 <- nlminb(
        start = theta_init,
        objective = obj_fn,
        gradient = grad_fn,
        lower = lower_theta,
        upper = upper_theta
      )
      list(theta = ref2$par, conv_flag = 0)
    }
  )

  theta_hat <- opt_res$theta
  conv_LN <- opt_res$conv_flag
  gamma_hat <- if (free_dim == 0) numeric(0) else theta_hat[seq_len(free_dim)]
  subVar <- theta_hat[length(theta_hat) - 1]
  cellVar <- theta_hat[length(theta_hat)]
  vare <- c(subVar, cellVar)
  fit_code <- 1L

  ord <- if ((posv$mct * ctx$mfs) < 3) 3L else 1L
  betas <- compose_beta(gamma_hat)

  hl_fit <- tryCatch({
    bobyqa(
      vare,
      pql_ll,
      reml = ctx$reml,
      eps = ctx$eps,
      ord = ord,
      betas = betas,
      intcol = ctx$intcol,
      posindy = posv$posindy,
      X = ctx$pred,
      offset = ctx$offset,
      Y = posv$Y,
      n_one = posv$n_onetwo,
      ytwo = posv$ytwo,
      fid = ctx$fid,
      cumsumy = ctx$cumsumy[gene_index, ],
      posind = ctx$posind[[gene_index]],
      nb = ctx$nb,
      k = ctx$k,
      nind = ctx$nind,
      lower = ctx$min,
      upper = ctx$max
    )
  }, error = function(e) {
    bobyqa(
      vare,
      pql_ll,
      reml = ctx$reml,
      eps = ctx$eps,
      ord = 1L,
      betas = betas,
      intcol = ctx$intcol,
      posindy = posv$posindy,
      X = ctx$pred,
      offset = ctx$offset,
      Y = posv$Y,
      n_one = posv$n_onetwo,
      ytwo = posv$ytwo,
      fid = ctx$fid,
      cumsumy = ctx$cumsumy[gene_index, ],
      posind = ctx$posind[[gene_index]],
      nb = ctx$nb,
      k = ctx$k,
      nind = ctx$nind,
      lower = ctx$min,
      upper = ctx$max
    )
  })
  vare <- hl_fit$par[1:2]
  fit_code <- 2L

  beta_anchor <- beta_star
  beta_anchor[ctx$intcol] <- beta_anchor[ctx$intcol] - vare[1] / 2

  repml <- opt_pml_constrained(
    ctx$pred,
    ctx$offset,
    posv$Y,
    ctx$fid - 1,
    ctx$cumsumy[gene_index, ],
    ctx$posind[[gene_index]] - 1,
    posv$posindy,
    beta_anchor,
    K,
    gamma_hat,
    vare,
    reml = ctx$reml,
    ctx$eps,
    1
  )

  hl_state <- hl_collect_state(
    repml = repml,
    sigma = vare,
    nind = ctx$nind,
    k = ctx$k,
    posv = posv
  )
  ps <- hl_state$per_subject_stats
  if (!is.null(ps)) {
    score <- hl_score_from_stats(ps, ctx$intcol, ctx$nb, ctx$k)
    hl_state$score <- score$U
    hl_state$score_subject <- score$per_subject
    hl_state$meat <- hl_meat_from_scores(score$per_subject)
    hess_blocks <- hl_hessian_blocks_from_stats(ps)
    eff <- hl_efficient_components(hess_blocks$H_bb, hess_blocks$H_bl, hess_blocks$H_ll, score$U, ctx$nb)
    hl_state$H_blocks <- hess_blocks
    hl_state$efficient_score <- eff$score_beta_eff
    hl_state$H_schur <- eff$H_beta_beta_schur
    hl_state$H_obs <- hl_hessian_from_stats(ps)
  }

  conv_code <- check_conv(repml, conv_LN, ctx$nb, vare, ctx$min, ctx$max)

  fccov <- matrix(NA_real_, ctx$nb, ctx$nb)
  if (conv_code != -25) {
    fccov <- Rfast::spdinv(repml$var)
  }

  loglik_val <- repml$loglik
  if (is.na(loglik_val)) {
    loglik_val <- repml$loglikp
  }

  res <- list(
    beta = repml$beta,
    dispersion = vare,
    hl_state = hl_state,
    repml = repml,
    conv_code = conv_code,
    fit_code = fit_code,
    loglik = loglik_val,
    posv = posv,
    nind = ctx$nind,
    k = ctx$k,
    fccov = fccov,
    constraint_residual = as.numeric(C %*% repml$beta)
  )

  if (compute_adj) {
    res <- hl_refresh_constrained_state(
      res = res,
      ctx = ctx,
      posv = posv,
      gene_index = gene_index
    )

    if (is.null(unconstrained_fit)) {
      stop("unconstrained_fit must be provided when compute_adj = TRUE.")
    }

    ps_con <- res$hl_state$per_subject_stats

    adj_res <- tryCatch(
      hl_adjusted_tests(
        unconstrained = unconstrained_fit,
        constrained = res,
        L = C,
        b = rhs,
        ridge_factor = ridge_factor,
        cond_limit = cond_limit
      ),
      error = function(e) {
        list(
          error = conditionMessage(e),
          ok = FALSE
        )
      }
    )
    res$hl_adj <- adj_res
  } else {
    res$hl_adj <- NULL
  }

  res
}

hl_recompute_aphl <- function(fit_res) {
  if (is.null(fit_res$repml) || is.null(fit_res$dispersion) ||
      is.null(fit_res$posv) || is.null(fit_res$nind) || is.null(fit_res$k)) {
    stop("fit_res must include repml, dispersion, posv, nind, and k.")
  }
  hl_compute_aphl(
    repml = fit_res$repml,
    sigma = fit_res$dispersion,
    nind = fit_res$nind,
    k = fit_res$k,
    posv = fit_res$posv
  )$aphl
}
