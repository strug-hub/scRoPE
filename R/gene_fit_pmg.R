fit_gene_pmg_unconstrained <- function(gene_index, posv, ctx) {
  with(ctx, {
    start <- c(rep(0, nb), 1)
    start[intcol] <- log(posv$mct) - moffset
    lower <- c(rep(-100, nb), min[1])
    upper <- c(rep(100, nb), max[1])

    fit <- tryCatch(
      {
        ref <- nlminb(
          start = start,
          objective = pmg_ll,
          gradient = pmg_der,
          posindy = posv$posindy,
          X = pred,
          offset = offset,
          Y = posv$Y,
          fid = fid,
          cumsumy = cumsumy[gene_index, ],
          posind = posind[[gene_index]],
          nb = nb,
          nind = nind,
          k = k,
          lower = lower,
          upper = upper
        )
        list(
          par = ref$par,
          optimizer = "nlminb",
          raw = ref
        )
      },
      error = function(e) {
        ref <- bobyqa(
          par = start,
          fn = pmg_ll,
          posindy = posv$posindy,
          X = pred,
          offset = offset,
          Y = posv$Y,
          fid = fid,
          cumsumy = cumsumy[gene_index, ],
          posind = posind[[gene_index]],
          nb = nb,
          nind = nind,
          k = k,
          lower = lower,
          upper = upper
        )
        list(
          par = ref$par,
          optimizer = "bobyqa",
          raw = ref
        )
      }
    )

    eval_at_opt <- pmg_ll_der_hes(
      para = fit$par,
      X = pred,
      offset = offset,
      Y = posv$Y,
      fid = fid,
      cumsumy = cumsumy[gene_index, ],
      posind = posind[[gene_index]],
      posindy = posv$posindy,
      nb = nb,
      nind = nind,
      k = k
    )

    H_obs <- eval_at_opt$observed_info
    npar <- nb + 1L
    naive_cov_full <- matrix(NA_real_, nrow = npar, ncol = npar)
    active <- is.finite(diag(eval_at_opt$hessian)) &
      (fit$par != lower) &
      (fit$par != upper)
    if (any(active)) {
      naive_active <- tryCatch(
        {
          solve(-eval_at_opt$hessian[active, active, drop = FALSE])
        },
        error = function(e) {
          tryCatch(
            Rfast::spdinv(-eval_at_opt$hessian[active, active, drop = FALSE]),
            error = function(e2) NULL
          )
        }
      )
      if (!is.null(naive_active)) {
        naive_cov_full[active, active] <- naive_active
      }
    }

    robust_out <- tryCatch(
      {
        compute_sandwich_variance2(
          eval_at_opt,
          nb = nb,
          compute_full = TRUE
        )
      },
      error = function(e) NULL
    )

    convergence <- 1L
    if (identical(fit$optimizer, "nlminb")) {
      conv_raw <- fit$raw$convergence
      if (!is.null(conv_raw) && conv_raw != 0) {
        convergence <- -50L
      }
    } else {
      conv_raw <- fit$raw$convergence
      if (!is.null(conv_raw) && conv_raw != 0) {
        convergence <- -50L
      }
    }
    if (all(is.na(naive_cov_full))) {
      convergence <- -25L
    }
    if (is.null(robust_out) && convergence > 0L) {
      convergence <- -27L
    }

    naive_cov_beta <- if (all(is.na(naive_cov_full))) {
      matrix(NA_real_, nrow = nb, ncol = nb)
    } else {
      naive_cov_full[seq_len(nb), seq_len(nb), drop = FALSE]
    }
    robust_cov_beta <- if (is.null(robust_out)) {
      matrix(NA_real_, nrow = nb, ncol = nb)
    } else {
      robust_out$Var_beta_adjusted
    }
    robust_cov_full <- if (is.null(robust_out)) {
      matrix(NA_real_, nrow = nb + 1L, ncol = nb + 1L)
    } else {
      robust_out$Var_full
    }

    godambe_input <- list(
      gradient = eval_at_opt$gradient,
      hessian = eval_at_opt$hessian,
      observed_info = eval_at_opt$observed_info,
      per_subject_gradients = eval_at_opt$per_subject_gradients
    )

    list(
      theta = list(
        beta = fit$par[seq_len(nb)],
        subVar = fit$par[nb + 1L]
      ),
      loglik = eval_at_opt$loglik,
      score = eval_at_opt$gradient,
      hessian = eval_at_opt$hessian,
      observed_info = eval_at_opt$observed_info,
      per_subject_gradients = eval_at_opt$per_subject_gradients,
      naive_cov = naive_cov_beta,
      naive_cov_full = naive_cov_full,
      robust_cov = robust_cov_beta,
      robust_cov_full = robust_cov_full,
      godambe = extract_godambe_components(godambe_input, nb),
      convergence = convergence,
      optimizer = fit$optimizer,
      raw_fit = fit$raw
    )
  })
}
