hl_fit_gene_core <- function(gene_index, posv, ctx) {
  with(ctx, {
    lmct <- log(posv$mct)
    betae <- numeric(nb)
    betae[intcol] <- lmct - moffset
    ord <- if ((posv$mct * mfs) < 3) 3L else 1L

    re_t <- tryCatch({
      ref <- switch(opt,
        "lbfgs" = lbfgs(
          c(betae, 1, cell_init),
          ptmg_ll_der,
          posindy = posv$posindy,
          X = pred,
          offset = offset,
          Y = posv$Y,
          n_one = posv$n_onetwo,
          ytwo = posv$ytwo,
          fid = fid,
          cumsumy = cumsumy[gene_index, ],
          posind = posind[[gene_index]],
          nb = nb,
          k = k,
          nind = nind,
          lower = c(rep(-100, nb), min[1], min[2]),
          upper = c(rep(100, nb), max[1], max[2]),
          control = list(ftol_abs = eps)
        ),
        "trust" = {
          out1 <- trust(
            objfun = ptmg_ll_der_hes3,
            parinit = c(betae, 0, 0),
            rinit = 2,
            rmax = 100,
            iterlim = 100,
            posindy = posv$posindy,
            X = pred,
            offset = offset,
            Y = posv$Y,
            n_one = posv$n_onetwo,
            ytwo = posv$ytwo,
            fid = fid,
            cumsumy = cumsumy[gene_index, ],
            posind = posind[[gene_index]],
            nb = nb,
            k = k,
            nind = nind
          )
          outp <- out1$argument
          outp[nb + 1] <- median(c(exp(outp[nb + 1]), min[1], max[1]))
          outp[nb + 2] <- median(c(exp(outp[nb + 2]), min[2], max[2]))
          list(par = outp)
        },
        stop("opt must be lbfgs or trust.")
      )
      c(ref$par, 1)
    }, error = function(e) {
      ref2 <- nlminb(
        start = c(betae, 1, cell_init),
        objective = ptmg_ll,
        gradient = ptmg_der,
        posindy = posv$posindy,
        X = pred,
        offset = offset,
        Y = posv$Y,
        n_one = posv$n_onetwo,
        ytwo = posv$ytwo,
        fam = id,
        fid = fid,
        cumsumy = cumsumy[gene_index, ],
        posind = posind[[gene_index]],
        nb = nb,
        k = k,
        nind = nind,
        lower = c(rep(-100, nb), min[1], min[2]),
        upper = c(rep(100, nb), max[1], max[2])
      )
      c(ref2$par, 0)
    })

    conv_LN <- re_t[length(re_t)]
    beta_ln <- re_t[seq_len(nb)]
    vare <- re_t[(nb + 1):(nb + 2)]

    if (ncell > 0) {
      invisible(get_cv(offset, pred, beta_ln, cell_ind, ncell, nind))
    } else {
      invisible(cv2)
    }
    fit_code <- 1L

    hl_fit <- tryCatch({
      bobyqa(
        vare,
        pql_ll,
        reml = reml,
        eps = eps,
        ord = ord,
        betas = betae,
        intcol = intcol,
        posindy = posv$posindy,
        X = pred,
        offset = offset,
        Y = posv$Y,
        n_one = posv$n_onetwo,
        ytwo = posv$ytwo,
        fid = fid,
        cumsumy = cumsumy[gene_index, ],
        posind = posind[[gene_index]],
        nb = nb,
        k = k,
        nind = nind,
        lower = min,
        upper = max
      )
    }, error = function(e) {
      bobyqa(
        vare,
        pql_ll,
        reml = reml,
        eps = eps,
        ord = 1L,
        betas = betae,
        intcol = intcol,
        posindy = posv$posindy,
        X = pred,
        offset = offset,
        Y = posv$Y,
        n_one = posv$n_onetwo,
        ytwo = posv$ytwo,
        fid = fid,
        cumsumy = cumsumy[gene_index, ],
        posind = posind[[gene_index]],
        nb = nb,
        k = k,
        nind = nind,
        lower = min,
        upper = max
      )
    })
    vare <- hl_fit$par[1:2]
    fit_code <- 2L

    betae[intcol] <- betae[intcol] - vare[1] / 2

    repml <- opt_pml_hl(
      pred,
      offset,
      posv$Y,
      fid - 1,
      cumsumy[gene_index, ],
      posind[[gene_index]] - 1,
      posv$posindy,
      nb,
      nind,
      k,
      betae,
      vare,
      reml,
      eps,
      1
    )

    hl_state <- hl_collect_state(
      repml = repml,
      sigma = vare,
      nind = nind,
      k = k,
      posv = posv
    )

    conv_code <- check_conv(repml, conv_LN, nb, vare, min, max)

    fccov <- matrix(NA_real_, nb, nb)
    if (conv_code != -25) {
      fccov <- Rfast::spdinv(repml$var)
    }

    loglik_val <- repml$loglik
    if (is.na(loglik_val)) {
      loglik_val <- repml$loglikp
    }

    list(
      repml = repml,
      vare = vare,
      hl_state = hl_state,
      conv_code = conv_code,
      fit_code = fit_code,
      loglik = loglik_val,
      posv = posv,
      nind = nind,
      k = k,
      fccov = fccov
    )
  })
}

fit_gene_hl <- function(gene_index, posv, ctx) {
  core <- hl_fit_gene_core(gene_index, posv, ctx)
  with(ctx, {
    hl_state <- core$hl_state
    hl_state$theta <- list(beta = core$repml$beta, subVar = core$vare[1], cellVar = core$vare[2])

    ps <- hl_state$per_subject_stats
    if (!is.null(ps)) {
      score <- hl_score_from_stats(ps, intcol, nb, k)
      hl_state$score <- score$U
      hl_state$score_subject <- score$per_subject
      hl_state$meat <- hl_meat_from_scores(score$per_subject)
      hess_blocks <- hl_hessian_blocks_from_stats(ps)
      hl_state$H_blocks <- hess_blocks
      eff <- hl_efficient_components(
        hess_blocks$H_bb,
        hess_blocks$H_bl,
        hess_blocks$H_ll,
        score$U,
        nb
      )
      hl_state$efficient_score <- eff$score_beta_eff
      hl_state$H_schur <- eff$H_beta_beta_schur
      hl_state$H_obs <- hl_hessian_from_stats(ps)
    }

    out_at_op <- ptmg_ll_der_hes4(
      c(core$repml$beta, log(core$vare[1]), log(core$vare[2])),
      pred,
      offset,
      posv$Y,
      posv$n_onetwo,
      posv$ytwo,
      fid,
      cumsumy[gene_index, ],
      posind[[gene_index]],
      posv$posindy,
      nb,
      nind,
      k
    )
    godambe_summary <- extract_godambe_components(out_at_op, nb)

    restemp <- c(
      core$repml$beta,
      core$vare[1],
      1 / core$vare[2],
      diag(core$fccov),
      core$conv_code,
      core$fit_code,
      core$loglik
    )
    if (output_re) {
      restemp <- c(restemp, core$repml$logw)
    }

    unconstrained_fit <- list(
      beta = core$repml$beta,
      dispersion = core$vare,
      hl_state = hl_state,
      repml = core$repml,
      conv_code = core$conv_code,
      fit_code = core$fit_code,
      loglik = core$loglik,
      posv = posv,
      nind = core$nind,
      k = core$k,
      fccov = core$fccov
    )

    list(
      vec = restemp,
      godambe = godambe_summary,
      loglik = core$loglik,
      loglik_ln = NA_real_,
      theta = list(beta = core$repml$beta, subVar = core$vare[1], cellVar = core$vare[2]),
      score = out_at_op$gradient,
      hl = hl_state,
      hl_unconstrained = unconstrained_fit,
      hl_adj = NULL
    )
  })
}
