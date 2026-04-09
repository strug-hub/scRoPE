fit_gene_unconstrained <- function(gene_index, posv, ctx, extra_starts = NULL) {
  # Nebula reference flow (LN branch):
  # 1. Primary optimisation via lbfgs/trust on (beta, log sigma^2, log phi).
  # 2. Compute gni = (ncell / k) * cellVar and cv2p via get_cv().
  # 3. Rescue branches:
  #    a) If gni < cutoff_cell -> bobyqa over (sigma^2, phi) (fit code 2).
  #    b) Else if kappa_obs = gni / (1 + cv2p) < kappa -> nlminb on sigma^2
  #       holding phi fixed (fit code 3).
  # 4. After any update, shift intercept by -sigma^2/2, call opt_pml(), and
  #    record fit code and diagnostics. We mirror this structure below.
  grad_tol <- 5e-4

  with(ctx, {
    lmct <- log(posv$mct)
    betae <- rep(0, nb)
    betae[intcol] <- lmct - moffset

    re_t <- tryCatch(
      {
        ref <- switch(opt,
          "lbfgs" = {
            lbfgs(
              c(betae, 1, cell_init),
              ptmg_ll_der,
              posindy = posv$posindy,
              X       = pred,
              offset  = offset,
              Y       = posv$Y,
              n_one   = posv$n_onetwo,
              ytwo    = posv$ytwo,
              fid     = fid,
              cumsumy = cumsumy[gene_index, ],
              posind  = posind[[gene_index]],
              nb      = nb,
              k       = k,
              nind    = nind,
              lower   = c(rep(-100, nb), min[1], min[2]),
              upper   = c(rep(100, nb), max[1], max[2]),
              control = list(ftol_abs = eps)
            )
          },
          "trust" = {
            out1 <- trust(
              objfun  = ptmg_ll_der_hes3,
              parinit = c(betae, 0, 0),
              rinit   = 2,
              rmax    = 100,
              iterlim = 100,
              posindy = posv$posindy,
              X       = pred,
              offset  = offset,
              Y       = posv$Y,
              n_one   = posv$n_onetwo,
              ytwo    = posv$ytwo,
              fid     = fid,
              cumsumy = cumsumy[gene_index, ],
              posind  = posind[[gene_index]],
              nb      = nb,
              k       = k,
              nind    = nind
            )
            outp <- out1$argument
            outp[nb + 1] <- median(c(exp(outp[nb + 1]), min[1], max[1]))
            outp[nb + 2] <- median(c(exp(outp[nb + 2]), min[2], max[2]))
            list(par = outp)
          },
          stop("opt must be lbfgs or trust.")
        )
        c(ref$par, 1)
      },
      error = function(e) {
        ref2 <- nlminb(
          start     = c(betae, 1, cell_init),
          objective = ptmg_ll,
          gradient  = ptmg_der,
          posindy   = posv$posindy,
          X         = pred,
          offset    = offset,
          Y         = posv$Y,
          n_one     = posv$n_onetwo,
          ytwo      = posv$ytwo,
          fam       = id,
          fid       = fid,
          cumsumy   = cumsumy[gene_index, ],
          posind    = posind[[gene_index]],
          nb        = nb,
          k         = k,
          nind      = nind,
          lower     = c(rep(-100, nb), min[1], min[2]),
          upper     = c(rep(100, nb), max[1], max[2])
        )
        c(ref2$par, 0)
      }
    )

    conv_LN <- re_t[length(re_t)]
    LN_betas <- re_t[1:nb]
    subVar <- re_t[nb + 1]
    cellVar <- re_t[nb + 2]
    fit <- 1
    gni <- (nind / k) * cellVar

    beta_ln <- LN_betas
    theta_ln_vec <- c(beta_ln, subVar, cellVar)

    evaluate_candidate <- function(theta_vec) {
      obj_val <- ptmg_ll(
        theta_vec,
        pred,
        offset,
        posv$Y,
        posv$n_onetwo,
        posv$ytwo,
        id,
        fid,
        cumsumy[gene_index, ],
        posind[[gene_index]],
        posv$posindy,
        nb,
        nind,
        k
      )
      grad_val <- ptmg_der(
        theta_vec,
        X = pred,
        offset = offset,
        Y = posv$Y,
        n_one = posv$n_onetwo,
        ytwo = posv$ytwo,
        fam = id,
        fid = fid,
        cumsumy = cumsumy[gene_index, ],
        posind = posind[[gene_index]],
        posindy = posv$posindy,
        nb = nb,
        nind = nind,
        k = k
      )
      list(
        theta = theta_vec,
        obj = obj_val,
        grad = grad_val,
        grad_norm = sqrt(sum(grad_val^2)),
        refined = FALSE
      )
    }

    candidates <- list(evaluate_candidate(theta_ln_vec))
    grad_vec_initial <- candidates[[1]]$grad
    grad_norm_initial <- candidates[[1]]$grad_norm

    lower_theta <- c(rep(-100, nb), min)
    upper_theta <- c(rep(100, nb), max)
    extra_attempts <- 0L
    extra_success <- 0L

    add_candidate_from_start <- function(start_theta, label = NULL) {
      start_theta <- as.numeric(start_theta)
      if (length(start_theta) != (nb + 2) || any(!is.finite(start_theta))) {
        return(FALSE)
      }
      best_obj <- min(vapply(candidates, function(c) c$obj, numeric(1)))
      tryCatch({
        opt <- nlminb(
          start     = start_theta,
          objective = ptmg_ll,
          gradient  = ptmg_der,
          lower     = lower_theta,
          upper     = upper_theta,
          X         = pred,
          offset    = offset,
          Y         = posv$Y,
          n_one     = posv$n_onetwo,
          ytwo      = posv$ytwo,
          fam       = id,
          fid       = fid,
          cumsumy   = cumsumy[gene_index, ],
          posind    = posind[[gene_index]],
          posindy   = posv$posindy,
          nb        = nb,
          nind      = nind,
          k         = k
        )
        if (is.finite(opt$objective) && (opt$convergence == 0 || opt$objective < best_obj - 1e-8)) {
          cand <- evaluate_candidate(opt$par)
          cand$refined <- TRUE
          candidates <<- c(candidates, list(cand))
          if (!is.null(label) && identical(label, "extra")) {
            extra_success <<- extra_success + 1L
          }
          return(TRUE)
        }
        FALSE
      }, error = function(e) FALSE)
    }

    if (grad_norm_initial > grad_tol) {
      add_candidate_from_start(theta_ln_vec)

      free_idx <- setdiff(seq_len(nb), intcol)
      if (length(free_idx) > 0) {
        theta_alt <- theta_ln_vec
        theta_alt[free_idx] <- 0
        add_candidate_from_start(theta_alt)
      }

      theta_zero <- c(rep(0, nb), 1, 1)
      add_candidate_from_start(theta_zero)
    }

    if (!is.null(extra_starts)) {
      extra_list <- lapply(extra_starts, as.numeric)
      extra_list <- extra_list[!vapply(extra_list, function(x) {
        is.null(x) || length(x) != (nb + 2) || any(!is.finite(x))
      }, logical(1))]
      if (length(extra_list) > 0) {
        seen <- character(0)
        for (start_vec in extra_list) {
          key <- paste0(round(start_vec, 12), collapse = "|")
          if (key %in% seen) {
            next
          }
          seen <- c(seen, key)
          extra_attempts <- extra_attempts + 1L
          add_candidate_from_start(start_vec, label = "extra")
        }
      }
    }

    best_idx <- which.min(vapply(candidates, function(c) c$obj, numeric(1)))
    best <- candidates[[best_idx]]
    theta_ln_vec <- best$theta
    beta_ln <- theta_ln_vec[seq_len(nb)]
    subVar <- theta_ln_vec[nb + 1]
    cellVar <- theta_ln_vec[nb + 2]
    theta_pre_rescue <- theta_ln_vec
    beta_pre_rescue <- beta_ln
    loglik_ln_max <- -best$obj
    conv_LN <- if (best$refined) 2L else conv_LN
    gni <- (nind / k) * cellVar
    debug_env <- ctx$debug_env
    if (is.null(debug_env)) {
      env_name <- getOption("scRoPE.debug_ln_env_name")
      if (is.character(env_name) && length(env_name) == 1 && nzchar(env_name) && exists(env_name, envir = .GlobalEnv)) {
        debug_env <- get(env_name, envir = .GlobalEnv)
      }
    }
    gene_labels <- ctx$gene_labels
    gene_label <- NA_character_
    if (!is.null(gene_labels) && length(gene_labels) >= gene_index) {
      gene_label <- gene_labels[gene_index]
    }
    initial_subVar <- subVar
    initial_cellVar <- cellVar
    rescue_label <- "none"
    variances_updated <- FALSE
    cv2p <- NA_real_
    kappa_obs <- NA_real_
    refined_flag <- best$refined
    grad_norm_final <- best$grad_norm

    ord <- if ((posv$mct * (nind / k)) < 3) 3 else 1
    base_betas <- rep(0, nb)
    base_betas[intcol] <- lmct - moffset
    loglik_ln_rescued <- NA_real_
    grad_norm_rescued <- NA_real_
    compute_loglik_ln <- function(beta_vec, sub_var, cell_var) {
      -ptmg_ll(
        c(beta_vec, sub_var, cell_var),
        pred,
        offset,
        posv$Y,
        posv$n_onetwo,
        posv$ytwo,
        id,
        fid,
        cumsumy[gene_index, ],
        posind[[gene_index]],
        posv$posindy,
        nb,
        nind,
        k
      )
    }

    if (gni < cutoff_cell) {
      re_t2 <- tryCatch(
        {
          bobyqa(
            c(subVar, cellVar),
            pql_ll,
            reml = 0,
            eps = eps,
            ord = ord,
            betas = base_betas,
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
        },
        error = function(e3) {
          alt_ord <- if (ord == 3) 1 else 3
          bobyqa(
            c(subVar, cellVar),
            pql_ll,
            reml = 0,
            eps = eps,
            ord = alt_ord,
            betas = base_betas,
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
        }
      )
      if (is.list(re_t2) && length(re_t2$par) >= 2) {
        subVar <- re_t2$par[1]
        cellVar <- re_t2$par[2]
        fit <- 2L
        rescue_label <- "bobyqa"
        variances_updated <- TRUE
      }
    } else {
      cv2p <- if (ncell > 0) {
        get_cv(offset, pred, beta_ln, cell_ind, ncell, nind)
      } else {
        cv2
      }
      kappa_obs <- gni / (1 + cv2p)
      if ((kappa_obs < 20) || ((kappa_obs < kappa) && (subVar < (8 / kappa_obs)))) {
        target_sigma <- min(max(8 / max(kappa_obs, 1e-08), min[1] * 1.5), max[1])
        start_sigma <- max(subVar, min[1] * (1 + 1e-03))
        if (start_sigma <= min[1] * (1 + 1e-03)) {
          start_sigma <- target_sigma
        }
        rescue_label <- "nlminb"
        varsig <- tryCatch(
          {
            nlminb(
              start = start_sigma,
              objective = pql_gamma_ll,
              gamma = cellVar,
              betas = base_betas,
              intcol = intcol,
              reml = 0,
              eps = eps,
              ord = ord,
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
              lower = min[1],
              upper = max[1]
            )
          },
          error = function(e4) {
            alt_ord <- if (ord == 3) 1 else 3
            nlminb(
              start = start_sigma,
              objective = pql_gamma_ll,
              gamma = cellVar,
              betas = base_betas,
              intcol = intcol,
              reml = 0,
              eps = eps,
              ord = alt_ord,
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
              lower = min[1],
              upper = max[1]
            )
          }
        )
        if (is.list(varsig) && !is.null(varsig$par)) {
          subVar <- max(varsig$par[1], min[1])
          fit <- 3L
          variances_updated <- TRUE
          if (subVar <= min[1] * (1 + 1e-04)) {
            start_pair <- c(max(subVar, min[1] * (1 + 1e-03), target_sigma), cellVar)
            re_t2 <- tryCatch(
              {
                bobyqa(
                  start_pair,
                  pql_ll,
                  reml = 0,
                  eps = eps,
                  ord = ord,
                  betas = base_betas,
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
              },
              error = function(e5) {
                alt_ord <- if (ord == 3) 1 else 3
                bobyqa(
                  start_pair,
                  pql_ll,
                  reml = 0,
                  eps = eps,
                  ord = alt_ord,
                  betas = base_betas,
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
              }
            )
            if (is.list(re_t2) && length(re_t2$par) >= 2) {
              subVar <- re_t2$par[1]
              cellVar <- re_t2$par[2]
              fit <- 2L
              variances_updated <- TRUE
            }
          }
        }
      }
    }

    entry <- list(
      method = "scrope",
      gene_index = gene_index,
      gene = gene_label,
      subVar_lbfgs = initial_subVar,
      cellVar_lbfgs = initial_cellVar,
      subVar_final = subVar,
      cellVar_final = cellVar,
      grad_norm_start = grad_norm_initial,
      grad_norm_final = grad_norm_final,
      extra_attempts = extra_attempts,
      extra_success = extra_success,
      rescue = rescue_label,
      kappa_obs = kappa_obs,
      cv2p = cv2p,
      gni = gni,
      fit_code = fit,
      conv_LN = conv_LN,
      loglik_ln_pre = loglik_ln_max,
      loglik_ln_rescued = NA_real_,
      grad_norm_rescued = NA_real_
    )
    godambe_summary <- NULL

    if (!allow_per_gene_switch) {
      final_betas <- if (use_betas == "betae") {
        temp <- rep(0, nb)
        temp[intcol] <- lmct - moffset
        temp[intcol] <- temp[intcol] - subVar / 2
        temp
      } else {
        temp <- beta_ln
        temp[intcol] <- temp[intcol] - subVar / 2
        temp
      }
      repml <- opt_pml(
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
        final_betas,
        c(subVar, cellVar),
        reml = 0,
        eps,
        1
      )
      conv_code <- check_conv(repml, conv_LN, nb, c(subVar, cellVar), min, max)

      fccov <- matrix(NA, nb, nb)
      if (conv_code != -25) {
        fccov <- Rfast::spdinv(repml$var)
      }

      out_at_op <- ptmg_ll_der_hes4(
        c(repml$beta, log(subVar), log(cellVar)),
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
      robust_try <- tryCatch(
        {
          compute_sandwich_variance2(out_at_op, nb, compute_full = FALSE)
        },
        error = function(e2) NULL
      )
      robust_var <- matrix(NA, nb, nb)
      if (is.null(robust_try)) {
        if (conv_code >= 0) {
          conv_code <- -27
        }
      } else {
        robust_var <- robust_try$Var_beta_adjusted
      }

      beta_rescued <- repml$beta
      theta_rescued_vec <- c(beta_rescued, subVar, cellVar)
      loglik_ln_rescued <- compute_loglik_ln(beta_rescued, subVar, cellVar)
      grad_norm_rescued <- NA_real_
      if (variances_updated) {
        post_eval <- evaluate_candidate(theta_rescued_vec)
        if (is.list(post_eval) && is.finite(post_eval$grad_norm)) {
          grad_norm_rescued <- post_eval$grad_norm
        }
      }
      loglik_ln <- loglik_ln_max

      loglik_val <- repml$loglik
      if (is.na(loglik_val)) {
        loglik_val <- repml$loglikp
      }

      entry$loglik_ln_rescued <- loglik_ln_rescued
      entry$grad_norm_rescued <- grad_norm_rescued
      debug_info <- entry
      if (!is.null(debug_env)) {
        entries <- get0("scrope", envir = debug_env, inherits = FALSE)
        if (is.null(entries)) {
          entries <- list()
        }
        entries[[length(entries) + 1]] <- entry
        assign("scrope", entries, envir = debug_env)
      }

      restemp <- c(
        beta_rescued,
        subVar,
        1 / cellVar,
        diag(fccov),
        diag(robust_var),
        conv_code,
        fit,
        loglik_val
      )
      if (output_re) {
        restemp <- c(restemp, repml$logw)
      }
      list(
        vec = restemp,
        godambe = godambe_summary,
        loglik = loglik_val,
        loglik_ln = loglik_ln,
        theta = list(beta = beta_rescued, subVar = subVar, cellVar = cellVar),
        theta_ln = list(
          beta = beta_pre_rescue,
          subVar = theta_pre_rescue[nb + 1],
          cellVar = theta_pre_rescue[nb + 2]
        ),
        theta_ln_rescued = list(beta = beta_rescued, subVar = subVar, cellVar = cellVar),
        theta_unrescued = list(
          beta = beta_pre_rescue,
          subVar = theta_pre_rescue[nb + 1],
          cellVar = theta_pre_rescue[nb + 2]
        ),
        score = out_at_op$gradient,
        ln_refinement = list(
          refined = refined_flag,
          grad_norm_start = grad_norm_initial,
          grad_norm_final = grad_norm_final,
          extra_attempts = extra_attempts,
          extra_success = extra_success
        ),
        debug = debug_info
      )
    } else {
      next_betas <- if (use_betas == "betae") {
        temp <- rep(0, nb)
        temp[intcol] <- lmct - moffset
        temp
      } else {
        beta_ln
      }
      vare <- c(subVar, cellVar)

      if (gni < cutoff_cell) {
        ord <- if ((posv$mct * (nind / k)) < 3) 3 else 1
        re_t2 <- tryCatch(
          {
            bobyqa(
              vare,
              pql_ll,
              reml = 0,
              eps = eps,
              ord = ord,
              betas = next_betas,
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
          },
          error = function(e3) {
            alt <- if (ord == 3) 1 else 3
            bobyqa(
              vare,
              pql_ll,
              reml = 0,
              eps = eps,
              ord = alt,
              betas = next_betas,
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
          }
        )
        vare <- re_t2$par[1:2]
        fit <- 2
        variances_updated <- TRUE
      } else {
        if (ncell > 0) {
          cv2p <- get_cv(offset, pred, LN_betas, cell_ind, ncell, nind)
        } else {
          cv2p <- cv2
        }
        kappa_obs <- gni / (1 + cv2p)
        if ((kappa_obs < 20) || ((kappa_obs < kappa) && (vare[1] < (8 / kappa_obs)))) {
          ord <- if ((posv$mct * (nind / k)) < 3) 3 else 1
          varsig <- tryCatch(
            {
              nlminb(
                start     = vare[1],
                objective = pql_gamma_ll,
                gamma     = vare[2],
                betas     = next_betas,
                intcol    = intcol,
                reml      = 0,
                eps       = eps,
                ord       = ord,
                posindy   = posv$posindy,
                X         = pred,
                offset    = offset,
                Y         = posv$Y,
                n_one     = posv$n_onetwo,
                ytwo      = posv$ytwo,
                fid       = fid,
                cumsumy   = cumsumy[gene_index, ],
                posind    = posind[[gene_index]],
                nb        = nb,
                k         = k,
                nind      = nind,
                lower     = min[1],
                upper     = max[1]
              )
            },
            error = function(e4) {
              alt <- if (ord == 3) 1 else 3
              nlminb(
                start     = vare[1],
                objective = pql_gamma_ll,
                gamma     = vare[2],
                betas     = next_betas,
                intcol    = intcol,
                reml      = 0,
                eps       = eps,
                ord       = alt,
                posindy   = posv$posindy,
                X         = pred,
                offset    = offset,
                Y         = posv$Y,
                n_one     = posv$n_onetwo,
                ytwo      = posv$ytwo,
                fid       = fid,
                cumsumy   = cumsumy[gene_index, ],
                posind    = posind[[gene_index]],
                nb        = nb,
                k         = k,
                nind      = nind,
                lower     = min[1],
                upper     = max[1]
              )
            }
          )
          vare[1] <- varsig$par[1]
          fit <- 3
          variances_updated <- TRUE
        }
      }

      next_betas[intcol] <- next_betas[intcol] - vare[1] / 2
      repml <- opt_pml(
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
        next_betas,
        vare,
        reml = 0,
        eps,
        1
      )
      conv_code <- check_conv(repml, conv_LN, nb, vare, min, max)

      fccov <- matrix(NA, nb, nb)
      if (conv_code != -25) {
        fccov <- Rfast::spdinv(repml$var)
      }

      out_at_op <- ptmg_ll_der_hes4(
        c(repml$beta, log(vare[1]), log(vare[2])),
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
      robust_try <- tryCatch(
        {
          compute_sandwich_variance2(out_at_op, nb, compute_full = FALSE)
        },
        error = function(e5) NULL
      )
      robust_var <- matrix(NA, nb, nb)
      if (is.null(robust_try)) {
        if (conv_code >= 0) {
          conv_code <- -27
        }
      } else {
        robust_var <- robust_try$Var_beta_adjusted
      }

      beta_rescued <- repml$beta
      theta_rescued_vec <- c(beta_rescued, vare[1], vare[2])
      loglik_ln_rescued <- compute_loglik_ln(beta_rescued, vare[1], vare[2])
      if (variances_updated) {
        post_eval <- evaluate_candidate(theta_rescued_vec)
        if (is.list(post_eval) && is.finite(post_eval$grad_norm)) {
          grad_norm_rescued <- post_eval$grad_norm
        }
      }
      loglik_ln <- loglik_ln_max

      loglik_val <- repml$loglik
      if (is.na(loglik_val)) {
        loglik_val <- repml$loglikp
      }

      entry$loglik_ln_rescued <- loglik_ln_rescued
      entry$grad_norm_rescued <- grad_norm_rescued
      debug_info <- entry
      if (!is.null(debug_env)) {
        entries <- get0("scrope", envir = debug_env, inherits = FALSE)
        if (is.null(entries)) {
          entries <- list()
        }
        entries[[length(entries) + 1]] <- entry
        assign("scrope", entries, envir = debug_env)
      }

      restemp <- c(
        beta_rescued,
        vare[1],
        1 / vare[2],
        diag(fccov),
        diag(robust_var),
        conv_code,
        fit,
        loglik_val
      )
      if (output_re) {
        restemp <- c(restemp, repml$logw)
      }
      list(
        vec = restemp,
        godambe = godambe_summary,
        loglik = loglik_val,
        loglik_ln = loglik_ln,
        theta = list(beta = beta_rescued, subVar = vare[1], cellVar = vare[2]),
        theta_ln = list(
          beta = beta_pre_rescue,
          subVar = theta_pre_rescue[nb + 1],
          cellVar = theta_pre_rescue[nb + 2]
        ),
        theta_ln_rescued = list(beta = beta_rescued, subVar = vare[1], cellVar = vare[2]),
        theta_unrescued = list(
          beta = beta_pre_rescue,
          subVar = theta_pre_rescue[nb + 1],
          cellVar = theta_pre_rescue[nb + 2]
        ),
        score = out_at_op$gradient,
        ln_refinement = list(
          refined = refined_flag,
          grad_norm_start = grad_norm_initial,
          grad_norm_final = grad_norm_final,
          extra_attempts = extra_attempts,
          extra_success = extra_success
        ),
        debug = debug_info
      )
    }
  })
}

fit_gene_constrained <- function(gene_index, posv, ctx, L, b) {
  if (ctx$allow_per_gene_switch) {
    stop("Constrained fits do not yet support per-gene LN/HL switching.")
  }
  if (ctx$opt != "lbfgs") {
    stop("Constrained fits currently require opt = 'lbfgs'.")
  }
  if (length(ctx$intcol) != 1L) {
    stop("Design matrix must have exactly one intercept column.")
  }
  intercept_idx <- ctx$intcol

  lmct <- log(posv$mct)

  constr <- constraint_reparameterization(L, b)
  beta_star <- constr$beta_star
  K <- constr$null_basis
  if (is.null(K)) {
    K <- matrix(0, nrow = ctx$nb, ncol = 0)
  }
  free_dim <- ncol(K)

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
    beta_full <- compose_beta(gamma)
    c(beta_full, disp)
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
    if (length(ctx$intcol) != 1L) {
      stop("Design matrix must have exactly one intercept column for constrained fit.")
    }
    sigma_idx <- ctx$nb + 1L
    grad_beta <- grad_full[seq_len(ctx$nb)]
    grad_gamma <- if (free_dim == 0) numeric(0) else as.vector(t(K) %*% grad_beta)
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
        grad_gamma <- if (free_dim == 0) numeric(0) else as.vector(t(K) %*% grad_beta)
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
  beta_hat <- compose_beta(gamma_hat)
  fit <- 1

  beta_anchor <- beta_star
  if (length(ctx$intcol) > 0) {
    beta_anchor[ctx$intcol] <- beta_anchor[ctx$intcol] - subVar / 2
  }

  theta_pre <- c(beta_hat, subVar, cellVar)
  loglik_ln_pre <- -ptmg_ll(
    c(beta_hat, subVar, cellVar),
    ctx$pred,
    ctx$offset,
    posv$Y,
    posv$n_onetwo,
    posv$ytwo,
    ctx$id,
    ctx$fid,
    ctx$cumsumy[gene_index, ],
    ctx$posind[[gene_index]],
    posv$posindy,
    ctx$nb,
    ctx$nind,
    ctx$k
  )

  loglik_ln_rescued <- loglik_ln_pre
  if (free_dim > 0) {
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
      c(subVar, cellVar),
      reml = 0,
      ctx$eps,
      1
    )
    conv_code <- check_conv(repml, conv_LN, ctx$nb, c(subVar, cellVar), ctx$min, ctx$max)
    final_betas <- repml$beta
    loglik_val <- repml$loglik
    if (is.na(loglik_val)) {
      loglik_val <- repml$loglikp
    }
    loglik_ln_rescued <- -ptmg_ll(
      c(final_betas, subVar, cellVar),
      ctx$pred,
      ctx$offset,
      posv$Y,
      posv$n_onetwo,
      posv$ytwo,
      ctx$id,
      ctx$fid,
      ctx$cumsumy[gene_index, ],
      ctx$posind[[gene_index]],
      posv$posindy,
      ctx$nb,
      ctx$nind,
      ctx$k
    )
    repml_logw <- repml$logw
    fccov <- matrix(NA_real_, ctx$nb, ctx$nb)
    if (conv_code != -25) {
      fccov <- tryCatch({
        Rfast::spdinv(repml$var)
      }, error = function(e) {
        solve(repml$var)
      })
    }
  } else {
    final_betas <- beta_anchor
    conv_code <- conv_LN
    loglik_val <- loglik_ln_pre
    repml_logw <- rep(0, ctx$k)
    fccov <- matrix(NA_real_, ctx$nb, ctx$nb)
    loglik_ln_rescued <- loglik_ln_pre
  }

  out_at_op <- ptmg_ll_der_hes4(
    c(final_betas, log(subVar), log(cellVar)),
    ctx$pred,
    ctx$offset,
    posv$Y,
    posv$n_onetwo,
    posv$ytwo,
    ctx$fid,
    ctx$cumsumy[gene_index, ],
    ctx$posind[[gene_index]],
    posv$posindy,
    ctx$nb,
    ctx$nind,
    ctx$k
  )
  chain_rule <- apply_chain_rule_transform(
    U = out_at_op$gradient,
    H = out_at_op$hessian,
    S = NULL,
    intercept_idx = intercept_idx,
    p = ctx$nb,
    U_per_subject = out_at_op$per_subject_gradients,
    sigma2 = subVar
  )
  out_at_op$gradient <- chain_rule$U
  out_at_op$hessian <- chain_rule$H
  if (!is.null(chain_rule$U_per_subject)) {
    out_at_op$per_subject_gradients <- chain_rule$U_per_subject
  }
  if (!is.null(chain_rule$S)) {
    out_at_op$S_override <- chain_rule$S
  }
  godambe_summary <- extract_godambe_components(out_at_op, ctx$nb)

  robust_try <- tryCatch(
    {
      compute_sandwich_variance2(out_at_op, ctx$nb, compute_full = FALSE)
    },
    error = function(e2) NULL
  )
  robust_var <- matrix(NA, ctx$nb, ctx$nb)
  if (is.null(robust_try)) {
    if (conv_code >= 0) {
      conv_code <- -27
    }
  } else {
    robust_var <- robust_try$Var_beta_adjusted
  }

  restemp <- c(
    final_betas,
    subVar,
    1 / cellVar,
    diag(fccov),
    diag(robust_var),
    conv_code,
    fit,
    loglik_val
  )
  if (ctx$output_re) {
    restemp <- c(restemp, repml_logw)
  }

  list(
    vec = restemp,
    godambe = godambe_summary,
    beta_constrained = final_betas,
    loglik = loglik_val,
    loglik_ln = loglik_ln_pre,
    loglik_ln_pre = loglik_ln_pre,
    loglik_ln_rescued = loglik_ln_rescued,
    theta = list(
      beta = final_betas,
      subVar = subVar,
      cellVar = cellVar
    ),
    theta_ln = list(
      beta = beta_hat,
      subVar = theta_pre[length(theta_pre) - 1],
      cellVar = theta_pre[length(theta_pre)]
    ),
    theta_ln_rescued = list(
      beta = final_betas,
      subVar = subVar,
      cellVar = cellVar
    ),
    theta_unrescued = list(
      beta = beta_hat,
      subVar = theta_pre[length(theta_pre) - 1],
      cellVar = theta_pre[length(theta_pre)]
    ),
    score = out_at_op$gradient,
    lrt_gate = list(
      converged = conv_code >= 0,
      interior = isTRUE(
        subVar > ctx$min[1] && subVar < ctx$max[1] &&
          cellVar > ctx$min[2] && cellVar < ctx$max[2]
      ),
      nuisance_interior = isTRUE(
        subVar > ctx$min[1] && subVar < ctx$max[1] &&
          cellVar > ctx$min[2] && cellVar < ctx$max[2]
      )
    ),
    debug = list(
      method = "scrope_constrained",
      gene_index = gene_index,
      fit_code = fit,
      conv_LN = conv_LN,
      conv_code = conv_code,
      free_dim = free_dim,
      loglik_ln_pre = loglik_ln_pre,
      loglik_ln_rescued = loglik_ln_rescued
    )
  )
}
