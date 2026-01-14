.hl_symmetrize_matrix <- function(mat) {
  0.5 * (mat + t(mat))
}

hl_state_godambe <- function(hl_state, nb) {
  if (is.null(hl_state)) {
    stop("hl_state must not be NULL.")
  }
  H_obs <- hl_state$H_obs
  if (is.null(H_obs)) {
    stop("hl_state$H_obs is missing; recompute Stage 4 caches before Stage 5.")
  }

  H_bb <- H_obs$H_bb
  H_bphi <- H_obs$H_bphi
  H_btau <- H_obs$H_btau
  H_phiphi <- H_obs$H_phiphi
  H_phitau <- H_obs$H_phitau
  H_tautau <- H_obs$H_tautau

  if (any(vapply(
    list(H_bb, H_bphi, H_btau, H_phiphi, H_phitau, H_tautau),
    is.null,
    logical(1)
  ))) {
    stop("Observed Hessian blocks missing from hl_state; verify Stage 4 implementation.")
  }

  if (nrow(H_bb) != nb) {
    stop("Dimension mismatch between hl_state$H_obs$H_bb and supplied nb.")
  }

  H_beta <- .hl_symmetrize_matrix(H_bb)
  H_beta_lambda <- cbind(H_bphi, H_btau)
  H_lambda <- rbind(
    cbind(H_phiphi, H_phitau),
    cbind(t(H_phitau), H_tautau)
  )
  H_lambda <- .hl_symmetrize_matrix(H_lambda)

  score <- hl_state$score
  if (is.null(score)) {
    stop("hl_state$score is missing; cannot build adjusted statistics.")
  }

  meat <- hl_state$meat
  if (is.null(meat)) {
    score_subject <- hl_state$score_subject
    if (is.null(score_subject)) {
      stop("Neither hl_state$meat nor hl_state$score_subject available for Godambe computation.")
    }
    meat <- tcrossprod(score_subject)
  }
  meat <- .hl_symmetrize_matrix(meat)

  list(
    H_beta = H_beta,
    H_beta_lambda = H_beta_lambda,
    H_lambda = H_lambda,
    S = meat,
    score = score
  )
}

hl_adjusted_tests <- function(unconstrained,
                              constrained,
                              L,
                              b = NULL,
                              ridge_factor = 1e-8,
                              cond_limit = 1e8) {
  if (is.null(constrained) || is.null(constrained$hl_state)) {
    stop("Constrained fit must include hl_state.")
  }
  if (is.null(unconstrained) || is.null(unconstrained$hl_state)) {
    stop("Unconstrained fit must include hl_state.")
  }

  beta_hat <- unconstrained$beta
  beta_con <- constrained$beta
  if (is.null(beta_hat) || is.null(beta_con)) {
    stop("Both fits must include beta coefficients.")
  }

  nb <- length(beta_hat)
  if (is.null(L)) {
    stop("Contrast matrix L must be provided.")
  }
  L <- as.matrix(L)
  if (ncol(L) != nb) {
    stop("Contrast matrix L must have ncol equal to number of fixed effects.")
  }

  d <- nrow(L)
  if (is.null(b)) {
    b <- rep(0, d)
  } else {
    if (length(b) != d) {
      stop("Vector b must have length equal to the number of rows of L.")
    }
  }

  godambe_con <- hl_state_godambe(constrained$hl_state, nb)
  godambe_un <- hl_state_godambe(unconstrained$hl_state, nb)

  H_beta <- godambe_con$H_beta
  H_beta_lambda <- godambe_con$H_beta_lambda
  H_lambda <- godambe_con$H_lambda
  S_full <- godambe_con$S
  score_full <- godambe_con$score

  idx_beta <- seq_len(nb)
  idx_lambda <- nb + seq_len(2L)

  S_bb <- S_full[idx_beta, idx_beta, drop = FALSE]
  S_bl <- S_full[idx_beta, idx_lambda, drop = FALSE]
  S_lb <- S_full[idx_lambda, idx_beta, drop = FALSE]
  S_ll <- S_full[idx_lambda, idx_lambda, drop = FALSE]

  H_lambda_inv_info <- invert_with_ridge(H_lambda,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )
  if (!H_lambda_inv_info$success) {
    stop("Failed to invert nuisance block H_lambda even after ridge adjustment.")
  }
  H_lambda_inv <- H_lambda_inv_info$inverse

  H_beta_dot_lambda <- H_beta - H_beta_lambda %*% H_lambda_inv %*% t(H_beta_lambda)
  H_beta_dot_lambda <- .hl_symmetrize_matrix(H_beta_dot_lambda)
  H_beta_dot_lambda_inv_info <- invert_with_ridge(H_beta_dot_lambda,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )
  if (!H_beta_dot_lambda_inv_info$success) {
    stop("Failed to invert Schur complement H_beta dot lambda even after ridge adjustment.")
  }
  H_beta_dot_lambda_inv <- H_beta_dot_lambda_inv_info$inverse

  score_beta <- score_full[idx_beta]
  score_lambda <- score_full[idx_lambda]
  score_eff <- score_beta - H_beta_lambda %*% H_lambda_inv %*% score_lambda

  S_eff <- S_bb -
    H_beta_lambda %*% H_lambda_inv %*% S_lb -
    S_bl %*% H_lambda_inv %*% t(H_beta_lambda) +
    H_beta_lambda %*% H_lambda_inv %*% S_ll %*% H_lambda_inv %*% t(H_beta_lambda)
  S_eff <- .hl_symmetrize_matrix(S_eff)

  Hpsi <- .hl_symmetrize_matrix(L %*% H_beta_dot_lambda_inv %*% t(L))
  Spsi <- .hl_symmetrize_matrix(L %*% S_eff %*% t(L))
  G_beta <- H_beta_dot_lambda_inv %*% S_eff %*% H_beta_dot_lambda_inv
  G_beta <- .hl_symmetrize_matrix(G_beta)
  Gpsi <- .hl_symmetrize_matrix(L %*% G_beta %*% t(L))

  Up <- as.numeric(L %*% H_beta_dot_lambda_inv %*% score_eff)

  H_beta_un <- godambe_un$H_beta
  H_beta_lambda_un <- godambe_un$H_beta_lambda
  H_lambda_un <- godambe_un$H_lambda
  S_full_un <- godambe_un$S

  S_bb_un <- S_full_un[idx_beta, idx_beta, drop = FALSE]
  S_bl_un <- S_full_un[idx_beta, idx_lambda, drop = FALSE]
  S_lb_un <- S_full_un[idx_lambda, idx_beta, drop = FALSE]
  S_ll_un <- S_full_un[idx_lambda, idx_lambda, drop = FALSE]

  H_lambda_un_inv_info <- invert_with_ridge(H_lambda_un,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )
  H_lambda_un_inv <- H_lambda_un_inv_info$inverse

  H_beta_dot_lambda_un <- H_beta_un - H_beta_lambda_un %*% H_lambda_un_inv %*% t(H_beta_lambda_un)
  H_beta_dot_lambda_un <- .hl_symmetrize_matrix(H_beta_dot_lambda_un)

  H_beta_dot_lambda_un_inv_info <- invert_with_ridge(H_beta_dot_lambda_un,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )
  H_beta_dot_lambda_un_inv <- H_beta_dot_lambda_un_inv_info$inverse

  S_eff_un <- S_bb_un -
    H_beta_lambda_un %*% H_lambda_un_inv %*% S_lb_un -
    S_bl_un %*% H_lambda_un_inv %*% t(H_beta_lambda_un) +
    H_beta_lambda_un %*% H_lambda_un_inv %*% S_ll_un %*% H_lambda_un_inv %*% t(H_beta_lambda_un)
  S_eff_un <- .hl_symmetrize_matrix(S_eff_un)

  Hpsi_un <- .hl_symmetrize_matrix(L %*% H_beta_dot_lambda_un_inv %*% t(L))
  G_beta_un <- H_beta_dot_lambda_un_inv %*% S_eff_un %*% H_beta_dot_lambda_un_inv
  G_beta_un <- .hl_symmetrize_matrix(G_beta_un)
  Gpsi_un <- .hl_symmetrize_matrix(L %*% G_beta_un %*% t(L))

  Spsi_inv_info <- invert_with_ridge(Spsi,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )
  Gpsi_inv_info <- invert_with_ridge(Gpsi,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )
  Hpsi_inv_info <- invert_with_ridge(Hpsi,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )

  psi_hat <- as.numeric(L %*% beta_hat)
  Delta <- psi_hat - b

  aphl_uncon <- unconstrained$hl_state$aphl
  aphl_con <- constrained$hl_state$aphl
  w_raw <- 2 * (aphl_uncon - aphl_con)
  if (!is.na(w_raw) && w_raw < 0 && abs(w_raw) <= 1e-6) {
    w_raw <- 0
  }
  flagged_negative_w_raw <- FALSE
  if (!is.na(w_raw) && w_raw < 0) {
    flagged_negative_w_raw <- TRUE
    w_raw <- NA_real_
  }

  scale <- NA_real_
  w_adj <- NA_real_
  method <- "lrt"
  lr_failure <- NA_character_

  wU_val <- if (Gpsi_inv_info$success) {
    val <- as.numeric(t(Up) %*% Gpsi_inv_info$inverse %*% Up)
    if (!is.na(val) && val < 0 && abs(val) <= 1e-8) 0 else val
  } else {
    NA_real_
  }

  q_prof_val <- if (Hpsi_inv_info$success) {
    val <- as.numeric(t(Up) %*% Hpsi_inv_info$inverse %*% Up)
    if (!is.na(val) && val < 0 && abs(val) <= 1e-8) 0 else val
  } else {
    NA_real_
  }

  fallback <- FALSE

  if (d == 1) {
    if (is.finite(wU_val) && is.finite(q_prof_val) && abs(q_prof_val) > .Machine$double.eps) {
      scale <- wU_val / q_prof_val
      w_adj <- scale * w_raw
      if ((isTRUE(Gpsi_inv_info$ridge > 0) || isTRUE(Hpsi_inv_info$ridge > 0)) && method == "lrt") {
        method <- "lrt_ridge"
      }
    } else {
      h_val <- as.numeric(Hpsi)
      g_val <- as.numeric(Gpsi)
      if (is.finite(h_val) && h_val > 0 && is.finite(g_val) && g_val > 0) {
        scale <- h_val / g_val
        w_adj <- scale * w_raw
        if ((isTRUE(Gpsi_inv_info$ridge > 0) || isTRUE(Hpsi_inv_info$ridge > 0)) && method == "lrt") {
          method <- "lrt_ridge"
        }
      } else {
        method <- "lrt_fail"
      }
    }
    wU <- NA_real_
    q_prof <- NA_real_
  } else {
    if (is.finite(wU_val) && is.finite(q_prof_val) && abs(q_prof_val) > .Machine$double.eps) {
      scale <- wU_val / q_prof_val
      w_adj <- scale * w_raw
      if ((isTRUE(Gpsi_inv_info$ridge > 0) || isTRUE(Hpsi_inv_info$ridge > 0)) && method == "lrt") {
        method <- "lrt_ridge"
      }
    } else {
      method <- "lrt_fail"
      if (!Gpsi_inv_info$success) {
        lr_failure <- "singular_Gpsi"
      } else if (!Hpsi_inv_info$success) {
        lr_failure <- "singular_Hpsi"
      } else {
        lr_failure <- "undefined_scale"
      }
    }
    wU <- wU_val
    q_prof <- q_prof_val
  }

  cov_psi_un <- .hl_symmetrize_matrix(Gpsi_un)

  if (!is.finite(w_adj)) {
    wald_inv_info <- invert_with_ridge(cov_psi_un,
      cond_limit = cond_limit,
      base_factor = ridge_factor
    )
    if (!wald_inv_info$success) {
      wald_inv_info <- invert_with_ridge(cov_psi_un,
        cond_limit = cond_limit,
        base_factor = max(ridge_factor, 1e-6)
      )
    }
    if (wald_inv_info$success) {
      method <- if (wald_inv_info$ridge > 0) "wald_ridge" else "wald"
      fallback <- TRUE
      cov_inv <- wald_inv_info$inverse
      w_adj <- as.numeric(t(Delta) %*% cov_inv %*% Delta)
      scale <- NA_real_
      if (is.na(w_adj) || !is.finite(w_adj)) {
        lr_failure <- "wald_fallback_failed"
      }
    } else {
      method <- "lrt_raw"
      fallback <- TRUE
      w_adj <- max(w_raw, 0)
      scale <- NA_real_
      if (!is.finite(w_raw)) {
        lr_failure <- "missing_profile_lr"
      }
    }
  }
  if (!is.finite(w_adj)) {
    if (is.na(lr_failure)) {
      lr_failure <- if (!is.finite(w_raw)) {
        "missing_profile_lr"
      } else {
        "undefined_adjustment"
      }
    }
  }
  if (is.finite(w_adj)) {
    lr_failure <- NA_character_
  } else if (flagged_negative_w_raw && is.na(lr_failure)) {
    lr_failure <- "negative_raw_lr"
  }
  lr_success <- is.finite(w_adj)

  Q_star <- if (Gpsi_inv_info$success) {
    val <- as.numeric(t(Up) %*% Gpsi_inv_info$inverse %*% Up)
    if (!is.na(val) && val < 0 && abs(val) <= 1e-8) 0 else val
  } else {
    NA_real_
  }

  Wald_star <- NA_real_
  Var_psi <- .hl_symmetrize_matrix(cov_psi_un)
  Gpsi_un_inv_info <- invert_with_ridge(cov_psi_un,
    cond_limit = cond_limit,
    base_factor = ridge_factor
  )
  if (Gpsi_un_inv_info$success) {
    Wald_mat <- Gpsi_un_inv_info$inverse
    Wald_star <- as.numeric(t(Delta) %*% Wald_mat %*% Delta)
    if (!is.na(Wald_star) && Wald_star < 0 && abs(Wald_star) <= 1e-8) {
      Wald_star <- 0
    }
  }
  wald_success <- is.finite(Wald_star)
  wald_failure <- NA_character_
  if (!wald_success) {
    if (!Gpsi_un_inv_info$success) {
      wald_failure <- "singular_cov"
    } else {
      wald_failure <- "undefined_stat"
    }
  }

  score_success <- isTRUE(Spsi_inv_info$success) && is.finite(Q_star)
  score_failure <- NA_character_
  if (!isTRUE(Spsi_inv_info$success)) {
    score_failure <- "singular_Spsi"
  } else if (!is.finite(Q_star)) {
    score_failure <- "undefined_stat"
  }

  list(
    dim = d,
    lr = list(
      raw = w_raw,
      adj = w_adj,
      scale = scale,
      method = method,
      fallback = fallback,
      success = lr_success,
      failure_reason = lr_failure
    ),
    wald = list(
      stat = Wald_star,
      delta = Delta,
      variance = Var_psi,
      ridge = list(
        H = if (Hpsi_inv_info$success) Hpsi_inv_info$ridge else NA_real_,
        G = if (Gpsi_inv_info$success) Gpsi_inv_info$ridge else NA_real_
      ),
      success = wald_success,
      failure_reason = wald_failure
    ),
    score = list(
      stat = Q_star,
      ridge = Spsi_inv_info$ridge,
      success = score_success,
      failure_reason = score_failure
    ),
    blocks = list(
      Hpsi = Hpsi,
      Spsi = Spsi,
      Gpsi = Gpsi,
      Up = Up,
      S_eff = S_eff,
      H_beta_dot_lambda = H_beta_dot_lambda
    ),
    nuisances = list(
      H_lambda = H_lambda,
      H_lambda_ridge = H_lambda_inv_info$ridge
    )
  )
}

hl_adjusted_contrasts <- function(unconstrained,
                                  ctx,
                                  gene_index,
                                  count,
                                  contrasts,
                                  posv = NULL,
                                  ridge_factor = 1e-8,
                                  cond_limit = 1e8) {
  if (is.null(contrasts) || length(contrasts) == 0) {
    return(vector("list", 0))
  }
  if (is.null(posv)) {
    posv <- call_posindy(count, gene_index - 1L, ctx$nind)
  }

  build_def <- function(def, idx) {
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

  defs <- vector("list", length(contrasts))
  for (idx in seq_along(contrasts)) {
    defs[[idx]] <- build_def(contrasts[[idx]], idx)
  }

  results <- vector("list", length(defs))
  names_vec <- character(length(defs))

  for (idx in seq_along(defs)) {
    def <- defs[[idx]]
    fit_con <- hl_fit_constrained(
      ctx = ctx,
      gene_index = gene_index,
      count = count,
      C = def$L,
      posv = posv,
      unconstrained_fit = unconstrained,
      compute_adj = TRUE,
      ridge_factor = ridge_factor,
      cond_limit = cond_limit,
      rhs = def$b
    )

    adj <- fit_con$hl_adj
    df <- nrow(def$L)

    lr_stat <- adj$lr$adj
    if (length(lr_stat) == 0L || all(is.na(lr_stat))) {
      lr_stat <- NA_real_
    }
    lr_p <- if (is.finite(lr_stat)) {
      stats::pchisq(max(lr_stat, 0), df = df, lower.tail = FALSE)
    } else {
      NA_real_
    }

    wald_stat <- adj$wald$stat
    if (length(wald_stat) == 0L || all(is.na(wald_stat))) {
      wald_stat <- NA_real_
    }
    wald_p <- if (is.finite(wald_stat)) {
      stats::pchisq(max(wald_stat, 0), df = df, lower.tail = FALSE)
    } else {
      NA_real_
    }

    score_stat <- adj$score$stat
    if (length(score_stat) == 0L || all(is.na(score_stat))) {
      score_stat <- NA_real_
    }
    score_p <- if (is.finite(score_stat)) {
      stats::pchisq(max(score_stat, 0), df = df, lower.tail = FALSE)
    } else {
      NA_real_
    }

    results[[idx]] <- list(
      label = def$label,
      df = df,
      contrast = def$L,
      rhs = def$b,
      constrained = fit_con,
      adj = adj,
      p_values = list(
        lr = lr_p,
        wald = wald_p,
        score = score_p
      )
    )
    names_vec[[idx]] <- def$label
  }

  names(results) <- names_vec
  results
}
