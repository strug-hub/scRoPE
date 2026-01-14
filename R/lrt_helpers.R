# Utilities for adjusted profile likelihood ratio tests

invert_posdef <- function(mat) {
  tryCatch({
    Rfast::spdinv(mat)
  }, error = function(e) {
    warning("Rfast::spdinv failed, using solve instead")
    solve(mat)
  })
}

invert_with_ridge <- function(mat, cond_limit = 1e8, base_factor = 1e-6, max_iter = 6) {
  d <- nrow(mat)
  if (d == 0) {
    return(list(success = FALSE, ridge = 0, inverse = NULL))
  }

  diag_scale <- mean(abs(diag(mat)))
  if (!is.finite(diag_scale) || diag_scale == 0) {
    diag_scale <- 1
  }
  base_ridge <- diag_scale * base_factor
  ridge <- 0

  cond0 <- tryCatch(kappa(mat), error = function(e) Inf)
  if (!is.finite(cond0) || cond0 > cond_limit) {
    ridge <- base_ridge
  }

  for (iter in seq_len(max_iter + 1)) {
    mat_adj <- if (ridge > 0) mat + diag(ridge, d) else mat
    inv_try <- tryCatch(invert_posdef(mat_adj), error = function(e) NULL)
    if (!is.null(inv_try) && all(is.finite(inv_try))) {
      return(list(success = TRUE, ridge = ridge, inverse = inv_try))
    }
    if (iter > max_iter) {
      break
    }
    ridge <- if (ridge == 0) base_ridge else ridge * 10
  }

  list(success = FALSE, ridge = ridge, inverse = NULL)
}

solve_sym <- function(mat, rhs) {
  chol_dec <- tryCatch({
    chol(mat)
  }, error = function(e) NULL)
  if (is.null(chol_dec)) {
    solve(mat, rhs)
  } else {
    backsolve(chol_dec, forwardsolve(t(chol_dec), rhs))
  }
}

compute_profile_lrt <- function(unconstrained, constrained, L, b) {
  beta_hat <- unconstrained$theta$beta
  beta_con <- constrained$theta$beta
  if (is.null(beta_hat) || is.null(beta_con)) {
    stop("Both unconstrained and constrained fits must provide theta$beta.")
  }

  pick_loglik <- function(fit) {
    if (!is.null(fit$loglik_ln) && is.finite(fit$loglik_ln)) {
      return(fit$loglik_ln)
    }
    fit$loglik
  }

  loglik_uncon <- pick_loglik(unconstrained)
  loglik_con <- pick_loglik(constrained)
  wP <- 2 * (loglik_uncon - loglik_con)
  if (!is.na(wP) && wP < 0 && abs(wP) <= 1e-6) {
    wP <- 0
  }

  godambe_con <- constrained$godambe
  if (is.null(godambe_con)) {
    stop("Constrained fit lacks Godambe summary; rerun to compute diagnostics.")
  }

  nb <- length(beta_hat)
  if (ncol(L) != nb) {
    stop("L must have the same number of columns as beta.")
  }

  H_beta <- godambe_con$H_beta
  H_beta_lambda <- godambe_con$H_beta_lambda
  H_lambda <- godambe_con$H_lambda
  score_full <- constrained$score
  if (is.null(score_full)) {
    stop("Constrained fit lacks gradient for score computation.")
  }
  score_beta <- score_full[seq_len(nb)]
  score_lambda <- score_full[(nb + 1):(nb + 2)]

  H_lambda_inv <- invert_posdef(H_lambda)
  score_eff <- score_beta - H_beta_lambda %*% H_lambda_inv %*% score_lambda

  H_beta_dot_lambda <- H_beta - H_beta_lambda %*% H_lambda_inv %*% t(H_beta_lambda)
  H_beta_dot_lambda_inv <- invert_posdef(H_beta_dot_lambda)
  Hpsi_raw <- L %*% H_beta_dot_lambda_inv %*% t(L)
  method <- "lrt"

  S_full <- godambe_con$S
  if (is.null(S_full)) {
    stop("Constrained fit lacks Godambe meat matrix S.")
  }
  idx_beta <- seq_len(nb)
  idx_lambda <- nb + seq_len(nrow(S_full) - nb)
  S_bb <- S_full[idx_beta, idx_beta, drop = FALSE]
  S_bl <- S_full[idx_beta, idx_lambda, drop = FALSE]
  S_lb <- S_full[idx_lambda, idx_beta, drop = FALSE]
  S_ll <- S_full[idx_lambda, idx_lambda, drop = FALSE]

  S_eff <- S_bb -
    H_beta_lambda %*% H_lambda_inv %*% S_lb -
    S_bl %*% H_lambda_inv %*% t(H_beta_lambda) +
    H_beta_lambda %*% H_lambda_inv %*% S_ll %*% H_lambda_inv %*% t(H_beta_lambda)
  S_eff <- 0.5 * (S_eff + t(S_eff))

  Spsi_raw <- L %*% S_eff %*% t(L)

  G_beta <- H_beta_dot_lambda_inv %*% S_eff %*% H_beta_dot_lambda_inv
  Gpsi_raw <- L %*% G_beta %*% t(L)

  Hpsi <- 0.5 * (Hpsi_raw + t(Hpsi_raw))
  Spsi <- 0.5 * (Spsi_raw + t(Spsi_raw))
  Gpsi <- 0.5 * (Gpsi_raw + t(Gpsi_raw))

  profile_score <- L %*% H_beta_dot_lambda_inv %*% score_eff

  d0 <- nrow(L)
  fallback <- FALSE

  Gpsi_inv_info <- invert_with_ridge(Gpsi)
  wU_val <- if (Gpsi_inv_info$success) {
    val <- as.numeric(t(profile_score) %*% Gpsi_inv_info$inverse %*% profile_score)
    if (!is.na(val) && val < 0 && abs(val) <= 1e-8) 0 else val
  } else {
    NA_real_
  }

  Hpsi_inv_info <- invert_with_ridge(Hpsi)
  qP_val <- if (Hpsi_inv_info$success) {
    val <- as.numeric(t(profile_score) %*% Hpsi_inv_info$inverse %*% profile_score)
    if (!is.na(val) && val < 0 && abs(val) <= 1e-8) 0 else val
  } else {
    NA_real_
  }

  scale <- NA_real_
  wP_adj <- NA_real_
  success <- FALSE
  failure_reason <- NA_character_
  set_failure <- function(code, overwrite = FALSE) {
    if (!overwrite && !is.na(failure_reason) && nzchar(failure_reason)) {
      return()
    }
    failure_reason <<- normalise_failure(code, "lrt")
  }
  if (!is.finite(wP)) {
    set_failure("missing_profile_lr", overwrite = TRUE)
  }
  negative_raw_lr <- is.finite(wP) && wP < -1e-6

  if (is.finite(wU_val) && is.finite(qP_val) && abs(qP_val) > .Machine$double.eps) {
    scale <- wU_val / qP_val
    wP_adj <- scale * wP
    if (Gpsi_inv_info$ridge > 0 && method == "lrt") {
      method <- "lrt_ridge"
    }
    if (is.finite(wP_adj)) {
      success <- TRUE
    }
  } else {
    h_val <- as.numeric(Hpsi)
    g_val <- as.numeric(Gpsi)
    if (is.finite(h_val) && h_val > 0 && is.finite(g_val) && g_val > 0) {
      scale <- h_val / g_val
      wP_adj <- scale * wP
      if (Gpsi_inv_info$ridge > 0 && method == "lrt") {
        method <- "lrt_ridge"
      }
      if (is.finite(wP_adj)) {
        success <- TRUE
      }
    } else {
      method <- "lrt_fail"
    }
  }

  if (!success) {
    if (!Hpsi_inv_info$success) {
      set_failure("singular_Hpsi")
    } else if (!Gpsi_inv_info$success) {
      set_failure("singular_Gpsi")
    } else if (!is.finite(qP_val) || abs(qP_val) <= .Machine$double.eps) {
      set_failure("undefined_scale")
    }
  } else if (!is.na(failure_reason)) {
    failure_reason <- NA_character_
  }

  wU <- if (d0 == 1) NA_real_ else wU_val
  qP <- if (d0 == 1) NA_real_ else qP_val

  ridge_info <- list(
    Hpsi = if (Hpsi_inv_info$success) Hpsi_inv_info$ridge else NA_real_,
    Gpsi = if (Gpsi_inv_info$success) Gpsi_inv_info$ridge else NA_real_,
    wald = 0
  )

  if (!is.finite(wP_adj)) {
    beta_diff <- as.matrix(L %*% beta_hat - b)
    wald_cov_raw <- L %*% G_beta %*% t(L)
    wald_inv_info <- invert_with_ridge(wald_cov_raw)
    if (wald_inv_info$success) {
      ridge_info$wald <- wald_inv_info$ridge
      method <- if (wald_inv_info$ridge > 0) "wald_ridge" else "wald"
      fallback <- TRUE
      cov_inv <- wald_inv_info$inverse
      wP_adj <- as.numeric(t(beta_diff) %*% cov_inv %*% beta_diff)
      scale <- NA_real_
      wU <- NA_real_
      qP <- NA_real_
      set_failure("undefined_adjustment")
    } else {
      method <- "lrt_raw"
      fallback <- TRUE
      ridge_info$wald <- NA_real_
      wP_adj <- max(wP, 0)
      scale <- NA_real_
      wU <- NA_real_
      qP <- NA_real_
      set_failure("wald_fallback_failed", overwrite = TRUE)
    }
  }

  if (!success) {
    if (is.na(failure_reason) && negative_raw_lr) {
      set_failure("negative_raw_lr", overwrite = TRUE)
    } else if (is.na(failure_reason)) {
      set_failure("undefined_adjustment", overwrite = TRUE)
    }
  }

  list(
    wP = wP,
    wU = wU,
    qP = qP,
    scale = scale,
    wP_adjusted = wP_adj,
    score = profile_score,
    Hpsi = Hpsi_raw,
    Gpsi = Gpsi_raw,
    Gbeta = G_beta,
    Jpsi = Spsi_raw,
    method = method,
    ridge = ridge_info,
    fallback = fallback,
    success = success,
    failure_reason = failure_reason
  )
}

compute_adj_score <- function(constrained, L, ridge_factor = 1e-8, check_scalar_identity = FALSE, tol = 1e-6) {
  godambe <- constrained$godambe
  if (is.null(godambe)) {
    stop("Constrained fit lacks Godambe summary; cannot compute adjusted score.")
  }
  score_full <- constrained$score
  if (is.null(score_full)) {
    stop("Constrained fit lacks gradient for score computation.")
  }

  H_beta <- godambe$H_beta
  H_beta_lambda <- godambe$H_beta_lambda
  H_lambda <- godambe$H_lambda

  nb <- nrow(H_beta)
  if (ncol(L) != nb) {
    stop("L must have the same number of columns as beta.")
  }

  n_lambda <- nrow(H_lambda)
  idx_beta <- seq_len(nb)
  idx_lambda <- nb + seq_len(n_lambda)

  H_lambda_inv <- invert_posdef(H_lambda)
  H_beta_dot_lambda <- H_beta - H_beta_lambda %*% H_lambda_inv %*% t(H_beta_lambda)
  H_beta_dot_lambda_inv <- invert_posdef(H_beta_dot_lambda)

  score_beta <- matrix(score_full[idx_beta], ncol = 1)
  score_lambda <- matrix(score_full[idx_lambda], ncol = 1)
  score_eff <- score_beta - H_beta_lambda %*% H_lambda_inv %*% score_lambda

  Upsi <- L %*% H_beta_dot_lambda_inv %*% score_eff

  S_full <- godambe$S
  if (is.null(S_full)) {
    stop("Constrained fit lacks Godambe meat matrix S.")
  }
  S_bb <- S_full[idx_beta, idx_beta, drop = FALSE]
  S_bl <- S_full[idx_beta, idx_lambda, drop = FALSE]
  S_lb <- S_full[idx_lambda, idx_beta, drop = FALSE]
  S_ll <- S_full[idx_lambda, idx_lambda, drop = FALSE]

  S_eff <- S_bb -
    H_beta_lambda %*% H_lambda_inv %*% S_lb -
    S_bl %*% H_lambda_inv %*% t(H_beta_lambda) +
    H_beta_lambda %*% H_lambda_inv %*% S_ll %*% H_lambda_inv %*% t(H_beta_lambda)
  S_eff <- 0.5 * (S_eff + t(S_eff))

  if (is.null(ridge_factor) || length(ridge_factor) == 0 || !is.finite(ridge_factor)) {
    ridge_factor <- 0
  }
  ridge_factor <- max(ridge_factor, 0)
  stat <- NA_real_
  pval <- NA_real_
  method <- "adj-score"

  base_factor <- if (ridge_factor > 0) ridge_factor else 0
  G_beta <- H_beta_dot_lambda_inv %*% S_eff %*% H_beta_dot_lambda_inv
  Gpsi <- L %*% G_beta %*% t(L)
  Gpsi <- 0.5 * (Gpsi + t(Gpsi))
  Gpsi_inv_info <- invert_with_ridge(Gpsi, base_factor = base_factor)

  if (Gpsi_inv_info$success) {
    Gpsi_inv <- Gpsi_inv_info$inverse
    stat <- as.numeric(t(Upsi) %*% Gpsi_inv %*% Upsi)
    if (!is.na(stat) && stat < 0 && abs(stat) <= 1e-8) {
      stat <- 0
    }
    pval <- pchisq(max(stat, 0), df = nrow(L), lower.tail = FALSE)
  }

  Hpsi <- L %*% H_beta_dot_lambda_inv %*% t(L)
  Hpsi <- 0.5 * (Hpsi + t(Hpsi))

  scalar_check <- NA_real_
  if (check_scalar_identity && nrow(L) == 1 && Gpsi_inv_info$success) {
    scalar_check <- 0
  }

  score_success <- isTRUE(Gpsi_inv_info$success) && is.finite(stat) && is.finite(pval)
  score_failure <- NA_character_
  if (!isTRUE(Gpsi_inv_info$success)) {
    score_failure <- normalise_failure("singular_Spsi", "score")
  } else if (!is.finite(stat) || !is.finite(pval)) {
    score_failure <- normalise_failure("undefined_stat", "score")
  }

  list(
    method = method,
    df = nrow(L),
    stat = stat,
    pval = pval,
    ridge = Gpsi_inv_info$ridge,
    success = score_success,
    score = Upsi,
    Spsi = L %*% S_eff %*% t(L),
    Gpsi = Gpsi,
    Hpsi = Hpsi,
    scalar_check = scalar_check,
    failure_reason = score_failure
  )
}
