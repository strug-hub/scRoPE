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

invert_general_with_ridge <- function(mat, cond_limit = 1e8, base_factor = 1e-6, max_iter = 6) {
  d <- nrow(mat)
  if (d == 0) {
    return(list(success = FALSE, ridge = 0, inverse = NULL))
  }

  mat <- 0.5 * (mat + t(mat))
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
    inv_try <- tryCatch(solve(mat_adj), error = function(e) NULL)
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

is_posdef_matrix <- function(mat, tol = 1e-10) {
  if (is.null(mat) || any(!is.finite(mat))) {
    return(FALSE)
  }
  mat <- 0.5 * (mat + t(mat))
  if (nrow(mat) == 1L) {
    return(as.numeric(mat) > tol)
  }
  !is.null(tryCatch(chol(mat), error = function(e) NULL))
}

extract_lrt_gate <- function(fit) {
  gate <- fit$lrt_gate
  if (!is.null(gate)) {
    return(gate)
  }

  converged <- NA
  if (!is.null(fit$convergence)) {
    converged <- isTRUE(fit$convergence > 0)
  } else if (!is.null(fit$debug$conv_code)) {
    converged <- isTRUE(fit$debug$conv_code >= 0)
  }

  list(
    converged = converged,
    interior = NA,
    nuisance_interior = NA
  )
}

scalar_profile_reparameterization <- function(L) {
  if (is.null(dim(L))) {
    L <- matrix(L, nrow = 1L)
  }
  if (!is.matrix(L) || nrow(L) != 1L) {
    stop("Scalar profile helper requires L to be a one-row matrix or vector.")
  }

  c_vec <- as.numeric(L[1, ])
  if (!all(is.finite(c_vec))) {
    stop("L must be finite for scalar profile calculations.")
  }
  if (sum(abs(c_vec)) <= 0) {
    stop("L must not be the zero vector.")
  }

  constr <- constraint_reparameterization(matrix(c_vec, nrow = 1L), 0)
  K <- constr$null_basis
  if (is.null(K)) {
    K <- matrix(0, nrow = length(c_vec), ncol = 0L)
  }

  list(
    L = matrix(c_vec, nrow = 1L),
    c = c_vec,
    a = c_vec / sum(c_vec^2),
    K = K,
    free_dim = ncol(K)
  )
}

scalar_profile_map <- function(target, n_lambda) {
  nb <- length(target$c)
  free_dim <- target$free_dim
  map <- matrix(0, nrow = nb + n_lambda, ncol = 1L + free_dim + n_lambda)
  map[seq_len(nb), 1L] <- target$a
  if (free_dim > 0L) {
    map[seq_len(nb), 1L + seq_len(free_dim)] <- target$K
  }
  if (n_lambda > 0L) {
    lambda_cols <- 1L + free_dim + seq_len(n_lambda)
    map[nb + seq_len(n_lambda), lambda_cols] <- diag(n_lambda)
  }
  map
}

compute_scalar_profile_quantities <- function(full_score,
                                              full_per_subject_scores,
                                              observed_info_full,
                                              L,
                                              ridge_factor = 1e-8,
                                              cond_limit = 1e8,
                                              tol = 1e-8) {
  target <- scalar_profile_reparameterization(L)
  full_score <- as.numeric(full_score)
  observed_info_full <- as.matrix(observed_info_full)
  full_per_subject_scores <- as.matrix(full_per_subject_scores)

  nb <- length(target$c)
  p <- length(full_score)
  n_lambda <- p - nb
  if (n_lambda < 0L) {
    stop("Full score dimension is smaller than the fixed-effect dimension.")
  }
  if (!all(dim(observed_info_full) == c(p, p))) {
    stop("observed_info_full must be a square matrix matching full_score.")
  }
  if (nrow(full_per_subject_scores) != p) {
    stop("full_per_subject_scores must have the same number of rows as full_score.")
  }

  map <- scalar_profile_map(target, n_lambda)
  score_red <- as.numeric(crossprod(map, full_score))
  per_subject_red <- crossprod(map, full_per_subject_scores)
  observed_info_red <- crossprod(map, observed_info_full %*% map)
  observed_info_red <- 0.5 * (observed_info_red + t(observed_info_red))

  eta_idx <- if (ncol(observed_info_red) > 1L) 2:ncol(observed_info_red) else integer(0)
  U_psi <- score_red[[1L]]
  U_eta <- score_red[eta_idx]
  U_psi_i <- as.numeric(per_subject_red[1L, ])
  U_eta_i <- per_subject_red[eta_idx, , drop = FALSE]

  H_psipsi <- matrix(observed_info_red[1L, 1L], nrow = 1L, ncol = 1L)
  H_psi_eta <- observed_info_red[1L, eta_idx, drop = FALSE]
  H_eta_psi <- observed_info_red[eta_idx, 1L, drop = FALSE]
  H_eta_eta <- observed_info_red[eta_idx, eta_idx, drop = FALSE]

  base_factor <- max(if (is.null(ridge_factor) || !is.finite(ridge_factor)) 0 else ridge_factor, 0)
  H_etaeta_inv_info <- if (length(eta_idx) == 0L) {
    list(success = TRUE, ridge = 0, inverse = matrix(0, nrow = 0L, ncol = 0L))
  } else {
    invert_general_with_ridge(
      H_eta_eta,
      cond_limit = cond_limit,
      base_factor = base_factor
    )
  }

  profile_subject_scores <- rep(NA_real_, ncol(per_subject_red))
  profile_score <- NA_real_
  Jpsi <- NA_real_
  Hpsi <- NA_real_
  Gpsi <- NA_real_
  wU <- NA_real_
  qP <- NA_real_
  scale <- NA_real_
  Vpsi_model <- NA_real_
  Vpsi_godambe <- NA_real_

  if (isTRUE(H_etaeta_inv_info$success)) {
    H_eta_eta_inv <- H_etaeta_inv_info$inverse
    if (length(eta_idx) == 0L) {
      profile_subject_scores <- U_psi_i
      profile_score <- U_psi
      Hpsi <- as.numeric(H_psipsi)
    } else {
      profile_subject_scores <- as.numeric(U_psi_i - H_psi_eta %*% H_eta_eta_inv %*% U_eta_i)
      profile_score <- as.numeric(U_psi - H_psi_eta %*% H_eta_eta_inv %*% U_eta)
      Hpsi <- as.numeric(H_psipsi - H_psi_eta %*% H_eta_eta_inv %*% H_eta_psi)
    }

    if (all(is.finite(profile_subject_scores))) {
      Jpsi <- sum(profile_subject_scores^2)
    }
    if (is.finite(profile_score) && is.finite(Jpsi) && Jpsi > tol) {
      wU <- (profile_score^2) / Jpsi
      if (wU < 0 && abs(wU) <= tol) {
        wU <- 0
      }
    }
    if (is.finite(profile_score) && is.finite(Hpsi) && Hpsi > tol) {
      qP <- (profile_score^2) / Hpsi
      if (qP < 0 && abs(qP) <= tol) {
        qP <- 0
      }
    }
    if (is.finite(Hpsi) && Hpsi > tol) {
      Vpsi_model <- 1 / Hpsi
    }
    if (is.finite(Hpsi) && Hpsi > tol && is.finite(Jpsi) && Jpsi > tol) {
      Gpsi <- (Hpsi^2) / Jpsi
      scale <- Hpsi / Jpsi
      if (is.finite(Gpsi) && Gpsi > tol) {
        Vpsi_godambe <- 1 / Gpsi
      }
    }
  }

  list(
    target = target,
    map = map,
    score_blocks = list(
      psi = U_psi,
      eta = U_eta,
      reduced = score_red
    ),
    per_subject_scores = list(
      psi = U_psi_i,
      eta = U_eta_i,
      reduced = per_subject_red
    ),
    observed_info_blocks = list(
      psi_psi = H_psipsi,
      psi_eta = H_psi_eta,
      eta_psi = H_eta_psi,
      eta_eta = H_eta_eta,
      reduced = observed_info_red
    ),
    profile = list(
      subject_scores = profile_subject_scores,
      score = profile_score,
      variability = Jpsi,
      sensitivity = Hpsi,
      godambe = Gpsi,
      wU = wU,
      qP = qP,
      scale = scale,
      Vpsi_model = Vpsi_model,
      Vpsi_godambe = Vpsi_godambe
    ),
    diagnostics = list(
      H_etaeta_invertible = isTRUE(H_etaeta_inv_info$success),
      H_etaeta_ridge = if (isTRUE(H_etaeta_inv_info$success)) H_etaeta_inv_info$ridge else NA_real_,
      Hpsi_positive = is.finite(Hpsi) && Hpsi > tol,
      Jpsi_positive = is.finite(Jpsi) && Jpsi > tol
    )
  )
}

compute_profile_lrt <- function(unconstrained,
                                constrained,
                                L,
                                b,
                                ridge_factor = 1e-8,
                                cond_limit = 1e8) {
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
  if (is.null(ridge_factor) || length(ridge_factor) == 0 || !is.finite(ridge_factor)) {
    ridge_factor <- 0
  }
  ridge_factor <- max(ridge_factor, 0)
  base_factor <- if (ridge_factor > 0) ridge_factor else 0

  if (nrow(L) == 1L) {
    profile_scalar <- compute_scalar_profile_quantities(
      full_score = constrained$score,
      full_per_subject_scores = constrained$per_subject_gradients,
      observed_info_full = constrained$observed_info,
      L = L,
      ridge_factor = ridge_factor,
      cond_limit = cond_limit
    )

    Hpsi_raw <- matrix(profile_scalar$profile$Vpsi_model, nrow = 1L, ncol = 1L)
    Gpsi_raw <- matrix(profile_scalar$profile$Vpsi_godambe, nrow = 1L, ncol = 1L)
    profile_score_val <- if (is.finite(profile_scalar$profile$score) &&
      is.finite(profile_scalar$profile$Vpsi_model)) {
      profile_scalar$profile$Vpsi_model * profile_scalar$profile$score
    } else {
      profile_scalar$profile$score
    }
    if (is.finite(profile_score_val) && abs(profile_score_val) <= 1e-8) {
      profile_score_val <- 0
    }
    profile_score <- matrix(profile_score_val, nrow = 1L, ncol = 1L)
    Hpsi <- matrix(profile_scalar$profile$sensitivity, nrow = 1L, ncol = 1L)
    Spsi <- matrix(profile_scalar$profile$variability, nrow = 1L, ncol = 1L)
    Gpsi <- matrix(profile_scalar$profile$godambe, nrow = 1L, ncol = 1L)
    wU_val <- profile_scalar$profile$wU
    qP_val <- profile_scalar$profile$qP
    scale <- profile_scalar$profile$scale

    method <- "lrt"
    fallback <- FALSE
    success <- FALSE
    failure_reason <- NA_character_
    set_failure <- function(code, overwrite = FALSE) {
      if (!overwrite && !is.na(failure_reason) && nzchar(failure_reason)) {
        return()
      }
      failure_reason <<- normalise_failure(code, "lrt")
    }
    adjusted_failure_reason <- NA_character_

    if (!is.finite(wP)) {
      set_failure("missing_profile_lr", overwrite = TRUE)
    }
    negative_raw_lr <- is.finite(wP) && wP < -1e-6

    wP_adj <- NA_real_
    if (is.finite(scale) && scale > 0 && is.finite(wP) && wP >= -1e-6) {
      wP_adj <- scale * max(wP, 0)
      if (!is.na(wP_adj) && wP_adj < 0 && abs(wP_adj) <= 1e-8) {
        wP_adj <- 0
      }
      if (is.finite(wP_adj) && wP_adj >= 0) {
        success <- TRUE
      } else {
        method <- "lrt_fail"
      }
    } else {
      method <- "lrt_fail"
    }

    fit_gate <- extract_lrt_gate(constrained)
    diagnostics <- list(
      converged = isTRUE(fit_gate$converged),
      nuisance_hessian_invertible = isTRUE(profile_scalar$diagnostics$H_etaeta_invertible),
      interior = isTRUE(fit_gate$interior),
      nuisance_interior = isTRUE(fit_gate$nuisance_interior),
      profiled_curvature_positive = isTRUE(profile_scalar$diagnostics$Hpsi_positive),
      scale_positive = is.finite(scale) && scale > 0
    )
    adjusted_available <- all(unlist(diagnostics))

    if (!adjusted_available) {
      if (!is.finite(wP)) {
        adjusted_failure_reason <- normalise_failure("missing_profile_lr", "lrt")
      } else if (negative_raw_lr) {
        adjusted_failure_reason <- normalise_failure("negative_raw_lr", "lrt")
      } else if (!diagnostics$converged) {
        adjusted_failure_reason <- normalise_failure("nonconverged_constrained", "lrt")
      } else if (!diagnostics$nuisance_hessian_invertible) {
        adjusted_failure_reason <- normalise_failure("singular_Hlambda", "lrt")
      } else if (!diagnostics$interior || !diagnostics$nuisance_interior) {
        adjusted_failure_reason <- normalise_failure("boundary_constrained", "lrt")
      } else if (!diagnostics$profiled_curvature_positive) {
        adjusted_failure_reason <- normalise_failure("nonpositive_curvature", "lrt")
      } else if (!diagnostics$scale_positive) {
        adjusted_failure_reason <- normalise_failure("nonpositive_scale", "lrt")
      } else {
        adjusted_failure_reason <- normalise_failure("undefined_adjustment", "lrt")
      }
    }

    success <- adjusted_available && is.finite(wP_adj) && wP_adj >= 0

    if (!success) {
      if (!is.na(adjusted_failure_reason)) {
        set_failure(adjusted_failure_reason, overwrite = TRUE)
      } else if (!isTRUE(profile_scalar$diagnostics$H_etaeta_invertible)) {
        set_failure("singular_Hpsi")
      } else if (!is.finite(qP_val) || abs(qP_val) <= .Machine$double.eps) {
        set_failure("undefined_scale")
      }
    } else if (!is.na(failure_reason)) {
      failure_reason <- NA_character_
    }

    wU <- NA_real_
    qP <- NA_real_
    ridge_info <- list(
      Hpsi = if (isTRUE(profile_scalar$diagnostics$Hpsi_positive)) 0 else NA_real_,
      Gpsi = if (is.finite(profile_scalar$profile$godambe) && profile_scalar$profile$godambe > 0) 0 else NA_real_,
      wald = 0
    )

    if (is.finite(wP_adj) && wP_adj < 0 && abs(wP_adj) <= 1e-8) {
      wP_adj <- 0
    }

    if (!adjusted_available || !is.finite(wP_adj) || wP_adj < 0) {
      beta_diff <- as.matrix(L %*% beta_hat - b)
      wald_cov_raw <- if (!is.null(unconstrained$robust_cov) &&
        all(dim(unconstrained$robust_cov) == c(nb, nb))) {
        L %*% unconstrained$robust_cov %*% t(L)
      } else {
        Gpsi_raw
      }
      wald_cov_raw <- 0.5 * (wald_cov_raw + t(wald_cov_raw))
      wald_inv_info <- invert_with_ridge(
        wald_cov_raw,
        cond_limit = cond_limit,
        base_factor = base_factor
      )
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

    return(list(
      wP = wP,
      wU = wU,
      qP = qP,
      scale = scale,
      wP_adjusted = wP_adj,
      score = profile_score,
      Hpsi = Hpsi,
      Gpsi = Gpsi,
      Gbeta = matrix(NA_real_, nrow = nb, ncol = nb),
      Jpsi = Spsi,
      Vpsi_model = Hpsi_raw,
      Vpsi_godambe = Gpsi_raw,
      method = method,
      ridge = ridge_info,
      fallback = fallback,
      success = success,
      adjusted_available = adjusted_available,
      adjusted_failure_reason = adjusted_failure_reason,
      gate = diagnostics,
      failure_reason = failure_reason
    ))
  }

  H_beta <- godambe_con$H_beta
  H_beta_lambda <- godambe_con$H_beta_lambda
  H_lambda <- godambe_con$H_lambda
  score_full <- constrained$score
  if (is.null(score_full)) {
    stop("Constrained fit lacks gradient for score computation.")
  }
  method <- "lrt"
  n_lambda <- nrow(H_lambda)
  idx_beta <- seq_len(nb)
  idx_lambda <- nb + seq_len(n_lambda)

  H_lambda_inv_info <- invert_general_with_ridge(
    H_lambda,
    cond_limit = cond_limit,
    base_factor = base_factor
  )

  score_beta <- matrix(score_full[idx_beta], ncol = 1)
  score_lambda <- matrix(score_full[idx_lambda], ncol = 1)
  profile_score <- matrix(NA_real_, nrow = nrow(L), ncol = 1)
  score_eff <- matrix(NA_real_, nrow = nb, ncol = 1)
  H_beta_dot_lambda <- matrix(NA_real_, nrow = nb, ncol = nb)
  Hpsi_raw <- matrix(NA_real_, nrow = nrow(L), ncol = nrow(L))
  Spsi_raw <- matrix(NA_real_, nrow = nrow(L), ncol = nrow(L))
  Gpsi_raw <- matrix(NA_real_, nrow = nrow(L), ncol = nrow(L))
  S_eff <- matrix(NA_real_, nrow = nb, ncol = nb)
  G_beta <- matrix(NA_real_, nrow = nb, ncol = nb)
  H_beta_dot_lambda_inv_info <- list(success = FALSE, ridge = NA_real_, inverse = NULL)

  if (H_lambda_inv_info$success) {
    H_lambda_inv <- H_lambda_inv_info$inverse
    score_eff <- score_beta - H_beta_lambda %*% H_lambda_inv %*% score_lambda

    H_beta_dot_lambda <- H_beta - H_beta_lambda %*% H_lambda_inv %*% t(H_beta_lambda)
    H_beta_dot_lambda <- 0.5 * (H_beta_dot_lambda + t(H_beta_dot_lambda))
    H_beta_dot_lambda_inv_info <- invert_general_with_ridge(
      H_beta_dot_lambda,
      cond_limit = cond_limit,
      base_factor = base_factor
    )

    if (H_beta_dot_lambda_inv_info$success) {
      H_beta_dot_lambda_inv <- H_beta_dot_lambda_inv_info$inverse
      Hpsi_raw <- L %*% H_beta_dot_lambda_inv %*% t(L)
      profile_score <- L %*% H_beta_dot_lambda_inv %*% score_eff

      S_full <- godambe_con$S
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

      Spsi_raw <- L %*% S_eff %*% t(L)
      G_beta <- H_beta_dot_lambda_inv %*% S_eff %*% H_beta_dot_lambda_inv
      G_beta <- 0.5 * (G_beta + t(G_beta))
      Gpsi_raw <- L %*% G_beta %*% t(L)
    }
  }

  Hpsi <- 0.5 * (Hpsi_raw + t(Hpsi_raw))
  Spsi <- 0.5 * (Spsi_raw + t(Spsi_raw))
  Gpsi <- 0.5 * (Gpsi_raw + t(Gpsi_raw))

  d0 <- nrow(L)
  fallback <- FALSE

  Gpsi_inv_info <- invert_with_ridge(
    Gpsi,
    cond_limit = cond_limit,
    base_factor = base_factor
  )
  wU_val <- if (Gpsi_inv_info$success) {
    val <- as.numeric(t(profile_score) %*% Gpsi_inv_info$inverse %*% profile_score)
    if (!is.na(val) && val < 0 && abs(val) <= 1e-8) 0 else val
  } else {
    NA_real_
  }

  Hpsi_inv_info <- invert_general_with_ridge(
    Hpsi,
    cond_limit = cond_limit,
    base_factor = base_factor
  )
  qP_val <- if (Hpsi_inv_info$success) {
    val <- as.numeric(t(profile_score) %*% Hpsi_inv_info$inverse %*% profile_score)
    if (!is.na(val) && val < 0 && abs(val) <= 1e-8) 0 else val
  } else {
    NA_real_
  }
  Hpsi_info <- if (Hpsi_inv_info$success) {
    0.5 * (Hpsi_inv_info$inverse + t(Hpsi_inv_info$inverse))
  } else {
    matrix(NA_real_, nrow = nrow(L), ncol = nrow(L))
  }
  Gpsi_info <- if (Gpsi_inv_info$success) {
    0.5 * (Gpsi_inv_info$inverse + t(Gpsi_inv_info$inverse))
  } else {
    matrix(NA_real_, nrow = nrow(L), ncol = nrow(L))
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
  adjusted_failure_reason <- NA_character_
  if (!is.finite(wP)) {
    set_failure("missing_profile_lr", overwrite = TRUE)
  }
  negative_raw_lr <- is.finite(wP) && wP < -1e-6

  if (is.finite(wU_val) && is.finite(qP_val) && abs(qP_val) > .Machine$double.eps) {
    scale_candidate <- wU_val / qP_val
    if (is.finite(scale_candidate) && scale_candidate > 0) {
      scale <- scale_candidate
      wP_adj <- scale * wP
      if (!is.na(wP_adj) && wP_adj < 0 && abs(wP_adj) <= 1e-8) {
        wP_adj <- 0
      }
      if (Gpsi_inv_info$ridge > 0 && method == "lrt") {
        method <- "lrt_ridge"
      }
      if (is.finite(wP_adj) && wP_adj >= 0) {
        success <- TRUE
      }
    } else {
      method <- "lrt_fail"
    }
  } else {
    h_val <- as.numeric(Hpsi)
    g_val <- as.numeric(Gpsi)
    if (is.finite(h_val) && h_val > 0 && is.finite(g_val) && g_val > 0) {
      scale <- h_val / g_val
      wP_adj <- scale * wP
      if (!is.na(wP_adj) && wP_adj < 0 && abs(wP_adj) <= 1e-8) {
        wP_adj <- 0
      }
      if (Gpsi_inv_info$ridge > 0 && method == "lrt") {
        method <- "lrt_ridge"
      }
      if (is.finite(wP_adj) && wP_adj >= 0) {
        success <- TRUE
      }
    } else {
      method <- "lrt_fail"
    }
  }

  fit_gate <- extract_lrt_gate(constrained)
  diagnostics <- list(
    converged = isTRUE(fit_gate$converged),
    nuisance_hessian_invertible = isTRUE(H_lambda_inv_info$success),
    interior = isTRUE(fit_gate$interior),
    nuisance_interior = isTRUE(fit_gate$nuisance_interior),
    profiled_curvature_positive = is_posdef_matrix(Hpsi_info),
    scale_positive = is.finite(scale) && scale > 0
  )
  adjusted_available <- all(unlist(diagnostics))

  if (!adjusted_available) {
    if (!is.finite(wP)) {
      adjusted_failure_reason <- normalise_failure("missing_profile_lr", "lrt")
    } else if (negative_raw_lr) {
      adjusted_failure_reason <- normalise_failure("negative_raw_lr", "lrt")
    } else if (!diagnostics$converged) {
      adjusted_failure_reason <- normalise_failure("nonconverged_constrained", "lrt")
    } else if (!diagnostics$nuisance_hessian_invertible) {
      adjusted_failure_reason <- normalise_failure("singular_Hlambda", "lrt")
    } else if (!diagnostics$interior || !diagnostics$nuisance_interior) {
      adjusted_failure_reason <- normalise_failure("boundary_constrained", "lrt")
    } else if (!diagnostics$profiled_curvature_positive) {
      adjusted_failure_reason <- normalise_failure("nonpositive_curvature", "lrt")
    } else if (!diagnostics$scale_positive) {
      adjusted_failure_reason <- normalise_failure("nonpositive_scale", "lrt")
    } else {
      adjusted_failure_reason <- normalise_failure("undefined_adjustment", "lrt")
    }
  }

  success <- adjusted_available && is.finite(wP_adj) && wP_adj >= 0

  if (!success) {
    if (!is.na(adjusted_failure_reason)) {
      set_failure(adjusted_failure_reason, overwrite = TRUE)
    } else if (!H_lambda_inv_info$success || !H_beta_dot_lambda_inv_info$success || !Hpsi_inv_info$success) {
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

  if (is.finite(wP_adj) && wP_adj < 0 && abs(wP_adj) <= 1e-8) {
    wP_adj <- 0
  }

  if (!adjusted_available || !is.finite(wP_adj) || wP_adj < 0) {
    beta_diff <- as.matrix(L %*% beta_hat - b)
    wald_cov_raw <- if (!is.null(unconstrained$robust_cov) &&
      all(dim(unconstrained$robust_cov) == c(nb, nb))) {
      L %*% unconstrained$robust_cov %*% t(L)
    } else {
      L %*% G_beta %*% t(L)
    }
    wald_cov_raw <- 0.5 * (wald_cov_raw + t(wald_cov_raw))
    wald_inv_info <- invert_with_ridge(
      wald_cov_raw,
      cond_limit = cond_limit,
      base_factor = base_factor
    )
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
    Hpsi = Hpsi_info,
    Gpsi = Gpsi_info,
    Gbeta = G_beta,
    Jpsi = Spsi_raw,
    Vpsi_model = Hpsi_raw,
    Vpsi_godambe = Gpsi_raw,
    method = method,
    ridge = ridge_info,
    fallback = fallback,
    success = success,
    adjusted_available = adjusted_available,
    adjusted_failure_reason = adjusted_failure_reason,
    gate = diagnostics,
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

  if (is.null(ridge_factor) || length(ridge_factor) == 0 || !is.finite(ridge_factor)) {
    ridge_factor <- 0
  }
  ridge_factor <- max(ridge_factor, 0)
  base_factor <- if (ridge_factor > 0) ridge_factor else 0

  H_lambda_inv_info <- invert_general_with_ridge(H_lambda, base_factor = base_factor)
  H_beta_dot_lambda <- matrix(NA_real_, nrow = nb, ncol = nb)
  H_beta_dot_lambda_inv_info <- list(success = FALSE, ridge = NA_real_, inverse = NULL)
  score_eff <- matrix(NA_real_, nrow = nb, ncol = 1)
  Upsi <- matrix(NA_real_, nrow = nrow(L), ncol = 1)
  Hpsi <- matrix(NA_real_, nrow = nrow(L), ncol = nrow(L))
  Gpsi <- matrix(NA_real_, nrow = nrow(L), ncol = nrow(L))
  S_eff <- matrix(NA_real_, nrow = nb, ncol = nb)

  score_beta <- matrix(score_full[idx_beta], ncol = 1)
  score_lambda <- matrix(score_full[idx_lambda], ncol = 1)

  if (H_lambda_inv_info$success) {
    H_lambda_inv <- H_lambda_inv_info$inverse
    H_beta_dot_lambda <- H_beta - H_beta_lambda %*% H_lambda_inv %*% t(H_beta_lambda)
    H_beta_dot_lambda_inv_info <- invert_general_with_ridge(H_beta_dot_lambda, base_factor = base_factor)

    score_eff <- score_beta - H_beta_lambda %*% H_lambda_inv %*% score_lambda

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

    if (H_beta_dot_lambda_inv_info$success) {
      H_beta_dot_lambda_inv <- H_beta_dot_lambda_inv_info$inverse
      Upsi <- L %*% H_beta_dot_lambda_inv %*% score_eff
      Hpsi <- L %*% H_beta_dot_lambda_inv %*% t(L)
      Hpsi <- 0.5 * (Hpsi + t(Hpsi))

      G_beta <- H_beta_dot_lambda_inv %*% S_eff %*% H_beta_dot_lambda_inv
      Gpsi <- L %*% G_beta %*% t(L)
      Gpsi <- 0.5 * (Gpsi + t(Gpsi))
    }
  } else {
    S_full <- godambe$S
    if (is.null(S_full)) {
      stop("Constrained fit lacks Godambe meat matrix S.")
    }
  }

  stat <- NA_real_
  pval <- NA_real_
  method <- "adj-score"
  Gpsi_inv_info <- invert_with_ridge(Gpsi, base_factor = base_factor)

  if (Gpsi_inv_info$success) {
    Gpsi_inv <- Gpsi_inv_info$inverse
    stat <- as.numeric(t(Upsi) %*% Gpsi_inv %*% Upsi)
    if (!is.na(stat) && stat < 0 && abs(stat) <= 1e-8) {
      stat <- 0
    }
    pval <- pchisq(max(stat, 0), df = nrow(L), lower.tail = FALSE)
  }

  scalar_check <- NA_real_
  if (check_scalar_identity && nrow(L) == 1 && Gpsi_inv_info$success) {
    scalar_check <- 0
  }

  score_success <- isTRUE(Gpsi_inv_info$success) && is.finite(stat) && is.finite(pval)
  score_failure <- NA_character_
  if (!isTRUE(H_lambda_inv_info$success) || !isTRUE(H_beta_dot_lambda_inv_info$success)) {
    score_failure <- normalise_failure("undefined_stat", "score")
  } else if (!isTRUE(Gpsi_inv_info$success)) {
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
