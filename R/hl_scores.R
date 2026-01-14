hl_build_counts <- function(posv, nind) {
  y <- numeric(nind)
  if (!is.null(posv$posindy) && length(posv$posindy) > 0) {
    idx <- as.integer(posv$posindy) + 1L
    y[idx] <- posv$Y
  }
  y
}

hl_score_from_stats <- function(per_stats, intcol, nb, k) {
  has_all <- all(c("h_beta", "h_eta_beta", "h_etaeta_beta",
                   "h_phi", "h_eta_phi", "h_etaeta_phi",
                   "h_tau", "h_eta_tau", "h_etaeta_tau",
                   "h_etaeta", "h_etaetaeta") %in% names(per_stats))
  if (!has_all) {
    stop("per_subject_stats missing required components.")
  }

  h_beta <- per_stats$h_beta
  h_eta_beta <- per_stats$h_eta_beta
  h_etaeta_beta <- per_stats$h_etaeta_beta
  h_phi <- per_stats$h_phi
  h_eta_phi <- per_stats$h_eta_phi
  h_etaeta_phi <- per_stats$h_etaeta_phi
  h_tau <- per_stats$h_tau
  h_eta_tau <- per_stats$h_eta_tau
  h_etaeta_tau <- per_stats$h_etaeta_tau
  h_etaeta <- per_stats$h_etaeta
  h_etaetaeta <- per_stats$h_etaetaeta

  per_subject <- matrix(0, nrow = nb + 2, ncol = k)

  for (i in seq_len(k)) {
    denom_i <- h_etaeta[i]
    if (!is.finite(denom_i) || abs(denom_i) < 1e-10) {
      denom_i <- sign(denom_i) * 1e-10
    }

    term_beta <- h_etaeta_beta[, i] / denom_i -
      (h_etaetaeta[i] * h_eta_beta[, i]) / (denom_i * denom_i)
    U_beta_i <- h_beta[, i] - 0.5 * term_beta

    term_phi <- h_etaeta_phi[i] / denom_i -
      (h_etaetaeta[i] * h_eta_phi[i]) / (denom_i * denom_i)
    U_phi_i <- h_phi[i] - 0.5 * term_phi

    term_tau <- h_etaeta_tau[i] / denom_i -
      (h_etaetaeta[i] * h_eta_tau[i]) / (denom_i * denom_i)
    U_tau_i <- h_tau[i] - 0.5 * term_tau

    per_subject[1:nb, i] <- U_beta_i
    per_subject[nb + 1, i] <- U_phi_i
    per_subject[nb + 2, i] <- U_tau_i
  }

  U_total <- rowSums(per_subject)
  list(U = U_total, per_subject = per_subject)
}

hl_meat_from_scores <- function(score_subject) {
  if (is.null(score_subject)) {
    return(NULL)
  }
  if (!is.matrix(score_subject)) {
    stop("score_subject must be a matrix of per-subject scores.")
  }
  base::tcrossprod(score_subject)
}

hl_hessian_blocks_from_stats <- function(per_stats) {
  if (!("h_beta" %in% names(per_stats)) && "hl_state" %in% names(per_stats)) {
    per_stats <- per_stats$hl_state$per_subject_stats
  }
  required <- c("h_beta", "h_eta_beta", "h_etaeta_beta",
                "h_phi", "h_eta_phi", "h_etaeta_phi",
                "h_tau", "h_eta_tau", "h_etaeta_tau",
                "h_etaeta", "h_etaetaeta")
  if (!all(required %in% names(per_stats))) {
    stop("per_subject_stats missing required components for Hessian blocks.")
  }

  h_beta <- per_stats$h_beta
  h_eta_beta <- per_stats$h_eta_beta
  h_etaeta_beta <- per_stats$h_etaeta_beta
  h_phi <- per_stats$h_phi
  h_eta_phi <- per_stats$h_eta_phi
  h_etaeta_phi <- per_stats$h_etaeta_phi
  h_tau <- per_stats$h_tau
  h_eta_tau <- per_stats$h_eta_tau
  h_etaeta_tau <- per_stats$h_etaeta_tau
  h_etaeta <- per_stats$h_etaeta
  h_etaetaeta <- per_stats$h_etaetaeta

  nb <- nrow(h_beta)
  k <- ncol(h_beta)

  H_bb <- matrix(0, nb, nb)
  H_bl <- matrix(0, nb, 2)
  H_ll <- matrix(0, 2, 2)

  for (i in seq_len(k)) {
    denom <- h_etaeta[i]
    if (!is.finite(denom) || abs(denom) < 1e-10) {
      denom <- sign(denom) * 1e-10
    }
    inv_denom <- 1 / denom
    inv_denom2 <- inv_denom * inv_denom

    H_bb <- H_bb - (h_etaeta_beta[, i, drop = FALSE] %*% t(h_eta_beta[, i, drop = FALSE])) * inv_denom2
    H_bb <- H_bb - h_etaeta_beta[, i, drop = FALSE] %*% t(h_etaeta_beta[, i, drop = FALSE]) * inv_denom
    H_bb <- H_bb - h_etaetaeta[i] * (h_eta_beta[, i, drop = FALSE] %*% t(h_eta_beta[, i, drop = FALSE])) * inv_denom2

    H_bl[, 1] <- H_bl[, 1] - (h_etaeta_beta[, i] * h_eta_phi[i]) * inv_denom2
    H_bl[, 1] <- H_bl[, 1] - h_etaeta_beta[, i] * h_etaeta_phi[i] * inv_denom
    H_bl[, 1] <- H_bl[, 1] - h_etaetaeta[i] * h_eta_beta[, i] * h_eta_phi[i] * inv_denom2

    H_bl[, 2] <- H_bl[, 2] - (h_etaeta_beta[, i] * h_eta_tau[i]) * inv_denom2
    H_bl[, 2] <- H_bl[, 2] - h_etaeta_beta[, i] * h_etaeta_tau[i] * inv_denom
    H_bl[, 2] <- H_bl[, 2] - h_etaetaeta[i] * h_eta_beta[, i] * h_eta_tau[i] * inv_denom2

    H_ll[1, 1] <- H_ll[1, 1] - (h_etaeta_phi[i] * h_eta_phi[i]) * inv_denom2 - h_etaeta_phi[i] * inv_denom - h_etaetaeta[i] * h_eta_phi[i]^2 * inv_denom2
    H_ll[2, 2] <- H_ll[2, 2] - (h_etaeta_tau[i] * h_eta_tau[i]) * inv_denom2 - h_etaeta_tau[i] * inv_denom - h_etaetaeta[i] * h_eta_tau[i]^2 * inv_denom2
    cross_phi_tau <- - (h_etaeta_phi[i] * h_eta_tau[i]) * inv_denom2 - h_etaeta_tau[i] * h_eta_phi[i] * inv_denom2 - h_etaetaeta[i] * h_eta_phi[i] * h_eta_tau[i] * inv_denom2
    H_ll[1, 2] <- H_ll[1, 2] + cross_phi_tau
    H_ll[2, 1] <- H_ll[2, 1] + cross_phi_tau
  }

  list(H_bb = H_bb, H_bl = H_bl, H_ll = H_ll)
}
hl_efficient_components <- function(H_bb, H_bl, H_ll, score, nb) {
  if (length(score) < nb + 2) {
    stop("Score length is incompatible with nb + dispersion components.")
  }
  if (nrow(H_bl) != nb || ncol(H_bl) != 2) {
    stop("H_bl must have dimensions nb x 2.")
  }
  if (nrow(H_ll) != 2 || ncol(H_ll) != 2) {
    stop("H_ll must be 2 x 2 (for log sigma^2 and phi).")
  }

  H_ll_inv <- tryCatch(
    solve(H_ll),
    error = function(e) {
      solve(H_ll + diag(1e-8, 2))
    }
  )

  H_beta_lambda <- H_bl
  H_lambda_beta <- t(H_bl)

  schur <- H_bb - H_beta_lambda %*% H_ll_inv %*% H_lambda_beta

  U_beta <- score[seq_len(nb)]
  U_lambda <- score[nb + seq_len(2)]
  U_eff <- U_beta - H_beta_lambda %*% H_ll_inv %*% U_lambda

  list(
    H_beta_beta = H_bb,
    H_beta_lambda = H_bl,
    H_lambda_lambda = H_ll,
    H_beta_beta_schur = schur,
    score_beta_eff = as.numeric(U_eff)
  )
}
