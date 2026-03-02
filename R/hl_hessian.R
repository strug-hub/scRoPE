hl_profile_hessian_from_stats <- function(per_stats) {
  if (is.null(per_stats$h_beta_beta) || is.null(per_stats$h_beta_phi) ||
      is.null(per_stats$h_phi_phi) || is.null(per_stats$h_tau_tau)) {
    stop("per_subject_stats is missing second-order blocks required for profile Hessian.")
  }

  nb <- nrow(per_stats$h_beta)
  k <- length(per_stats$h_etaeta)

  H_bb_total <- matrix(0, nb, nb)
  H_bphi_total <- numeric(nb)
  H_btau_total <- numeric(nb)
  H_phiphi_total <- 0
  H_phitau_total <- 0
  H_tautau_total <- 0

  per_subject <- vector("list", k)

  for (i in seq_len(k)) {
    a <- per_stats$h_etaeta[i]
    if (!is.finite(a) || abs(a) < 1e-12) {
      a <- sign(a) * 1e-12
    }

    d_beta <- per_stats$h_eta_beta[, i]
    d_phi <- per_stats$h_eta_phi[i]
    d_tau <- per_stats$h_eta_tau[i]

    h_beta_beta_i <- per_stats$h_beta_beta[[i]]
    h_beta_phi_i <- per_stats$h_beta_phi[, i]
    h_phi_phi_i <- per_stats$h_phi_phi[i]
    h_tau_tau_i <- per_stats$h_tau_tau[i]

    H_bb_i <- -(h_beta_beta_i - (tcrossprod(d_beta) / a))
    H_bphi_i <- -(h_beta_phi_i - (d_beta * d_phi / a))
    H_btau_i <- d_beta * (d_tau / a)
    H_phiphi_i <- -(h_phi_phi_i - (d_phi * d_phi) / a)
    H_phitau_i <- d_phi * (d_tau / a)
    H_tautau_i <- -(h_tau_tau_i - (d_tau * d_tau) / a)

    per_subject[[i]] <- list(
      H_bb = H_bb_i,
      H_bphi = H_bphi_i,
      H_btau = H_btau_i,
      H_phiphi = H_phiphi_i,
      H_phitau = H_phitau_i,
      H_tautau = H_tautau_i
    )

    H_bb_total <- H_bb_total + H_bb_i
    H_bphi_total <- H_bphi_total + H_bphi_i
    H_btau_total <- H_btau_total + H_btau_i
    H_phiphi_total <- H_phiphi_total + H_phiphi_i
    H_phitau_total <- H_phitau_total + H_phitau_i
    H_tautau_total <- H_tautau_total + H_tautau_i
  }

  list(
    H_bb = H_bb_total,
    H_bphi = matrix(H_bphi_total, ncol = 1),
    H_btau = matrix(H_btau_total, ncol = 1),
    H_phiphi = matrix(H_phiphi_total, nrow = 1, ncol = 1),
    H_phitau = matrix(H_phitau_total, nrow = 1, ncol = 1),
    H_tautau = matrix(H_tautau_total, nrow = 1, ncol = 1),
    per_subject = per_subject
  )
}

hl_hessian_from_stats <- function(per_stats) {
  profile <- hl_profile_hessian_from_stats(per_stats)
  laplace <- hl_laplace_hessian_from_stats(per_stats)

  list(
    H_bb = profile$H_bb + laplace$H_bb,
    H_bphi = profile$H_bphi + laplace$H_bphi,
    H_btau = profile$H_btau + laplace$H_btau,
    H_phiphi = profile$H_phiphi + laplace$H_phiphi,
    H_phitau = profile$H_phitau + laplace$H_phitau,
    H_tautau = profile$H_tautau + laplace$H_tautau,
    profile = profile,
    laplace = laplace
  )
}

hl_hessian_phi_tau_to_sigma2_phi <- function(H_obs, sigma2) {
  if (sigma2 <= 0) {
    stop("sigma2 must be positive for the Hessian transformation.")
  }

  nb <- nrow(H_obs$H_bb)
  H_old <- rbind(
    cbind(H_obs$H_bb, H_obs$H_bphi, H_obs$H_btau),
    cbind(t(H_obs$H_bphi), H_obs$H_phiphi, H_obs$H_phitau),
    cbind(t(H_obs$H_btau), t(H_obs$H_phitau), H_obs$H_tautau)
  )

  J <- diag(nb + 2)
  J[nb + 1, nb + 1] <- 0
  J[nb + 1, nb + 2] <- 1
  J[nb + 2, nb + 1] <- 1 / sigma2
  J[nb + 2, nb + 2] <- 0

  H_new <- t(J) %*% H_old %*% J

  list(
    matrix = H_new,
    blocks = list(
      H_bb = H_new[1:nb, 1:nb, drop = FALSE],
      H_bsig = matrix(H_new[1:nb, nb + 1], ncol = 1),
      H_bphi = matrix(H_new[1:nb, nb + 2], ncol = 1),
      H_sigsig = matrix(H_new[nb + 1, nb + 1], nrow = 1, ncol = 1),
      H_sigphi = matrix(H_new[nb + 1, nb + 2], nrow = 1, ncol = 1),
      H_phiphi = matrix(H_new[nb + 2, nb + 2], nrow = 1, ncol = 1)
    )
  )
}

hl_laplace_hessian_from_stats <- function(per_stats) {
  required <- c(
    "h_etaeta", "h_etaetaeta", "h_etaetaetaeta",
    "h_etaeta_beta", "h_etaeta_phi", "h_etaeta_tau",
    "h_eta_beta", "h_eta_phi", "h_eta_tau",
    "h_etaeta_beta_beta", "g_beta_beta",
    "h_etaeta_beta_phi", "f_beta_phi", "g_beta_phi",
    "h_etaeta_phi_phi", "g_phi_phi",
    "b_phi", "b_beta", "g_tau_tau"
  )
  missing <- setdiff(required, names(per_stats))
  if (length(missing) > 0) {
    stop("per_subject_stats missing required components for Laplace Hessian: ",
      paste(missing, collapse = ", "))
  }

  nb <- nrow(per_stats$h_beta)
  k <- length(per_stats$h_etaeta)

  H_bb_total <- matrix(0, nb, nb)
  H_bphi_total <- numeric(nb)
  H_btau_total <- numeric(nb)
  H_phiphi_total <- 0
  H_phitau_total <- 0
  H_tautau_total <- 0

  per_subject <- vector("list", k)

  for (i in seq_len(k)) {
    a <- per_stats$h_etaeta[i]
    b <- per_stats$h_etaetaeta[i]
    c <- per_stats$h_etaetaetaeta[i]

    if (!is.finite(a) || abs(a) < 1e-12) {
      a <- sign(a) * 1e-12
    }

    d_beta <- per_stats$h_eta_beta[, i, drop = TRUE]
    d_phi <- per_stats$h_eta_phi[i]
    d_tau <- per_stats$h_eta_tau[i]

    e_beta <- per_stats$h_etaeta_beta[, i, drop = TRUE]
    e_phi <- per_stats$h_etaeta_phi[i]
    e_tau <- per_stats$h_etaeta_tau[i]

    b_beta_vec <- per_stats$b_beta[, i, drop = TRUE]
    b_phi_val <- per_stats$b_phi[i]
    b_tau_val <- e_tau

    g_bb <- per_stats$g_beta_beta[[i]]
    h_etaeta_bb <- per_stats$h_etaeta_beta_beta[[i]]
    h_etaeta_bphi <- per_stats$h_etaeta_beta_phi[[i]]
    f_bphi <- per_stats$f_beta_phi[[i]]
    g_bphi <- per_stats$g_beta_phi[[i]]
    h_etaeta_phiphi <- per_stats$h_etaeta_phi_phi[i]
    g_phiphi <- per_stats$g_phi_phi[i]

    g_btau <- matrix(0, nb, 1)
    h_etaeta_btau <- matrix(0, nb, 1)
    f_btau <- rep(0, nb)

    g_phitau <- 0
    h_etaeta_phitau <- 0
    f_phitau <- 0

    g_tautau <- per_stats$g_tau_tau[i]
    # In the current model, h_{eta,tau,tau} comes from the random-effect term.
    # If not explicitly provided, use the closed form f_{tau,tau} = -d_tau.
    f_tautau <- if (!is.null(per_stats$h_eta_tau_tau)) {
      per_stats$h_eta_tau_tau[i]
    } else {
      -d_tau
    }

    term_bb <- g_bb / a -
      (tcrossprod(e_beta) + b * h_etaeta_bb +
        tcrossprod(b_beta_vec, d_beta) + tcrossprod(d_beta, b_beta_vec)) / (a * a) +
      (2 * b * tcrossprod(e_beta, d_beta) + 2 * b * tcrossprod(d_beta, e_beta) + c * tcrossprod(d_beta)) / (a^3) -
      (2 * b * b * tcrossprod(d_beta)) / (a^4)

    term_bphi <- (g_bphi / a) -
      ((e_beta * e_phi) + b * f_bphi + b_beta_vec * d_phi + b_phi_val * d_beta) / (a * a) +
      (2 * b * (e_beta * d_phi + d_beta * e_phi) + c * d_beta * d_phi) / (a^3) -
      (2 * b * b * d_beta * d_phi) / (a^4)

    term_btau <- (g_btau / a) -
      ((e_beta * e_tau) + b * f_btau + b_beta_vec * d_tau + b_tau_val * d_beta) / (a * a) +
      (2 * b * (e_beta * d_tau + d_beta * e_tau) + c * d_beta * d_tau) / (a^3) -
      (2 * b * b * d_beta * d_tau) / (a^4)

    term_phiphi <- (g_phiphi / a) -
      ((e_phi * e_phi) + b * h_etaeta_phiphi + 2 * b_phi_val * d_phi) / (a * a) +
      (4 * b * e_phi * d_phi + c * d_phi * d_phi) / (a^3) -
      (2 * b * b * d_phi * d_phi) / (a^4)

    term_phitau <- (g_phitau / a) -
      ((e_phi * e_tau) + b_phi_val * d_tau + b_tau_val * d_phi) / (a * a) +
      (2 * b * (e_phi * d_tau + d_phi * e_tau) + c * d_phi * d_tau) / (a^3) -
      (2 * b * b * d_phi * d_tau) / (a^4)

    term_tautau <- (g_tautau / a) -
      ((e_tau * e_tau) + b * f_tautau + 2 * b_tau_val * d_tau) / (a * a) +
      (4 * b * e_tau * d_tau + c * d_tau * d_tau) / (a^3) -
      (2 * b * b * d_tau * d_tau) / (a^4)

    per_subject[[i]] <- list(
      H_bb = term_bb,
      H_bphi = matrix(term_bphi, ncol = 1),
      H_btau = matrix(term_btau, ncol = 1),
      H_phiphi = matrix(term_phiphi, nrow = 1, ncol = 1),
      H_phitau = matrix(term_phitau, nrow = 1, ncol = 1),
      H_tautau = matrix(term_tautau, nrow = 1, ncol = 1)
    )

    H_bb_total <- H_bb_total + term_bb
    H_bphi_total <- H_bphi_total + term_bphi
    H_btau_total <- H_btau_total + term_btau
    H_phiphi_total <- H_phiphi_total + term_phiphi
    H_phitau_total <- H_phitau_total + term_phitau
    H_tautau_total <- H_tautau_total + term_tautau
  }

  list(
    H_bb = H_bb_total,
    H_bphi = matrix(H_bphi_total, ncol = 1),
    H_btau = matrix(H_btau_total, ncol = 1),
    H_phiphi = matrix(H_phiphi_total, nrow = 1, ncol = 1),
    H_phitau = matrix(H_phitau_total, nrow = 1, ncol = 1),
    H_tautau = matrix(H_tautau_total, nrow = 1, ncol = 1),
    per_subject = per_subject
  )
}
