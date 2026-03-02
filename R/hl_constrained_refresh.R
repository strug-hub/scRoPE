hl_refresh_constrained_state <- function(res,
                                         ctx,
                                         posv,
                                         gene_index) {
  nb <- ctx$nb
  sigma2 <- res$dispersion[1]
  phi <- res$dispersion[2]
  if (sigma2 <= 0 || phi <= 0) {
    stop("dispersion parameters must be positive to refresh constrained state.")
  }

  theta_full <- c(res$beta, log(sigma2), log(phi))

  ptmg_out <- ptmg_ll_der_hes4(
    theta_full,
    ctx$pred,
    ctx$offset,
    posv$Y,
    posv$n_onetwo,
    posv$ytwo,
    ctx$fid,
    ctx$cumsumy[gene_index, ],
    ctx$posind[[gene_index]],
    posv$posindy,
    nb,
    ctx$nind,
    ctx$k
  )

  chain <- apply_chain_rule_transform(
    U = ptmg_out$gradient,
    H = ptmg_out$hessian,
    S = NULL,
    intercept_idx = ctx$intcol,
    p = nb,
    U_per_subject = ptmg_out$per_subject_gradients,
    sigma2 = sigma2
  )

  score_full <- chain$U
  per_subject_full <- chain$U_per_subject
  if (is.null(per_subject_full)) {
    per_subject_full <- ptmg_out$per_subject_gradients
  }
  meat_full <- if (!is.null(chain$S)) {
    chain$S
  } else if (!is.null(per_subject_full)) {
    tcrossprod(per_subject_full)
  } else {
    matrix(0, nb + 2, nb + 2)
  }

  idx_beta <- seq_len(nb)
  # ptmg_ll_der_hes4 parameter order is: beta, log(sigma2), log(phi)
  idx_tau <- nb + 1L
  idx_log_phi <- nb + 2L

  score_nat <- c(
    score_full[idx_beta],
    score_full[idx_log_phi] / phi,
    score_full[idx_tau]
  )

  if (!is.null(per_subject_full)) {
    per_subject_nat <- rbind(
      per_subject_full[idx_beta, , drop = FALSE],
      per_subject_full[idx_log_phi, , drop = FALSE] / phi,
      per_subject_full[idx_tau, , drop = FALSE]
    )
  } else {
    per_subject_nat <- matrix(0, nrow = nb + 2, ncol = ctx$k)
  }

  meat_full[idx_log_phi, ] <- meat_full[idx_log_phi, ] / phi
  meat_full[, idx_log_phi] <- meat_full[, idx_log_phi] / phi
  reorder_idx <- c(idx_beta, idx_log_phi, idx_tau)
  meat_nat <- meat_full[reorder_idx, reorder_idx, drop = FALSE]
  meat_nat <- 0.5 * (meat_nat + t(meat_nat))

  hessian_full <- chain$H
  hessian_full[idx_log_phi, ] <- hessian_full[idx_log_phi, ] / phi
  hessian_full[, idx_log_phi] <- hessian_full[, idx_log_phi] / phi
  # Nonlinear term for y = log(phi): d2/dphi2 = (d2/dy2 - d/dy) / phi^2.
  hessian_full[idx_log_phi, idx_log_phi] <- hessian_full[idx_log_phi, idx_log_phi] -
    score_full[idx_log_phi] / (phi * phi)
  hessian_nat <- hessian_full[reorder_idx, reorder_idx, drop = FALSE]
  hessian_nat <- 0.5 * (hessian_nat + t(hessian_nat))

  H_bb <- hessian_nat[idx_beta, idx_beta, drop = FALSE]
  H_bphi <- matrix(hessian_nat[idx_beta, nb + 1L], ncol = 1)
  H_btau <- matrix(hessian_nat[idx_beta, nb + 2L], ncol = 1)
  H_phiphi <- matrix(hessian_nat[nb + 1L, nb + 1L], ncol = 1)
  H_phitau <- matrix(hessian_nat[nb + 1L, nb + 2L], ncol = 1)
  H_tautau <- matrix(hessian_nat[nb + 2L, nb + 2L], ncol = 1)

  H_bl <- cbind(H_bphi, H_btau)
  H_ll <- rbind(
    c(H_phiphi, H_phitau),
    c(H_phitau, H_tautau)
  )

  eff_con <- hl_efficient_components(
    H_bb,
    H_bl,
    H_ll,
    score_nat,
    nb
  )

  res$hl_state$score <- score_nat
  res$hl_state$score_subject <- per_subject_nat
  res$hl_state$meat <- meat_nat
  res$hl_state$H_blocks <- list(
    H_bb = H_bb,
    H_bl = H_bl,
    H_ll = H_ll
  )
  res$hl_state$H_schur <- eff_con$H_beta_beta_schur
  res$hl_state$efficient_score <- eff_con$score_beta_eff
  res$hl_state$H_obs <- list(
    H_bb = H_bb,
    H_bphi = H_bphi,
    H_btau = H_btau,
    H_phiphi = H_phiphi,
    H_phitau = H_phitau,
    H_tautau = H_tautau
  )

  res
}
