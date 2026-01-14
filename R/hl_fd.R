hl_score_at_theta <- function(theta, ctx, gene_index, count, posv = NULL) {
  nb <- ctx$nb
  if (length(theta) != nb + 2) {
    stop("theta must have length nb + 2 (beta coefficients, phi, log sigma^2).")
  }

  if (is.null(posv)) {
    posv <- call_posindy(count, gene_index - 1L, ctx$nind)
  }

  beta <- theta[seq_len(nb)]
  phi <- theta[nb + 1L]
  tau <- theta[nb + 2L]
  sigma2 <- exp(tau)
  if (phi <= 0) {
    stop("phi must be positive.")
  }
  para <- c(beta, tau, log(phi))

  out <- ptmg_ll_der_hes4(
    para,
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

  grad <- out$gradient
  grad_beta <- grad[seq_len(nb)]
  grad_tau <- grad[nb + 1L]
  grad_log_phi <- grad[nb + 2L]
  grad_phi <- grad_log_phi / phi

  c(grad_beta, grad_phi, grad_tau)
}

hl_fd_hessian <- function(theta, ctx, gene_index, count, posv = NULL, step = 1e-4) {
  nb <- ctx$nb
  p <- nb + 2L
  if (length(theta) != p) {
    stop("theta must have length nb + 2 (beta coefficients, phi, log sigma^2).")
  }
  if (!is.numeric(step) || length(step) != 1L || step <= 0) {
    stop("step must be a positive numeric scalar.")
  }

  if (is.null(posv)) {
    posv <- call_posindy(count, gene_index - 1L, ctx$nind)
  }

  hessian <- matrix(0, p, p)

  for (j in seq_len(p)) {
    scale_j <- max(1, abs(theta[j]))
    step_j <- step * scale_j

    ei <- rep(0, p)
    ei[j] <- step_j

    s_plus <- hl_score_at_theta(theta + ei, ctx, gene_index, count, posv)
    s_minus <- hl_score_at_theta(theta - ei, ctx, gene_index, count, posv)

    hessian[, j] <- (s_plus - s_minus) / (2 * step_j)
  }

  hessian <- 0.5 * (hessian + t(hessian))

  list(hessian = hessian)
}
