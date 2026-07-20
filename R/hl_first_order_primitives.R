# Estimator-aligned primitives for the smooth first-order nested HL route.
#
# All scientific derivatives in this file use the canonical order
# (beta, tau = log(sigma2), xi = log(NB size)).  They are deliberately
# independent of the historical PTMG derivative cache.

hl_fo_logspace_add <- function(a, b) {
  z <- pmax(a, b)
  z + log1p(exp(-abs(a - b)))
}

hl_fo_softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

hl_fo_lgamma_increment <- function(phi, y) {
  if (phi <= 1e5) {
    return(lgamma(phi + y) - lgamma(phi))
  }
  vapply(y, function(count) {
    if (count == 0) return(0)
    if (count == floor(count) && count <= 1e5) {
      return(count * log(phi) +
        sum(log1p(seq.int(0, count - 1) / phi)))
    }
    lgamma(phi + count) - lgamma(phi)
  }, numeric(1))
}

hl_fo_r_minus_log1p_odds <- function(log_odds, r) {
  out <- r - hl_fo_softplus(log_odds)
  small <- log_odds < log(1e-3)
  if (any(small)) {
    t <- exp(log_odds[small])
    total <- numeric(length(t))
    power <- t
    for (k in 2:14) {
      power <- power * t
      total <- total + (-1)^(k + 1L) * (k - 1) / k * power
    }
    out[small] <- total
  }
  out
}

hl_fo_asymmetry <- function(x) {
  if (!length(x)) {
    return(0)
  }
  sqrt(sum((x - t(x))^2)) / max(1, sqrt(sum(x^2)))
}

hl_fo_d_minus_expm1 <- function(d) {
  if (abs(d) >= 1e-3) {
    return(d - expm1(d))
  }
  # -sum_{k=2}^infinity d^k/k!, evaluated without subtracting the common
  # first-order term from expm1(d).
  term <- d * d / 2
  total <- term
  for (k in 3:12) {
    term <- term * d / k
    total <- total + term
  }
  -total
}

hl_fo_gamma_normalizer_derivatives <- function(alpha, log_alpha) {
  if (alpha <= 1e5) {
    return(list(
      value = alpha * log_alpha - lgamma(alpha) - alpha,
      first = log_alpha - digamma(alpha),
      second = 1 / alpha - trigamma(alpha),
      third = -1 / alpha^2 - psigamma(alpha, deriv = 2L)
    ))
  }
  inv <- 1 / alpha
  inv2 <- inv * inv
  remainder <- inv / 12 - inv^3 / 360 + inv^5 / 1260 -
    inv^7 / 1680 + inv^9 / 1188
  list(
    value = 0.5 * (log_alpha - log(2 * pi)) - remainder,
    first = inv / 2 + inv2 / 12 - inv^4 / 120 +
      inv^6 / 252 - inv^8 / 240,
    second = -inv2 / 2 - inv^3 / 6 + inv^5 / 30 -
      inv^7 / 42 + inv^9 / 30,
    third = inv^3 + inv^4 / 2 - inv^6 / 6 +
      inv^8 / 6 - 3 * inv^10 / 10
  )
}

hl_fo_nb_primitives <- function(y, log_m, xi) {
  y <- as.numeric(y)
  log_m <- as.numeric(log_m)
  phi <- exp(xi)
  if (!is.finite(phi) || phi <= 0 || length(y) != length(log_m)) {
    return(list(ok = FALSE, reason = "nonfinite_nb_domain"))
  }

  r <- stats::plogis(log_m - xi)
  one_minus_r <- stats::plogis(xi - log_m)
  log_d <- hl_fo_logspace_add(log_m, xi)
  inv_d <- exp(-log_d)

  log_odds <- log_m - xi
  loglik <- hl_fo_lgamma_increment(phi, y) - lgamma(y + 1) -
    phi * hl_fo_softplus(log_odds) -
    y * hl_fo_softplus(-log_odds)

  s <- y * one_minus_r - phi * r
  t <- -(y + phi) * r * one_minus_r
  u <- -(y + phi) * r * one_minus_r * (1 - 2 * r)
  w <- -(y + phi) * r * one_minus_r * (1 - 6 * r + 6 * r^2)

  l_phi <- digamma(y + phi) - digamma(phi) + xi + 1 - log_d -
    (y + phi) * inv_d
  l_phiphi <- trigamma(y + phi) - trigamma(phi) + 1 / phi - inv_d +
    y * inv_d^2 - r * inv_d
  l_phiphiphi <- psigamma(y + phi, deriv = 2L) -
    psigamma(phi, deriv = 2L) - 1 / phi^2 + inv_d^2 -
    2 * (y - exp(log_m)) * inv_d^3

  s_phi <- r * y * one_minus_r / phi - r^2
  s_phiphi <- -2 * s_phi * inv_d
  t_phi <- -r * one_minus_r +
    (y + phi) * r * one_minus_r * (1 - 2 * r) / phi
  u_phi <- -r * one_minus_r * (1 - 2 * r) +
    (y + phi) * r * one_minus_r * (1 - 6 * r + 6 * r^2) / phi
  t_phiphi <- -r * (
    2 * inv_d - 4 * (y + 2 * phi) * inv_d^2 +
      6 * phi * (y + phi) * inv_d^3
  )

  stable_large_size <- phi > 1e5 && all(y == floor(y)) && all(y <= 1e5)
  if (stable_large_size) {
    q_sums <- vapply(y, function(count) {
      if (count == 0) return(0)
      sum(seq.int(0, count - 1) / (phi + seq.int(0, count - 1)))
    }, numeric(1))
    q1_sums <- vapply(y, function(count) {
      if (count == 0) return(0)
      q <- seq.int(0, count - 1) / (phi + seq.int(0, count - 1))
      sum(q * (1 - q))
    }, numeric(1))
    q2_sums <- vapply(y, function(count) {
      if (count == 0) return(0)
      q <- seq.int(0, count - 1) / (phi + seq.int(0, count - 1))
      sum(q * (1 - q) * (1 - 2 * q))
    }, numeric(1))
    r_minus_f <- hl_fo_r_minus_log1p_odds(log_odds, r)
    l_xi_stable <- phi * r_minus_f + y * r - q_sums
    l_xixi_stable <- phi * (r_minus_f + r^2) -
      y * r * one_minus_r + q1_sums
    l_xixixi_stable <- phi * (r_minus_f + 2 * r^3) +
      y * r * one_minus_r * (1 - 2 * r) - q2_sums
    s_xi_stable <- r * s
    s_xixi_stable <- r * s * (2 * r - 1)
    t_xi_stable <- r * one_minus_r * (y - 2 * r * (y + phi))
    t_xixi_stable <- r * one_minus_r * (
      -(1 - 2 * r) * (y - 2 * r * (y + phi)) +
        2 * r * one_minus_r * (y + phi) - 2 * r * phi
    )
    u_xi_stable <- t_xi_stable * (1 - 2 * r) +
      2 * t * r * one_minus_r

    l_phi <- l_xi_stable / phi
    l_phiphi <- (l_xixi_stable - l_xi_stable) / phi^2
    l_phiphiphi <- (
      l_xixixi_stable - l_xi_stable - 3 * phi^2 * l_phiphi
    ) / phi^3
    s_phi <- s_xi_stable / phi
    s_phiphi <- (s_xixi_stable - s_xi_stable) / phi^2
    t_phi <- t_xi_stable / phi
    t_phiphi <- (t_xixi_stable - t_xi_stable) / phi^2
    u_phi <- u_xi_stable / phi
  }

  out <- list(
    ok = TRUE,
    phi = phi,
    r = r,
    loglik = loglik,
    s = s,
    t = t,
    u = u,
    w = w,
    l_phi = l_phi,
    l_phiphi = l_phiphi,
    l_phiphiphi = l_phiphiphi,
    s_phi = s_phi,
    s_phiphi = s_phiphi,
    t_phi = t_phi,
    t_phiphi = t_phiphi,
    u_phi = u_phi
  )
  if (stable_large_size) {
    out$l_xi <- l_xi_stable
    out$l_xixi <- l_xixi_stable
    out$l_xixixi <- l_xixixi_stable
    out$s_xi <- s_xi_stable
    out$s_xixi <- s_xixi_stable
    out$t_xi <- t_xi_stable
    out$t_xixi <- t_xixi_stable
    out$u_xi <- u_xi_stable
  } else {
    out$l_xi <- phi * l_phi
    out$l_xixi <- phi * l_phi + phi^2 * l_phiphi
    out$l_xixixi <- phi * l_phi + 3 * phi^2 * l_phiphi +
      phi^3 * l_phiphiphi
    out$s_xi <- phi * s_phi
    out$s_xixi <- phi * s_phi + phi^2 * s_phiphi
    out$t_xi <- phi * t_phi
    out$t_xixi <- phi * t_phi + phi^2 * t_phiphi
    out$u_xi <- phi * u_phi
  }

  vals <- unlist(out[setdiff(names(out), c("ok", "reason"))], use.names = FALSE)
  if (any(!is.finite(vals))) {
    return(list(ok = FALSE, reason = "nonfinite_nb_arithmetic"))
  }
  out
}

hl_fo_gamma_primitives <- function(eta, tau) {
  v <- exp(tau)
  if (!is.finite(v) || v <= 0 || !is.finite(eta)) {
    return(list(ok = FALSE, reason = "nonfinite_gamma_domain"))
  }
  evm1 <- expm1(v)
  if (!is.finite(evm1) || evm1 <= 0) {
    return(list(ok = FALSE, reason = "nonfinite_gamma_domain"))
  }

  log_alpha <- -log(evm1)
  alpha <- exp(log_alpha)
  log_lambda <- log_alpha - v / 2
  lambda <- exp(log_lambda)
  lambda_z <- exp(log_lambda + eta)
  z <- exp(eta)
  if (any(!is.finite(c(alpha, lambda, lambda_z, z)))) {
    return(list(ok = FALSE, reason = "nonfinite_gamma_arithmetic"))
  }

  dot_alpha <- -v * alpha * (1 + alpha)
  ddot_alpha <- v * alpha * (1 + alpha) * (v * (1 + 2 * alpha) - 1)
  dddot_alpha <- v * alpha * (1 + alpha) * (
    -1 + 3 * v * (1 + 2 * alpha) -
      v^2 * (1 + 6 * alpha + 6 * alpha^2)
  )

  d <- eta - v / 2
  exp_d <- exp(d)
  f <- hl_fo_d_minus_expm1(d)
  f1 <- -expm1(d)
  f2 <- -exp_d
  f3 <- -exp_d
  d1 <- -v / 2
  d2 <- d1
  d3 <- d1
  normalizer <- hl_fo_gamma_normalizer_derivatives(alpha, log_alpha)

  g <- alpha * f + normalizer$value
  g_eta <- alpha * f1
  g_etaeta <- alpha * f2
  g_tau <- dot_alpha * (f + normalizer$first) + alpha * f1 * d1
  g_eta_tau <- dot_alpha * f1 + alpha * f2 * d1
  g_etaeta_tau <- -(dot_alpha * exp_d + alpha * exp_d * d1)
  g_tautau <- normalizer$second * dot_alpha^2 +
    2 * f1 * dot_alpha * d1 + alpha * f2 * d1^2 +
    (f + normalizer$first) * ddot_alpha + alpha * f1 * d2
  g_tautautau <- normalizer$third * dot_alpha^3 +
    3 * f2 * dot_alpha * d1^2 + alpha * f3 * d1^3 +
    3 * (
      normalizer$second * dot_alpha * ddot_alpha +
        f1 * (ddot_alpha * d1 + dot_alpha * d2) +
        alpha * f2 * d1 * d2
    ) + (f + normalizer$first) * dddot_alpha + alpha * f1 * d3
  g_eta_tautau <- ddot_alpha * f1 +
    2 * dot_alpha * f2 * d1 + alpha * (f3 * d1^2 + f2 * d2)
  g_etaeta_tautau <- -exp_d * (
    ddot_alpha + 2 * dot_alpha * d1 + alpha * (d1^2 + d2)
  )

  out <- list(
    ok = TRUE,
    v = v,
    alpha = alpha,
    lambda = lambda,
    log_alpha = log_alpha,
    log_lambda = log_lambda,
    g = g,
    g_eta = g_eta,
    g_etaeta = g_etaeta,
    g_etaetaeta = g_etaeta,
    g_etaetaetaeta = g_etaeta,
    g_tau = g_tau,
    g_eta_tau = g_eta_tau,
    g_etaeta_tau = g_etaeta_tau,
    g_etaetaeta_tau = g_etaeta_tau,
    g_tautau = g_tautau,
    g_tautautau = g_tautautau,
    g_eta_tautau = g_eta_tautau,
    g_etaeta_tautau = g_etaeta_tautau
  )
  vals <- unlist(out[setdiff(names(out), c("ok", "reason"))], use.names = FALSE)
  if (any(!is.finite(vals))) {
    return(list(ok = FALSE, reason = "nonfinite_gamma_arithmetic"))
  }
  out
}

hl_fo_mode_components <- function(y, log_m_base, eta, tau, xi) {
  nb <- hl_fo_nb_primitives(y, log_m_base + eta, xi)
  if (!isTRUE(nb$ok)) {
    return(nb)
  }
  gp <- hl_fo_gamma_primitives(eta, tau)
  if (!isTRUE(gp$ok)) {
    return(gp)
  }
  list(
    ok = TRUE,
    h_eta = sum(nb$s) + gp$g_eta,
    h_etaeta = sum(nb$t) + gp$g_etaeta,
    nb = nb,
    gamma = gp
  )
}

hl_fo_subject_raw <- function(y, X, log_offset, beta, eta, rho_full,
                              include_laplace_constant = TRUE) {
  p <- length(beta)
  d <- p + 2L
  tau_idx <- p + 1L
  xi_idx <- p + 2L
  tau <- rho_full[1]
  xi <- rho_full[2]
  log_m <- as.numeric(log_offset + X %*% beta + eta)

  nb <- hl_fo_nb_primitives(y, log_m, xi)
  if (!isTRUE(nb$ok)) {
    return(nb)
  }
  gp <- hl_fo_gamma_primitives(eta, tau)
  if (!isTRUE(gp$ok)) {
    return(gp)
  }

  h_grad <- numeric(d)
  h_hess <- matrix(0, d, d)
  h_third <- array(0, dim = c(d, d, d))
  d_eta <- numeric(d)
  e_etaeta <- numeric(d)
  k_eta <- matrix(0, d, d)
  f_etaetaeta <- numeric(d)
  m_etaeta <- matrix(0, d, d)

  h_grad[seq_len(p)] <- as.numeric(crossprod(X, nb$s))
  h_grad[tau_idx] <- gp$g_tau
  h_grad[xi_idx] <- sum(nb$l_xi)

  h_hess[seq_len(p), seq_len(p)] <- crossprod(X, X * nb$t)
  h_hess[seq_len(p), xi_idx] <- as.numeric(crossprod(X, nb$s_xi))
  h_hess[xi_idx, seq_len(p)] <- h_hess[seq_len(p), xi_idx]
  h_hess[tau_idx, tau_idx] <- gp$g_tautau
  h_hess[xi_idx, xi_idx] <- sum(nb$l_xixi)
  h_third[tau_idx, tau_idx, tau_idx] <- gp$g_tautautau
  h_third[xi_idx, xi_idx, xi_idx] <- sum(nb$l_xixixi)

  for (a in seq_len(p)) {
    for (b in seq_len(p)) {
      for (cc in seq_len(p)) {
        h_third[a, b, cc] <- sum(X[, a] * X[, b] * X[, cc] * nb$u)
      }
      bb_xi <- sum(X[, a] * X[, b] * nb$t_xi)
      h_third[a, b, xi_idx] <- bb_xi
      h_third[a, xi_idx, b] <- bb_xi
      h_third[xi_idx, a, b] <- bb_xi
    }
    b_xixi <- sum(X[, a] * nb$s_xixi)
    h_third[a, xi_idx, xi_idx] <- b_xixi
    h_third[xi_idx, a, xi_idx] <- b_xixi
    h_third[xi_idx, xi_idx, a] <- b_xixi
  }

  d_eta[seq_len(p)] <- as.numeric(crossprod(X, nb$t))
  d_eta[tau_idx] <- gp$g_eta_tau
  d_eta[xi_idx] <- sum(nb$s_xi)

  e_etaeta[seq_len(p)] <- as.numeric(crossprod(X, nb$u))
  e_etaeta[tau_idx] <- gp$g_etaeta_tau
  e_etaeta[xi_idx] <- sum(nb$t_xi)

  k_eta[seq_len(p), seq_len(p)] <- crossprod(X, X * nb$u)
  k_eta[seq_len(p), xi_idx] <- as.numeric(crossprod(X, nb$t_xi))
  k_eta[xi_idx, seq_len(p)] <- k_eta[seq_len(p), xi_idx]
  k_eta[tau_idx, tau_idx] <- gp$g_eta_tautau
  k_eta[xi_idx, xi_idx] <- sum(nb$s_xixi)

  f_etaetaeta[seq_len(p)] <- as.numeric(crossprod(X, nb$w))
  f_etaetaeta[tau_idx] <- gp$g_etaetaeta_tau
  f_etaetaeta[xi_idx] <- sum(nb$u_xi)

  m_etaeta[seq_len(p), seq_len(p)] <- crossprod(X, X * nb$w)
  m_etaeta[seq_len(p), xi_idx] <- as.numeric(crossprod(X, nb$u_xi))
  m_etaeta[xi_idx, seq_len(p)] <- m_etaeta[seq_len(p), xi_idx]
  m_etaeta[tau_idx, tau_idx] <- gp$g_etaeta_tautau
  m_etaeta[xi_idx, xi_idx] <- sum(nb$t_xixi)

  a_mode <- sum(nb$t) + gp$g_etaeta
  b_mode <- sum(nb$u) + gp$g_etaetaeta
  c_mode <- sum(nb$w) + gp$g_etaetaetaeta
  h_eta <- sum(nb$s) + gp$g_eta
  h_value <- sum(nb$loglik) + gp$g

  raw <- list(
    ok = TRUE,
    h = h_value,
    h_eta = h_eta,
    a = a_mode,
    b = b_mode,
    c = c_mode,
    h_grad = h_grad,
    h_hess = h_hess,
    h_third = h_third,
    d = d_eta,
    e = e_etaeta,
    k = k_eta,
    f = f_etaetaeta,
    M = m_etaeta,
    nb = nb,
    gamma = gp,
    coordinate_order = c(paste0("beta", seq_len(p)), "tau", "xi")
  )
  hl_fo_profile_subject(raw, include_laplace_constant = include_laplace_constant)
}

hl_fo_profile_subject <- function(raw, include_laplace_constant = TRUE) {
  if (!isTRUE(raw$ok)) {
    return(raw)
  }
  a <- raw$a
  if (!is.finite(a) || a >= 0) {
    return(list(ok = FALSE, reason = "nonnegative_mode_curvature", raw = raw))
  }
  d <- raw$d
  e <- raw$e
  b <- raw$b
  cc <- raw$c

  p_grad <- raw$h_grad
  p_hess <- raw$h_hess - tcrossprod(d) / a
  ntheta <- length(d)
  p_third <- raw$h_third
  for (r in seq_len(ntheta)) {
    for (s in seq_len(ntheta)) {
      for (tt in seq_len(ntheta)) {
        p_third[r, s, tt] <- raw$h_third[r, s, tt] -
          (raw$k[r, s] * d[tt] + raw$k[r, tt] * d[s] +
             raw$k[s, tt] * d[r]) / a +
          (e[r] * d[s] * d[tt] + e[s] * d[r] * d[tt] +
             e[tt] * d[r] * d[s]) / a^2 -
          b * d[r] * d[s] * d[tt] / a^3
      }
    }
  }

  abar_grad <- e - b * d / a
  abar_hess <- raw$M -
    (tcrossprod(raw$f, d) + tcrossprod(d, raw$f) + b * raw$k) / a +
    (cc * tcrossprod(d) + b * (tcrossprod(e, d) + tcrossprod(d, e))) / a^2 -
    b^2 * tcrossprod(d) / a^3

  ell_grad <- raw$h_grad - 0.5 * (e / a - b * d / a^2)
  ell_hess <- p_hess - 0.5 * (abar_hess / a - tcrossprod(abar_grad) / a^2)

  ell_hess_expanded <- raw$h_hess - tcrossprod(d) / a - raw$M / (2 * a) +
    (tcrossprod(raw$f, d) + tcrossprod(d, raw$f) + b * raw$k +
       tcrossprod(e)) / (2 * a^2) -
    (cc * tcrossprod(d) + 2 * b * (tcrossprod(e, d) + tcrossprod(d, e))) /
      (2 * a^3) +
    b^2 * tcrossprod(d) / a^4

  laplace_constant <- if (isTRUE(include_laplace_constant)) 0.5 * log(2 * pi) else 0
  list(
    ok = TRUE,
    raw = raw,
    p = raw$h,
    p_grad = p_grad,
    p_hess = p_hess,
    p_third = p_third,
    abar_grad = abar_grad,
    abar_hess = abar_hess,
    ell = raw$h + laplace_constant - 0.5 * log(-a),
    ell_grad = ell_grad,
    ell_hess = ell_hess,
    ell_hess_expanded = ell_hess_expanded,
    mode_first = -d / a,
    mode_second = -raw$k / a +
      (tcrossprod(e, d) + tcrossprod(d, e)) / a^2 -
      b * tcrossprod(d) / a^3,
    coordinate_order = raw$coordinate_order
  )
}

hl_fo_contract_subject <- function(subject, transform) {
  if (!isTRUE(subject$ok)) {
    return(subject)
  }
  transform <- as.matrix(transform)
  d_out <- ncol(transform)
  if (nrow(transform) == d_out && d_out > 0L &&
      max(abs(transform - diag(d_out))) == 0) {
    return(list(
      ok = TRUE,
      raw = subject$raw,
      p = subject$p,
      ell = subject$ell,
      p_grad = subject$p_grad,
      p_hess = subject$p_hess,
      p_third = subject$p_third,
      ell_grad = subject$ell_grad,
      ell_hess = subject$ell_hess,
      ell_hess_expanded = subject$ell_hess_expanded,
      mode_first = subject$mode_first,
      mode_second = subject$mode_second
    ))
  }
  p3 <- array(0, dim = c(d_out, d_out, d_out))
  if (d_out > 0) {
    for (a in seq_len(d_out)) {
      for (b in seq_len(d_out)) {
        for (cc in seq_len(d_out)) {
          acc <- 0
          for (i in seq_len(nrow(transform))) {
            if (transform[i, a] == 0) next
            for (j in seq_len(nrow(transform))) {
              if (transform[j, b] == 0) next
              for (kk in seq_len(nrow(transform))) {
                if (transform[kk, cc] == 0) next
                acc <- acc + transform[i, a] * transform[j, b] *
                  transform[kk, cc] * subject$p_third[i, j, kk]
              }
            }
          }
          p3[a, b, cc] <- acc
        }
      }
    }
  }
  list(
    ok = TRUE,
    raw = subject$raw,
    p = subject$p,
    ell = subject$ell,
    p_grad = as.numeric(crossprod(transform, subject$p_grad)),
    p_hess = crossprod(transform, subject$p_hess %*% transform),
    p_third = p3,
    ell_grad = as.numeric(crossprod(transform, subject$ell_grad)),
    ell_hess = crossprod(transform, subject$ell_hess %*% transform),
    ell_hess_expanded = crossprod(
      transform,
      subject$ell_hess_expanded %*% transform
    ),
    mode_first = as.numeric(crossprod(transform, subject$mode_first)),
    mode_second = crossprod(transform, subject$mode_second %*% transform)
  )
}
