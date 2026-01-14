.hl_aphl_value <- function(repml, sigma, nind, k, n_onetwo, ytwo, posindy_len) {
  if (is.null(repml)) {
    stop("repml must be a list returned from opt_pml().")
  }
  if (!is.numeric(sigma) || length(sigma) < 2) {
    stop("sigma must contain at least two dispersion parameters.")
  }

  sigma <- as.numeric(sigma)
  loglik_inner <- repml$loglik
  if (is.nan(loglik_inner)) {
    loglik_inner <- repml$loglikp
  }

  exps <- exp(sigma[1])
  alpha <- 1 / (exps - 1)
  lambda <- 1 / (sqrt(exps) * (exps - 1))
  gamma <- sigma[2]

  adj_term <- nind * gamma * log(gamma) +
    k * alpha * log(lambda) -
    k * lgamma(alpha)

  pos_count <- posindy_len - sum(n_onetwo)
  gamma_term <- sum(lgamma(ytwo + gamma)) -
    pos_count * lgamma(gamma) +
    n_onetwo[1] * log(gamma) +
    n_onetwo[2] * log(gamma + 1)

  loglik <- loglik_inner + adj_term + gamma_term
  loglik <- loglik - 0.5 * repml$logdet + log1p(repml$second)

 per_stats <- repml$per_subject_stats
 aphl_subject <- if (!is.null(per_stats)) {
   ps_loglik <- per_stats$loglik
    ps_logdet <- per_stats$log_vw
    if (is.null(ps_logdet) && !is.null(per_stats$log_H_etaeta)) {
      ps_logdet <- per_stats$log_H_etaeta
    }
    if (!is.null(ps_loglik) && !is.null(ps_logdet)) {
      ps_loglik - 0.5 * ps_logdet
    } else {
      rep(NA_real_, k)
    }
 } else {
   rep(NA_real_, k)
  }

  list(
    aphl = loglik,
    loglik_inner = loglik_inner,
    aphl_subject = aphl_subject
  )
}

hl_collect_aphl_state <- function(repml, sigma, nind, k, posv) {
  if (is.null(posv$posindy) || is.null(posv$n_onetwo) || is.null(posv$ytwo)) {
    stop("posv must include posindy, n_onetwo, and ytwo components.")
  }

  stats <- .hl_aphl_value(
    repml = repml,
    sigma = sigma,
    nind = nind,
    k = k,
    n_onetwo = posv$n_onetwo,
    ytwo = posv$ytwo,
    posindy_len = length(posv$posindy)
  )

  list(
    aphl = stats$aphl,
    loglik_inner = stats$loglik_inner,
    aphl_subject = stats$aphl_subject,
    eta_hat = as.numeric(repml$logw),
    h_etaeta = as.numeric(repml$h_etaeta),
    logdet = repml$logdet,
    second = repml$second,
    per_subject_stats = repml$per_subject_stats
  )
}
