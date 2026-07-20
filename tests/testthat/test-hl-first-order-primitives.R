.hl_fo_primitive <- function(name) {
  getFromNamespace(name, "scRoPE")
}

.hl_fo_trace_info <- function(label, trace) {
  paste(
    sprintf(
      "%s: selector=%s selected_ell=%d discrepancy=%.17g stable=%s",
      label,
      trace$selector,
      trace$selected_ell,
      trace$analytic_discrepancy,
      trace$stable_minimum
    ),
    paste(capture.output(dput(trace)), collapse = "\n"),
    sep = "\nfull trace:\n"
  )
}

.hl_fo_expect_trace <- function(trace, tolerance, label) {
  expect_true(is.finite(trace$analytic_discrepancy),
    info = .hl_fo_trace_info(label, trace))
  expect_true(trace$analytic_discrepancy <= tolerance,
    info = .hl_fo_trace_info(label, trace))
  invisible(trace)
}

.hl_fo_nb_trace <- function(y, log_m, phi, variable, lower_name, analytic_name) {
  nb_fn <- .hl_fo_primitive("hl_fo_nb_primitives")
  if (identical(variable, "log_m")) {
    point <- log_m
    lower <- function(x) nb_fn(y, x[[1L]], log(phi))[[lower_name]]
  } else if (identical(variable, "phi")) {
    point <- phi
    lower <- function(x) nb_fn(y, log_m, log(x[[1L]]))[[lower_name]]
  } else if (identical(variable, "xi")) {
    point <- log(phi)
    lower <- function(x) nb_fn(y, log_m, x[[1L]])[[lower_name]]
  } else {
    stop("Unknown NB trace variable.")
  }
  analytic <- nb_fn(y, log_m, log(phi))[[analytic_name]]
  hl_fo_fd_trace(lower, point = point, coordinate = 1L, analytic = analytic)
}

.hl_fo_gamma_trace <- function(eta, tau, variable, lower_name, analytic_name) {
  gamma_fn <- .hl_fo_primitive("hl_fo_gamma_primitives")
  if (identical(variable, "eta")) {
    point <- eta
    lower <- function(x) gamma_fn(x[[1L]], tau)[[lower_name]]
  } else if (identical(variable, "tau")) {
    point <- tau
    lower <- function(x) gamma_fn(eta, x[[1L]])[[lower_name]]
  } else {
    stop("Unknown Gamma trace variable.")
  }
  analytic <- gamma_fn(eta, tau)[[analytic_name]]
  hl_fo_fd_trace(lower, point = point, coordinate = 1L, analytic = analytic)
}

.hl_fo_subject_data <- function() {
  list(
    y = c(0, 1, 4, 2, 7),
    X = cbind(`(Intercept)` = 1, x = c(-1.2, -0.4, 0.1, 0.8, 1.5)),
    log_offset = log(c(0.8, 1.1, 0.9, 1.4, 1.05)),
    theta = c(beta0 = -0.65, beta1 = 0.28, tau = log(0.55), xi = log(6.5))
  )
}

.hl_fo_mode_root <- function(theta, data = .hl_fo_subject_data()) {
  mode_fn <- .hl_fo_primitive("hl_fo_mode_components")
  p <- ncol(data$X)
  beta <- theta[seq_len(p)]
  tau <- theta[[p + 1L]]
  xi <- theta[[p + 2L]]
  log_m_base <- as.numeric(data$log_offset + data$X %*% beta)
  score <- function(eta) {
    out <- mode_fn(data$y, log_m_base, eta, tau, xi)
    if (!isTRUE(out$ok)) {
      return(NA_real_)
    }
    out$h_eta
  }

  center <- exp(tau) / 2
  width <- max(1, 0.25 * (1 + abs(center)))
  bracket <- c(center - width, center + width)
  for (ii in 0:60) {
    values <- c(score(bracket[[1L]]), score(bracket[[2L]]))
    if (all(is.finite(values)) && values[[1L]] > 0 && values[[2L]] < 0) {
      return(uniroot(score,
        interval = bracket,
        tol = .Machine$double.eps^0.75,
        maxiter = 1000L
      )$root)
    }
    width <- width * 2
    bracket <- c(center - width, center + width)
  }
  stop("Independent mode solver did not find a bracket.")
}

.hl_fo_profile_at <- function(theta, data = .hl_fo_subject_data()) {
  subject_fn <- .hl_fo_primitive("hl_fo_subject_raw")
  p <- ncol(data$X)
  eta <- .hl_fo_mode_root(theta, data)
  out <- subject_fn(
    y = data$y,
    X = data$X,
    log_offset = data$log_offset,
    beta = theta[seq_len(p)],
    eta = eta,
    rho_full = theta[p + 1:2],
    include_laplace_constant = TRUE
  )
  if (!isTRUE(out$ok)) {
    stop("Subject evaluator failed: ", out$reason)
  }
  out$eta_hat <- eta
  out
}

test_that("complete NB values and log-mean derivatives match independent calculations", {
  nb_fn <- .hl_fo_primitive("hl_fo_nb_primitives")
  cases <- list(
    list(y = 0, log_m = -2.1, phi = 0.35),
    list(y = 3, log_m = 0.2, phi = 4.7),
    list(y = 18, log_m = 2.4, phi = 85)
  )
  chain <- list(
    s = c(lower = "loglik", analytic = "s"),
    t = c(lower = "s", analytic = "t"),
    u = c(lower = "t", analytic = "u"),
    w = c(lower = "u", analytic = "w")
  )

  for (case_id in seq_along(cases)) {
    case <- cases[[case_id]]
    out <- nb_fn(case$y, case$log_m, log(case$phi))
    expect_true(out$ok)
    expect_equal(
      out$loglik,
      dnbinom(case$y, mu = exp(case$log_m), size = case$phi, log = TRUE),
      tolerance = 1e-12
    )
    expect_equal(out$r, exp(case$log_m) / (case$phi + exp(case$log_m)), tolerance = 1e-14)

    for (derivative in names(chain)) {
      trace <- .hl_fo_nb_trace(
        case$y,
        case$log_m,
        case$phi,
        variable = "log_m",
        lower_name = chain[[derivative]][["lower"]],
        analytic_name = chain[[derivative]][["analytic"]]
      )
      .hl_fo_expect_trace(
        trace,
        1e-6,
        sprintf("NB case %d d(%s)/d(log_m) -> %s", case_id, chain[[derivative]][["lower"]], derivative)
      )
    }
  }
})

test_that("every displayed NB natural-phi and log-size mixed derivative passes lower-order differences", {
  cases <- list(
    list(y = 0, log_m = -1.6, phi = 0.45),
    list(y = 5, log_m = 0.7, phi = 3.2),
    list(y = 14, log_m = 2.0, phi = 40)
  )
  phi_chain <- list(
    l_phi = c(lower = "loglik", analytic = "l_phi"),
    l_phiphi = c(lower = "l_phi", analytic = "l_phiphi"),
    s_phi = c(lower = "s", analytic = "s_phi"),
    s_phiphi = c(lower = "s_phi", analytic = "s_phiphi"),
    t_phi = c(lower = "t", analytic = "t_phi"),
    t_phiphi = c(lower = "t_phi", analytic = "t_phiphi"),
    u_phi = c(lower = "u", analytic = "u_phi")
  )
  xi_chain <- list(
    l_xi = c(lower = "loglik", analytic = "l_xi"),
    l_xixi = c(lower = "l_xi", analytic = "l_xixi"),
    l_xixixi = c(lower = "l_xixi", analytic = "l_xixixi"),
    s_xi = c(lower = "s", analytic = "s_xi"),
    s_xixi = c(lower = "s_xi", analytic = "s_xixi"),
    t_xi = c(lower = "t", analytic = "t_xi"),
    t_xixi = c(lower = "t_xi", analytic = "t_xixi"),
    u_xi = c(lower = "u", analytic = "u_xi")
  )

  for (case_id in seq_along(cases)) {
    case <- cases[[case_id]]
    for (derivative in names(phi_chain)) {
      trace <- .hl_fo_nb_trace(
        case$y, case$log_m, case$phi,
        variable = "phi",
        lower_name = phi_chain[[derivative]][["lower"]],
        analytic_name = phi_chain[[derivative]][["analytic"]]
      )
      .hl_fo_expect_trace(
        trace,
        1e-6,
        sprintf("NB case %d natural-phi %s", case_id, derivative)
      )
    }
    for (derivative in names(xi_chain)) {
      trace <- .hl_fo_nb_trace(
        case$y, case$log_m, case$phi,
        variable = "xi",
        lower_name = xi_chain[[derivative]][["lower"]],
        analytic_name = xi_chain[[derivative]][["analytic"]]
      )
      .hl_fo_expect_trace(
        trace,
        1e-6,
        sprintf("NB case %d log-size %s", case_id, derivative)
      )
    }
  }
})

test_that("complete Gamma values and all derivatives through eta-eta-tau-tau pass lower-order differences", {
  gamma_fn <- .hl_fo_primitive("hl_fo_gamma_primitives")
  cases <- list(
    list(eta = -1.1, tau = log(0.08)),
    list(eta = 0.35, tau = log(0.7)),
    list(eta = 1.4, tau = log(2.2))
  )
  eta_chain <- list(
    g_eta = c(lower = "g", analytic = "g_eta"),
    g_etaeta = c(lower = "g_eta", analytic = "g_etaeta"),
    g_etaetaeta = c(lower = "g_etaeta", analytic = "g_etaetaeta"),
    g_etaetaetaeta = c(lower = "g_etaetaeta", analytic = "g_etaetaetaeta")
  )
  tau_chain <- list(
    g_tau = c(lower = "g", analytic = "g_tau"),
    g_tautau = c(lower = "g_tau", analytic = "g_tautau"),
    g_tautautau = c(lower = "g_tautau", analytic = "g_tautautau"),
    g_eta_tau = c(lower = "g_eta", analytic = "g_eta_tau"),
    g_etaeta_tau = c(lower = "g_etaeta", analytic = "g_etaeta_tau"),
    g_etaetaeta_tau = c(lower = "g_etaetaeta", analytic = "g_etaetaeta_tau"),
    g_eta_tautau = c(lower = "g_eta_tau", analytic = "g_eta_tautau"),
    g_etaeta_tautau = c(lower = "g_etaeta_tau", analytic = "g_etaeta_tautau")
  )

  for (case_id in seq_along(cases)) {
    case <- cases[[case_id]]
    out <- gamma_fn(case$eta, case$tau)
    expect_true(out$ok)
    expect_equal(
      out$g,
      dgamma(exp(case$eta), shape = out$alpha, rate = out$lambda, log = TRUE) + case$eta,
      tolerance = 1e-12
    )

    for (derivative in names(eta_chain)) {
      trace <- .hl_fo_gamma_trace(
        case$eta, case$tau,
        variable = "eta",
        lower_name = eta_chain[[derivative]][["lower"]],
        analytic_name = eta_chain[[derivative]][["analytic"]]
      )
      .hl_fo_expect_trace(
        trace,
        1e-6,
        sprintf("Gamma case %d eta derivative %s", case_id, derivative)
      )
    }
    for (derivative in names(tau_chain)) {
      trace <- .hl_fo_gamma_trace(
        case$eta, case$tau,
        variable = "tau",
        lower_name = tau_chain[[derivative]][["lower"]],
        analytic_name = tau_chain[[derivative]][["analytic"]]
      )
      .hl_fo_expect_trace(
        trace,
        1e-6,
        sprintf("Gamma case %d tau derivative %s", case_id, derivative)
      )
    }
  }
})

test_that("subject mode implicit first and second derivatives use independently re-solved modes", {
  data <- .hl_fo_subject_data()
  theta <- data$theta
  subject <- .hl_fo_profile_at(theta, data)

  expect_lt(abs(subject$raw$h_eta), 1e-10)
  expect_lt(subject$raw$a, 0)

  for (rr in seq_along(theta)) {
    first_trace <- hl_fo_fd_trace(
      function(x) .hl_fo_mode_root(x, data),
      point = theta,
      coordinate = rr,
      analytic = subject$mode_first[[rr]]
    )
    .hl_fo_expect_trace(first_trace, 1e-6, sprintf("mode first coordinate %d", rr))

    second_trace <- hl_fo_fd_trace(
      function(x) .hl_fo_profile_at(x, data)$mode_first,
      point = theta,
      coordinate = rr,
      analytic = subject$mode_second[, rr]
    )
    .hl_fo_expect_trace(second_trace, 1e-5, sprintf("mode second column %d", rr))
  }
})

test_that("profiled p derivatives through third order use fully re-solved modes", {
  data <- .hl_fo_subject_data()
  theta <- data$theta
  subject <- .hl_fo_profile_at(theta, data)
  gp <- subject$raw$gamma
  independent_h <- sum(dnbinom(
    data$y,
    mu = exp(data$log_offset + data$X %*% theta[1:2] + subject$eta_hat),
    size = exp(theta[[4L]]),
    log = TRUE
  )) + dgamma(exp(subject$eta_hat), shape = gp$alpha, rate = gp$lambda, log = TRUE) + subject$eta_hat
  expect_equal(unname(subject$p), unname(independent_h), tolerance = 1e-10)

  for (rr in seq_along(theta)) {
    gradient_trace <- hl_fo_fd_trace(
      function(x) .hl_fo_profile_at(x, data)$p,
      point = theta,
      coordinate = rr,
      analytic = subject$p_grad[[rr]]
    )
    .hl_fo_expect_trace(gradient_trace, 1e-7, sprintf("p gradient coordinate %d", rr))

    hessian_trace <- hl_fo_fd_trace(
      function(x) .hl_fo_profile_at(x, data)$p_grad,
      point = theta,
      coordinate = rr,
      analytic = subject$p_hess[, rr]
    )
    .hl_fo_expect_trace(hessian_trace, 1e-6, sprintf("p Hessian column %d", rr))

    third_trace <- hl_fo_fd_trace(
      function(x) as.numeric(.hl_fo_profile_at(x, data)$p_hess),
      point = theta,
      coordinate = rr,
      analytic = as.numeric(subject$p_third[, , rr])
    )
    .hl_fo_expect_trace(third_trace, 1e-5, sprintf("complete p third-derivative slice %d", rr))
  }
})

test_that("Laplace gradient and Hessian pass lower-order differences and both Hessian formulas agree", {
  data <- .hl_fo_subject_data()
  theta <- data$theta
  subject <- .hl_fo_profile_at(theta, data)

  expect_equal(
    subject$ell,
    subject$p + 0.5 * log(2 * pi) - 0.5 * log(-subject$raw$a),
    tolerance = 1e-12
  )
  expect_lte(hl_fo_sdisc(subject$ell_hess, subject$ell_hess_expanded), 1e-10)

  for (rr in seq_along(theta)) {
    gradient_trace <- hl_fo_fd_trace(
      function(x) .hl_fo_profile_at(x, data)$ell,
      point = theta,
      coordinate = rr,
      analytic = subject$ell_grad[[rr]]
    )
    .hl_fo_expect_trace(gradient_trace, 1e-6, sprintf("Laplace gradient coordinate %d", rr))

    hessian_trace <- hl_fo_fd_trace(
      function(x) .hl_fo_profile_at(x, data)$ell_grad,
      point = theta,
      coordinate = rr,
      analytic = subject$ell_hess[, rr]
    )
    .hl_fo_expect_trace(hessian_trace, 1e-5, sprintf("Laplace Hessian column %d", rr))
  }
})

test_that("contracted-coordinate profile and Laplace derivatives retain full-coordinate algebra", {
  contract_fn <- .hl_fo_primitive("hl_fo_contract_subject")
  data <- .hl_fo_subject_data()
  theta <- data$theta
  K <- c(1, -0.45)
  K <- K / sqrt(sum(K^2))
  transform <- matrix(0, nrow = 4L, ncol = 3L)
  transform[1:2, 1L] <- K
  transform[3L, 2L] <- 1
  transform[4L, 3L] <- 1
  reduced <- c(gamma = 0.2, tau = theta[[3L]], xi = theta[[4L]])
  anchor <- theta - drop(transform %*% reduced)
  at_reduced <- function(value) {
    full_theta <- anchor + drop(transform %*% value)
    full <- .hl_fo_profile_at(full_theta, data)
    contract_fn(full, transform)
  }

  full <- .hl_fo_profile_at(theta, data)
  contracted <- at_reduced(reduced)
  expect_equal(contracted$p_grad, as.numeric(crossprod(transform, full$p_grad)), tolerance = 1e-12)
  expect_equal(contracted$p_hess, crossprod(transform, full$p_hess %*% transform), tolerance = 1e-12)
  expect_equal(contracted$ell_grad, as.numeric(crossprod(transform, full$ell_grad)), tolerance = 1e-12)
  expect_equal(contracted$ell_hess, crossprod(transform, full$ell_hess %*% transform), tolerance = 1e-12)
  expect_lte(hl_fo_sdisc(contracted$ell_hess, contracted$ell_hess_expanded), 1e-10)

  for (rr in seq_along(reduced)) {
    p_gradient_trace <- hl_fo_fd_trace(
      function(x) at_reduced(x)$p,
      point = reduced,
      coordinate = rr,
      analytic = contracted$p_grad[[rr]]
    )
    .hl_fo_expect_trace(p_gradient_trace, 1e-7, sprintf("contracted p gradient coordinate %d", rr))

    p_hessian_trace <- hl_fo_fd_trace(
      function(x) at_reduced(x)$p_grad,
      point = reduced,
      coordinate = rr,
      analytic = contracted$p_hess[, rr]
    )
    .hl_fo_expect_trace(p_hessian_trace, 1e-6, sprintf("contracted p Hessian column %d", rr))

    # Reduced coordinate one is the contracted beta direction; as above,
    # exclude only the unavailable pure active-dispersion third block.
    required <- if (rr == 1L) {
      matrix(TRUE, nrow = length(reduced), ncol = length(reduced))
    } else {
      outer(seq_along(reduced), seq_along(reduced), function(a, b) {
        a == 1L | b == 1L
      })
    }
    p_third_trace <- hl_fo_fd_trace(
      function(x) as.numeric(at_reduced(x)$p_hess[required]),
      point = reduced,
      coordinate = rr,
      analytic = as.numeric(contracted$p_third[, , rr][required])
    )
    .hl_fo_expect_trace(p_third_trace, 1e-5, sprintf("required contracted p third slice %d", rr))

    ell_gradient_trace <- hl_fo_fd_trace(
      function(x) at_reduced(x)$ell,
      point = reduced,
      coordinate = rr,
      analytic = contracted$ell_grad[[rr]]
    )
    .hl_fo_expect_trace(ell_gradient_trace, 1e-6, sprintf("contracted ell gradient coordinate %d", rr))

    ell_hessian_trace <- hl_fo_fd_trace(
      function(x) at_reduced(x)$ell_grad,
      point = reduced,
      coordinate = rr,
      analytic = contracted$ell_hess[, rr]
    )
    .hl_fo_expect_trace(ell_hessian_trace, 1e-5, sprintf("contracted ell Hessian column %d", rr))
  }
})

test_that("B1-scale Gamma and NB primitives remain finite and differentiable", {
  gamma_fn <- .hl_fo_primitive("hl_fo_gamma_primitives")
  nb_fn <- .hl_fo_primitive("hl_fo_nb_primitives")
  eta <- 2e-5
  tau <- log(1e-8)
  gamma <- gamma_fn(eta, tau)
  expect_true(gamma$ok)
  expect_true(all(is.finite(unlist(gamma[setdiff(names(gamma), "ok")]))) )
  h_tau <- 1e-4
  gamma_plus <- gamma_fn(eta, tau + h_tau)
  gamma_minus <- gamma_fn(eta, tau - h_tau)
  expect_lte(hl_fo_sdisc(
    (gamma_plus$g - gamma_minus$g) / (2 * h_tau), gamma$g_tau
  ), 1e-8)
  expect_lte(hl_fo_sdisc(
    (gamma_plus$g_tau - gamma_minus$g_tau) / (2 * h_tau),
    gamma$g_tautau
  ), 1e-8)
  expect_lte(hl_fo_sdisc(
    (gamma_plus$g_tautau - gamma_minus$g_tautau) / (2 * h_tau),
    gamma$g_tautautau
  ), 1e-8)

  y <- c(0, 1, 2, 5)
  log_m <- log(c(0.2, 1, 3, 7))
  xi <- log(1e8)
  nb <- nb_fn(y, log_m, xi)
  expect_true(nb$ok)
  expect_true(all(is.finite(unlist(nb[setdiff(names(nb), "ok")]))) )
  expect_equal(
    nb$loglik,
    dnbinom(y, size = exp(xi), mu = exp(log_m), log = TRUE),
    tolerance = 1e-7
  )
  h_xi <- 1e-2
  nb_plus <- nb_fn(y, log_m, xi + h_xi)
  nb_minus <- nb_fn(y, log_m, xi - h_xi)
  expect_lte(hl_fo_sdisc(
    (nb_plus$loglik - nb_minus$loglik) / (2 * h_xi), nb$l_xi
  ), 1e-8)
  expect_lte(hl_fo_sdisc(
    (nb_plus$l_xi - nb_minus$l_xi) / (2 * h_xi), nb$l_xixi
  ), 1e-8)
  expect_lte(hl_fo_sdisc(
    (nb_plus$l_xixi - nb_minus$l_xixi) / (2 * h_xi), nb$l_xixixi
  ), 1e-8)
})
