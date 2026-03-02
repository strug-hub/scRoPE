test_that("hl_laplace_hessian_from_stats includes b_i f_{i,tau,tau} with fallback", {
  fn <- getFromNamespace("hl_laplace_hessian_from_stats", "scRoPE")

  per_stats <- list(
    h_beta = matrix(0, nrow = 1, ncol = 1),
    h_etaeta = c(2),
    h_etaetaeta = c(3),
    h_etaetaetaeta = c(5),
    h_etaeta_beta = matrix(0, nrow = 1, ncol = 1),
    h_etaeta_phi = c(0),
    h_etaeta_tau = c(7),
    h_eta_beta = matrix(0, nrow = 1, ncol = 1),
    h_eta_phi = c(0),
    h_eta_tau = c(11),
    h_etaeta_beta_beta = list(matrix(0, nrow = 1, ncol = 1)),
    g_beta_beta = list(matrix(0, nrow = 1, ncol = 1)),
    h_etaeta_beta_phi = list(0),
    f_beta_phi = list(0),
    g_beta_phi = list(0),
    h_etaeta_phi_phi = c(0),
    g_phi_phi = c(0),
    b_phi = c(0),
    b_beta = matrix(0, nrow = 1, ncol = 1),
    g_tau_tau = c(13)
  )

  out <- fn(per_stats)

  a <- per_stats$h_etaeta[1]
  b <- per_stats$h_etaetaeta[1]
  c4 <- per_stats$h_etaetaetaeta[1]
  d_tau <- per_stats$h_eta_tau[1]
  e_tau <- per_stats$h_etaeta_tau[1]
  g_tautau <- per_stats$g_tau_tau[1]
  f_tautau <- -d_tau

  expected <- (g_tautau / a) -
    ((e_tau * e_tau) + b * f_tautau + 2 * e_tau * d_tau) / (a * a) +
    (4 * b * e_tau * d_tau + c4 * d_tau * d_tau) / (a^3) -
    (2 * b * b * d_tau * d_tau) / (a^4)

  expect_equal(out$H_tautau[1, 1], expected)
})

test_that("hl_laplace_hessian_from_stats uses explicit h_eta_tau_tau when provided", {
  fn <- getFromNamespace("hl_laplace_hessian_from_stats", "scRoPE")

  per_stats <- list(
    h_beta = matrix(0, nrow = 1, ncol = 1),
    h_etaeta = c(2),
    h_etaetaeta = c(3),
    h_etaetaetaeta = c(5),
    h_etaeta_beta = matrix(0, nrow = 1, ncol = 1),
    h_etaeta_phi = c(0),
    h_etaeta_tau = c(7),
    h_eta_beta = matrix(0, nrow = 1, ncol = 1),
    h_eta_phi = c(0),
    h_eta_tau = c(11),
    h_eta_tau_tau = c(17),
    h_etaeta_beta_beta = list(matrix(0, nrow = 1, ncol = 1)),
    g_beta_beta = list(matrix(0, nrow = 1, ncol = 1)),
    h_etaeta_beta_phi = list(0),
    f_beta_phi = list(0),
    g_beta_phi = list(0),
    h_etaeta_phi_phi = c(0),
    g_phi_phi = c(0),
    b_phi = c(0),
    b_beta = matrix(0, nrow = 1, ncol = 1),
    g_tau_tau = c(13)
  )

  out <- fn(per_stats)

  a <- per_stats$h_etaeta[1]
  b <- per_stats$h_etaetaeta[1]
  c4 <- per_stats$h_etaetaetaeta[1]
  d_tau <- per_stats$h_eta_tau[1]
  e_tau <- per_stats$h_etaeta_tau[1]
  g_tautau <- per_stats$g_tau_tau[1]
  f_tautau <- per_stats$h_eta_tau_tau[1]

  expected <- (g_tautau / a) -
    ((e_tau * e_tau) + b * f_tautau + 2 * e_tau * d_tau) / (a * a) +
    (4 * b * e_tau * d_tau + c4 * d_tau * d_tau) / (a^3) -
    (2 * b * b * d_tau * d_tau) / (a^4)

  expect_equal(out$H_tautau[1, 1], expected)
})
