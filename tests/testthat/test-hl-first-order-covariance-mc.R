.hl_fo_covariance_get <- function(name) {
  getFromNamespace(name, "scRoPE")
}

test_that("canonical covariance Jacobians match numerical lower-order maps", {
  beta <- c(-0.8, 0.25, -0.15)
  theta <- c(beta, tau = log(0.45), xi = log(7.5))
  beta_map <- matrix(
    c(
      0.5, -0.2, 0.1,
      0, 1.4, 0,
      0, 0, 0.8
    ),
    3,
    3,
    byrow = TRUE
  )
  public_map <- .hl_fo_covariance_get("hl_fo_canonical_to_public")
  natural_map <- .hl_fo_covariance_get("hl_fo_canonical_to_natural")
  public <- public_map(beta, theta[4L], theta[5L], beta_map = beta_map)
  natural <- natural_map(beta, theta[4L], theta[5L])

  for (jj in seq_along(theta)) {
    public_trace <- hl_fo_fd_trace(
      function(x) public_map(x[1:3], x[4L], x[5L], beta_map = beta_map)$value,
      point = theta,
      coordinate = jj,
      analytic = public$jacobian[, jj]
    )
    natural_trace <- hl_fo_fd_trace(
      function(x) natural_map(x[1:3], x[4L], x[5L])$value,
      point = theta,
      coordinate = jj,
      analytic = natural$jacobian[, jj]
    )
    expect_lte(public_trace$analytic_discrepancy, 1e-8)
    expect_lte(natural_trace$analytic_discrepancy, 1e-8)
  }
})

test_that("seeded 1e5-draw covariance-map check satisfies the locked gate", {
  set.seed(8675309, kind = "L'Ecuyer-CMRG")
  beta <- c(-0.8, 0.25, -0.15)
  tau <- log(0.45)
  xi <- log(7.5)
  beta_map <- matrix(
    c(
      0.5, -0.2, 0.1,
      0, 1.4, 0,
      0, 0, 0.8
    ),
    3,
    3,
    byrow = TRUE
  )
  base_factor <- matrix(
    c(
      0.20, 0, 0, 0, 0,
      0.03, 0.16, 0, 0, 0,
      -0.02, 0.01, 0.13, 0, 0,
      0.015, -0.01, 0.005, 0.11, 0,
      -0.01, 0.02, -0.005, 0.01, 0.09
    ),
    5,
    5,
    byrow = TRUE
  )
  V <- tcrossprod(base_factor)
  map <- .hl_fo_covariance_get("hl_fo_canonical_to_public")(
    beta, tau, xi, beta_map = beta_map
  )$jacobian
  mapped_direct <- map %*% V %*% t(map)
  mapped_helper <- .hl_fo_covariance_get("hl_fo_cov_congruence")(V, map)
  expect_lte(hl_fo_sdisc(mapped_direct, mapped_helper), 1e-12)

  normal <- matrix(stats::rnorm(1e5 * nrow(V)), ncol = nrow(V))
  perturbation <- normal %*% chol(V)
  mapped_perturbation <- perturbation %*% t(map)
  V_mc <- stats::cov(mapped_perturbation)
  relative_frobenius <- sqrt(sum((V_mc - mapped_direct)^2)) /
    max(1, sqrt(sum(mapped_direct^2)))
  expect_lte(relative_frobenius, 0.02)
})
