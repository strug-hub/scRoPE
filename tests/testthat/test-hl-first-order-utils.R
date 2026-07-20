hl_fo_test_fd_jacobian <- function(fn, x, step = 1e-6) {
  f0 <- fn(x)
  out <- matrix(0, nrow = length(f0), ncol = length(x))
  for (jj in seq_along(x)) {
    plus <- minus <- x
    h <- step * (1 + abs(x[[jj]]))
    plus[[jj]] <- plus[[jj]] + h
    minus[[jj]] <- minus[[jj]] - h
    out[, jj] <- (fn(plus) - fn(minus)) / (2 * h)
  }
  out
}

test_that("first-order numerical conventions use the declared formulas", {
  a <- matrix(c(1, -3, 5, 2), 2, 2)
  b <- matrix(c(1.5, -2, 4, 2), 2, 2)
  expected <- max(abs(a - b) / (1 + pmax(abs(a), abs(b))))
  expect_equal(hl_fo_sdisc(a, b), expected)
  expect_equal(hl_fo_sdisc(numeric(0), numeric(0)), 0)
  expect_error(hl_fo_sdisc(a, as.vector(b)), "equally shaped")
  expect_equal(hl_fo_sdisc(c(1, Inf), c(1, Inf)), Inf)

  M <- matrix(c(1, 2, 2, 4), 2, 2)
  info <- hl_fo_svd_info(M)
  expect_equal(info$rank, 1L)
  expect_equal(
    info$threshold,
    max(dim(M)) * .Machine$double.eps * info$singular_values[[1L]]
  )
  expect_equal(M %*% info$pinv %*% M, M, tolerance = 1e-12)
  expect_equal(hl_fo_num_rank(diag(c(1, 1e-16)))$rank, 1L)

  zero <- matrix(0, 3, 2)
  expect_equal(hl_fo_num_rank(zero)$rank, 0L)
  expect_equal(hl_fo_pinv(zero), matrix(0, 2, 3))

  symmetric <- matrix(c(2, 1, 1, 3), 2, 2)
  expect_equal(hl_fo_symmetry_defect(symmetric), 0)
  nonsymmetric <- matrix(c(2, 0, 1, 3), 2, 2)
  expect_equal(
    hl_fo_symmetry_defect(nonsymmetric),
    sqrt(sum((nonsymmetric - t(nonsymmetric))^2)) /
      max(1, sqrt(sum(nonsymmetric^2)))
  )
})

test_that("RHS consistency follows the common SVD projection rule", {
  M <- rbind(c(1, 0), c(2, 0), c(0, 1))
  good <- hl_fo_rhs_consistency(M, c(3, 6, -1))
  bad <- hl_fo_rhs_consistency(M, c(3, 7, -1))

  expect_true(good$consistent)
  expect_lte(good$residual_norm, good$tolerance)
  expect_false(bad$consistent)
  expect_gt(bad$residual_norm, bad$tolerance)

  full_row_rank <- matrix(c(-1.98935170, 0.61982575, -0.05612874,
                            -0.1557955, -1.4707524, -0.4781501,
                            0.4179416, 1.3586796, -0.1027877), 3, 3)
  expect_true(
    hl_fo_rhs_consistency(full_row_rank, c(0.4, -0.05, -1.4))$consistent
  )
})

test_that("affine constraints are rank-reduced deterministically", {
  C <- rbind(
    c(1, 0, 0),
    c(1, 0, 0),
    c(0, 2, 0)
  )
  rhs <- c(2, 2, 6)
  fit <- hl_fo_affine_constraint(C, rhs)
  replay <- hl_fo_affine_constraint(C, rhs)

  expect_true(fit$success)
  expect_equal(fit$rank, 2L)
  expect_equal(fit$selected_rows, c(1L, 3L))
  expect_equal(fit$discarded_rows, 2L)
  expect_equal(fit$beta_star, c(2, 3, 0), tolerance = 1e-14)
  expect_equal(fit$K, matrix(c(0, 0, 1), 3, 1), tolerance = 1e-14)
  expect_equal(fit$K, replay$K, tolerance = 0)
  expect_equal(crossprod(fit$K), diag(1), tolerance = 1e-14)
  expect_equal(C %*% fit$K, matrix(0, 3, 1), tolerance = 1e-14)
  expect_equal(as.vector(C %*% fit$beta_star), rhs, tolerance = 1e-14)
  expect_lte(fit$residual_scaled, 100 * .Machine$double.eps)
  expect_true(fit$residual_satisfied)
  trace_residual <- hl_fo_constraint_residual(
    C, fit$beta_star + as.vector(fit$K %*% 1.25), rhs
  )
  expect_true(trace_residual$satisfied)
  expect_equal(trace_residual$residual, numeric(nrow(C)), tolerance = 1e-14)

  infeasible <- hl_fo_affine_constraint(C, c(2, 3, 6))
  expect_false(infeasible$success)
  expect_equal(infeasible$status, "constraint_infeasible")
})

test_that("affine constraints cover K equals I and k equals zero", {
  unrestricted <- hl_fo_affine_constraint(C = NULL, p = 3L)
  expect_true(unrestricted$success)
  expect_equal(unrestricted$rank, 0L)
  expect_equal(unrestricted$beta_star, numeric(3))
  expect_equal(unrestricted$K, diag(3))
  expect_equal(unrestricted$free_dim, 3L)

  fixed <- hl_fo_affine_constraint(diag(3), c(1, -2, 0.5))
  expect_true(fixed$success)
  expect_equal(fixed$beta_star, c(1, -2, 0.5), tolerance = 1e-14)
  expect_equal(dim(fixed$K), c(3L, 0L))
  expect_equal(fixed$free_dim, 0L)

  C <- matrix(c(1, 1, 0), 1, 3)
  one_row <- hl_fo_affine_constraint(C, 2)
  expect_equal(one_row$beta_star, c(1, 1, 0), tolerance = 1e-14)
  expect_equal(C %*% one_row$K, matrix(0, 1, 2), tolerance = 1e-14)
  expect_equal(crossprod(one_row$K), diag(2), tolerance = 1e-14)
  for (jj in seq_len(ncol(one_row$K))) {
    pivot <- which(abs(one_row$K[, jj]) == max(abs(one_row$K[, jj])))[1L]
    expect_gte(one_row$K[pivot, jj], 0)
  }
})

test_that("canonical natural and public transformations have exact Jacobians", {
  beta <- c(0.3, -0.2, 1.1)
  theta <- c(beta, log(0.4), log(8))

  natural <- hl_fo_canonical_to_natural(
    theta[1:3], theta[4], theta[5]
  )
  expect_equal(natural$value, c(beta, 0.4, 8), tolerance = 1e-14)
  numerical_natural <- hl_fo_test_fd_jacobian(function(z) {
    hl_fo_canonical_to_natural(z[1:3], z[4], z[5])$value
  }, theta)
  expect_lte(hl_fo_sdisc(natural$jacobian, numerical_natural), 1e-8)

  canonical <- hl_fo_natural_to_canonical(
    natural$beta, natural$sigma2, natural$phi
  )
  expect_equal(canonical$value, theta, tolerance = 1e-14)
  expect_equal(canonical$jacobian %*% natural$jacobian, diag(5),
               tolerance = 1e-14)

  public <- hl_fo_canonical_to_public(beta, theta[4], theta[5])
  expect_equal(public$value, c(beta, 0.4, 1 / 8), tolerance = 1e-14)
  expect_equal(public$jacobian[5, 5], -1 / 8, tolerance = 1e-14)
  numerical_public <- hl_fo_test_fd_jacobian(function(z) {
    hl_fo_canonical_to_public(z[1:3], z[4], z[5])$value
  }, theta)
  expect_lte(hl_fo_sdisc(public$jacobian, numerical_public), 1e-8)

  public_inverse <- hl_fo_public_to_canonical(
    public$beta, public$sigma2, public$delta
  )
  expect_equal(public_inverse$value, theta, tolerance = 1e-14)
  expect_equal(public_inverse$jacobian %*% public$jacobian, diag(5),
               tolerance = 1e-14)
})

test_that("mean-one seed transformation includes the intercept-variance term", {
  beta <- c(0.4, -1.2, 0.2)
  tau <- log(0.6)
  xi <- log(5)
  mapped <- hl_fo_canonical_to_mean_one(beta, tau, xi, intercept_idx = 2L)

  expected_beta <- beta
  expected_beta[2] <- expected_beta[2] + 0.3
  expect_equal(mapped$beta, expected_beta, tolerance = 1e-14)
  expect_equal(mapped$jacobian[2, 4], 0.3, tolerance = 1e-14)

  theta <- c(beta, tau, xi)
  numerical <- hl_fo_test_fd_jacobian(function(z) {
    hl_fo_canonical_to_mean_one(z[1:3], z[4], z[5], 2L)$value
  }, theta)
  expect_lte(hl_fo_sdisc(mapped$jacobian, numerical), 1e-8)

  inverse <- hl_fo_mean_one_to_canonical(
    mapped$beta, mapped$sigma2, mapped$phi, intercept_idx = 2L
  )
  expect_equal(inverse$value, theta, tolerance = 1e-14)
  expect_equal(inverse$jacobian %*% mapped$jacobian, diag(5),
               tolerance = 1e-14)

  V <- crossprod(matrix(seq_len(25) / 50, 5, 5)) + diag(5)
  V_mapped <- hl_fo_cov_congruence(V, mapped$jacobian)
  V_roundtrip <- hl_fo_cov_congruence(V_mapped, inverse$jacobian)
  expect_lte(hl_fo_sdisc(V, V_roundtrip), 1e-10)
})

test_that("design map handles arbitrary intcol and raw constant value", {
  Z <- cbind(
    varying_a = c(0, 2, 4, 6),
    raw_intercept = rep(-3, 4),
    varying_b = c(-4, -2, 0, 2)
  )
  means <- colMeans(Z)
  scales <- sqrt(colMeans(sweep(Z, 2, means)^2))
  # The constant column scale is not part of the varying-column map.
  scales[2] <- 0
  map <- hl_fo_design_map(means, scales, intercept_idx = 2L,
                          raw_constant = -3)
  X <- Z
  X[, 1] <- (Z[, 1] - means[1]) / scales[1]
  X[, 2] <- 1
  X[, 3] <- (Z[, 3] - means[3]) / scales[3]

  beta_internal <- c(0.7, -0.5, 1.2)
  original <- hl_fo_internal_to_original_beta(beta_internal, map)
  expect_equal(as.vector(Z %*% original$beta),
               as.vector(X %*% beta_internal), tolerance = 1e-13)
  roundtrip <- hl_fo_original_to_internal_beta(original$beta, map)
  expect_equal(roundtrip$beta, beta_internal, tolerance = 1e-14)
  expect_equal(map$to_internal %*% map$to_original, diag(3),
               tolerance = 1e-14)

  C_report <- matrix(c(1, -0.5, 2), 1, 3)
  rhs <- as.numeric(C_report %*% original$beta)
  transformed <- hl_fo_transform_constraint_to_internal(
    C_report, rhs, map$to_original
  )
  expect_equal(as.numeric(transformed$C %*% beta_internal), rhs,
               tolerance = 1e-14)
  expect_equal(transformed$rhs, rhs)

  V <- crossprod(matrix(c(1, 0.2, -0.1, 0.2, 2, 0.3,
                          -0.1, 0.3, 1.5), 3, 3))
  V_original <- hl_fo_cov_congruence(V, map$to_original)
  V_back <- hl_fo_cov_congruence(V_original, map$to_internal)
  expect_lte(hl_fo_sdisc(V, V_back), 1e-10)

  public <- hl_fo_canonical_to_public(
    beta_internal, log(0.4), log(8), beta_map = map$to_original
  )
  canonical <- hl_fo_public_to_canonical(
    public$beta, public$sigma2, public$delta,
    beta_inverse = map$to_internal
  )
  expect_equal(canonical$value,
               c(beta_internal, log(0.4), log(8)), tolerance = 1e-13)
  expect_equal(canonical$jacobian %*% public$jacobian, diag(5),
               tolerance = 1e-13)
})

test_that("active dispersion helpers preserve canonical order including q zero", {
  both <- hl_fo_active_coordinates(c("xi", "tau"))
  expect_equal(both$names, c("tau", "xi"))
  expect_equal(both$P_A, diag(2))
  expect_equal(hl_fo_select_active_rho(1, 2, c("xi", "tau")),
               c(tau = 1, xi = 2))

  xi_only <- hl_fo_active_coordinates("xi")
  expect_equal(xi_only$P_A, matrix(c(0, 1), 1, 2))
  insertion <- hl_fo_active_insertion(2L, "xi")
  expect_equal(dim(insertion), c(4L, 3L))
  expect_equal(insertion[4, 3], 1)
  expect_equal(insertion[3, ], c(0, 0, 0))

  expanded <- hl_fo_expand_active_covariance(diag(c(2, 3, 4)), 2L, "xi")
  expect_equal(diag(expanded$covariance), c(2, 3, 0, 4))

  q0 <- hl_fo_active_coordinates(character(0))
  expect_equal(q0$q, 0L)
  q0_insertion <- hl_fo_active_insertion(2L, character(0))
  expect_equal(dim(q0_insertion), c(4L, 2L))
  q0_cov <- hl_fo_expand_active_covariance(diag(c(2, 3)), 2L,
                                           character(0))$covariance
  expect_equal(diag(q0_cov), c(2, 3, 0, 0))
})

test_that("untagged beta seeds are rejected and tagged seeds map exactly", {
  expect_error(hl_fo_canonicalize_beta_seed(c(1, 2)), "Untagged")
  expect_error(
    hl_fo_canonicalize_beta_seed(list(beta = c(1, 2))),
    "Untagged"
  )

  canonical_seed <- hl_fo_make_beta_seed(
    c(1, 2, 3), "canonical_conditional", sigma2 = 0.4, phi = 8,
    source = "core$repml$beta"
  )
  canonical <- hl_fo_canonicalize_beta_seed(canonical_seed, expected_p = 3L)
  expect_equal(canonical$beta, c(1, 2, 3))
  expect_equal(canonical$tau, log(0.4))
  expect_equal(canonical$xi, log(8))

  mean_one_seed <- hl_fo_make_beta_seed(
    c(1, 2, 3), "mean_one", sigma2 = 0.4, phi = 8,
    source = "raw betae"
  )
  shifted <- hl_fo_canonicalize_beta_seed(mean_one_seed, intercept_idx = 3L)
  expect_equal(shifted$beta, c(1, 2, 2.8), tolerance = 1e-14)
  expect_error(
    hl_fo_canonicalize_beta_seed(
      hl_fo_make_beta_seed(c(1, 2), "mean_one"), intercept_idx = 1L
    ),
    "requires sigma2"
  )
})

test_that("robust Wald rank-reduces contrasts without ridge promotion", {
  beta <- c(1, -0.5, 0.25)
  V <- diag(c(0.04, 0.09, 0.16))
  R <- rbind(c(1, 0, 0), c(2, 0, 0), c(0, 1, 0))
  fit <- hl_fo_robust_wald(beta, V, R, rhs = c(0, 0, 0),
                           n_subjects = 240L)

  expect_true(fit$success)
  expect_equal(fit$status, "ok")
  expect_equal(fit$df, 2L)
  expect_equal(fit$contrast_rank, 2L)
  expect_equal(fit$covariance_rank, 2L)
  expect_equal(fit$n_subjects, 240L)
  expect_equal(fit$ridge, 0)
  expect_equal(fit$statistic, 1^2 / 0.04 + (-0.5)^2 / 0.09,
               tolerance = 1e-12)
  expect_equal(fit$p_value,
               stats::pchisq(fit$statistic, 2, lower.tail = FALSE))

  singular <- hl_fo_robust_wald(beta, diag(c(0, 0.09, 0.16)),
                                matrix(c(1, 0, 0), 1, 3))
  expect_false(singular$success)
  expect_equal(singular$status, "wald_covariance_rank_deficient")
  expect_equal(singular$ridge, 0)

  indefinite <- hl_fo_robust_wald(
    beta, diag(c(-0.04, 0.09, 0.16)), matrix(c(0, 1, 0), 1, 3)
  )
  expect_false(indefinite$success)
  expect_equal(indefinite$status, "wald_covariance_invalid")

  inconsistent <- hl_fo_robust_wald(
    beta, V, rbind(c(1, 0, 0), c(2, 0, 0)), rhs = c(0, 1)
  )
  expect_false(inconsistent$success)
  expect_equal(inconsistent$status, "wald_contrast_inconsistent")
})

test_that("Wald inference is invariant under a full beta reporting map", {
  beta_internal <- c(0.7, -0.5, 1.2)
  V_internal <- matrix(c(0.2, 0.01, -0.02,
                         0.01, 0.1, 0.03,
                         -0.02, 0.03, 0.3), 3, 3)
  map <- hl_fo_design_map(
    column_means = c(2, 0, -1),
    column_scales = c(4, 0, 2),
    intercept_idx = 2L,
    raw_constant = -3
  )
  beta_report <- as.vector(map$to_original %*% beta_internal)
  V_report <- hl_fo_cov_congruence(V_internal, map$to_original)
  R_report <- matrix(c(1, -0.5, 2), 1, 3)
  rhs <- 0.1

  report_fit <- hl_fo_robust_wald(
    beta_report, V_report, R_report, rhs = rhs, n_subjects = 300L
  )
  internal_contrast <- hl_fo_transform_constraint_to_internal(
    R_report, rhs, map$to_original
  )
  internal_fit <- hl_fo_robust_wald(
    beta_internal, V_internal, internal_contrast$C,
    rhs = internal_contrast$rhs, n_subjects = 300L
  )

  expect_true(report_fit$success)
  expect_true(internal_fit$success)
  expect_equal(report_fit$statistic, internal_fit$statistic,
               tolerance = 1e-12)
  expect_equal(report_fit$p_value, internal_fit$p_value,
               tolerance = 1e-12)
})
