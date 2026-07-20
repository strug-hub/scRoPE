test_that("first-order robust Wald rank-reduces contrasts without a ridge", {
  wald_fn <- getFromNamespace("hl_fo_unrestricted_wald", "scRoPE")
  affine_fn <- getFromNamespace("hl_fo_affine_constraint", "scRoPE")
  fit <- list(
    beta = c(0.5, -0.25),
    covariance = list(V_beta = matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)),
    data = list(m = 60L),
    point_success = TRUE,
    covariance_success = TRUE,
    constraint = affine_fn(C = NULL, p = 2L),
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  R <- rbind(c(1, 0), c(2, 0), c(0, 1))

  out <- wald_fn(fit, R, rhs = c(0, 0, 0))
  expected <- as.numeric(crossprod(fit$beta, solve(fit$covariance$V_beta, fit$beta)))

  expect_true(out$success)
  expect_identical(out$method, "wald")
  expect_identical(out$status, "ok")
  expect_identical(out$df, 2L)
  expect_identical(out$contrast_rank, 2L)
  expect_identical(out$covariance_rank, 2L)
  expect_identical(out$m, 60L)
  expect_equal(out$stat, expected, tolerance = 1e-12)
  expect_equal(out$p, pchisq(expected, 2, lower.tail = FALSE))
  expect_null(out$ridge)
})

test_that("first-order robust Wald rejects inconsistent or rank-deficient inputs", {
  wald_fn <- getFromNamespace("hl_fo_unrestricted_wald", "scRoPE")
  affine_fn <- getFromNamespace("hl_fo_affine_constraint", "scRoPE")
  fit <- list(
    beta = c(0.5, -0.25),
    covariance = list(V_beta = diag(c(0.04, 0))),
    data = list(m = 30L),
    point_success = TRUE,
    covariance_success = TRUE,
    constraint = affine_fn(C = NULL, p = 2L),
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )

  rank_fail <- wald_fn(fit, diag(2))
  expect_false(rank_fail$success)
  expect_identical(rank_fail$status, "wald_covariance_rank_deficient")
  expect_true(is.na(rank_fail$stat))

  inconsistent <- wald_fn(
    fit,
    rbind(c(1, 0), c(2, 0)),
    rhs = c(0, 1)
  )
  expect_false(inconsistent$success)
  expect_identical(inconsistent$status, "wald_contrast_inconsistent")

  untagged <- fit
  untagged$coordinate_tag <- NULL
  tag_fail <- wald_fn(untagged, matrix(c(1, 0), nrow = 1))
  expect_false(tag_fail$success)
  expect_identical(tag_fail$status, "wald_coordinate_mismatch")

  constrained <- fit
  constrained$constraint <- affine_fn(
    C = matrix(c(1, 0), nrow = 1L), rhs = 0
  )
  constrained_fail <- wald_fn(constrained, matrix(c(1, 0), nrow = 1L))
  expect_false(constrained_fail$success)
  expect_identical(constrained_fail$status, "wald_requires_unrestricted_fit")

  uncertified <- fit
  uncertified$covariance_success <- FALSE
  uncertified_fail <- wald_fn(uncertified, matrix(c(1, 0), nrow = 1L))
  expect_false(uncertified_fail$success)
  expect_identical(
    uncertified_fail$status,
    "wald_unrestricted_covariance_not_certified"
  )
})

hl_score_synthetic_state <- function(R = matrix(c(1, 0), nrow = 1)) {
  p <- 2L
  q <- 1L
  d <- 2L * p + q
  C_R <- t(R) %*% solve(tcrossprod(R))
  K <- matrix(c(0, 1), nrow = 2)
  Q <- cbind(C_R, K)
  T <- matrix(0, d, d)
  T[1:2, 1:2] <- Q
  T[3:4, 3:4] <- Q
  T[5, 5] <- 1

  A_transformed <- rbind(
    c(2.0, 0.5, -0.2, 0.3, 0.1),
    c(-0.1, 2.0, 0.2, 0.0, 0.1),
    c(0.4, -0.3, 1.5, 0.2, 0.0),
    c(0.0, 0.1, -0.2, 1.2, 0.3),
    c(0.2, 0.0, 0.1, -0.1, 1.8)
  )
  u_subject <- rbind(
    c(0.6, -0.2, 0.4, 0.1),
    c(0.3, -0.1, -0.4, 0.2),
    c(-0.2, 0.5, 0.1, -0.4),
    c(0.4, -0.3, 0.2, -0.3),
    c(-0.1, 0.2, -0.2, 0.1)
  )
  # The nuisance equations are certified roots; the tested row is unsolved.
  u_subject[2:5, 4] <- -rowSums(u_subject[2:5, 1:3, drop = FALSE])

  list(
    psi_subject = solve(t(T), u_subject),
    A_zeta = solve(t(T), A_transformed) %*% solve(T),
    beta = c(0, 0.4),
    rhs = 0,
    q = q,
    transformed_A = A_transformed,
    transformed_u = u_subject
  )
}

test_that("full-chi Score post-root algebra uses adjusted subject meat directly", {
  score_fn <- getFromNamespace("hl_fo_score_postroot", "scRoPE")
  R <- matrix(c(1, 0), nrow = 1)
  fx <- hl_score_synthetic_state(R)

  out <- score_fn(
    psi_subject = fx$psi_subject,
    A_zeta = fx$A_zeta,
    R = R,
    beta = fx$beta,
    rhs = fx$rhs,
    q = fx$q
  )

  nuisance <- 2:5
  H_expected <- fx$transformed_A[1, nuisance, drop = FALSE] %*%
    solve(fx$transformed_A[nuisance, nuisance, drop = FALSE])
  q_expected <- fx$transformed_u[1, , drop = FALSE] -
    H_expected %*% fx$transformed_u[nuisance, , drop = FALSE]
  meat_expected <- tcrossprod(q_expected)
  stat_expected <- as.numeric(rowSums(q_expected)^2 / meat_expected)

  expect_true(out$success)
  expect_identical(out$status, "ok")
  expect_identical(out$method, "score_ee")
  expect_equal(out$H, H_expected, tolerance = 1e-12)
  expect_equal(out$q_subject, q_expected, tolerance = 1e-12)
  expect_equal(out$adjusted_meat, meat_expected, tolerance = 1e-12)
  expect_equal(out$adjusted_meat, out$adjusted_meat_block, tolerance = 1e-12)
  expect_lte(out$meat_crosscheck, 1e-12)
  expect_equal(out$stat, stat_expected, tolerance = 1e-12)
  expect_identical(out$contrast_rank, 1L)
  expect_identical(out$covariance_rank, 1L)
  expect_lt(out$J_rank, nrow(fx$A_zeta))
})

test_that("full-chi Score statistic is invariant to contrast scaling", {
  score_fn <- getFromNamespace("hl_fo_score_postroot", "scRoPE")
  R1 <- matrix(c(1, 0), nrow = 1)
  R2 <- 2 * R1
  fx1 <- hl_score_synthetic_state(R1)
  fx2 <- hl_score_synthetic_state(R2)

  out1 <- score_fn(fx1$psi_subject, fx1$A_zeta, R1, fx1$beta, 0, fx1$q)
  out2 <- score_fn(fx2$psi_subject, fx2$A_zeta, R2, fx2$beta, 0, fx2$q)

  expect_true(out1$success)
  expect_true(out2$success)
  expect_equal(out1$stat, out2$stat, tolerance = 1e-12)
  expect_equal(out1$p, out2$p, tolerance = 1e-12)
})

test_that("full-chi Score post-root algebra handles q = 0 and k = 0", {
  score_fn <- getFromNamespace("hl_fo_score_postroot", "scRoPE")
  psi_subject <- rbind(
    c(0.4, -0.1, 0.2),
    c(0.3, -0.2, -0.1)
  )
  A_zeta <- matrix(c(2, 0.3, 0.2, 1.5), nrow = 2, byrow = TRUE)

  out <- score_fn(
    psi_subject = psi_subject,
    A_zeta = A_zeta,
    R = 1,
    beta = 0,
    rhs = 0,
    q = 0
  )
  H_expected <- matrix(0.3 / 1.5, nrow = 1)
  q_expected <- psi_subject[1, , drop = FALSE] -
    H_expected %*% psi_subject[2, , drop = FALSE]

  expect_true(out$success)
  expect_identical(out$status, "ok")
  expect_identical(dim(out$K), c(1L, 0L))
  expect_identical(out$nuisance_indices, 2L)
  expect_equal(out$H, H_expected, tolerance = 1e-12)
  expect_equal(out$q_subject, q_expected, tolerance = 1e-12)
  expect_identical(out$adjusted_sensitivity_rank, 1L)
  expect_identical(out$adjusted_meat_rank, 1L)
})

test_that("full-chi Score reports singular nuisance and adjusted meat failures", {
  score_fn <- getFromNamespace("hl_fo_score_postroot", "scRoPE")
  R <- matrix(c(1, 0), nrow = 1)
  fx <- hl_score_synthetic_state(R)

  A_bad <- fx$A_zeta
  A_bad[5, ] <- 0
  nuisance_fail <- score_fn(fx$psi_subject, A_bad, R, fx$beta, 0, fx$q)
  expect_false(nuisance_fail$success)
  expect_identical(nuisance_fail$status, "score_nuisance_singular")

  psi_bad <- matrix(0, nrow = nrow(fx$psi_subject), ncol = ncol(fx$psi_subject))
  meat_fail <- score_fn(psi_bad, fx$A_zeta, R, fx$beta, 0, fx$q)
  expect_false(meat_fail$success)
  expect_identical(meat_fail$status, "score_meat_rank_deficient")
})

test_that("objective-ratio audit remains unsupported when Wald is finite", {
  unsupported_fn <- getFromNamespace(
    "hl_fo_objective_ratio_unsupported",
    "scRoPE"
  )
  wald <- list(method = "fallback", stat = 3.2, p = 0.07, success = TRUE)
  out <- unsupported_fn(
    raw_difference = -2.5,
    failure_reason = "negative_raw_objective_difference",
    wald = wald
  )

  expect_false(out$supported)
  expect_false(out$success)
  expect_identical(out$status, "objective_ratio_unsupported")
  expect_identical(out$failure_reason, "negative_raw_objective_difference")
  expect_equal(out$audit$raw_difference, -2.5)
  expect_identical(out$wald$method, "wald")
  expect_true(out$wald$success)
})
