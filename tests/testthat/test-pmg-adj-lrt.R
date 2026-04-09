test_that("PMG adjusted LRT is null at the exact unconstrained estimate", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(1, 0), nrow = 1)
  res <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = list(list(L = L, b = fit_u$theta$beta[1], label = "intercept")),
    unconstrained = fit_u
  )[[1]]

  expect_equal(res$adj$lr$wP, 0, tolerance = 1e-8)
  expect_equal(res$adj$lr$wP_adjusted, 0, tolerance = 1e-6)
  expect_equal(as.numeric(res$adj$lr$score), 0, tolerance = 1e-6)
  expect_true(isTRUE(res$adj$lr$success))
  expect_true(isTRUE(res$adj$lr$adjusted_available))
  expect_false(isTRUE(res$adj$lr$fallback))
  expect_equal(dim(res$adj$lr$Hpsi), c(1, 1))
  expect_equal(dim(res$adj$lr$Gpsi), c(1, 1))
  expect_equal(dim(res$adj$lr$Jpsi), c(1, 1))
  expect_true(as.numeric(res$adj$lr$Hpsi) > 0)
  expect_true(as.numeric(res$adj$lr$Gpsi) > 0)
})

test_that("PMG adjusted LRT is available off-null when the scale gate passes", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")
  direct_lrt_fn <- getFromNamespace("compute_profile_lrt", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(0, 1), nrow = 1)
  b <- 0

  res <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = list(list(L = L, b = b, label = "slope")),
    unconstrained = fit_u
  )[[1]]

  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = b,
    unconstrained_fit = fit_u
  )
  direct_lrt <- direct_lrt_fn(fit_u, fit_c, L, b)

  expect_gt(res$adj$lr$wP, 0)
  expect_gt(res$adj$lr$wP_adjusted, 0)
  expect_true(isTRUE(res$adj$lr$adjusted_available))
  expect_false(isTRUE(res$adj$lr$fallback))
  expect_true(is.finite(res$p_values$lr))
  expect_equal(res$adj$lr$wP, direct_lrt$wP, tolerance = 1e-8)
  expect_equal(res$adj$lr$wP_adjusted, direct_lrt$wP_adjusted, tolerance = 1e-8)
  expect_equal(res$adj$lr$Gpsi, direct_lrt$Gpsi, tolerance = 1e-8)
  expect_equal(res$adj$lr$Hpsi, direct_lrt$Hpsi, tolerance = 1e-8)
  expect_equal(res$adj$lr$Jpsi, direct_lrt$Jpsi, tolerance = 1e-8)
  expect_equal(res$adj$lr$Vpsi_model, direct_lrt$Vpsi_model, tolerance = 1e-8)
  expect_equal(res$adj$lr$Vpsi_godambe, direct_lrt$Vpsi_godambe, tolerance = 1e-8)
})

test_that("PMG adjusted LRT gates out nonpositive curvature and falls back to Wald", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")
  direct_lrt_fn <- getFromNamespace("compute_profile_lrt", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(1, 0), nrow = 1)
  b <- 0

  res <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = list(list(L = L, b = b, label = "intercept")),
    unconstrained = fit_u
  )[[1]]

  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = b,
    unconstrained_fit = fit_u
  )
  direct_lrt <- direct_lrt_fn(fit_u, fit_c, L, b)

  expect_false(isTRUE(res$adj$lr$adjusted_available))
  expect_true(isTRUE(res$adj$lr$fallback))
  expect_match(res$adj$lr$method, "^wald")
  expect_equal(res$adj$lr$adjusted_failure_reason, "nonpositive_curvature")
  expect_false(res$adj$lr$gate$profiled_curvature_positive)
  expect_false(res$adj$lr$gate$scale_positive)
  expect_gt(res$adj$lr$wP_adjusted, 0)
  expect_equal(res$adj$lr$wP_adjusted, direct_lrt$wP_adjusted, tolerance = 1e-8)
})
