test_that("PMG adjusted score is null at the exact unconstrained estimate", {
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

  expect_equal(res$adj$wald$stat, 0, tolerance = 1e-8)
  expect_equal(res$adj$score$stat, 0, tolerance = 1e-6)
  expect_equal(as.numeric(res$adj$score$score), 0, tolerance = 1e-6)
  expect_true(isTRUE(res$adj$score$success))
  expect_equal(dim(res$adj$score$score), c(1, 1))
  expect_equal(dim(res$adj$score$Hpsi), c(1, 1))
  expect_equal(dim(res$adj$score$Gpsi), c(1, 1))
  expect_equal(dim(res$adj$score$Spsi), c(1, 1))
})

test_that("PMG adjusted score is positive off-null and matches the direct helper", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")
  direct_score_fn <- getFromNamespace("compute_adj_score", "scRoPE")

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
  direct_score <- direct_score_fn(fit_c, L)

  expect_gt(res$adj$score$stat, 0)
  expect_gt(res$adj$wald$stat, 0)
  expect_true(is.finite(res$p_values$score))
  expect_true(is.finite(res$p_values$wald))
  expect_equal(res$adj$score$stat, direct_score$stat, tolerance = 1e-8)
  expect_equal(res$adj$score$Gpsi, direct_score$Gpsi, tolerance = 1e-8)
  expect_equal(res$adj$score$Hpsi, direct_score$Hpsi, tolerance = 1e-8)
  expect_equal(res$adj$score$Spsi, direct_score$Spsi, tolerance = 1e-8)
})
