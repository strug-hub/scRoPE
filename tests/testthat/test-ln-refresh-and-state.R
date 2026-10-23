test_that("lower-level LN dispatch fails closed outside an explicit joint context", {
  unconstrained <- getFromNamespace("fit_gene_unconstrained", "scRoPE")
  constrained <- getFromNamespace("fit_gene_constrained", "scRoPE")
  for (ctx in list(list(), list(ln_estimator = NULL),
                   list(ln_estimator = "legacy_post_pml"),
                   list(ln_estimator = "unknown"))) {
    expect_error(unconstrained(1L, NULL, ctx), "Legacy-LN dispatch is retired")
    expect_error(constrained(1L, NULL, ctx, NULL, NULL),
                 "Legacy-LN dispatch is retired")
  }
  expect_error(
    getFromNamespace("ln_operational_attempt_gene", "scRoPE")(1L, NULL, list()),
    "legacy-LN operational fallback is retired"
  )
})

test_that("maintained generic dispatch forwards only to joint LN fitters", {
  local_mocked_bindings(
    fit_gene_joint_ln_unconstrained = function(gene_index, posv, ctx, extra_starts) {
      list(kind = "joint", gene = gene_index, extra_starts = extra_starts)
    },
    fit_gene_joint_ln_constrained = function(gene_index, posv, ctx, L, b) {
      list(kind = "joint_constrained", gene = gene_index, L = L, b = b)
    },
    .package = "scRoPE"
  )
  ctx <- list(ln_estimator = "joint_ln")
  extra <- list(c(1, 2, 3))
  expect_identical(
    getFromNamespace("fit_gene_unconstrained", "scRoPE")(2L, NULL, ctx, extra),
    list(kind = "joint", gene = 2L, extra_starts = extra)
  )
  L <- matrix(c(0, 1), nrow = 1L)
  expect_identical(
    getFromNamespace("fit_gene_constrained", "scRoPE")(3L, NULL, ctx, L, 0.5),
    list(kind = "joint_constrained", gene = 3L, L = L, b = 0.5)
  )
})
