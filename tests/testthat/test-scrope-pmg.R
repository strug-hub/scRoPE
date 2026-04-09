test_that("scrope PGMM integrates exact PMG fits and adjusted tests", {
  fx <- pmg_subject_level_test_fixture()
  prepare_fn <- getFromNamespace("hl_prepare_data", "scRoPE")
  posindy_fn <- getFromNamespace("call_posindy", "scRoPE")
  fit_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  adj_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")

  res <- scrope(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    model = "PGMM",
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    additional_tests = c("score", "lrt"),
    lrt_details = TRUE
  )

  expect_equal(res$algorithm, "PGMM")
  expect_null(res$random_effect)
  expect_equal(colnames(res$overdispersion), "Subject")
  expect_false("Cell" %in% colnames(res$overdispersion))

  sm <- res$summary
  expect_true(all(c(
    "logFC_z_subject",
    "se_z_subject",
    "se_robust_z_subject",
    "p_robust_z_subject",
    "score_adj_z_subject",
    "lrt_adj_z_subject",
    "lrt_method_z_subject"
  ) %in% colnames(sm)))

  prep <- prepare_fn(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
  posv <- posindy_fn(prep$count, prep$gid[[1]] - 1L, prep$ctx$nind)
  fit_u <- fit_fn(
    gene_index = prep$gid[[1]],
    posv = posv,
    ctx = prep$ctx
  )
  adj <- adj_fn(
    ctx = prep$ctx,
    gene_index = prep$gid[[1]],
    posv = posv,
    contrasts = list(list(
      L = matrix(c(0, 1), nrow = 1),
      b = 0,
      label = "z_subject"
    )),
    unconstrained = fit_u
  )[[1]]

  slope_scale <- prep$design_sds[2]
  expect_equal(sm$logFC_z_subject[1], fit_u$theta$beta[2] / slope_scale, tolerance = 1e-8)
  expect_equal(sm$se_z_subject[1], sqrt(fit_u$naive_cov[2, 2]) / slope_scale, tolerance = 1e-8)
  expect_equal(sm$se_robust_z_subject[1], sqrt(fit_u$robust_cov[2, 2]) / slope_scale, tolerance = 1e-8)
  expect_equal(sm$p_robust_z_subject[1], adj$p_values$wald, tolerance = 1e-8)
  expect_equal(sm$score_adj_z_subject[1], adj$adj$score$stat, tolerance = 1e-8)
  expect_equal(sm$lrt_adj_z_subject[1], adj$adj$lr$wP_adjusted, tolerance = 1e-8)
  expect_equal(sm$lrt_method_z_subject[1], adj$adj$lr$method)
  expect_true(nrow(res$lrt_details) >= 1)
})

test_that("scrope PGMM rejects cell-level predictors", {
  fx <- pmg_interior_test_fixture()

  expect_error(
    scrope(
      count = fx$count,
      id = fx$id,
      pred = fx$pred,
      model = "PGMM",
      verbose = FALSE,
      cpc = 0,
      mincp = 1
    ),
    "subject-level predictors"
  )
})
