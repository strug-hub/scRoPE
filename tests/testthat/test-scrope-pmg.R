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
    pgmm_control = list(iter.max = 150L, eval.max = 200L),
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
  expect_true(all(c(
    "convergence_basis",
    "gate_optimizer_converged",
    "gate_information_scaled_convergence_ok",
    "gate_H_red_positive_definite",
    "max_abs_newton_step_scaled",
    "gradient_quadratic_lr_bound",
    "gate_constraint_satisfied",
    "retry_triggered",
    "retry_used",
    "optimizer_attempts"
  ) %in% names(res$lrt_details)))
  expect_true(all(res$lrt_details$convergence_basis %in% c(
    "optimizer", "raw_gradient", "information_scaled"
  )))
})

test_that("top-level scrope PGMM matches frozen legacy PMM values for cell-level predictors", {
  fx <- pmg_interior_test_fixture()

  current <- scrope(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    model = "PGMM",
    pgmm_control = list(iter.max = 150L, eval.max = 200L),
    verbose = FALSE,
    cpc = 0,
    mincp = 1
  )
  legacy <- pmg_frozen_legacy("cell")
  receipt <- pmg_cell_endpoint_receipt(current$diagnostics$raw[[1L]], fx, legacy)
  pmg_write_cell_endpoint_receipt(receipt, current$summary)

  expect_equal(current$summary$logFC_x, legacy$summary$logFC_x, tolerance = 1e-6)
  expect_equal(current$summary$se_x, legacy$summary$se_x, tolerance = 1e-6)
  expect_equal(current$summary$p_naive_x, legacy$summary$p_x, tolerance = 1e-6)
  expect_identical(as.numeric(current$overdispersion$Subject), receipt$theta[3L])
  # A frozen optimizer endpoint is not an exact MLE: its sigma differs from
  # the independent profile root by 5.16e-8. Require both endpoints' measured
  # accuracy instead of exact cross-platform reproduction of that residual.
  expect_true(pmg_cell_endpoint_is_precise(receipt), info = paste(
    "PGMM endpoint failed the independent root/curvature accuracy certificate;",
    "inspect _diagnostics/pgmm-cell-endpoint.txt (or SCROPE_TEST_DIAGNOSTICS_DIR)."
  ))
  expect_lt(abs(receipt$reference$shape_score), 1e-11)
  expect_lt(abs(receipt$reference$slope_score), 1e-12)
  expect_equal(current$summary$se_x, receipt$reference$model_se, tolerance = 1e-6)
  expect_equal(current$summary$se_robust_x, receipt$reference$robust_se, tolerance = 1e-6)
  expect_equal(current$convergence, legacy$convergence)
  expect_equal(current$algorithm, legacy$algorithm)
})

test_that("PGMM endpoint accuracy checks reject material numerical errors", {
  fx <- pmg_interior_test_fixture()
  fit_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  eval_fn <- getFromNamespace("pmg_ll_der_hes", "scRoPE")
  raw <- fit_fn(fx$gene_index, fx$posv, fx$prep$ctx,
                control = list(iter.max = 150L, eval.max = 200L))
  legacy <- pmg_frozen_legacy("cell")
  expect_true(pmg_cell_endpoint_is_precise(pmg_cell_endpoint_receipt(raw, fx, legacy)))
  # A real 1e-4 change in sigma must not pass by inflating its Newton budget.
  perturbed <- raw
  perturbed$theta$subVar <- perturbed$theta$subVar + 1e-4
  ctx <- fx$prep$ctx
  evaluated <- eval_fn(
    c(perturbed$theta$beta, perturbed$theta$subVar), ctx$pred, ctx$offset,
    fx$posv$Y, ctx$fid, ctx$cumsumy[fx$gene_index, ],
    ctx$posind[[fx$gene_index]], fx$posv$posindy, ctx$nb, ctx$nind, ctx$k
  )
  perturbed$score <- evaluated$gradient
  perturbed$observed_info <- evaluated$observed_info
  expect_false(pmg_cell_endpoint_is_precise(pmg_cell_endpoint_receipt(perturbed, fx, legacy)))
  wrong_status <- pmg_cell_endpoint_receipt(raw, fx, legacy)
  wrong_status$optimizer_diagnostics$optimizer_success <- FALSE
  expect_false(pmg_cell_endpoint_is_precise(wrong_status))
})

test_that("top-level scrope PGMM matches frozen legacy PMM values for subject-level predictors", {
  fx <- pmg_subject_level_test_fixture()

  current <- scrope(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    model = "PGMM",
    pgmm_control = list(iter.max = 150L, eval.max = 200L),
    verbose = FALSE,
    cpc = 0,
    mincp = 1
  )
  legacy <- pmg_frozen_legacy("subject")

  expect_equal(current$summary$logFC_z_subject, legacy$summary$logFC_z_subject, tolerance = 1e-6)
  expect_equal(current$summary$se_z_subject, legacy$summary$se_z_subject, tolerance = 1e-6)
  expect_equal(current$summary$p_naive_z_subject, legacy$summary$p_z_subject, tolerance = 1e-6)
  expect_equal(current$overdispersion$Subject, as.numeric(legacy$overdispersion), tolerance = 1e-6)
  expect_equal(current$convergence, legacy$convergence)
  expect_equal(current$algorithm, legacy$algorithm)
})

test_that("scrope PGMM allows cell-level predictors and returns results", {
  fx <- pmg_interior_test_fixture()

  res <- scrope(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    model = "PGMM",
    verbose = FALSE,
    cpc = 0,
    mincp = 1
  )

  expect_equal(res$algorithm, "PGMM")
  expect_true(all(c("logFC_x", "se_x", "se_robust_x", "p_robust_x") %in% colnames(res$summary)))
  expect_equal(nrow(res$summary), 1)
})

test_that("scrope PGMM messages when predictors vary within subject", {
  fx <- pmg_interior_test_fixture()

  expect_message(
    scrope(
      count = fx$count,
      id = fx$id,
      pred = fx$pred,
      model = "PGMM",
      verbose = TRUE,
      cpc = 0,
      mincp = 1
    ),
    "cell-level predictor\\(s\\) detected"
  )
})
