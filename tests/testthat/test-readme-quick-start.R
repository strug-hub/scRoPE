test_that("README quick-start PGMM example completes with adjusted tests", {
  data(sample_data)

  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

  expect_no_warning(
    fit <- scrope(
      sample_data$count,
      sample_data$sid,
      pred = design,
      model = "PGMM",
      ncore = 1,
      additional_tests = c("score", "lrt"),
      lrt_details = TRUE
    )
  )

  expect_equal(nrow(fit$summary), nrow(sample_data$count))
  expect_equal(nrow(fit$lrt_details), nrow(sample_data$count) * 3L)
  expect_false(any(fit$lrt_details$fallback))
  expect_true(all(fit$lrt_details$method %in% c("lrt", "lrt_ridge")))
  expect_true(all(c(
    "score_adj_X1",
    "lrt_adj_X1",
    "score_adj_X2",
    "lrt_adj_X2",
    "score_adj_cccontrol",
    "lrt_adj_cccontrol"
  ) %in% names(fit$summary)))
})

test_that("README explicit joint-LN exposes Wald and calibration-pending tests", {
  data(sample_data)
  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

  fit <- scrope(
    sample_data$count,
    sample_data$sid,
    pred = design,
    ncore = 1,
    verbose = FALSE,
    additional_tests = c("score", "lrt"),
    lrt_details = TRUE,
    ln_estimator = "joint_ln"
  )

  expect_equal(nrow(fit$summary), nrow(sample_data$count))
  expect_true(all(fit$provenance$per_gene$route %in% c("R0", "R2")))
  r0 <- fit$provenance$per_gene$route == "R0"
  robust_success <- fit$summary$wald_success_X1
  expect_true(all(robust_success[r0]))
  expect_true(all(!robust_success[!r0]))
  expect_true(all(is.na(fit$lrt_details$stat)))
  expect_true(all(is.na(fit$lrt_details$p_value)))
  expect_true(all(is.na(fit$lrt_details$score_stat)))
  expect_true(all(is.na(fit$lrt_details$score_p_value)))
  expect_true(all(!fit$lrt_details$success))
  expect_true(all(!fit$lrt_details$score_success))
  expect_true(all(fit$lrt_details$method == "unavailable_joint_ln"))
  expect_true(all(
    fit$lrt_details$failure_reason == "joint_ln_calibration_pending"
  ))
  expect_false(any(fit$lrt_details$fallback))
})

test_that("README PGMM examples support subject- and cell-level predictors", {
  data(sample_data)

  design_subject <- model.matrix(~ X2, data = sample_data$pred)
  expect_no_warning(
    fit_pmg_subject <- scrope(
      sample_data$count,
      sample_data$sid,
      pred = design_subject,
      model = "PGMM",
      ncore = 1
    )
  )

  design_cell <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
  expect_message(
    fit_pmg_cell <- scrope(
      sample_data$count,
      sample_data$sid,
      pred = design_cell,
      model = "PGMM",
      ncore = 1
    ),
    "cell-level predictor\\(s\\) detected"
  )

  expect_equal(nrow(fit_pmg_subject$summary), nrow(sample_data$count))
  expect_equal(nrow(fit_pmg_cell$summary), nrow(sample_data$count))
  expect_true(all(c(
    "logFC_X1",
    "se_robust_X1",
    "p_robust_X1"
  ) %in% names(fit_pmg_cell$summary)))
})
