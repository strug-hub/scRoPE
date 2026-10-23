test_that("README PGMM example reports adjusted-test availability and fallbacks", {
  data(sample_data)

  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
  warnings <- character(0)

  fit <- withCallingHandlers(
    scrope(
      sample_data$count,
      sample_data$sid,
      pred = design,
      offset = sample_data$offset,
      model = "PGMM",
      ncore = 1,
      additional_tests = c("score", "lrt"),
      lrt_details = TRUE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(nrow(fit$summary), nrow(sample_data$count))
  expect_equal(nrow(fit$lrt_details), nrow(sample_data$count) * 3L)
  details <- fit$lrt_details
  expected_keys <- outer(rownames(sample_data$count),
                         c("X1", "X2", "cccontrol"), paste, sep = "/")
  expect_setequal(paste(details$gene, details$contrast, sep = "/"),
                  as.vector(expected_keys))
  expect_false(anyDuplicated(details[c("gene", "contrast")]) > 0L)
  expect_identical(details$success, details$adjusted_available)
  accepted <- details$adjusted_available
  # This frozen, unselected teaching realization documents one curvature
  # rejection. Keep its fallback branch covered rather than permitting an
  # all-success result to bypass the fallback assertions below.
  expect_identical(sum(accepted), 23L)
  expect_identical(sum(details$fallback), 1L)
  curvature_rejection <- details[
    details$gene == "G08" & details$contrast == "X2", , drop = FALSE
  ]
  expect_identical(nrow(curvature_rejection), 1L)
  expect_false(curvature_rejection$adjusted_available)
  expect_false(curvature_rejection$success)
  expect_false(curvature_rejection$gate_profiled_curvature_positive)
  expect_identical(curvature_rejection$adjusted_failure_reason,
                   "nonpositive_curvature")
  expect_identical(curvature_rejection$failure_reason, "nonpositive_curvature")
  expect_identical(curvature_rejection$method, "wald")
  expect_true(curvature_rejection$fallback)
  expect_true(all(details$method[accepted] %in% c("lrt", "lrt_ridge")))
  expect_false(any(details$fallback[accepted]))
  expect_true(all(is.finite(details$stat[accepted])))
  expect_true(all(details$p_value[accepted] >= 0 &
                    details$p_value[accepted] <= 1))
  expect_true(all(is.na(details$adjusted_failure_reason[accepted])))
  expect_true(all(!is.na(details$adjusted_failure_reason[!accepted]) &
                    nzchar(details$adjusted_failure_reason[!accepted])))
  expect_true(all(details$score_success))
  expect_true(all(is.finite(details$score_stat)))

  # All requests remain in the output. A failed adjusted LRT can carry an
  # explicitly labeled Wald fallback; that value is not an accepted LRT.
  for (i in seq_len(nrow(details))) {
    detail <- details[i, ]
    row <- fit$summary[fit$summary$gene == detail$gene, , drop = FALSE]
    name <- detail$contrast
    expect_identical(row[[paste0("lrt_method_", name)]], detail$method)
    expect_identical(row[[paste0("lrt_success_", name)]], detail$success)
    expect_equal(row[[paste0("lrt_adj_", name)]], detail$stat)
    expect_equal(row[[paste0("p_lrt_adj_", name)]], detail$p_value)
    expect_identical(row[[paste0("score_success_", name)]],
                     detail$score_success)
    if (!detail$adjusted_available) {
      expect_identical(row[[paste0("lrt_adjusted_failure_", name)]],
                       detail$adjusted_failure_reason)
      if (detail$fallback) {
        expect_true(detail$method %in% c("wald", "wald_ridge"))
        expect_true(is.finite(detail$stat))
        expect_true(is.finite(detail$p_value))
        if (detail$method == "wald") {
          expected_wald <- (row[[paste0("logFC_", name)]] /
                              row[[paste0("se_robust_", name)]])^2
          expect_equal(detail$stat, expected_wald, tolerance = 1e-8)
          expect_equal(detail$p_value,
                       pchisq(expected_wald, df = 1L, lower.tail = FALSE),
                       tolerance = 1e-8)
        }
      } else {
        expect_true(is.na(detail$stat))
        expect_true(is.na(detail$p_value))
      }
    }
  }

  # Assert the complete warning stream rather than suppressing numerical
  # failures or treating a disclosed fallback as a data-generation failure.
  expected_warnings <- character(0)
  n_fallback <- sum(details$method %in% c("wald", "wald_ridge", "lrt_fail", "lrt_raw"))
  n_ridge <- sum(details$method %in% c("lrt_ridge", "wald_ridge"))
  if (n_fallback > 0L) expected_warnings <- c(expected_warnings, sprintf(
    "Adjusted LRT fell back or was gated out for %d gene-contrast(s); see lrt_method_* columns.",
    n_fallback
  ))
  if (n_ridge > 0L) expected_warnings <- c(expected_warnings, sprintf(
    "Adjusted LRT applied ridge regularisation for %d gene-contrast(s).", n_ridge
  ))
  expect_identical(warnings, expected_warnings)
  expect_true(all(c(
    "score_adj_X1",
    "lrt_adj_X1",
    "score_adj_X2",
    "lrt_adj_X2",
    "score_adj_cccontrol",
    "lrt_adj_cccontrol"
  ) %in% names(fit$summary)))
})

test_that("README explicit joint-LN exposes Wald and experimental optional tests", {
  data(sample_data)
  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

  fit <- scrope(
    sample_data$count,
    sample_data$sid,
    pred = design,
    offset = sample_data$offset,
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
  details <- fit$lrt_details
  expect_true(any(details$success))
  expect_true(any(details$score_success))
  for (test in c("lrt", "score")) {
    success <- if (test == "lrt") details$success else details$score_success
    statistic <- if (test == "lrt") details$stat else details$score_stat
    p_value <- if (test == "lrt") details$p_value else details$score_p_value
    failure <- if (test == "lrt") details$failure_reason else details$score_failure_reason
    expect_true(all(is.finite(statistic[success])))
    expect_equal(p_value[success],
                 pchisq(statistic[success], df = 1L, lower.tail = FALSE),
                 tolerance = 1e-14)
    expect_true(all(is.na(statistic[!success])))
    expect_true(all(is.na(p_value[!success])))
    expect_true(all(!is.na(failure[!success]) & nzchar(failure[!success])))
  }
  expect_true(all(details$method[details$success] ==
                    "joint_ln_pace_adjusted_ratio"))
  optional_status <- subset(fit$test_status, test %in% c("score", "lrt"))
  expect_true(all(optional_status$status[optional_status$success] ==
                    "experimental_asymptotic"))
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
      offset = sample_data$offset,
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
      offset = sample_data$offset,
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
