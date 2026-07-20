test_that("public HL sample-data replay uses only the corrected first-order path", {
  data(sample_data)
  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

  fit <- scrope_hl(
    sample_data$count[1L, , drop = FALSE],
    sample_data$sid,
    pred = design,
    ncore = 1L,
    verbose = FALSE,
    additional_tests = c("score", "lrt"),
    lrt_details = TRUE,
    keep_diagnostics = TRUE
  )

  expect_identical(nrow(fit$summary), 1L)
  expect_identical(fit$summary$gene, "A")
  expect_true(fit$summary$fit_success)
  expect_true(fit$summary$covariance_success)
  expect_identical(fit$summary$fit_status, "ok")
  expect_identical(fit$summary$covariance_status, "ok")
  expect_equal(fit$summary$`logFC_(Intercept)`, -1.940801, tolerance = 2e-6)
  expect_equal(fit$summary$logFC_X1, -0.01672225, tolerance = 1e-6)
  expect_equal(fit$summary$se_robust_X1, 0.0315635, tolerance = 1e-6)
  expect_equal(fit$overdispersion$Subject, 0.08432375, tolerance = 2e-7)
  expect_equal(fit$overdispersion$Cell, 0.9284692, tolerance = 2e-7)

  score_columns <- grep("^hl_score_success_", names(fit$summary), value = TRUE)
  lrt_columns <- grep("^hl_lrt_success_", names(fit$summary), value = TRUE)
  expect_length(score_columns, 3L)
  expect_length(lrt_columns, 3L)
  expect_true(all(vapply(fit$summary[score_columns], isTRUE, logical(1))))
  expect_false(any(vapply(fit$summary[lrt_columns], isTRUE, logical(1))))
  expect_true(all(grepl(
    "objective_ratio_unsupported",
    unlist(fit$summary[grep("^hl_lrt_failure_", names(fit$summary))])
  )))
  expect_true(all(fit$hl_adj_details$test[fit$hl_adj_details$test ==
    "objective_ratio"] == "objective_ratio"))
  expect_false(any(fit$hl_adj_details$success[
    fit$hl_adj_details$test == "objective_ratio"
  ]))

  state <- fit$diagnostics$first_order[[1L]]
  expect_identical(state$coordinate_tag,
                   "canonical_conditional_beta_tau_xi")
  expect_identical(state$active, 1:2)
  expect_true(state$point_success)
  expect_true(state$covariance_success)
  expect_true(is.matrix(state$augmented$A))
  expect_true(is.matrix(state$augmented$J))
  expect_false(isTRUE(state$covariance$sandwich$ridge > 0))
  expect_false("ord" %in% names(state))
  expect_true(all(vapply(
    fit$diagnostics$hl_adj[[1L]],
    function(test) {
      identical(test$objective_ratio$method, "unsupported_objective_ratio") &&
        !isTRUE(test$objective_ratio$success) &&
        identical(test$score$null_fit$object_type, "full_chi_score_null") &&
        identical(test$score$constrained_fit$object_type,
                  "reduced_constrained_maximizer")
    },
    logical(1)
  )))

  expect_match(paste(deparse(body(scrope_hl)), collapse = "\n"),
               "hl_fo_public_scrope_hl")
  expect_false(grepl(
    "ptmg|ord[[:space:]]*=[[:space:]]*3|fit_gene_hl|hl_adjusted_tests",
    paste(deparse(body(scrope_hl)), collapse = "\n"),
    ignore.case = TRUE
  ))
})
