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

  expect_identical(fit$provenance$schema_version, 2L)
  expect_identical(fit$provenance$route, "HL")
  expect_identical(fit$provenance$requested_model, "HL")
  expect_identical(fit$provenance$effective_model, "HL")
  expect_identical(fit$provenance$estimator, "first_order_nested_hl")
  expect_identical(fit$provenance$model_covariance, "not_reported")
  expect_identical(
    fit$provenance$robust_covariance,
    "hl_first_order_augmented_donor_sandwich"
  )
  expect_identical(
    fit$provenance$test_capabilities,
    data.frame(
      test = c("wald", "score", "lrt"),
      public_available = c(TRUE, TRUE, FALSE),
      status = c("available", "available", "objective_ratio_unsupported"),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(fit$provenance$per_gene$gene_id, fit$summary$gene_id)
  expect_identical(fit$provenance$per_gene$gene, fit$summary$gene)
  expect_identical(
    fit$provenance$per_gene$fit_success,
    fit$summary$fit_success
  )
  expect_identical(
    fit$provenance$per_gene$covariance_success,
    fit$summary$covariance_success
  )

  expect_match(paste(deparse(body(scrope_hl)), collapse = "\n"),
               "hl_fo_public_scrope_hl")
  expect_false(grepl(
    "ptmg|ord[[:space:]]*=[[:space:]]*3|fit_gene_hl|hl_adjusted_tests",
    paste(deparse(body(scrope_hl)), collapse = "\n"),
    ignore.case = TRUE
  ))
})

test_that("standalone HL compatibility arguments warn when ignored", {
  data(sample_data)
  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
  warnings <- character(0)

  fit <- withCallingHandlers(
    scrope_hl(
      sample_data$count[1L, , drop = FALSE],
      sample_data$sid,
      pred = design,
      cpc = 0,
      mincp = 1,
      opt = "trust",
      ncore = 2L,
      fmaxsize = 1e8,
      use_betas = "final_betas",
      adj_ridge_factor = 1e-6,
      adj_cond_limit = 1e6,
      keep_diagnostics = FALSE,
      verbose = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("opt='trust' is ignored", warnings, fixed = TRUE)))
  expect_true(any(grepl("use_betas is ignored", warnings, fixed = TRUE)))
  expect_true(any(grepl("ncore > 1 is ignored", warnings, fixed = TRUE)))
  expect_true(any(grepl("fmaxsize is ignored", warnings, fixed = TRUE)))
  expect_true(any(grepl("adj_ridge_factor is ignored", warnings, fixed = TRUE)))
  expect_true(any(grepl("adj_cond_limit is ignored", warnings, fixed = TRUE)))
  expect_null(fit$diagnostics)
  expect_identical(fit$provenance$route, "HL")
  expect_identical(
    fit$provenance$test_capabilities$public_available,
    c(TRUE, TRUE, FALSE)
  )

  expect_error(
    scrope_hl(
      sample_data$count[1L, , drop = FALSE],
      sample_data$sid,
      pred = design,
      opt = "unsupported",
      verbose = FALSE
    ),
    "opt must be 'lbfgs' or 'trust'",
    fixed = TRUE
  )
})
