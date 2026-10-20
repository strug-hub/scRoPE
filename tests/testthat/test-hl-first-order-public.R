test_that("public HL synthetic example uses the certified first-order path", {
  data(sample_data)
  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

  fit <- scrope_hl(
    sample_data$count[1L, , drop = FALSE],
    sample_data$sid,
    pred = design,
    offset = sample_data$offset,
    ncore = 1L,
    verbose = FALSE,
    additional_tests = c("score", "lrt"),
    lrt_details = TRUE,
    keep_diagnostics = TRUE
  )

  expect_identical(nrow(fit$summary), 1L)
  expect_identical(fit$summary$gene, rownames(sample_data$count)[1L])
  expect_true(fit$summary$fit_success)
  expect_true(fit$summary$covariance_success)
  expect_identical(fit$summary$fit_status, "ok")
  expect_identical(fit$summary$covariance_status, "ok")

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

  # Reconstruct the original-design map directly from the supplied matrix,
  # without using the package's centering or public-output helpers.
  predictor_names <- colnames(design)
  p <- ncol(design)
  intercept <- which(predictor_names == "(Intercept)")
  varying <- setdiff(seq_len(p), intercept)
  means <- colMeans(design)
  scales <- sqrt(colMeans(sweep(design, 2L, means, "-")^2))
  map <- diag(p)
  map[cbind(varying, varying)] <- 1 / scales[varying]
  map[intercept, varying] <- -means[varying] / scales[varying]
  expect_equal(unname(state$design_transform$M_X), map, tolerance = 1e-12)
  expected_beta <- as.numeric(map %*% state$beta)

  # A direct solve of the donor estimating contributions supplies an
  # independent covariance reconstruction; no sandwich helper is called.
  influence <- solve(state$augmented$A, state$augmented$psi_subject)
  expected_augmented_covariance <- tcrossprod(influence)
  free_beta <- seq_len(ncol(state$K))
  expected_internal_covariance <- state$K %*%
    expected_augmented_covariance[free_beta, free_beta, drop = FALSE] %*%
    t(state$K)
  expected_covariance <- map %*% expected_internal_covariance %*% t(map)
  expected_se <- sqrt(diag(expected_covariance))
  expected_wald <- expected_beta^2 / diag(expected_covariance)
  expected_p <- pchisq(expected_wald, df = 1L, lower.tail = FALSE)
  expect_equal(state$augmented$J,
               tcrossprod(state$augmented$psi_subject), tolerance = 1e-10)
  expect_equal(unname(state$covariance$V_beta),
               unname(expected_internal_covariance), tolerance = 1e-8)
  for (j in seq_len(p)) {
    nm <- predictor_names[j]
    expect_equal(fit$summary[[paste0("logFC_", nm)]], expected_beta[j],
                 tolerance = 1e-10)
    expect_equal(fit$summary[[paste0("se_robust_", nm)]], expected_se[j],
                 tolerance = 1e-8)
    expect_equal(fit$summary[[paste0("wald_stat_", nm)]], expected_wald[j],
                 tolerance = 1e-8)
    expect_equal(fit$summary[[paste0("p_robust_", nm)]], expected_p[j],
                 tolerance = 1e-8)
    expect_true(fit$summary[[paste0("wald_success_", nm)]])
    expect_true(is.na(fit$summary[[paste0("se_", nm)]]))
    expect_true(is.na(fit$summary[[paste0("p_naive_", nm)]]))
  }
  expect_equal(fit$overdispersion$Subject,
               unname(exp(state$rho_full[1L])), tolerance = 1e-10)
  expect_equal(fit$overdispersion$Cell,
               unname(exp(-state$rho_full[2L])), tolerance = 1e-10)
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
      offset = sample_data$offset,
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
      offset = sample_data$offset,
      opt = "unsupported",
      verbose = FALSE
    ),
    "opt must be 'lbfgs' or 'trust'",
    fixed = TRUE
  )
})
