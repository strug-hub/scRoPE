test_that("automatic covariance coordinates agree across native and recovered rows", {
  fx <- readRDS(test_path("fixtures", "jln-research.rds"))
  x <- fx$ordinary$input
  X <- x$X
  X[, 2L] <- 100 * X[, 2L]
  X[, 3L] <- 0.02 * X[, 3L]
  count <- rbind(native = x$y, recovered = fx$zero$input$y)
  fit <- suppressWarnings(scrope(count, x$id, X, offset = x$offset,
    cpc = 0, mincp = 1, verbose = FALSE, keep_diagnostics = TRUE))
  strict <- suppressWarnings(scrope(count, x$id, X, offset = x$offset,
    cpc = 0, mincp = 1, verbose = FALSE, ln_recovery = "none", keep_diagnostics = TRUE))
  expect_identical(fit$recovery$routes$route, c("ordinary_jln", "fixed_v"))
  # Fresh native fits have different objective-instance IDs and RNG receipts.
  # Numerical results can have roundoff; all other public metadata stays exact.
  fields <- setdiff(names(strict), c("diagnostics", "summary", "overdispersion"))
  expect_identical(fit$recovery$primary[fields], strict[fields])
  expect_recovery_refit_table(fit$recovery$primary$summary, strict$summary)
  expect_recovery_refit_table(fit$recovery$primary$overdispersion, strict$overdispersion,
    numerical = c("Subject", "Cell"))
  expect_recovery_refit_matrix(fit$recovery$primary$diagnostics$raw[[1L]]$robust_cov,
    strict$diagnostics$raw[[1L]]$robust_cov)
  expect_null(fit$recovery$primary$diagnostics$raw[[1L]]$robust_cov_public)
  expect_recovery_refit_table(fit$summary[1L, ], strict$summary[1L, ])
  # Preservation within the same fit still requires bitwise identity.
  expect_identical(fit$summary[1L, ], fit$recovery$primary$summary[1L, ])
  expect_identical(fit$overdispersion[1L, ], fit$recovery$primary$overdispersion[1L, ])
  expect_identical(fit$diagnostics$raw[[1L]]$robust_cov,
    fit$recovery$primary$diagnostics$raw[[1L]]$robust_cov)
  M <- fit$provenance$coefficient_map_internal_to_public
  for (j in seq_len(nrow(count))) {
    raw <- fit$diagnostics$raw[[j]]
    reported <- unlist(fit$summary[j, paste0("se_robust_", colnames(X))])^2
    expect_identical(raw$robust_cov_coordinates, "standardized_design")
    expect_identical(raw$robust_cov_public_coordinates, "centered_public_design")
    expect_equal(raw$robust_cov_public, M %*% raw$robust_cov %*% t(M),
      tolerance = 1e-12)
    expect_equal(unname(diag(raw$robust_cov_public)), unname(reported),
      tolerance = 1e-12)
  }
})

test_that("recovery provenance describes the active fit and preserves native failures", {
  x <- readRDS(test_path("fixtures", "jln-research.rds"))$zero$input
  fit <- suppressWarnings(scrope(x$y, x$id, x$X, offset = x$offset,
    cpc = 0, mincp = 1, verbose = FALSE, keep_diagnostics = FALSE,
    additional_tests = c("score", "lrt"), lrt_details = TRUE, output_re = TRUE))
  active <- fit$provenance$per_gene
  strict <- suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,
    cpc=0,mincp=1,verbose=FALSE,recovery="none"))
  native <- strict$provenance$per_gene
  expect_identical(fit$recovery$native_status$enabled_failure_reason,native$enabled_failure_reason)
  expect_identical(fit$availability$native_convergence,strict$convergence)
  expect_identical(active$enabled_route_code, "fixed_v")
  expect_true(active$inference_available)
  expect_true(is.na(active$enabled_failure_reason))
  expect_true(is.na(active$resolution_veto_classes))
  expect_identical(active$resolution_mode, "automatic_recovery")
  expect_false(is.na(native$enabled_failure_reason))
  expect_false(is.na(native$resolution_veto_classes))
  strict_fields <- grep("^strict_", names(active), value = TRUE)
  expect_identical(active[strict_fields], native[strict_fields])
  expect_null(fit$diagnostics)
  expect_null(fit$recovery$diagnostics)
  expect_null(fit$recovery$primary)
  expect_true(all(fit$recovery$coefficients$available))
  wald <- fit$test_status$test %in% c("wald","wald_model")
  expect_true(all(is.na(fit$test_status$failure_stage[wald])))
  expect_true(all(fit$test_status$failure_stage[!wald] == "recovery_test_scope"))
  expect_true(all(!fit$test_status$success[!wald]))
  expect_true(all(is.na(fit$summary[grep("^p_(score|lrt)_", names(fit$summary))])))
  expect_true(all(is.na(fit$random_effect)))
  expect_true(all(is.na(fit$lrt_details$p_value)))
})

test_that("partial recovery keeps covariance and status on the estimable coordinates", {
  x <- readRDS(test_path("fixtures", "jln-face.rds"))
  X <- x$X
  X[, 2L] <- 100 * X[, 2L]
  X[, 4L] <- 0.02 * X[, 4L]
  run <- function(keep) suppressWarnings(scrope(x$y, x$id, X, offset = x$offset,
    cpc = 0, mincp = 1, verbose = FALSE, keep_diagnostics = keep))
  fit <- run(TRUE)
  compact <- run(FALSE)
  expect_identical(fit$recovery$routes$route, "identifiable_nuisance_limit")
  expect_recovery_refit_table(fit$summary, compact$summary)
  expect_identical(fit$provenance$per_gene, compact$provenance$per_gene)
  expect_identical(compact$recovery$coefficients, fit$recovery$coefficients)
  expect_identical(compact$provenance$per_gene$enabled_failure_reason,
    "some_coefficients_nonestimable")
  expect_true(is.na(compact$provenance$per_gene$resolution_veto_classes))
  expect_identical(compact$availability$availability, "partial_inference")
  raw <- fit$diagnostics$raw[[1L]]
  keep <- fit$recovery$coefficients$available
  expect_identical(keep, c(FALSE, TRUE, FALSE, TRUE))
  expect_true(all(is.na(raw$robust_cov[!keep, ])))
  expect_true(all(is.na(raw$robust_cov[, !keep])))
  expect_true(all(is.finite(raw$robust_cov[keep, keep])))
  M <- fit$provenance$coefficient_map_internal_to_public[keep, keep, drop = FALSE]
  expect_equal(raw$robust_cov_public[keep, keep],
    M %*% raw$robust_cov[keep, keep] %*% t(M), tolerance = 1e-12)
  reported <- unlist(fit$summary[paste0("se_robust_", colnames(X)[keep])])^2
  expect_equal(unname(diag(raw$robust_cov_public)[keep]), unname(reported),
    tolerance = 1e-12)
  expect_true(all(fit$test_status$failure_stage[!keep] %in% c("coefficient_estimability","model_coefficient_inference")))
  expect_true(all(grepl("nonestimable",fit$test_status$failure_reason[!keep])))
  expect_true(all(is.na(fit$test_status$failure_stage[keep])))
})

test_that("refit comparison tolerates only fitted numerical roundoff", {
  baseline <- data.frame(gene_id = 1L, gene = "synthetic", route = "ordinary_jln",
    wald_success_x = TRUE, logFC_x = 0.2, se_robust_x = 0.1, p_robust_x = 0.05,
    logFC_region = NA_real_, fixed_donor_variance = 1e-4,
    finite_cluster_multiplier = 1, stringsAsFactors = FALSE)
  close <- baseline
  for (name in c("logFC_x", "se_robust_x", "p_robust_x"))
    close[[name]] <- close[[name]] + 1e-13
  expect_length(recovery_refit_table_differences(close, baseline), 0L)
  for (name in c("logFC_x", "se_robust_x", "p_robust_x")) {
    changed <- baseline
    changed[[name]] <- changed[[name]] + 1e-6
    expect_gt(length(recovery_refit_table_differences(changed, baseline)), 0L)
  }
  for (name in c("gene_id", "gene", "route", "wald_success_x",
                 "fixed_donor_variance", "finite_cluster_multiplier")) {
    changed <- baseline
    changed[[name]] <- switch(name, gene_id = 2L, gene = "other", route = "withheld",
      wald_success_x = FALSE, baseline[[name]] + 1e-13)
    expect_gt(length(recovery_refit_table_differences(changed, baseline)), 0L)
  }
  changed <- baseline
  changed$gene_id <- as.double(changed$gene_id)
  expect_gt(length(recovery_refit_table_differences(changed, baseline)), 0L)
  changed <- baseline
  changed$logFC_region <- 0
  expect_gt(length(recovery_refit_table_differences(changed, baseline)), 0L)
  changed$logFC_region <- NaN
  expect_gt(length(recovery_refit_table_differences(changed, baseline)), 0L)
  changed <- baseline
  rownames(changed) <- "other"
  expect_gt(length(recovery_refit_table_differences(changed, baseline)), 0L)
  expect_gt(length(recovery_refit_table_differences(baseline[rev(names(baseline))],
    baseline)), 0L)
  expect_gt(length(recovery_refit_table_differences(baseline[-1L], baseline)), 0L)

  research <- data.frame(gene_id = 1L, route = "fixed_v", coefficient = 0.2,
    se_robust = 0.1, p_robust = 0.05, ci_lower = 0, ci_upper = 0.4,
    fixed_donor_variance = 1e-4, finite_cluster_multiplier = 1)
  changed <- research
  changed$ci_lower <- 1e-13
  expect_length(recovery_refit_table_differences(changed, research), 0L)
  changed$ci_lower <- 1e-6
  expect_gt(length(recovery_refit_table_differences(changed, research)), 0L)

  covariance <- matrix(c(1, 0.1, 0.1, 2), 2,
    dimnames = list(c("x", "z"), c("x", "z")))
  expect_length(recovery_refit_numeric_differences(covariance + 1e-13, covariance), 0L)
  expect_gt(length(recovery_refit_numeric_differences(covariance + 1e-6,
    covariance)), 0L)
  changed <- covariance
  dimnames(changed) <- NULL
  expect_gt(length(recovery_refit_numeric_differences(changed, covariance)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(0, 0L)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(Inf, -Inf)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(NA_real_, NaN)), 0L)
})
