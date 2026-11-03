test_that("automatic covariance coordinates agree across native and recovered rows", {
  fx <- readRDS(test_path("fixtures", "jln-research.rds"))
  x <- fx$ordinary$input
  X <- x$X
  X[, 2L] <- 100 * X[, 2L]
  X[, 3L] <- 0.02 * X[, 3L]
  count <- rbind(native = x$y, recovered = fx$zero$input$y)
  fit <- suppressWarnings(scrope(count, x$id, X, offset = x$offset,
    cpc = 0, mincp = 1, verbose = FALSE))
  strict <- suppressWarnings(scrope(count, x$id, X, offset = x$offset,
    cpc = 0, mincp = 1, verbose = FALSE, ln_recovery = "none"))
  expect_identical(fit$recovery$routes$route, c("ordinary_jln", "fixed_v"))
  # Fresh native fits have different objective-instance IDs and RNG receipts.
  # Their public results/provenance and native covariance must be preserved.
  fields <- setdiff(names(strict), "diagnostics")
  expect_identical(fit$recovery$primary[fields], strict[fields])
  expect_identical(fit$recovery$primary$diagnostics$raw[[1L]]$robust_cov,
    strict$diagnostics$raw[[1L]]$robust_cov)
  expect_null(fit$recovery$primary$diagnostics$raw[[1L]]$robust_cov_public)
  expect_identical(fit$summary[1L, ], strict$summary[1L, ])
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
  native <- fit$recovery$primary$provenance$per_gene
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
  expect_null(fit$recovery$primary$diagnostics)
  expect_true(all(fit$recovery$coefficients$available))
  wald <- fit$test_status$test == "wald"
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
  expect_identical(fit$summary, compact$summary)
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
  expect_true(all(fit$test_status$failure_stage[!keep] == "coefficient_estimability"))
  expect_true(all(is.na(fit$test_status$failure_stage[keep])))
})
