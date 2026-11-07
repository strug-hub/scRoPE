test_that("refit comparisons use each element's scale", {
  baseline <- c(large = 1e6, ordinary = 0.1, tiny = 1e-16, zero = 0)
  expect_length(recovery_refit_numeric_differences(
    baseline * (1 + 1e-12), baseline), 0L)
  for (i in seq_len(3L)) {
    changed <- baseline
    changed[i] <- changed[i] * 1.001
    expect_gt(length(recovery_refit_numeric_differences(changed, baseline)), 0L)
  }
  changed <- baseline
  changed[4L] <- 1e-30
  expect_gt(length(recovery_refit_numeric_differences(changed, baseline)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(1e308, -1e308)), 0L)
  expect_length(recovery_refit_numeric_differences(numeric(0), numeric(0)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(1, c(1, 1))), 0L)
})

test_that("only signed coefficient and interval values have a zero-scale floor", {
  baseline <- data.frame(beta = 0, ci_lower = 0, ci_upper = 0,
    robust_se = 1e-16, p_value = 1e-100, statistic = 1e-20,
    status = "available", fixed_donor_variance = 1e-4)
  numerical <- c("beta", "ci_lower", "ci_upper", "robust_se", "p_value", "statistic")
  differences <- function(x) recovery_refit_table_differences(x, baseline, numerical)
  changed <- baseline
  changed[c("beta", "ci_lower", "ci_upper")] <- 1e-13
  expect_length(differences(changed), 0L)
  for (name in c("beta", "ci_lower", "ci_upper")) {
    changed <- baseline
    changed[[name]] <- 1e-9
    expect_gt(length(differences(changed)), 0L)
  }
  for (name in c("robust_se", "p_value", "statistic")) {
    changed <- baseline
    changed[[name]] <- changed[[name]] * (1 + 1e-12)
    expect_length(differences(changed), 0L)
    changed[[name]] <- baseline[[name]] * 1.01
    expect_gt(length(differences(changed)), 0L)
    changed[[name]] <- 0
    expect_gt(length(differences(changed)), 0L)
  }
  changed <- baseline
  changed$fixed_donor_variance <- changed$fixed_donor_variance + 1e-16
  expect_gt(length(differences(changed)), 0L)
  changed <- baseline
  changed$status <- "withheld"
  expect_gt(length(differences(changed)), 0L)
})

test_that("refit comparison preserves schema and all nonfinite states", {
  baseline <- c(a = 0.1, b = NA_real_, c = NaN, d = Inf, e = -Inf)
  expect_length(recovery_refit_numeric_differences(baseline, baseline), 0L)
  for (i in seq_along(baseline)) {
    changed <- baseline
    changed[i] <- c(NA_real_, NaN, NA_real_, -Inf, Inf)[i]
    expect_gt(length(recovery_refit_numeric_differences(changed, baseline)), 0L)
  }
  expect_gt(length(recovery_refit_numeric_differences(unname(baseline), baseline)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(1, 1L)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(matrix(1, 1L), 1)), 0L)
  expect_gt(length(recovery_refit_numeric_differences(structure(1, class = "other"), 1)), 0L)
})

test_that("matrix coordinate scales preserve tiny variances and near-zero covariances", {
  baseline <- diag(c(4e-16, 2))
  changed <- baseline
  changed[1L, 1L] <- baseline[1L, 1L] * (1 + 1e-12)
  expect_length(recovery_refit_matrix_differences(changed, baseline), 0L)
  changed[1L, 1L] <- baseline[1L, 1L] * 1.01506
  expect_gt(length(recovery_refit_matrix_differences(changed, baseline)), 0L)
  changed <- baseline
  changed[1L, 2L] <- changed[2L, 1L] <- 1e-30
  expect_length(recovery_refit_matrix_differences(changed, baseline), 0L)
  # A generic vector has no covariance scale or license for an absolute floor.
  expect_gt(length(recovery_refit_numeric_differences(changed, baseline)), 0L)
  changed[1L, 2L] <- sqrt(baseline[1L, 1L]) * sqrt(baseline[2L, 2L]) * 1e-6
  expect_gt(length(recovery_refit_matrix_differences(changed, baseline)), 0L)
  changed <- baseline
  dimnames(changed) <- list(c("x", "z"), c("x", "z"))
  expect_gt(length(recovery_refit_matrix_differences(changed, baseline)), 0L)
  changed <- baseline
  changed[1L, 1L] <- NA_real_
  expect_gt(length(recovery_refit_matrix_differences(changed, baseline)), 0L)
  baseline[1L, 1L] <- 0
  changed <- baseline
  changed[1L, 1L] <- 1e-40
  expect_gt(length(recovery_refit_matrix_differences(changed, baseline)), 0L)
  expect_gt(length(recovery_refit_matrix_differences(matrix(1, 1, 2), matrix(1, 1, 2))), 0L)
})
