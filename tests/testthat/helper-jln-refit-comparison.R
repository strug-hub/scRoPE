# Independent numerical fits need not be bitwise identical across platforms.
# 1e-10 relative error is well below the 1e-6 root/Newton certificates. Only
# explicitly named fitted values receive tolerance; metadata, types and
# nonfinite states are exact. See tests/NUMERICAL_TESTING.md for the policy.
recovery_refit_attributes_identical <- function(actual, expected) {
  # Attribute storage order is immaterial, as in identical()'s default behavior.
  a <- attributes(actual)
  b <- attributes(expected)
  if (is.null(a) || is.null(b)) return(identical(a, b))
  identical(a[order(names(a))], b[order(names(b))])
}

recovery_refit_numeric_differences <- function(actual, expected, absolute = 0,
    comparison_scale = NULL) {
  if (!is.double(actual) || !is.double(expected) ||
      length(actual) != length(expected) ||
      !recovery_refit_attributes_identical(actual, expected))
    return("numerical type or attributes differ")
  if (!identical(is.na(actual), is.na(expected)) ||
      !identical(is.nan(actual), is.nan(expected)) ||
      !identical(is.finite(actual), is.finite(expected)))
    return("missing or nonfinite states differ")
  finite <- is.finite(actual)
  if (!identical(actual[!finite], expected[!finite]))
    return("nonfinite values differ")
  # Compare each element separately: an unrelated large entry must not mask
  # drift in a small SE, p-value or covariance entry. Scaling before subtracting
  # also avoids overflow when finite values have opposite signs.
  a <- actual[finite]
  b <- expected[finite]
  scale <- pmax(abs(a), abs(b))
  if (!is.null(comparison_scale)) {
    if (length(comparison_scale) != length(actual) ||
        any(!is.finite(comparison_scale[finite]) | comparison_scale[finite] < 0))
      return("comparison scale is unavailable")
    scale <- pmax(scale, comparison_scale[finite])
  }
  nonzero <- scale > 0
  error <- abs(a[nonzero] / scale[nonzero] - b[nonzero] / scale[nonzero])
  budget <- 1e-10 + absolute / scale[nonzero]
  bad <- which(nonzero)[error > budget]
  if (!length(bad)) return(character(0))
  indices <- which(finite)[bad]
  sprintf("fitted values differ at element(s) %s (relative budget 1e-10, absolute budget %g)",
    paste(head(indices, 8L), collapse = ", "), absolute)
}

recovery_refit_matrix_differences <- function(actual, expected) {
  if (!is.matrix(actual) || !is.matrix(expected) ||
      !is.double(actual) || !is.double(expected) ||
      nrow(actual) != ncol(actual) ||
      !identical(dim(actual), dim(expected))) return("square numerical matrices required")
  # Covariance/information off-diagonals may vanish through cancellation.
  # Measure them in their coordinate scales, while each diagonal still has
  # its own relative budget. No unit-sized or global-matrix floor is used.
  diagonal_scale <- pmax(sqrt(abs(diag(actual))), sqrt(abs(diag(expected))))
  scale <- outer(diagonal_scale, diagonal_scale)
  recovery_refit_numeric_differences(actual, expected, comparison_scale = scale)
}

recovery_refit_absolute_budget <- function(name) {
  # Signed coefficients and interval endpoints may cross zero under roundoff.
  # SEs, variances, statistics and probabilities receive no absolute floor.
  if (grepl("^logFC_", name) ||
      name %in% c("coefficient", "beta", "ci_lower", "ci_upper")) 1e-12 else 0
}

recovery_refit_summary_fields <- function(x) {
  unique(c(grep("^(logFC_|se_|p_|wald_stat_|model_wald_stat_)", names(x), value = TRUE),
    intersect(c("coefficient", "ci_lower", "ci_upper"), names(x))))
}

recovery_refit_table_differences <- function(actual, expected,
    numerical = recovery_refit_summary_fields(actual)) {
  if (!is.data.frame(actual) || !is.data.frame(expected) ||
      !recovery_refit_attributes_identical(actual, expected))
    return("table structure differs")
  if (!length(numerical) || !all(numerical %in% names(actual)))
    return("fitted numerical columns are missing")
  unlist(lapply(names(actual), function(name) {
    difference <- if (name %in% numerical) {
      recovery_refit_numeric_differences(actual[[name]], expected[[name]],
        absolute = recovery_refit_absolute_budget(name))
    } else if (!identical(actual[[name]], expected[[name]])) {
      "exact metadata differs"
    } else character(0)
    if (length(difference)) paste0(name, ": ", difference) else character(0)
  }), use.names = FALSE)
}

expect_recovery_refit_table <- function(actual, expected,
    numerical = recovery_refit_summary_fields(actual)) {
  differences <- recovery_refit_table_differences(actual, expected, numerical)
  expect_identical(length(differences), 0L, info = paste(differences, collapse = "\n"))
}

expect_recovery_refit_numeric <- function(actual, expected, absolute = 0) {
  differences <- recovery_refit_numeric_differences(actual, expected, absolute)
  expect_identical(length(differences), 0L, info = paste(differences, collapse = "\n"))
}

expect_recovery_refit_matrix <- function(actual, expected) {
  differences <- recovery_refit_matrix_differences(actual, expected)
  expect_identical(length(differences), 0L, info = paste(differences, collapse = "\n"))
}
