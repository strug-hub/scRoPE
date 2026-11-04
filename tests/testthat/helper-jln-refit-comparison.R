# Independent numerical fits need not be bitwise identical across platforms.
# 1e-10 is well below the 1e-6 root/Newton certificates. Only explicitly named
# fitted values receive tolerance; metadata, types and nonfinite states are exact.
recovery_refit_attributes_identical <- function(actual, expected) {
  # Attribute storage order is immaterial, as in identical()'s default behavior.
  a <- attributes(actual)
  b <- attributes(expected)
  if (is.null(a) || is.null(b)) return(identical(a, b))
  identical(a[order(names(a))], b[order(names(b))])
}

recovery_refit_numeric_differences <- function(actual, expected) {
  if (!is.double(actual) || !is.double(expected) ||
      !recovery_refit_attributes_identical(actual, expected))
    return("numerical type or attributes differ")
  if (!identical(is.na(actual), is.na(expected)) ||
      !identical(is.nan(actual), is.nan(expected)) ||
      !identical(is.finite(actual), is.finite(expected)))
    return("missing or nonfinite states differ")
  finite <- is.finite(actual)
  if (!identical(actual[!finite], expected[!finite]))
    return("nonfinite values differ")
  comparison <- all.equal(actual, expected, tolerance = 1e-10)
  if (isTRUE(comparison)) character(0) else as.character(comparison)
}

recovery_refit_summary_fields <- function(x) {
  unique(c(grep("^(logFC_|se_|p_|wald_stat_)", names(x), value = TRUE),
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
      recovery_refit_numeric_differences(actual[[name]], expected[[name]])
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

expect_recovery_refit_numeric <- function(actual, expected) {
  differences <- recovery_refit_numeric_differences(actual, expected)
  expect_identical(length(differences), 0L, info = paste(differences, collapse = "\n"))
}
