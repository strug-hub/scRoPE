# Test-only helpers for the mandatory installed-package parallel gate.
# Receipts contain explicit fields, never raw testthat objects or session paths.
scrope_gate_sanitizer <- function(private_values = character()) {
  private_values <- unique(as.character(private_values))
  private_values <- private_values[!is.na(private_values) & nzchar(private_values) &
                                     !private_values %in% c(".", "..", "/", "\\")]
  private_values <- private_values[order(nchar(private_values), decreasing = TRUE)]
  function(message) {
    message <- paste(as.character(message), collapse = " ")
    # Quoted paths can contain spaces. Preserve surrounding explanatory text
    # only when an explicit matching quote marks the end of the whole path.
    path_start <- "(?:[A-Za-z]:[\\\\/]|\\\\\\\\|~/|/)"
    message <- gsub(paste0("([\"'`])", path_start, "[^\\r\\n]*?\\1"),
                    "\\1<path>\\1", message, perl = TRUE)
    # An unquoted path's end is ambiguous when spaces are allowed. Redact
    # conservatively through the end of its line, not just its first token.
    message <- gsub(paste0("(?<![[:alnum:]_])", path_start, "[^\\r\\n]*"),
                    "<path>", message, perl = TRUE)
    for (value in private_values) {
      message <- gsub(value, "<redacted>", message, fixed = TRUE)
    }
    message <- gsub("[\r\n\t]+", " ", message)
    if (nchar(message) > 500L) message <- paste0(substr(message, 1L, 497L), "...")
    unname(message)
  }
}

scrope_gate_test_receipt <- function(summary, expected, sanitize) {
  columns <- c("test", "failed", "error", "warning", "skipped", "passed")
  valid_schema <- is.data.frame(summary) && all(columns %in% names(summary))
  receipt <- list(test = expected, status = "failed", rows = if (is.data.frame(summary)) {
    nrow(summary)
  } else 0L)
  if (!valid_schema || nrow(summary) != 1L ||
      !identical(summary$test, expected) || anyNA(summary[columns])) {
    receipt$messages <- "Required regression did not return exactly its expected test row."
    return(receipt)
  }
  counts <- summary[setdiff(columns, "test")]
  valid_counts <- all(vapply(counts, function(x) {
    (is.numeric(x) || is.logical(x)) && length(x) == 1L &&
      is.finite(x) && x >= 0 && x == floor(x)
  }, logical(1)))
  if (!valid_counts) {
    receipt$messages <- "Required regression returned invalid result counts."
    return(receipt)
  }
  receipt$counts <- lapply(counts, unname)
  receipt$status <- if (summary$failed == 0L && !summary$error &&
      summary$warning == 0L && !summary$skipped && summary$passed >= 1L) {
    "passed"
  } else "failed"
  messages <- character()
  if ("result" %in% names(summary)) {
    for (condition in summary$result[[1L]]) {
      if (inherits(condition, c("expectation_failure", "expectation_error",
                                "expectation_warning", "expectation_skip"))) {
        messages <- c(messages, sanitize(conditionMessage(condition)))
        if (length(messages) >= 5L) break
      }
    }
  }
  if (receipt$status == "failed" && !length(messages)) {
    messages <- "Required regression failed, warned, was skipped, or had no passing assertions."
  }
  receipt$messages <- messages
  receipt
}

scrope_gate_run_required_tests <- function(required_tests, run_test,
                                          save_checkpoint, sanitize) {
  state <- list(status = "running", tests = list(), checkpoint = NULL)
  for (file in names(required_tests)) {
    expected <- unname(required_tests[[file]])
    state$tests[[file]] <- list(test = expected, status = "running")
    state$checkpoint <- list(stage = "regression_started", test = file)
    save_checkpoint(state)
    # A test-file error is a failed result, not permission to omit the next test.
    state$tests[[file]] <- tryCatch({
      summary <- run_test(file, expected)
      scrope_gate_test_receipt(summary, expected, sanitize)
    }, error = function(condition) {
      list(test = expected, status = "failed", execution_error = TRUE,
           messages = sanitize(conditionMessage(condition)))
    })
    state$checkpoint <- list(stage = "regression_completed", test = file)
    save_checkpoint(state)
  }
  passed <- length(state$tests) == length(required_tests) &&
    length(required_tests) > 0L && all(vapply(state$tests, function(test) {
      identical(test$status, "passed")
    }, logical(1)))
  state$status <- if (passed) "passed" else "failed"
  state$checkpoint <- list(stage = "regressions_completed")
  save_checkpoint(state)
  state
}

scrope_gate_assert_passed <- function(state) {
  if (!identical(state$status, "passed")) {
    stop("Installed-package two-worker gate failed; inspect the bounded test receipts.",
         call. = FALSE)
  }
  invisible(state)
}
