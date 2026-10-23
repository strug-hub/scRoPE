installed_gate_helpers <- local({
  functions <- new.env(parent = baseenv())
  sys.source(testthat::test_path("..", "validation", "parallel_gate_helpers.R"),
             envir = functions)
  functions
})

installed_gate_row <- function(description, failed = 0L, error = FALSE,
                               warning = 0L, skipped = FALSE, passed = 1L) {
  data.frame(test = description, failed = failed, error = error, warning = warning,
             skipped = skipped, passed = passed, stringsAsFactors = FALSE)
}

test_that("the installed gate runs the second regression after a first failure", {
  required <- c("first.R" = "first required test", "second.R" = "second required test")
  for (first_outcome in c("assertion", "error")) {
    calls <- character()
    checkpoints <- list()
    state <- installed_gate_helpers$scrope_gate_run_required_tests(
      required,
      run_test = function(file, expected) {
        calls <<- c(calls, file)
        if (file == "first.R" && first_outcome == "error") {
          stop("controlled failure in /private/example-user/test.R")
        }
        installed_gate_row(expected, failed = as.integer(file == "first.R"))
      },
      save_checkpoint = function(state) {
        checkpoints[[length(checkpoints) + 1L]] <<- state
      },
      sanitize = installed_gate_helpers$scrope_gate_sanitizer("example-user")
    )
    expect_identical(calls, names(required))
    expect_identical(state$status, "failed")
    expect_identical(state$tests[[1L]]$status, "failed")
    expect_identical(state$tests[[2L]]$status, "passed")
    expect_length(checkpoints, 5L)
    expect_identical(checkpoints[[2L]]$tests[[1L]]$status, "failed")
    expect_identical(checkpoints[[3L]]$checkpoint$test, "second.R")
    expect_identical(checkpoints[[4L]]$tests[[2L]]$status, "passed")
    expect_error(installed_gate_helpers$scrope_gate_assert_passed(state),
                 "two-worker gate failed")
    expect_false(any(grepl("example-user|/private", unlist(state), fixed = FALSE)))
  }
})

test_that("the installed gate rejects skips warnings missing and wrong test rows", {
  run <- function(summary) {
    installed_gate_helpers$scrope_gate_run_required_tests(
      c("only.R" = "required test"), function(...) summary, function(...) NULL,
      installed_gate_helpers$scrope_gate_sanitizer()
    )
  }
  good <- installed_gate_row("required test")
  expect_identical(run(good)$status, "passed")
  expect_invisible(installed_gate_helpers$scrope_gate_assert_passed(run(good)))
  bad <- list(installed_gate_row("required test", skipped = TRUE),
              installed_gate_row("required test", warning = 1L),
              installed_gate_row("required test", error = TRUE),
              installed_gate_row("required test", passed = 0L),
              installed_gate_row("another test"), good[FALSE, ], rbind(good, good))
  for (summary in bad) {
    state <- run(summary)
    expect_identical(state$status, "failed")
    expect_error(installed_gate_helpers$scrope_gate_assert_passed(state),
                 "two-worker gate failed")
  }
})

test_that("gate receipts redact paths and bound saved messages", {
  sanitize <- installed_gate_helpers$scrope_gate_sanitizer(c("example-user", "example-host"))
  messages <- c("/Users/example-user/package/file.R:12 failed",
                "C:\\Users\\example-user\\file.R failed", "\\\\example-host\\share\\file.R failed",
                "~/private/file.R failed", "example-user on example-host", strrep("a", 1000L))
  clean <- vapply(messages, sanitize, character(1))
  expect_false(any(grepl("example-user|example-host|/Users|C:|~/|\\\\", clean)))
  expect_true(all(nchar(clean) <= 500L))
  summary <- installed_gate_row("required test", failed = 7L)
  summary$result <- list(rep(list(structure(
    list(message = messages[[1L]], call = NULL),
    class = c("expectation_failure", "condition")
  )), 7L))
  receipt <- installed_gate_helpers$scrope_gate_test_receipt(summary, "required test", sanitize)
  expect_length(receipt$messages, 5L)
  expect_false("result" %in% names(receipt))
  expect_identical(receipt$counts$failed, 7L)
  expect_false(any(grepl("example-user|/Users", receipt$messages)))
})


test_that("gate receipts redact complete paths with spaces without changing fit labels", {
  sanitize <- installed_gate_helpers$scrope_gate_sanitizer()
  quoted <- c("failed at '/unlisted-root/Research Group/file.R' while loading",
              'failed at "C:\\Research Group\\file.R" while loading',
              "failed at `\\\\unlisted-host\\Research Group\\file.R` while loading",
              "failed at '~/Research Group/file.R' while loading")
  for (message in quoted) {
    clean <- sanitize(message)
    expect_false(grepl("Research|Group|file.R|unlisted", clean))
    expect_match(clean, "while loading", fixed = TRUE)
  }
  unquoted <- c("failed at /unlisted-root/Research Group/file.R: detail",
                "failed at C:\\Research Group\\file.R: detail",
                "failed at \\\\unlisted-host\\Research Group\\file.R: detail")
  for (message in unquoted) {
    expect_identical(sanitize(message), "failed at <path>")
  }
  scientific <- "raw[[1L]]$diagnostics$curvature: rejected; score = NA; status = joint_ln_fit_rejected"
  expect_identical(sanitize(scientific), scientific)
  expect_identical(sanitize(paste0(unquoted[[1L]], "\n", scientific)),
                   paste("failed at <path>", scientific))
})
