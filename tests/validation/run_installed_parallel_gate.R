#!/usr/bin/env Rscript

# Execute both installed-package PSOCK regressions; a skip is a failed gate.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_installed_parallel_gate.R OUTPUT_RDS EXPECTED_PACKAGE_PATH",
       call. = FALSE)
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
source(file.path(dirname(script), "parallel_gate_helpers.R"), local = TRUE)
output <- args[[1L]]
private_values <- c(args, dirname(args), getwd(), tempdir(), path.expand("~"),
                    Sys.getenv(c("USER", "LOGNAME", "USERNAME", "HOSTNAME",
                                 "USERPROFILE", "COMPUTERNAME")),
                    unname(Sys.info()[c("nodename", "login", "user", "effective_user")]))
sanitize <- scrope_gate_sanitizer(private_values)
record <- list(
  purpose = "installed_package_parallel_software_gate_not_calibration",
  status = "started", checkpoint = list(stage = "installation_validation"),
  r_version = R.version.string, platform = R.version$platform,
  tests = list()
)
# Both regression helpers write their paired receipts alongside this gate receipt,
# where the existing CI artifact upload retains them. Paths are not saved in it.
tryCatch({
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
  Sys.setenv(SCROPE_TEST_DIAGNOSTICS_DIR = dirname(output))
  saveRDS(record, output)
  expected_package_path <- normalizePath(args[[2L]], mustWork = TRUE)
  if (!file.exists(file.path(expected_package_path, "Meta", "package.rds"))) {
    stop("EXPECTED_PACKAGE_PATH must be an installed package, not a source tree.",
         call. = FALSE)
  }
  suppressPackageStartupMessages(
    library(scRoPE, lib.loc = dirname(expected_package_path))
  )
  loaded_package_path <- normalizePath(
    getNamespaceInfo(asNamespace("scRoPE"), "path"), mustWork = TRUE
  )
  record$namespace_matches_expected <- identical(loaded_package_path, expected_package_path)
  if (!record$namespace_matches_expected) {
    stop("Loaded scRoPE does not match the expected installation.", call. = FALSE)
  }
  dll <- getLoadedDLLs()[["scRoPE"]]
  if (is.null(dll)) stop("The installed scRoPE DLL is not loaded.", call. = FALSE)
  dll_path <- normalizePath(dll[["path"]], mustWork = TRUE)
  record$dll_inside_expected_install <- startsWith(
    dll_path, paste0(expected_package_path, .Platform$file.sep)
  )
  if (!record$dll_inside_expected_install) {
    stop("scRoPE DLL is outside the expected installation.", call. = FALSE)
  }
  record$dll_md5 <- unname(tools::md5sum(dll_path))
  record$package_version <- as.character(utils::packageVersion(
    "scRoPE", lib.loc = dirname(expected_package_path)
  ))
  record$testthat_version <- as.character(utils::packageVersion("testthat"))
  record$parallelly_version <- as.character(utils::packageVersion("parallelly"))
  record$available_cores <- unname(parallelly::availableCores())
  record$core_details <- parallelly::availableCores(which = "all")
  names(record$core_details) <- vapply(names(record$core_details), sanitize, character(1))
  record$checkpoint <- list(stage = "resource_validation")
  saveRDS(record, output)
  if (!is.finite(record$available_cores) || record$available_cores < 2L) {
    stop(paste0("Two allocated cores are required; parallelly detected ",
                record$available_cores, ". This gate cannot pass by skipping."),
         call. = FALSE)
  }

  # Verify two real workers can load this exact installation and DLL. Resource
  # detection is never overridden to make this mandatory gate pass.
  probe_workers <- function() {
    cluster <- parallelly::makeClusterPSOCK(
      2L, rscript_libs = unique(c(dirname(expected_package_path), .libPaths()))
    )
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(expected) {
      suppressPackageStartupMessages(library(scRoPE, lib.loc = dirname(expected)))
      namespace_path <- normalizePath(
        getNamespaceInfo(asNamespace("scRoPE"), "path"), mustWork = TRUE
      )
      dll <- getLoadedDLLs()[["scRoPE"]]
      if (is.null(dll)) stop("Worker did not load the scRoPE DLL.")
      dll_path <- normalizePath(dll[["path"]], mustWork = TRUE)
      namespace_matches <- identical(namespace_path, expected)
      dll_inside <- startsWith(dll_path, paste0(expected, .Platform$file.sep))
      stopifnot(namespace_matches, dll_inside)
      list(pid = Sys.getpid(), namespace_matches_expected = namespace_matches,
           dll_inside_expected_install = dll_inside,
           dll_md5 = unname(tools::md5sum(dll_path)),
           package_version = as.character(utils::packageVersion("scRoPE")),
           r_version = R.version.string, platform = R.version$platform)
    }, expected_package_path)
  }
  record$checkpoint <- list(stage = "worker_probe")
  saveRDS(record, output)
  record$workers <- probe_workers()
  worker_pids <- vapply(record$workers, `[[`, integer(1), "pid")
  record$two_distinct_worker_processes <- length(unique(worker_pids)) == 2L &&
    !Sys.getpid() %in% worker_pids
  record$worker_dlls_match_parent <- all(vapply(record$workers, function(worker) {
    identical(worker$dll_md5, record$dll_md5)
  }, logical(1)))
  if (!record$two_distinct_worker_processes || !record$worker_dlls_match_parent) {
    stop("The two-worker installation/DLL identity probe failed.", call. = FALSE)
  }

  test_dir <- file.path(expected_package_path, "tests", "testthat")
  required_tests <- c(
    "test-ln-joint-psock-namespace.R" =
      "installed joint-LN PSOCK fits use the running package",
    "test-ln-joint-routes-modes.R" =
      "joint LN single-core and multicore results agree"
  )
  record$checkpoint <- list(stage = "installed_test_validation")
  saveRDS(record, output)
  if (!file.exists(file.path(test_dir, "helper-ln-joint.R"))) {
    stop("Required installed tests are missing; install with --install-tests.",
         call. = FALSE)
  }
  state <- scrope_gate_run_required_tests(
    required_tests,
    run_test = function(file, expected) {
      if (!file.exists(file.path(test_dir, file))) {
        stop("Required installed regression is missing; install with --install-tests.",
             call. = FALSE)
      }
      result <- testthat::test_file(
        file.path(test_dir, file), desc = expected, package = "scRoPE",
        load_package = "installed", load_helpers = TRUE, reporter = "silent",
        stop_on_failure = FALSE, stop_on_warning = FALSE
      )
      as.data.frame(result)
    },
    save_checkpoint = function(state) {
      record$status <<- state$status
      record$tests <<- state$tests
      record$checkpoint <<- state$checkpoint
      saveRDS(record, output)
    },
    sanitize = sanitize
  )
  scrope_gate_assert_passed(state)
  cat("Installed-package two-worker gate: PASS (both required tests ran; no skips)\n")
  cat("Saved gate and paired regression receipts in the requested output directory.\n")
}, error = function(condition) {
  record$status <- "failed"
  record$failure <- sanitize(conditionMessage(condition))
  # If the requested receipt directory itself is unwritable, preserve the safe
  # failure message without exposing a second unsanitized filesystem condition.
  tryCatch(saveRDS(record, output), error = function(...) invisible(NULL))
  stop(record$failure, call. = FALSE)
})
