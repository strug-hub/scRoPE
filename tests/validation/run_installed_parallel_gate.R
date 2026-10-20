#!/usr/bin/env Rscript

# Execute the two installed-package PSOCK regressions; a skip is a failed gate.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_installed_parallel_gate.R OUTPUT_RDS EXPECTED_PACKAGE_PATH",
       call. = FALSE)
}
output <- args[[1L]]
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
if (!identical(loaded_package_path, expected_package_path)) {
  stop("Loaded scRoPE from ", loaded_package_path, "; expected ",
       expected_package_path, ".", call. = FALSE)
}
dll <- getLoadedDLLs()[["scRoPE"]]
if (is.null(dll)) stop("The installed scRoPE DLL is not loaded.", call. = FALSE)
dll_path <- normalizePath(dll[["path"]], mustWork = TRUE)
if (!startsWith(dll_path, paste0(expected_package_path, .Platform$file.sep))) {
  stop("scRoPE DLL is outside the expected installation: ", dll_path,
       call. = FALSE)
}
record <- list(
  purpose = "installed_package_parallel_software_gate_not_calibration",
  status = "started", package_path = loaded_package_path, dll_path = dll_path,
  dll_md5 = unname(tools::md5sum(dll_path)),
  package_version = as.character(utils::packageVersion(
    "scRoPE", lib.loc = dirname(expected_package_path)
  )),
  available_cores = parallelly::availableCores(),
  core_details = parallelly::availableCores(which = "all"),
  session = utils::sessionInfo()
)
saveRDS(record, output)
fail_gate <- function(message) {
  record$status <- "failed"
  record$failure <- message
  saveRDS(record, output)
  stop(message, call. = FALSE)
}
if (!is.finite(record$available_cores) || record$available_cores < 2L) {
  fail_gate(paste0("Two allocated cores are required; parallelly detected ",
                   record$available_cores, ". This gate cannot pass by skipping."))
}

# Verify two real worker processes can load this exact installation and DLL.
# Do not override parallelly's resource detection to make the gate pass.
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
    stopifnot(identical(namespace_path, expected),
              startsWith(dll_path, paste0(expected, .Platform$file.sep)))
    list(pid = Sys.getpid(), package_path = namespace_path, dll_path = dll_path,
         dll_md5 = unname(tools::md5sum(dll_path)))
  }, expected_package_path)
}
workers <- tryCatch(probe_workers(), error = identity)
if (inherits(workers, "error")) fail_gate(conditionMessage(workers))
record$workers <- workers
worker_pids <- vapply(workers, `[[`, integer(1), "pid")
if (length(unique(worker_pids)) != 2L || Sys.getpid() %in% worker_pids ||
    !all(vapply(workers, function(x) identical(x$dll_md5, record$dll_md5), logical(1)))) {
  fail_gate("The two-worker installation/DLL identity probe failed.")
}

test_dir <- file.path(expected_package_path, "tests", "testthat")
required_tests <- c(
  "test-ln-joint-psock-namespace.R" =
    "installed joint-LN PSOCK fits use the running package",
  "test-ln-joint-routes-modes.R" =
    "joint LN single-core and multicore results agree"
)
if (!all(file.exists(file.path(test_dir, names(required_tests)))) ||
    !file.exists(file.path(test_dir, "helper-ln-joint.R"))) {
  fail_gate("Required installed tests are missing; install with --install-tests.")
}
record$tests <- list()
for (file in names(required_tests)) {
  result <- tryCatch(testthat::test_file(
    file.path(test_dir, file), desc = unname(required_tests[[file]]),
    package = "scRoPE", load_package = "installed", load_helpers = TRUE,
    reporter = "summary", stop_on_failure = FALSE, stop_on_warning = FALSE
  ), error = identity)
  if (inherits(result, "error")) fail_gate(conditionMessage(result))
  summary <- as.data.frame(result)
  record$tests[[file]] <- summary
  expected <- unname(required_tests[[file]])
  columns <- c("test", "failed", "error", "warning", "skipped", "passed")
  if (!all(columns %in% names(summary)) || nrow(summary) != 1L ||
      !identical(summary$test, expected) || anyNA(summary[columns]) ||
      any(summary$failed > 0L | summary$error | summary$warning > 0L |
            summary$skipped | summary$passed < 1L)) {
    fail_gate(paste0("Required parallel regression failed, warned, was skipped, ",
                     "or was not executed: ", expected))
  }
}
record$status <- "passed"
saveRDS(record, output)
cat("Installed-package two-worker gate: PASS (both required tests ran; no skips)\n")
cat("Saved:", normalizePath(output, mustWork = TRUE), "\n")
