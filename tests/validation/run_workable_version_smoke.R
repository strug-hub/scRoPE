#!/usr/bin/env Rscript

# Installed-package execution and schemas only; not inferential calibration.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_workable_version_smoke.R OUTPUT_RDS EXPECTED_PACKAGE_PATH",
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

data(sample_data, package = "scRoPE", lib.loc = dirname(expected_package_path))
design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
stopifnot(nrow(sample_data$count) >= 2L)
counts_two <- sample_data$count[1:2, , drop = FALSE]
counts_one <- sample_data$count[1L, , drop = FALSE]
common <- list(
  id = sample_data$sid, pred = design, offset = sample_data$offset,
  verbose = FALSE, cpc = 0, mincp = 1, ncore = 1L, keep_diagnostics = FALSE
)
timed_fit <- function(count, ...) {
  started <- proc.time()[["elapsed"]]
  fit <- do.call(scrope, c(list(count = count), common, list(...)))
  list(fit = fit, elapsed_seconds = proc.time()[["elapsed"]] - started)
}

# Valid inputs ensure an unrelated input error cannot pass a retirement check.
expect_retired <- function(label, thunk, selector) {
  condition <- tryCatch({ thunk(); NULL }, error = identity)
  if (!inherits(condition, "error")) {
    stop(label, " unexpectedly returned instead of rejecting a retired route.",
         call. = FALSE)
  }
  message <- conditionMessage(condition)
  stopifnot(
    grepl(selector, message, fixed = TRUE),
    grepl("retired", message, fixed = TRUE),
    grepl("0.2.5.9005 or earlier", message, fixed = TRUE)
  )
  list(rejected = TRUE, message = message)
}
legacy_inputs <- list(
  count = counts_one, id = sample_data$sid, pred = design,
  offset = sample_data$offset, verbose = FALSE, ncore = 1L, cpc = 0, mincp = 1
)
retirement <- list(
  legacy_estimator = expect_retired("legacy estimator", function() {
    timed_fit(counts_one, model = "LN", ln_estimator = "legacy_post_pml")
  }, "ln_estimator='legacy_post_pml'"),
  legacy_fallback = expect_retired("legacy fallback", function() {
    timed_fit(counts_one, model = "LN",
              ln_operational_fallback = "legacy_post_pml")
  }, "ln_operational_fallback='legacy_post_pml'"),
  nebula_default = expect_retired("default nebula", function() {
    do.call(nebula, legacy_inputs)
  }, "nebula(method='LN')"),
  nebula_explicit_ln = expect_retired("explicit nebula LN", function() {
    do.call(nebula, c(legacy_inputs, list(method = "LN")))
  }, "nebula(method='LN')"),
  nebula_pmm_ln = expect_retired("nebula PMM LN", function() {
    do.call(nebula, c(legacy_inputs, list(model = "PMM", method = "LN")))
  }, "nebula(method='LN')")
)

results <- list(
  metadata = list(
    purpose = "software_plumbing_only_not_calibration",
    n_subjects = length(unique(sample_data$sid)),
    n_cells = ncol(sample_data$count),
    mean_cells_per_subject = ncol(sample_data$count) /
      length(unique(sample_data$sid)),
    package_version = as.character(utils::packageVersion(
      "scRoPE", lib.loc = dirname(expected_package_path)
    )),
    package_path = loaded_package_path, dll_path = dll_path,
    dll_md5 = unname(tools::md5sum(dll_path)), session = utils::sessionInfo()
  ),
  retirement = retirement,
  routes = list(
    ln_default = timed_fit(counts_two,
                          additional_tests = c("score", "lrt"), lrt_details = TRUE),
    ln_joint_explicit = timed_fit(
      counts_two, model = "LN", ln_estimator = "joint_ln",
      ln_operational_fallback = "none",
      additional_tests = c("score", "lrt"), lrt_details = TRUE
    ),
    hl = timed_fit(counts_one, model = "HL",
                   additional_tests = c("score", "lrt"), lrt_details = TRUE),
    pgmm = timed_fit(counts_two, model = "PGMM")
  )
)
expected_routes <- c(ln_default = "LN", ln_joint_explicit = "LN",
                     hl = "HL", pgmm = "PGMM")
expected_schemas <- c(ln_default = 3L, ln_joint_explicit = 3L, hl = 2L, pgmm = 2L)
expected_estimators <- c(
  ln_default = "ln_joint_root", ln_joint_explicit = "ln_joint_root",
  hl = "first_order_nested_hl", pgmm = "pgmm_exact_marginal"
)
for (name in names(expected_routes)) {
  fit <- results$routes[[name]]$fit
  stopifnot(
    is.data.frame(fit$summary), nrow(fit$summary) > 0L,
    is.data.frame(fit$overdispersion),
    length(fit$convergence) == nrow(fit$summary),
    length(fit$algorithm) == nrow(fit$summary),
    is.list(fit$provenance),
    identical(fit$provenance$schema_version, unname(expected_schemas[[name]])),
    identical(fit$provenance$route, unname(expected_routes[[name]])),
    identical(fit$provenance$estimator, unname(expected_estimators[[name]])),
    is.data.frame(fit$provenance$test_capabilities),
    identical(fit$provenance$test_capabilities$test, c("wald", "score", "lrt"))
  )
}
default <- results$routes$ln_default$fit
joint <- results$routes$ln_joint_explicit$fit
for (field in c("summary", "overdispersion", "convergence", "algorithm",
                "provenance", "lrt_details")) {
  stopifnot(isTRUE(all.equal(default[[field]], joint[[field]], tolerance = 1e-10)))
}
for (fit in list(default, joint)) {
  stopifnot(
    is.data.frame(fit$lrt_details), nrow(fit$lrt_details) > 0L,
    all(c("contrast", "p_value", "score_p_value", "success", "score_success") %in%
          names(fit$lrt_details)),
    identical(fit$provenance$test_capabilities$public_available,
              c(TRUE, TRUE, TRUE)),
    is.data.frame(fit$test_status),
    is.null(fit$diagnostics)
  )
  # Existing additional_tests opts into experimental asymptotic tails only.
  # Check accepted and rejected rows separately, never replacing failure by p=1.
  details <- fit$lrt_details
  stopifnot(any(details$success), any(details$score_success))
  for (test in c("score", "lrt")) {
    success <- if (test == "score") details$score_success else details$success
    stat <- if (test == "score") details$score_stat else details$stat
    p_value <- if (test == "score") details$score_p_value else details$p_value
    failure <- if (test == "score") details$score_failure_reason else details$failure_reason
    stopifnot(
      all(is.finite(stat[success])),
      isTRUE(all.equal(p_value[success],
                      pchisq(stat[success], 1L, lower.tail = FALSE),
                      tolerance = 1e-12)),
      all(is.na(stat[!success])), all(is.na(p_value[!success])),
      all(!is.na(failure[!success]) & nzchar(failure[!success]))
    )
  }
  stopifnot(all(details$method[details$success] ==
                  "joint_ln_pace_adjusted_ratio"))
  optional <- subset(fit$test_status, test %in% c("score", "lrt"))
  stopifnot(all(optional$status[optional$success] == "experimental_asymptotic"))
}
hl <- results$routes$hl$fit
stopifnot(is.data.frame(hl$hl_adj_details), nrow(hl$hl_adj_details) > 0L)
hl_objective_rows <- hl$hl_adj_details$test == "objective_ratio"
stopifnot(
  all(c("score", "objective_ratio") %in% hl$hl_adj_details$test),
  any(hl_objective_rows),
  identical(hl$provenance$test_capabilities$public_available, c(TRUE, TRUE, FALSE)),
  all(!hl$hl_adj_details$success[hl_objective_rows]),
  all(hl$hl_adj_details$failure_reason[hl_objective_rows] ==
        "objective_ratio_unsupported")
)
saveRDS(results, output)
stopifnot(identical(readRDS(output), results))
cat("Workable-version software smoke: PASS\n")
cat("Purpose: software plumbing only; not calibration\n")
cat("Saved:", normalizePath(output, mustWork = TRUE), "\n")
for (name in names(results$routes)) {
  cat(name, "elapsed_seconds=",
      format(results$routes[[name]]$elapsed_seconds, digits = 5), "\n")
}
