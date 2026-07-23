#!/usr/bin/env Rscript

# Small end-to-end software smoke test for the unified sc-RoPE routes.
# This checks execution and schemas only; it is not an inferential calibration.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: run_workable_version_smoke.R OUTPUT_RDS EXPECTED_PACKAGE_PATH",
    call. = FALSE
  )
}
output <- args[[1L]]
expected_package_path <- normalizePath(args[[2L]], mustWork = TRUE)

suppressPackageStartupMessages(library(scRoPE))
loaded_package_path <- normalizePath(find.package("scRoPE"), mustWork = TRUE)
if (!identical(loaded_package_path, expected_package_path)) {
  stop(
    "Loaded scRoPE from ", loaded_package_path,
    "; expected release candidate at ", expected_package_path, ".",
    call. = FALSE
  )
}
data(sample_data, package = "scRoPE")

design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
counts_two <- sample_data$count[seq_len(min(2L, nrow(sample_data$count))), , drop = FALSE]
counts_one <- sample_data$count[1L, , drop = FALSE]
common <- list(
  id = sample_data$sid,
  pred = design,
  offset = sample_data$offset,
  verbose = FALSE,
  cpc = 0,
  mincp = 1,
  ncore = 1L,
  keep_diagnostics = FALSE
)

timed_fit <- function(count, ...) {
  started <- proc.time()[["elapsed"]]
  fit <- do.call(scrope, c(list(count = count), common, list(...)))
  list(fit = fit, elapsed_seconds = proc.time()[["elapsed"]] - started)
}

results <- list(
  metadata = list(
    purpose = "software_plumbing_only_not_calibration",
    n_subjects = length(unique(sample_data$sid)),
    n_cells = ncol(sample_data$count),
    mean_cells_per_subject = ncol(sample_data$count) /
      length(unique(sample_data$sid)),
    package_version = as.character(utils::packageVersion("scRoPE")),
    package_path = loaded_package_path
  ),
  routes = list(
    ln_legacy = timed_fit(counts_two, model = "LN"),
    ln_joint = timed_fit(
      counts_two,
      model = "LN",
      ln_estimator = "joint_ln",
      additional_tests = c("score", "lrt"),
      lrt_details = TRUE
    ),
    hl = timed_fit(
      counts_one,
      model = "HL",
      additional_tests = c("score", "lrt"),
      lrt_details = TRUE
    ),
    pgmm = timed_fit(counts_two, model = "PGMM")
  )
)

expected_routes <- c(
  ln_legacy = "LN",
  ln_joint = "LN",
  hl = "HL",
  pgmm = "PGMM"
)
for (name in names(expected_routes)) {
  fit <- results$routes[[name]]$fit
  stopifnot(
    is.data.frame(fit$summary),
    nrow(fit$summary) > 0L,
    is.data.frame(fit$overdispersion),
    length(fit$convergence) == nrow(fit$summary),
    length(fit$algorithm) == nrow(fit$summary),
    is.list(fit$provenance),
    identical(fit$provenance$schema_version, 2L),
    identical(fit$provenance$route, unname(expected_routes[[name]])),
    is.data.frame(fit$provenance$test_capabilities),
    identical(fit$provenance$test_capabilities$test, c("wald", "score", "lrt"))
  )
}

joint <- results$routes$ln_joint$fit
stopifnot(
  is.data.frame(joint$lrt_details),
  nrow(joint$lrt_details) > 0L,
  all(c(
    "contrast", "p_value", "score_p_value", "success", "score_success"
  ) %in% names(joint$lrt_details)),
  identical(
    joint$provenance$test_capabilities$public_available,
    c(TRUE, FALSE, FALSE)
  ),
  all(is.na(joint$lrt_details$p_value)),
  all(is.na(joint$lrt_details$score_p_value)),
  all(!joint$lrt_details$success),
  all(!joint$lrt_details$score_success)
)

hl <- results$routes$hl$fit
hl_objective_rows <- hl$hl_adj_details$test == "objective_ratio"
stopifnot(
  is.data.frame(hl$hl_adj_details),
  nrow(hl$hl_adj_details) > 0L,
  all(c("score", "objective_ratio") %in% hl$hl_adj_details$test),
  any(hl_objective_rows),
  identical(hl$provenance$test_capabilities$public_available, c(TRUE, TRUE, FALSE)),
  all(!hl$hl_adj_details$success[hl_objective_rows]),
  all(
    hl$hl_adj_details$failure_reason[hl_objective_rows] ==
      "objective_ratio_unsupported"
  )
)

saveRDS(results, output)
roundtrip <- readRDS(output)
stopifnot(identical(names(roundtrip$routes), names(results$routes)))

cat("Workable-version software smoke: PASS\n")
cat("Purpose: software plumbing only; not calibration\n")
cat("Saved:", normalizePath(output, mustWork = TRUE), "\n")
for (name in names(results$routes)) {
  cat(name, "elapsed_seconds=", format(results$routes[[name]]$elapsed_seconds, digits = 5), "\n")
}
