#!/usr/bin/env Rscript

# Bounded synthetic timing harness for the shadow-only joint-LN candidate
# resolver. This measures classifier/probe plumbing, not estimator runtime and
# not inferential calibration.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop(
    "Usage: run_joint_ln_availability_runtime.R OUTPUT_RDS [REPETITIONS]",
    call. = FALSE
  )
}
output_path <- args[[1L]]
repetitions <- if (length(args) == 2L) as.integer(args[[2L]]) else 50L
if (is.na(repetitions) || repetitions < 20L || repetitions > 200L) {
  stop("REPETITIONS must be an integer from 20 through 200.", call. = FALSE)
}
output_parent <- normalizePath(dirname(output_path), mustWork = TRUE)
output_path <- file.path(output_parent, basename(output_path))

suppressPackageStartupMessages(library(scRoPE))

resolver <- getFromNamespace("ln_joint_resolve_candidates", "scRoPE")
resolution_thresholds <- getFromNamespace(
  "ln_joint_resolution_thresholds",
  "scRoPE"
)
check_names <- getFromNamespace("ln_joint_candidate_check_names", "scRoPE")
objective_provenance <- getFromNamespace(
  "ln_joint_objective_provenance",
  "scRoPE"
)
curvature_diagnostics <- getFromNamespace(
  "ln_joint_curvature_diagnostics",
  "scRoPE"
)
selected_thresholds <- getFromNamespace("ln_joint_thresholds", "scRoPE")
scaled_max <- getFromNamespace("ln_joint_scaled_max", "scRoPE")
root_metric <- getFromNamespace("ln_joint_root_metric", "scRoPE")
newton_metric <- getFromNamespace("ln_joint_probe_newton_metric", "scRoPE")
probe_full_log <- getFromNamespace("ln_joint_probe_full_log", "scRoPE")

installed_tree_manifest <- function() {
  package_root <- normalizePath(
    system.file(package = "scRoPE"),
    mustWork = TRUE
  )
  files <- list.files(
    package_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  files <- files[file.info(files)$isdir %in% FALSE]
  relative <- substring(files, nchar(package_root) + 2L)
  hashes <- as.character(tools::md5sum(files))
  order_index <- order(relative)
  stats::setNames(hashes[order_index], relative[order_index])
}

runtime_refresh_candidate <- function(candidate) {
  state <- candidate$state
  thresholds <- selected_thresholds()
  candidate$donor_additivity_scaled <- scaled_max(
    rowSums(state$U_subject_reduced),
    state$U_reduced
  )
  candidate$root_residual <- root_metric(
    state$U_subject_reduced,
    state$U_reduced
  )
  candidate$newton_displacement <- newton_metric(state)
  candidate$curvature <- curvature_diagnostics(state$A_reduced)
  dispersion <- exp(candidate$theta_full_log[
    candidate$n_beta + seq_len(2L)
  ])
  boundary_tolerance <- sqrt(.Machine$double.eps) *
    (1 + pmax(abs(candidate$bounds$lower), abs(candidate$bounds$upper)))
  candidate$interior <- all(
    dispersion > candidate$bounds$lower + boundary_tolerance &
      dispersion < candidate$bounds$upper - boundary_tolerance
  )
  candidate$checks <- c(
    donor_additivity = candidate$donor_additivity_scaled <=
      thresholds$donor_additivity,
    root = candidate$root_residual <= thresholds$root,
    newton = candidate$newton_displacement <= thresholds$newton,
    curvature = isTRUE(candidate$curvature$positive_definite) &&
      candidate$curvature$scaled_min_eigenvalue >=
        thresholds$curvature_min_eigen,
    condition = candidate$curvature$scaled_condition_number <=
      thresholds$condition_number,
    interior = candidate$interior,
    constraint = candidate$constraint_residual <= thresholds$constraint
  )
  candidate$failure <- if (all(candidate$checks)) {
    NA_character_
  } else {
    "synthetic_complete_check_vector_failure"
  }
  candidate
}

runtime_candidate <- function(
    label,
    objective = 10,
    theta_full_log = c(0.2, -0.1, log(0.8), log(1.2)),
    certified = TRUE,
    core = TRUE) {
  n_beta <- 2L
  score <- rep(0, length(theta_full_log))
  score_subject <- cbind(score / 2, score / 2)
  state <- list(
    theta_reduced_log = as.numeric(theta_full_log),
    full_log = as.numeric(theta_full_log),
    objective_provenance = objective_provenance(),
    U_reduced = score,
    U_subject_reduced = score_subject,
    A_reduced = diag(length(theta_full_log)),
    raw = list(value = as.numeric(objective)),
    beta = theta_full_log[seq_len(n_beta)],
    sigma2 = exp(theta_full_log[[n_beta + 1L]]),
    gamma = exp(theta_full_log[[n_beta + 2L]])
  )
  candidate <- list(
    label = label,
    start_provenance = "synthetic_runtime_fixture",
    required = isTRUE(core),
    core = isTRUE(core),
    finite_endpoint = TRUE,
    certified = isTRUE(certified),
    objective = as.numeric(objective),
    theta_reduced = c(
      theta_full_log[seq_len(n_beta)],
      exp(theta_full_log[n_beta + seq_len(2L)])
    ),
    theta_reduced_log = as.numeric(theta_full_log),
    theta_full_log = as.numeric(theta_full_log),
    n_beta = n_beta,
    parameterization = "beta, log(sigma2), log(gamma_nb_size)",
    objective_provenance = objective_provenance(),
    objective_instance = list(
      scheme = "synthetic_joint_ln_objective_context_v1",
      id = "synthetic-runtime",
      gene_index = 1L
    ),
    bounds = list(lower = c(0.1, 0.1), upper = c(10, 10)),
    probe_bounds = list(
      lower = c(-100, -100, log(0.1), log(0.1)),
      upper = c(100, 100, log(10), log(10))
    ),
    probe_map = list(
      schema_version = 1L,
      beta_star = c(0, 0),
      null_basis = diag(2),
      n_beta = 2L,
      free_dim = 2L
    ),
    optimizer = list(convergence = 0L, message = "synthetic"),
    state = state,
    donor_additivity_scaled = scaled_max(rowSums(score_subject), score),
    root_residual = root_metric(score_subject, score),
    newton_displacement = newton_metric(state),
    curvature = curvature_diagnostics(state$A_reduced),
    interior = TRUE,
    constraint_residual = 0,
    checks = setNames(rep(TRUE, length(check_names())), check_names()),
    failure = NA_character_
  )
  runtime_refresh_candidate(candidate)
}

runtime_update_state <- function(candidate, full_log, score = NULL) {
  n_beta <- as.integer(candidate$n_beta)
  candidate$theta_full_log <- as.numeric(full_log)
  candidate$theta_reduced_log <- as.numeric(full_log)
  candidate$theta_reduced <- c(
    full_log[seq_len(n_beta)],
    exp(full_log[n_beta + seq_len(2L)])
  )
  candidate$state$theta_reduced_log <- as.numeric(full_log)
  candidate$state$full_log <- as.numeric(full_log)
  candidate$state$beta <- full_log[seq_len(n_beta)]
  candidate$state$sigma2 <- exp(full_log[[n_beta + 1L]])
  candidate$state$gamma <- exp(full_log[[n_beta + 2L]])
  if (!is.null(score)) {
    score <- as.numeric(score)
    candidate$state$U_reduced <- score
    candidate$state$U_subject_reduced <- cbind(score / 2, score / 2)
  }
  candidate$state$raw$value <- as.numeric(candidate$objective)
  runtime_refresh_candidate(candidate)
}

runtime_uncertified <- function(selected, kind) {
  candidate <- selected
  candidate$label <- paste0("synthetic_", kind)
  candidate$required <- FALSE
  candidate$core <- FALSE
  candidate$certified <- FALSE
  candidate$objective <- selected$objective
  full_log <- selected$theta_full_log
  n_beta <- as.integer(selected$n_beta)

  if (identical(kind, "dominated")) {
    candidate$objective <- selected$objective + 0.01
    candidate <- runtime_update_state(candidate, full_log)
  } else if (identical(kind, "boundary")) {
    full_log[[n_beta + 1L]] <- log(candidate$bounds$lower[[1L]])
    candidate <- runtime_update_state(candidate, full_log)
  } else if (identical(kind, "same_basin")) {
    delta <- c(2e-5, rep(0, length(full_log) - 1L))
    full_log <- full_log + delta
    candidate$objective <- selected$objective + 0.5 * sum(delta^2)
    candidate <- runtime_update_state(candidate, full_log, score = -delta)
  } else {
    stop("Unknown runtime fixture kind: ", kind)
  }
  candidate
}

runtime_evaluator <- function(candidate, theta_reduced_log) {
  initial <- candidate$state
  theta_reduced_log <- as.numeric(theta_reduced_log)
  delta <- theta_reduced_log - initial$theta_reduced_log
  score <- as.numeric(initial$U_reduced - initial$A_reduced %*% delta)
  n_subject <- ncol(initial$U_subject_reduced)
  list(
    full_log = probe_full_log(theta_reduced_log, candidate),
    theta_reduced_log = theta_reduced_log,
    objective_provenance = objective_provenance(),
    U_reduced = score,
    U_subject_reduced = matrix(
      rep(score / n_subject, n_subject),
      nrow = length(score)
    ),
    A_reduced = initial$A_reduced,
    raw = list(
      value = candidate$objective -
        as.numeric(crossprod(initial$U_reduced, delta)) +
        0.5 * as.numeric(crossprod(
          delta,
          initial$A_reduced %*% delta
        ))
    )
  )
}

selected <- runtime_candidate("certified_a")
agreeing <- runtime_candidate(
  "certified_b",
  objective = 10 + 1e-10,
  theta_full_log = selected$theta_full_log + c(1e-10, 0, 0, 0)
)
scenarios <- list(
  fast_path = list(
    candidates = list(selected, agreeing, runtime_uncertified(selected, "dominated")),
    expected_strict_route = "R0",
    expected_class = "dominated_uncertified",
    review = FALSE
  ),
  immediate_veto = list(
    candidates = list(selected, agreeing, runtime_uncertified(selected, "boundary")),
    expected_strict_route = "R2",
    expected_class = "competitive_boundary",
    review = FALSE
  ),
  narrow_review = list(
    candidates = list(selected, agreeing, runtime_uncertified(selected, "same_basin")),
    expected_strict_route = "R2",
    expected_class = "same_basin_review_nonveto",
    review = TRUE
  )
)

quantile_summary <- function(x) {
  if (!length(x)) {
    return(c(median = NA_real_, p90 = NA_real_, p99 = NA_real_))
  }
  c(
    median = stats::median(x),
    p90 = unname(stats::quantile(x, 0.90, names = FALSE, type = 8)),
    p99 = unname(stats::quantile(x, 0.99, names = FALSE, type = 8))
  )
}

peak_heap_mb <- function(gc_result) {
  max_used <- match("max used", colnames(gc_result))
  mb_column <- max_used + 1L
  if (is.na(max_used) || mb_column > ncol(gc_result) ||
      !identical(colnames(gc_result)[[mb_column]], "(Mb)")) {
    NA_real_
  } else {
    sum(gc_result[, mb_column])
  }
}

wall_clock_seconds <- function() as.numeric(Sys.time())

run_scenario <- function(name, scenario) {
  run_once <- function(timed_probe = FALSE) {
    probe_elapsed <- 0
    probe_calls <- 0L
    probe <- function(candidate, theta_reduced_log) {
      started <- wall_clock_seconds()
      on.exit({
        probe_elapsed <<- probe_elapsed +
          (wall_clock_seconds() - started)
        probe_calls <<- probe_calls + 1L
      })
      if (!isTRUE(scenario$review)) {
        stop("Nonreview scenario unexpectedly entered the basin probe")
      }
      runtime_evaluator(candidate, theta_reduced_log)
    }
    started <- wall_clock_seconds()
    result <- resolver(
      scenario$candidates,
      resolution_thresholds = resolution_thresholds(),
      basin_probe = probe,
      enable_prospective = FALSE
    )
    elapsed <- wall_clock_seconds() - started
    list(
      result = result,
      elapsed = elapsed,
      probe_elapsed = if (timed_probe) probe_elapsed else 0,
      probe_calls = probe_calls
    )
  }

  warm_one <- run_once()
  warm_two <- run_once()
  stopifnot(identical(warm_one$result, warm_two$result))
  classes <- vapply(
    warm_one$result$resolution$candidate_resolutions,
    `[[`,
    character(1),
    "class"
  )
  stopifnot(
    identical(warm_one$result$resolution$strict_route,
              scenario$expected_strict_route),
    scenario$expected_class %in% classes,
    identical(warm_one$result$resolution$mode, "shadow_default_off")
  )

  invisible(gc(reset = TRUE))
  elapsed <- numeric(repetitions)
  review_elapsed <- numeric(repetitions)
  review_entered <- logical(repetitions)
  counts <- matrix(
    0L,
    nrow = repetitions,
    ncol = 4L,
    dimnames = list(NULL, c("state", "objective", "score", "hessian"))
  )
  last_result <- NULL
  for (index in seq_len(repetitions)) {
    current <- run_once(timed_probe = TRUE)
    last_result <- current$result
    elapsed[[index]] <- current$elapsed
    review_elapsed[[index]] <- current$probe_elapsed
    review_entered[[index]] <- current$probe_calls > 0L
    evaluation_counts <- current$result$resolution$evaluation_counts
    counts[index, ] <- c(
      evaluation_counts$probe_state_evaluations,
      evaluation_counts$probe_objective_evaluations,
      evaluation_counts$probe_score_evaluations,
      evaluation_counts$probe_hessian_evaluations
    )
  }
  heap <- peak_heap_mb(gc())
  list(
    scenario = name,
    repetitions = repetitions,
    resolver_elapsed_seconds = quantile_summary(elapsed),
    review_fraction = mean(review_entered),
    review_elapsed_seconds = quantile_summary(
      review_elapsed[review_entered]
    ),
    probe_evaluation_counts_total = colSums(counts),
    probe_evaluation_counts_per_resolution = colMeans(counts),
    total_engine_counts_available = isTRUE(
      last_result$resolution$evaluation_counts$total_engine_counts_available
    ),
    total_engine_counts_reason =
      last_result$resolution$evaluation_counts$total_engine_counts_reason,
    peak_r_heap_mb = heap,
    peak_memory_scope = "R_gc_heap_proxy_not_process_peak_rss",
    expected_class = scenario$expected_class,
    strict_route = last_result$resolution$strict_route,
    prospective_route = last_result$resolution$prospective_route,
    would_nonveto = last_result$resolution$would_nonveto,
    serial_replay_identical = TRUE,
    raw_elapsed_seconds = elapsed,
    raw_review_elapsed_seconds = review_elapsed
  )
}

started_at <- Sys.time()
scenario_results <- Map(run_scenario, names(scenarios), scenarios)
names(scenario_results) <- names(scenarios)
result <- list(
  metadata = list(
    purpose = "bounded_software_runtime_only_not_inferential_calibration",
    scope = "candidate_classifier_and_bounded_probe_only",
    package_version = as.character(utils::packageVersion("scRoPE")),
    installed_tree_md5_manifest = installed_tree_manifest(),
    resolution_policy = resolution_thresholds()$policy_version,
    resolution_thresholds = resolution_thresholds(),
    repetitions = repetitions,
    platform = R.version$platform,
    r_version = R.version.string,
    started_at = started_at,
    finished_at = Sys.time(),
    total_optimizer_evaluation_counts = "unavailable_by_design",
    peak_memory_measure = "R_gc_heap_proxy_not_process_peak_rss"
  ),
  scenarios = scenario_results
)
saveRDS(result, output_path, version = 3L)
roundtrip <- readRDS(output_path)
stopifnot(identical(roundtrip, result))

summary_table <- do.call(rbind, lapply(scenario_results, function(x) {
  data.frame(
    scenario = x$scenario,
    repetitions = x$repetitions,
    median_seconds = unname(x$resolver_elapsed_seconds[["median"]]),
    p90_seconds = unname(x$resolver_elapsed_seconds[["p90"]]),
    p99_seconds = unname(x$resolver_elapsed_seconds[["p99"]]),
    review_fraction = x$review_fraction,
    median_review_seconds = unname(x$review_elapsed_seconds[["median"]]),
    probe_state_evaluations =
      unname(x$probe_evaluation_counts_total[["state"]]),
    peak_r_heap_mb = x$peak_r_heap_mb,
    stringsAsFactors = FALSE
  )
}))
print(summary_table, row.names = FALSE)
cat("Saved:", normalizePath(output_path, mustWork = TRUE), "\n")
cat("Scope: candidate classifier and bounded probe only; not calibration.\n")
