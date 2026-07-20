#!/usr/bin/env Rscript

# Active estimator-level gate for the versioned difficult-case fixture
# contract.  This is deliberately separate from testthat because each fit runs
# the complete multistart nested solver and cold replay.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
output <- if (length(args) >= 2L) {
  normalizePath(args[[2L]], mustWork = FALSE)
} else {
  file.path(repo_root, "temp", "hl_first_order_validation", "fixture_contract_v2")
}
case_filter <- if (length(args) >= 3L && nzchar(args[[3L]])) {
  strsplit(args[[3L]], ",", fixed = TRUE)[[1L]]
} else {
  NULL
}
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Repository root does not contain DESCRIPTION: ", repo_root)
}

source(
  file.path(repo_root, "tests", "testthat", "helper-hl-first-order.R"),
  local = FALSE
)
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required for the active fixture contract gate.")
}
pkgload::load_all(repo_root, recompile = FALSE, helpers = FALSE, quiet = TRUE)
make_data <- getFromNamespace("hl_fo_make_data", "scRoPE")
fit_first_order <- getFromNamespace("hl_fo_fit", "scRoPE")
terminal_boundary_faces <- getFromNamespace(
  "hl_fo_terminal_boundary_faces", "scRoPE"
)

dir.create(output, recursive = TRUE, showWarnings = FALSE)
fit_output <- file.path(output, "fits")
dir.create(fit_output, recursive = TRUE, showWarnings = FALSE)

d1 <- hl_fo_fixture("D1", fresh = TRUE)
d2 <- hl_fo_fixture("D2", fresh = TRUE)
cases <- list(
  Z1 = hl_fo_fixture("Z1", fresh = TRUE),
  D1_low_v_high_phi = d1$cases$low_v_high_phi,
  D1_low_v_low_phi = d1$cases$low_v_low_phi,
  D1_high_v_high_phi = d1$cases$high_v_high_phi,
  D1_high_v_low_phi = d1$cases$high_v_low_phi,
  Z2 = hl_fo_fixture("Z2", fresh = TRUE),
  D2_low_v_high_phi = d2$cases$low_v_high_phi,
  D2_low_v_low_phi = d2$cases$low_v_low_phi
)
if (!is.null(case_filter)) {
  unknown <- setdiff(case_filter, names(cases))
  if (length(unknown)) {
    stop("Unknown fixture contract case: ", paste(unknown, collapse = ", "))
  }
  cases <- cases[case_filter]
}

.scalar_character <- function(x) {
  if (is.null(x) || !length(x)) NA_character_ else as.character(x[[1L]])
}

.scalar_numeric <- function(x) {
  if (is.null(x) || !length(x) || !is.finite(x[[1L]])) NA_real_ else
    as.numeric(x[[1L]])
}

.cold_max <- function(fit) {
  if (is.null(fit$cold_replay)) return(Inf)
  fields <- c("beta", "eta", "F", "g", "H_F", "A", "J")
  values <- unlist(fit$cold_replay[fields], use.names = FALSE)
  if (!length(values) || any(!is.finite(values))) Inf else max(values)
}

.relative_guard_distance <- function(fit) {
  if (is.null(fit$guard_bounds) || is.null(fit$rho)) return(NA_real_)
  lower <- fit$guard_bounds$lower
  upper <- fit$guard_bounds$upper
  if (length(lower) != length(fit$rho) || any(upper <= lower)) return(NA_real_)
  min(c((fit$rho - lower) / (upper - lower),
        (upper - fit$rho) / (upper - lower)))
}

.fit_case <- function(label, fixture) {
  data <- make_data(
    y = fixture$y,
    X = fixture$X,
    log_offset = fixture$log_pi,
    subject = fixture$subject_index,
    intcol = fixture$intcol
  )
  seed <- list(
    beta = fixture$truth$beta,
    rho_full = c(fixture$truth$tau, fixture$truth$xi),
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    fit_first_order(data = data, seed = seed),
    error = function(error) list(
      success = FALSE,
      point_success = FALSE,
      covariance_success = FALSE,
      fit_status = "unhandled_error",
      fit_failure_reason = conditionMessage(error)
    )
  )
  elapsed <- proc.time()[["elapsed"]] - started
  saveRDS(
    list(fixture = fixture, fit = fit, elapsed_seconds = elapsed),
    file.path(fit_output, paste0(label, ".rds")),
    version = 3L
  )

  filter <- hl_fo_public_filter_summary(fixture)
  point_success <- isTRUE(fit$point_success)
  covariance_success <- isTRUE(fit$covariance_success)
  fit_status <- .scalar_character(fit$fit_status)
  expected <- .scalar_character(fixture$expected)
  guard_distance <- .relative_guard_distance(fit)
  cold_max <- .cold_max(fit)
  hessian_eigenvalues <- if (!is.null(fit$H_F) && length(fit$H_F)) {
    eigen(0.5 * (fit$H_F + t(fit$H_F)), symmetric = TRUE,
          only.values = TRUE)$values
  } else {
    numeric(0)
  }
  hessian_negative <- length(hessian_eigenvalues) == 2L &&
    all(is.finite(hessian_eigenvalues)) && max(hessian_eigenvalues) < 0
  A_rank <- if (is.null(fit$augmented$A)) NA_integer_ else
    hl_fo_num_rank(fit$augmented$A)$rank
  J_rank <- if (is.null(fit$augmented$J)) NA_integer_ else
    hl_fo_num_rank(fit$augmented$J)$rank
  A_dim <- if (is.null(fit$augmented$A)) NA_integer_ else nrow(fit$augmented$A)
  J_dim <- if (is.null(fit$augmented$J)) NA_integer_ else nrow(fit$augmented$J)
  stationarity <- fit$cold_replay$stationarity
  boundary_faces <- if (is.null(fit$boundary_faces)) character(0) else
    sort(unique(fit$boundary_faces))
  boundary_round_count <- length(fit$boundary_rounds)
  all_final_starts_cross_outward <- FALSE
  if (boundary_round_count) {
    final_round <- fit$boundary_rounds[[boundary_round_count]]
    final_start_faces <- lapply(
      final_round$runs,
      terminal_boundary_faces,
      lower = final_round$bounds$lower,
      upper = final_round$bounds$upper
    )
    all_final_starts_cross_outward <-
      length(final_start_faces) > 0L && all(lengths(final_start_faces) > 0L)
  }

  if (identical(expected, "outer_boundary")) {
    boundary_contract <- fixture$boundary_contract
    contract_pass <- filter$passes && !point_success &&
      identical(fit_status, "outer_boundary") &&
      identical(boundary_faces, sort(boundary_contract$faces)) &&
      identical(length(fit$guard_expansions), boundary_contract$expansion_count) &&
      identical(boundary_round_count, boundary_contract$round_count) &&
      all_final_starts_cross_outward
  } else if (identical(expected, "certified_interior_success")) {
    contract_pass <- filter$passes && point_success && covariance_success &&
      fit_status %in% c("ok", "multiple_near_tied_maxima") &&
      isTRUE(stationarity$success) && is.finite(guard_distance) &&
      guard_distance > 1e-6 && hessian_negative &&
      identical(A_rank, A_dim) && identical(J_rank, J_dim) &&
      is.finite(cold_max) && cold_max <= 1e-8
  } else {
    contract_pass <- FALSE
  }

  data.frame(
    contract_version = fixture$contract_version,
    case = label,
    fixture_id = fixture$id,
    seed = fixture$seed,
    m = fixture$m,
    n_i = fixture$n_i,
    total_count = filter$total_count,
    nonzero_cells = filter$nonzero_cells,
    count_per_cell = filter$count_per_cell,
    public_filter_pass = filter$passes,
    spec_expected = .scalar_character(fixture$spec_expected),
    expected = expected,
    fit_status = fit_status,
    fit_reason = .scalar_character(
      if (is.null(fit$fit_failure_reason)) fit$reason else fit$fit_failure_reason
    ),
    point_success = point_success,
    covariance_success = covariance_success,
    sigma2_hat = if (point_success) exp(fit$rho_full[[1L]]) else NA_real_,
    phi_hat = if (point_success) exp(fit$rho_full[[2L]]) else NA_real_,
    gradient_max = if (point_success) max(abs(fit$g)) else NA_real_,
    hessian_eigen_max = if (length(hessian_eigenvalues))
      max(hessian_eigenvalues) else NA_real_,
    relative_guard_distance = guard_distance,
    boundary_faces = if (length(boundary_faces))
      paste(boundary_faces, collapse = "|") else NA_character_,
    boundary_rounds = boundary_round_count,
    all_final_starts_cross_outward = all_final_starts_cross_outward,
    guard_expansions = length(fit$guard_expansions),
    cold_max = if (is.finite(cold_max)) cold_max else NA_real_,
    A_rank = A_rank,
    A_dim = A_dim,
    J_rank = J_rank,
    J_dim = J_dim,
    contract_pass = contract_pass,
    elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )
}

rows <- vector("list", length(cases))
for (ii in seq_along(cases)) {
  label <- names(cases)[[ii]]
  cat("Running", label, "...\n")
  rows[[ii]] <- .fit_case(label, cases[[ii]])
  print(rows[[ii]], row.names = FALSE)
}
results <- do.call(rbind, rows)
write.csv(results, file.path(output, "contract_results.csv"), row.names = FALSE)
saveRDS(
  list(
    artifact = "first_order_hl_difficult_fixture_contract_v2",
    contract_version = hl_fo_fixture_contract_version(),
    package_version = unname(read.dcf(
      file.path(repo_root, "DESCRIPTION"), fields = "Version"
    )[[1L]]),
    results = results
  ),
  file.path(output, "contract_manifest.rds"),
  version = 3L
)

if (!all(results$contract_pass)) {
  failed <- results$case[!results$contract_pass]
  stop("Fixture contract failed: ", paste(failed, collapse = ", "))
}
cat("All versioned difficult-case fixture contracts passed.\n")
