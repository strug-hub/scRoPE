history_native <- function(x) {
  f <- suppressWarnings(scrope(x$y, x$id, x$X, offset = x$offset,
    cpc = 0, mincp = 1, verbose = FALSE, recovery = "none", keep_diagnostics = TRUE))
  list(success = FALSE, joint = f$diagnostics$raw[[1]]$ln_joint)
}

test_that("retained recovery endpoints are reevaluated and provenance checked", {
  x <- readRDS(test_path("fixtures", "jln-face-competitor.rds"))
  primary <- history_native(x)
  dat <- scRoPE:::ln_research_context(x$y, x$X, x$id, x$offset)
  records <- primary$joint$candidate_records
  best <- records[[which.min(vapply(records, `[[`, numeric(1), "objective"))]]
  original <- records[[1L]]
  primary$joint$candidate_records[[1L]]$recovery <- list(
    original = original, certificate = best,
    optimizer = list(par = best$theta_reduced, value = best$objective),
    polish = list(theta_reduced = best$theta_reduced, objective = best$objective),
    adopted = FALSE, status = "certificate_failed")
  audit <- scRoPE:::ln_recovery_known_points(primary, dat)
  expect_true(audit$success)
  expect_true(all(c("S_direct:recovery_original", "S_direct:recovery_certificate",
    "S_direct:recovery_optimizer", "S_direct:recovery_polish") %in% audit$observations$label))
  expect_equal(min(audit$observations$q), best$objective, tolerance = 1e-10)
  # Inconsistent stored scalars cannot conceal a better point or invent a competitor.
  primary$joint$candidate_records[[1L]]$recovery$certificate$objective <- 1e12
  reevaluated <- scRoPE:::ln_recovery_known_points(primary, dat)
  expect_false(reevaluated$success)
  expect_identical(reevaluated$reason, "retained_objective_audit_failed")
  primary$joint$candidate_records[[1L]]$recovery$certificate$objective <- best$objective
  primary$joint$candidate_records[[1L]]$recovery$certificate$objective_instance <- NULL
  invalid <- scRoPE:::ln_recovery_known_points(primary, dat)
  expect_false(invalid$success)
  expect_identical(invalid$reason, "retained_objective_audit_failed")
  primary$joint$candidate_records[[1L]]$recovery$certificate <- best
  primary$joint$candidate_records[[1L]]$recovery$optimizer$par <- 1
  expect_false(scRoPE:::ln_recovery_known_points(primary, dat)$success)
})

test_that("an unadopted recovery competitor vetoes a local limiting root", {
  x <- readRDS(test_path("fixtures", "jln-face-competitor.rds"))
  primary <- history_native(x)
  records <- primary$joint$candidate_records
  best <- records[[which.min(vapply(records, `[[`, numeric(1), "objective"))]]
  # Retain only the worse local trajectory as active; preserve an actual
  # same-dataset better endpoint in a failed native-recovery certificate.
  primary$joint$candidate_records <- records[1L]
  primary$joint$candidate_records[[1L]]$recovery <- list(
    original = records[[1L]], certificate = best,
    adopted = FALSE, status = "certificate_failed")
  f <- scRoPE:::ln_recovery_route(primary, x$y, x$X, x$id, x$offset)
  expect_identical(f$route, "withheld")
  expect_identical(f$result$reason, "reduced_root_worse_than_known_full_candidate")
  audit <- f$face_audit$competitor_audit$observations
  expect_true(any(!audit$passes & grepl("recovery_certificate", audit$label)))
})

test_that("warning anchors use actual native recovery history", {
  x <- readRDS(test_path("fixtures", "jln-research.rds"))$warning$input
  primary <- history_native(x)
  dat <- scRoPE:::ln_research_context(x$y, x$X, x$id, x$offset)
  expect_true(scRoPE:::ln_research_warning_fit(primary, dat)$success)
  records <- primary$joint$candidate_records
  j <- which(vapply(records, function(r) isTRUE(r$core) && isTRUE(r$certified), logical(1)))[1L]
  primary$joint$candidate_records[[j]]$recovery <- list(adopted = TRUE, status = "adopted")
  expect_identical(scRoPE:::ln_research_warning_fit(primary, dat)$reason,
    "insufficient_distinct_original_core_starts")
  primary$joint$candidate_records[[j]]$recovery <- list(adopted = FALSE,
    status = "certificate_failed", certificate = records[[j]])
  primary$joint$candidate_records[[j]]$recovery$certificate$objective_instance <- NULL
  expect_identical(scRoPE:::ln_research_warning_fit(primary, dat)$reason,
    "invalid_retained_observations")
})

test_that("all-free refinement refuses unverifiable historical evidence", {
  x <- readRDS(test_path("fixtures", "jln-research.rds"))$ordinary$input
  primary <- history_native(x)
  dat <- scRoPE:::ln_research_context(x$y, x$X, x$id, x$offset)
  r <- primary$joint$candidate_records[[1L]]
  r$theta_reduced_log[1L] <- r$theta_reduced_log[1L] + 1e-5
  r$certified <- FALSE; r$checks["newton"] <- FALSE
  primary$joint$candidate_records[[1L]] <- r
  expect_true(scRoPE:::ln_recovery_refine(primary, dat)$success)
  r$objective_instance <- NULL
  primary$joint$candidate_records[[1L]]$recovery <- list(certificate = r)
  f <- scRoPE:::ln_recovery_refine(primary, dat)
  expect_false(f$success)
  expect_identical(f$reason, "refinement_retained_competitor_not_resolved")
})
