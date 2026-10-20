# Small, generic software fixtures only.  The 64 fixed counts below are the
# existing solver fixture (seed 271828, L'Ecuyer-CMRG); embedding them avoids
# consuming RNG merely to construct this retention regression.
.hl_retention_ns <- function(name) getFromNamespace(name, "scRoPE")

.hl_retention_case <- function() {
  m <- 16L
  subject <- rep(seq_len(m), each = 4L)
  x <- rep(seq(-1, 1, length.out = m), each = 4L) +
    rep(c(-0.3, -0.1, 0.1, 0.3), m)
  X <- cbind(`(Intercept)` = 1, x = x)
  log_offset <- rep(log(c(0.8, 0.95, 1.1, 1.25)), m)
  y <- c(0, 2, 2, 2, 0, 2, 2, 1, 0, 0, 0, 2, 0, 0, 0, 3,
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 2, 2,
         2, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0,
         0, 0, 1, 1, 0, 0, 1, 0, 3, 1, 0, 11, 0, 1, 3, 4)
  list(y = y, X = X, subject = subject, log_offset = log_offset,
       data = .hl_retention_ns("hl_fo_make_data")(
         y, X, log_offset, subject, intcol = 1L
       ),
       seed = list(beta = c(4, -4), rho_full = log(c(0.45, 6)),
                   coordinate_tag = "canonical_conditional_beta_tau_xi"),
       fixed = c(sigma2 = 0.45, phi = 6))
}

# Only history-bearing fields are excluded here.  Numerical results, root and
# covariance certificates, reasons, constraints, per-subject derivatives and
# reference state remain available for bit-for-bit comparison.
.hl_retention_science <- function(x) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  optional <- c("trace", "outer_trace", "outer_runs", "candidates",
                "root_trace", "null_runs", "root_candidates",
                "trust_trace", "polish_trace", "evaluator_trace",
                "boundary_rounds")
  if (!is.null(names(x))) x <- x[!names(x) %in% optional]
  lapply(x, .hl_retention_science)
}

.hl_retention_call <- function(name, args, keep_diagnostics, legacy = FALSE) {
  if (!legacy || identical(name, "scrope_hl")) {
    args$keep_diagnostics <- keep_diagnostics
  }
  do.call(.hl_retention_ns(name), args)
}

.hl_retention_observe_public <- function(args, keep_diagnostics) {
  namespace <- asNamespace("scRoPE")
  snapshots <- list()
  completed_fits <- list()
  observer <- function(event, fit = NULL) {
    frames <- sys.frames()
    eligible <- vapply(frames, function(frame) {
      exists("gene_fits", frame, inherits = FALSE) &&
        exists("coefficient_wald", frame, inherits = FALSE) &&
        exists("gg", frame, inherits = FALSE)
    }, logical(1))
    if (!any(eligible)) return(invisible(NULL))
    frame <- frames[[tail(which(eligible), 1L)]]
    objects <- mget(c("gene_fits", "test_objects", "coefficient_wald"), frame)
    snapshots[[length(snapshots) + 1L]] <<- data.frame(
      event = event, gene_index = get("gg", frame),
      fits_length = length(objects$gene_fits),
      fits_completed = sum(!vapply(objects$gene_fits, is.null, logical(1))),
      tests_length = length(objects$test_objects),
      wald_length = length(objects$coefficient_wald),
      accumulator_bytes = as.numeric(object.size(objects)),
      stringsAsFactors = FALSE
    )
    if (identical(event, "fit_exit")) {
      completed_fits[[length(completed_fits) + 1L]] <<- fit
    }
    invisible(NULL)
  }
  old_option <- options(scrope.retention.test.observer = observer)
  on.exit(options(old_option), add = TRUE)
  invisible(capture.output(trace(
    "hl_fo_fit", where = namespace, print = FALSE,
    tracer = quote(getOption("scrope.retention.test.observer")("fit_entry")),
    exit = quote(getOption("scrope.retention.test.observer")(
      "fit_exit", returnValue()
    ))
  )))
  on.exit(invisible(capture.output(untrace("hl_fo_fit", where = namespace))),
          add = TRUE)
  invisible(capture.output(trace(
    "hl_fo_public_scrope_hl", where = namespace, print = FALSE,
    exit = quote(getOption("scrope.retention.test.observer")("public_exit"))
  )))
  on.exit(invisible(capture.output(untrace(
    "hl_fo_public_scrope_hl", where = namespace
  ))), add = TRUE)
  result <- .hl_retention_call("scrope_hl", args, keep_diagnostics)
  list(result = result, snapshots = do.call(rbind, snapshots),
       fits = completed_fits)
}

# One bounded battery; the evidence script reuses exactly these definitions.
.hl_retention_collect <- function(keep_diagnostics, legacy = FALSE) {
  case <- .hl_retention_case()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else if (
      exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    ) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(97231L)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  call <- function(name, args) {
    .hl_retention_call(name, args, keep_diagnostics, legacy)
  }

  fixed <- call("hl_fo_fit", list(
    data = case$data, active = character(0),
    fixed_natural = case$fixed, seed = case$seed
  ))
  score <- call("hl_fo_score_test", list(
    data = case$data, R = matrix(c(0, 1), 1L), rhs = 0.3,
    active = character(0), fixed_natural = case$fixed,
    unrestricted_fit = fixed
  ))
  make_data <- .hl_retention_ns("hl_fo_make_data")
  zero_data <- make_data(rep(0, 12), matrix(1, 12, 1L),
                         numeric(12), rep(1:4, each = 3L), intcol = 1L)
  zero <- call("hl_fo_fit", list(
    data = zero_data, active = character(0), fixed_natural = case$fixed,
    seed = list(beta = 0, rho_full = case$seed$rho_full,
                coordinate_tag = case$seed$coordinate_tag)
  ))
  singular_data <- make_data(rep(1, 12), matrix(1, 12, 1L),
                             numeric(12), rep(1:4, each = 3L), intcol = 1L)
  singular <- call("hl_fo_fit", list(
    data = singular_data, active = character(0), fixed_natural = case$fixed,
    seed = list(beta = 0, rho_full = case$seed$rho_full,
                coordinate_tag = case$seed$coordinate_tag)
  ))
  rank_data <- make_data(case$y, cbind(case$X, duplicate = case$X[, 2L]),
                         case$log_offset, case$subject, intcol = 1L)
  rank <- call("hl_fo_fit", list(
    data = rank_data, active = character(0), fixed_natural = case$fixed
  ))

  evaluate_twice <- function(name, args) {
    evaluator <- call(name, args)
    first <- evaluator$evaluate(case$seed$rho_full, cold = TRUE)
    second <- evaluator$evaluate(case$seed$rho_full, cold = TRUE)
    list(first = first, second = second, count = evaluator$count(),
         history = evaluator$trace())
  }
  outer <- evaluate_twice("hl_fo_make_outer_evaluator", list(
    data = case$data, rho_template = case$seed$rho_full, active = 1:2,
    beta_star = numeric(2), K = diag(2), seed_beta = case$seed$beta
  ))
  affine <- .hl_retention_ns("hl_fo_affine_constraint")(
    matrix(c(0, 1), 1L), rhs = 0.3, p = 2L
  )
  score_evaluator <- evaluate_twice("hl_fo_score_make_evaluator", list(
    data = case$data, rho_template = case$seed$rho_full, active = 1:2,
    beta_star = affine$beta_star, K = affine$K,
    seed_beta = case$seed$beta,
    constraint_R = matrix(c(0, 1), 1L), constraint_rhs = 0.3
  ))
  public <- .hl_retention_observe_public(list(
    count = rbind(generic_a = case$y, generic_b = case$y,
                   zero_adverse = rep(0, length(case$y))),
    id = case$subject, pred = case$X, offset = exp(case$log_offset),
    cpc = -1, mincp = 0, ncore = 1L, verbose = FALSE,
    additional_tests = c("score", "lrt"), lrt_details = TRUE, output_re = TRUE
  ), keep_diagnostics)
  list(fixed = fixed, score = score, zero = zero, singular = singular,
       rank = rank, outer = outer, score_evaluator = score_evaluator,
       public = public,
       rng_unchanged = identical(seed_before, .Random.seed) &&
         identical(kind_before, RNGkind()))
}

test_that("HL output retention preserves computations and avoids gene accumulation", {
  full <- .hl_retention_collect(TRUE)
  compact <- .hl_retention_collect(FALSE)

  expect_true(full$rng_unchanged)
  expect_true(compact$rng_unchanged)
  expect_true(full$fixed$point_success)
  expect_true(full$fixed$covariance_success)
  expect_true(full$score$success)
  expect_identical(full$zero$fit_status, "inner_no_finite_root")
  expect_identical(full$rank$fit_status, "inner_rank_deficient")
  expect_true(full$singular$point_success)
  expect_false(full$singular$covariance_success)

  for (name in c("fixed", "score", "zero", "singular", "rank")) {
    expect_identical(.hl_retention_science(compact[[name]]),
                     .hl_retention_science(full[[name]]), info = name)
  }
  expect_gt(length(full$fixed$inner$trace), 0L)
  expect_length(compact$fixed$inner$trace, 0L)
  expect_gt(length(full$fixed$outer_trace), 0L)
  expect_length(compact$fixed$outer_trace, 0L)
  expect_lt(as.numeric(object.size(compact$fixed)),
            as.numeric(object.size(full$fixed)))
  expect_identical(compact$fixed$certification, full$fixed$certification)
  expect_identical(compact$fixed$covariance, full$fixed$covariance)
  expect_identical(compact$fixed$inner$mode_diagnostics,
                   full$fixed$inner$mode_diagnostics)
  expect_identical(compact$score$root_certification,
                   full$score$root_certification)

  for (name in c("outer", "score_evaluator")) {
    expect_true(full[[name]]$first$success, info = name)
    expect_true(compact[[name]]$first$success, info = name)
    expect_equal(full[[name]]$count, 2L, info = name)
    expect_identical(compact[[name]]$count, full[[name]]$count, info = name)
    expect_length(full[[name]]$history, 2L)
    expect_length(compact[[name]]$history, 0L)
    expect_identical(.hl_retention_science(compact[[name]]$first),
                     .hl_retention_science(full[[name]]$first), info = name)
    expect_identical(.hl_retention_science(compact[[name]]$second),
                     .hl_retention_science(full[[name]]$second), info = name)
  }

  expect_null(compact$public$result$diagnostics)
  expected <- full$public$result
  expected$diagnostics <- NULL
  observed <- compact$public$result
  observed$diagnostics <- NULL
  expect_identical(observed, expected)
  expect_identical(nrow(compact$public$result$summary), 3L)
  expect_identical(compact$public$result$summary$fit_status[[3L]],
                   "inner_no_finite_root")
  expect_true(all(full$public$result$summary$hl_score_success_x[1:2]))
  expect_false(any(full$public$result$summary$hl_lrt_success_x[1:2]))
  expect_true(all(full$public$result$summary$hl_lrt_failure_x[1:2] ==
                    "objective_ratio_unsupported"))
  expect_lt(as.numeric(object.size(compact$public$result)),
            as.numeric(object.size(full$public$result)))
  # Three per-gene fits plus two real constrained fits for the active Score
  # calls; the all-zero gene fails before a Score call is attempted.
  expect_equal(length(full$public$fits), 5L)
  expect_equal(length(compact$public$fits), 5L)
  for (ii in seq_along(full$public$fits)) {
    expect_identical(.hl_retention_science(compact$public$fits[[ii]]),
                     .hl_retention_science(full$public$fits[[ii]]))
    if (length(full$public$fits[[ii]]$outer_runs)) {
      run_summary <- function(run) {
        fields <- c("success", "status", "reason", "start", "evaluations",
                    "stationarity", "invalid_guard_faces")
        run[intersect(fields, names(run))]
      }
      expect_identical(lapply(compact$public$fits[[ii]]$outer_runs, run_summary),
                       lapply(full$public$fits[[ii]]$outer_runs, run_summary))
      expect_true(all(vapply(compact$public$fits[[ii]]$outer_runs,
                             function(run) length(run$trust_trace) == 0L,
                             logical(1))))
      expect_true(all(vapply(compact$public$fits[[ii]]$outer_runs,
                             function(run) is.null(run$state$subjects) &&
                               is.null(run$state$augmented), logical(1))))
    }
  }
  expect_true(all(compact$public$snapshots$fits_length == 0L))
  expect_true(all(compact$public$snapshots$tests_length == 0L))
  expect_true(all(compact$public$snapshots$wald_length == 0L))
  expect_true(all(compact$public$snapshots$fits_completed == 0L))
  end_full <- tail(full$public$snapshots, 1L)
  end_compact <- tail(compact$public$snapshots, 1L)
  expect_identical(end_full$event, "public_exit")
  expect_equal(end_full$fits_completed, 3L)
  expect_lt(end_compact$accumulator_bytes, end_full$accumulator_bytes)
})
