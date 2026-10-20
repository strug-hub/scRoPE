# Controlled generic endpoints test software behavior, not scientific availability.
ln_recovery_test_case <- function() {
  prep <- ln_joint_prepare_fixture(ln_joint_test_fixture())
  map <- ln_joint_map(prep$ctx$nb)
  instance <- ln_joint_objective_instance_registry("register", 1L, prep$posv, prep$ctx)
  theta <- c(rep(0, prep$ctx$nb), 1, 1)
  candidate <- ln_joint_candidate_state(
    theta, ln_joint_objective(theta, 1L, prep$posv, prep$ctx, map),
    list(par = theta, convergence = 7L, message = "controlled initial endpoint"),
    "generic_core", "direct", TRUE, TRUE, 1L, prep$posv, prep$ctx, map, instance
  )
  candidate$start <- theta + c(rep(0.1, prep$ctx$nb), 0, 0)
  candidate$optimizer_endpoint <- theta
  candidate$optimizer_objective <- candidate$objective
  candidate$polished <- FALSE
  list(candidate = candidate, gene_index = 1L, posv = prep$posv, ctx = prep$ctx,
       map = map, lower = c(rep(-100, map$free_dim), prep$ctx$min),
       upper = c(rep(100, map$free_dim), prep$ctx$max), objective_instance = instance)
}

ln_recovery_effective <- function(x) x[setdiff(names(x), "recovery")]

test_that("one recovery uses the exact endpoint and a fresh complete certificate", {
  args <- ln_recovery_test_case()
  original_lbfgs <- nloptr::lbfgs
  original_polish <- ln_joint_polish
  original_certificate <- ln_joint_candidate_state
  calls <- c(optimizer = 0L, polish = 0L, certificate = 0L)
  local_mocked_bindings(
    lbfgs = function(x0, fn, lower, upper, control) {
      calls[["optimizer"]] <<- calls[["optimizer"]] + 1L
      expect_identical(x0, args$candidate$theta_reduced)
      expect_false(identical(x0, args$candidate$start))
      expect_identical(lower, args$lower)
      expect_identical(upper, args$upper)
      expect_identical(control, list(ftol_abs = args$ctx$eps, ftol_rel = 0,
                                    xtol_rel = 1e-6, maxeval = 1000L))
      original_lbfgs(x0, fn, lower = lower, upper = upper, control = control)
    },
    ln_joint_polish = function(..., max_iter) {
      calls[["polish"]] <<- calls[["polish"]] + 1L
      expect_identical(max_iter, 25L)
      original_polish(..., max_iter = max_iter)
    },
    ln_joint_candidate_state = function(..., objective_instance) {
      calls[["certificate"]] <<- calls[["certificate"]] + 1L
      expect_identical(objective_instance, args$objective_instance)
      original_certificate(..., objective_instance = objective_instance)
    }, .package = "scRoPE"
  )
  result <- do.call(ln_joint_recover_candidate, args)
  expect_identical(calls, c(optimizer = 1L, polish = 1L, certificate = 1L))
  expect_true(result$certified)
  expect_true(all(result$checks))
  expect_identical(result$recovery$status, "adopted")
  expect_identical(result$recovery$original, args$candidate[setdiff(names(args$candidate), "state")])
  expect_false(result$recovery$original$certified)
  expect_false(all(result$recovery$original$checks))
  expect_equal(result$state$theta_reduced_log, result$theta_reduced_log, tolerance = 0)
  expect_false(identical(result$theta_reduced, args$candidate$theta_reduced))
  fields <- c("label", "start_provenance", "required", "core", "start",
              "objective_provenance", "objective_instance")
  expect_identical(result[fields], args$candidate[fields])
  expect_lte(result$objective, args$candidate$objective + result$recovery$objective_allowance)
  args$candidate <- result
  expect_identical(do.call(ln_joint_recover_candidate, args), result)
  expect_identical(calls, c(optimizer = 1L, polish = 1L, certificate = 1L))
  record <- ln_joint_candidate_records(list(result))[[1L]]
  expect_identical(record$recovery, result$recovery)
  expect_identical(ln_joint_resolve_candidates(list(result))$failure_reason,
                   "insufficient_certified_core_starts")
  distinct <- args$candidate; distinct$label <- "distinct_failed_start"
  distinct$finite_endpoint <- FALSE; distinct$certified <- FALSE; distinct$objective <- NA_real_
  expect_identical(ln_joint_resolve_candidates(list(result, distinct))$route, "R2")
})

test_that("ineligible roles, maps, metadata and derivative failures bypass recovery", {
  args <- ln_recovery_test_case()
  local_mocked_bindings(lbfgs = function(...) stop("must not optimize"), .package = "scRoPE")
  variants <- list(
    list(required = FALSE), list(core = FALSE), list(start_provenance = "pql_start"),
    list(certified = TRUE), list(finite_endpoint = FALSE), list(state = NULL),
    list(objective = Inf), list(theta_reduced = c(NA_real_, 0, 0, 1, 1)),
    list(theta_reduced = c(101, 0, 0, 1, 1)), list(objective_provenance = list()),
    list(objective_instance = NULL), list(objective_instance = list(id = "wrong")),
    list(checks = c(root = FALSE)), list(failure = "derivative_evaluation_failed", state = NULL)
  )
  for (change in variants) {
    a <- args
    for (name in names(change)) a$candidate[name] <- change[name]
    expect_identical(do.call(ln_joint_recover_candidate, a), a$candidate)
  }
  a <- args; a$map$constrained <- TRUE
  expect_identical(do.call(ln_joint_recover_candidate, a), a$candidate)
  a <- args; a$upper[[1L]] <- 200
  expect_identical(do.call(ln_joint_recover_candidate, a), a$candidate)
  a <- args; a$candidate$state$U_subject_reduced[1, 1] <- 100
  expect_identical(do.call(ln_joint_recover_candidate, a), a$candidate)
  a <- args; a$candidate$constraint_residual <- 1
  expect_identical(do.call(ln_joint_recover_candidate, a), a$candidate)
})

test_that("a real boundary certificate is excluded", {
  a <- ln_recovery_test_case(); theta <- a$candidate$theta_reduced
  theta[[length(theta)]] <- a$ctx$max[[2L]]
  a$candidate <- ln_joint_candidate_state(
    theta, ln_joint_objective(theta, 1L, a$posv, a$ctx, a$map), list(convergence = 0L),
    "generic_boundary", "direct", TRUE, TRUE, 1L, a$posv, a$ctx, a$map, a$objective_instance
  )
  expect_false(a$candidate$checks[["interior"]])
  local_mocked_bindings(lbfgs = function(...) stop("must not optimize"), .package = "scRoPE")
  expect_identical(do.call(ln_joint_recover_candidate, a), a$candidate)
})

test_that("errors and unusable optimizer returns retain the original and cannot retry", {
  a <- ln_recovery_test_case()
  for (kind in c("error", "nonfinite", "out_of_bounds", "wrong_length")) {
    calls <- 0L
    result <- with_mocked_bindings({
      first <- do.call(ln_joint_recover_candidate, a)
      b <- a; b$candidate <- first
      expect_identical(do.call(ln_joint_recover_candidate, b), first)
      first
    }, lbfgs = function(...) {
      calls <<- calls + 1L
      if (kind == "error") stop("controlled failure")
      par <- a$candidate$theta_reduced
      if (kind == "nonfinite") par[[1L]] <- NA_real_
      if (kind == "out_of_bounds") par[[1L]] <- 101
      if (kind == "wrong_length") par <- par[-1L]
      list(par = par, convergence = -1L, message = "controlled return")
    }, .package = "scRoPE")
    expect_identical(calls, 1L)
    expect_identical(ln_recovery_effective(result), a$candidate)
    expect_true(result$recovery$attempted)
    expect_false(result$recovery$adopted)
    expect_identical(result$recovery$status, if (kind == "error") "optimizer_error" else "unusable_optimizer_endpoint")
  }
})

test_that("polish and certification errors, nonfinite values and failed checks retain failure", {
  a <- ln_recovery_test_case()
  for (stage in c("objective", "polish", "certificate")) {
    mocks <- list(lbfgs = function(...) list(par = a$candidate$theta_reduced))
    name <- switch(stage, objective = "ln_joint_objective", polish = "ln_joint_polish",
                   certificate = "ln_joint_candidate_state")
    mocks[[name]] <- function(...) stop("controlled stage error")
    result <- do.call(with_mocked_bindings, c(list(quote(do.call(ln_joint_recover_candidate, a))), mocks, list(.package = "scRoPE")))
    expect_identical(ln_recovery_effective(result), a$candidate)
    expect_false(result$recovery$adopted)
    expect_match(result$recovery$status, "_error$")
  }
  result <- with_mocked_bindings(do.call(ln_joint_recover_candidate, a),
    lbfgs = function(...) list(par = a$candidate$theta_reduced),
    ln_joint_polish = function(...) list(theta_reduced = a$candidate$theta_reduced, objective = NA_real_),
    .package = "scRoPE")
  expect_identical(result$recovery$status, "unusable_polished_endpoint")
  expect_identical(ln_recovery_effective(result), a$candidate)
  result <- with_mocked_bindings(do.call(ln_joint_recover_candidate, a),
    lbfgs = function(...) list(par = a$candidate$theta_reduced),
    ln_joint_polish = function(...) list(theta_reduced = a$candidate$theta_reduced, objective = a$candidate$objective),
    .package = "scRoPE")
  expect_identical(result$recovery$status, "certificate_failed")
  expect_false(result$recovery$certificate$certified)
  expect_identical(ln_recovery_effective(result), a$candidate)
  result <- with_mocked_bindings(do.call(ln_joint_recover_candidate, a),
    lbfgs = function(...) list(par = a$candidate$theta_reduced),
    ln_joint_objective = function(...) Inf, .package = "scRoPE")
  expect_identical(result$recovery$status, "nonfinite_optimizer_objective")
  expect_identical(ln_recovery_effective(result), a$candidate)
})

test_that("a certified recovery must respect the original objective allowance", {
  a <- ln_recovery_test_case()
  # Deliberately inconsistent objective values isolate the adoption guard only.
  for (factor in c(0.5, 2)) {
    ans <- with_mocked_bindings(do.call(ln_joint_recover_candidate, a),
      ln_joint_candidate_state = function(...) {
        result <- a$candidate
        result$certified <- TRUE; result$checks[] <- TRUE
        result$objective <- a$candidate$objective + factor * 1e-12 * (1 + abs(a$candidate$objective))
        result
      }, .package = "scRoPE")
    expect_identical(ans$recovery$adopted, factor < 1)
    expect_identical(ans$recovery$status, if (factor < 1) "adopted" else "objective_increase")
    if (factor > 1) expect_identical(ln_recovery_effective(ans), a$candidate)
  }
})

test_that("termination codes are diagnostic and warnings are retained", {
  a <- ln_recovery_test_case(); original_lbfgs <- nloptr::lbfgs
  result <- NULL
  expect_warning(result <- with_mocked_bindings(do.call(ln_joint_recover_candidate, a),
    lbfgs = function(...) {
      warning("controlled warning")
      ans <- original_lbfgs(...); ans$convergence <- -1L; ans
    }, .package = "scRoPE"), "controlled warning")
  expect_true(result$certified)
  expect_true(result$recovery$adopted)
  expect_identical(result$recovery$optimizer$convergence, -1L)
  expect_identical(result$recovery$warnings, "controlled warning")
})

test_that("the public wrapper resolves a recovered required start without a duplicate veto", {
  original_nlminb <- stats::nlminb
  original_polish <- ln_joint_polish
  optimizer_calls <- 0L; polish_calls <- 0L
  local_mocked_bindings(
    nlminb = function(start, objective, ...) {
      optimizer_calls <<- optimizer_calls + 1L
      if (optimizer_calls == 1L) return(list(par = start, objective = objective(start),
        convergence = 0L, message = "controlled failed core endpoint"))
      original_nlminb(start, objective, ...)
    },
    ln_joint_polish = function(theta_reduced, gene_index, posv, ctx, map, ...) {
      polish_calls <<- polish_calls + 1L
      if (polish_calls == 1L) return(list(theta_reduced = theta_reduced,
        objective = ln_joint_objective(theta_reduced, gene_index, posv, ctx, map)))
      original_polish(theta_reduced, gene_index, posv, ctx, map, ...)
    }, .package = "scRoPE"
  )
  fit <- ln_joint_fit_fixture()$fit
  raw <- fit$diagnostics$raw[[1L]]; records <- raw$ln_joint$candidate_records
  recovered <- Filter(function(x) isTRUE(x$recovery$adopted), records)
  expect_length(recovered, 1L)
  expect_identical(recovered[[1L]]$label, "S_direct")
  expect_true(recovered[[1L]]$required)
  expect_true(recovered[[1L]]$core)
  expect_false(recovered[[1L]]$recovery$original$certified)
  expect_identical(recovered[[1L]]$objective_instance, recovered[[1L]]$recovery$certificate$objective_instance)
  expect_equal(sum(vapply(records, function(x) identical(x$label, "S_direct"), logical(1))), 1)
  expect_identical(raw$route, "R0")
  expect_true(raw$inference_available)
  expect_true(raw$ln_joint$covariance$success)
  expect_true(all(vapply(Filter(function(x) !isTRUE(x$core), records), function(x) is.null(x$recovery), logical(1))))
})

test_that("the engine bypasses optional, PQL and constrained candidates", {
  a <- ln_recovery_test_case()
  local_mocked_bindings(
    nlminb = function(start, objective, ...) list(par = start, objective = objective(start), convergence = 0L),
    ln_joint_polish = function(theta_reduced, gene_index, posv, ctx, map, ...)
      list(theta_reduced = theta_reduced, objective = ln_joint_objective(theta_reduced, gene_index, posv, ctx, map)),
    ln_joint_pql_start = function(best, gene_index, posv, ctx, map)
      list(triggered = TRUE, status = "generated", start = c(rep(0, map$free_dim), 1, 1)),
    .package = "scRoPE"
  )
  fit <- ln_joint_fit_engine(1L, a$posv, a$ctx, a$map,
    extra_starts = list(a$candidate$theta_reduced), objective_instance = a$objective_instance)
  excluded <- Filter(function(x) !isTRUE(x$core), fit$candidates)
  expect_true(all(c("primary_optimizer", "pql_completed_start", "declared_extra_1") %in% vapply(excluded, `[[`, "", "label")))
  expect_true(all(vapply(excluded, function(x) is.null(x$recovery), logical(1))))
  pql <- Filter(function(x) identical(x$label, "pql_completed_start"), excluded)[[1L]]
  expect_true(pql$required)
  expect_false(pql$core)
  expect_false(pql$certified)
  constrained <- ln_joint_map(a$ctx$nb, matrix(c(1, 0, 0), nrow = 1L), 0)
  fit <- ln_joint_fit_engine(1L, a$posv, a$ctx, constrained, objective_instance = a$objective_instance)
  expect_true(all(vapply(fit$candidates, function(x) is.null(x$recovery), logical(1))))
})


test_that("certified endpoints are preserved and incomplete certificates cannot be adopted", {
  a <- ln_recovery_test_case()
  successful <- do.call(ln_joint_recover_candidate, a)
  b <- a; b$candidate <- ln_recovery_effective(successful)
  with_mocked_bindings(
    expect_identical(do.call(ln_joint_recover_candidate, b), b$candidate),
    lbfgs = function(...) stop("already certified"), .package = "scRoPE")
  result <- with_mocked_bindings(do.call(ln_joint_recover_candidate, a),
    ln_joint_candidate_state = function(...) {
      x <- successful; x$checks <- c(root = TRUE); x
    }, .package = "scRoPE")
  expect_identical(result$recovery$status, "certificate_failed")
  expect_identical(ln_recovery_effective(result), a$candidate)
  for (status in c(-1L, 0L, 1L, 7L)) {
    a$candidate$optimizer$convergence <- status
    result <- with_mocked_bindings(do.call(ln_joint_recover_candidate, a),
      lbfgs = function(...) stop("eligible regardless of status"), .package = "scRoPE")
    expect_true(result$recovery$attempted)
    expect_identical(result$recovery$status, "optimizer_error")
  }
})
