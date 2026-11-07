test_that("LN worker results are bounded by public output rather than donor history", {
  native <- function(donors) list(
    vec = c(0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.2, 0.3, 1, 1, -10),
    route = "R0", inference_available = TRUE, failure_reason = NA_character_,
    start_provenance = "direct", mode_method = "none", mode_status = "not_requested",
    mode = list(method = "none", status = "not_requested", success = FALSE,
      max_score_residual = NA_real_, max_newton_displacement = NA_real_,
      values = rep(NA_real_, donors), donor_status = rep("not_requested", donors)),
    godambe = list(per_subject_gradients = matrix(1, 4L, donors)),
    ln_joint = list(
      state = list(U_subject_full = matrix(1, 4L, donors)),
      candidate_records = list(list(history = matrix(1, donors, 8L))),
      resolution = list(strict_route = "R0", prospective_route = "R0",
        strict_inference_route = "R0", enabled_inference_route = "R0",
        selected_root_certified = TRUE, would_nonveto = FALSE,
        strict_failure_reason = NA_character_, prospective_failure_reason = NA_character_,
        global_veto_classes = character(0), policy_version = "native", mode = "enabled",
        candidate_resolutions = list(matrix(1, donors, 8L)))),
    lrt = list(x = list(stat = 1, pval = 0.3, success = TRUE,
      constrained = list(state = matrix(1, donors, 8L)),
      score_test = list(per_subject_gradients = matrix(1, donors, 8L)),
      raw = list(wP_minus_qP = 0, wP_over_qP = 1, adjusted_available = TRUE,
        adjusted_failure_reason = NA_character_, donor_scores = matrix(1, donors, 8L),
        gate = list(converged = TRUE, unrestricted_route = "R0",
          diagnostic_history = matrix(1, donors, 8L)))))
  )
  small <- native(20L)
  large <- native(20000L)
  compact <- ln_finalize_gene(large, FALSE, keep_lrt_details = TRUE)
  expect_identical(compact, ln_finalize_gene(small, FALSE, keep_lrt_details = TRUE))
  expect_identical(ln_finalize_gene(large, TRUE), large)
  expect_identical(ln_joint_status_record(compact), ln_joint_status_record(large))
  expect_identical(ln_availability_build(list(compact), 1L, "gene"),
    ln_availability_build(list(large), 1L, "gene"))
  expect_lt(as.numeric(object.size(compact)), as.numeric(object.size(large)) / 100)
  expect_null(compact$ln_joint$state)
  expect_null(compact$ln_joint$candidate_records)
  expect_null(compact$mode$values)
  expect_null(compact$mode$donor_status)
  expect_null(compact$lrt$x$constrained)
  expect_null(compact$lrt$x$score_test)
  expect_null(compact$lrt$x$raw$donor_scores)
  expect_identical(compact$lrt$x$raw$gate, list(converged = TRUE, unrestricted_route = "R0"))

  # The same native certification/failure inputs must survive for withheld rows.
  large$route <- "R2"
  large$inference_available <- FALSE
  large$failure_reason <- "competitive_boundary"
  large$ln_joint$resolution$global_veto_classes <- "competitive_boundary"
  compact <- ln_finalize_gene(large, FALSE)
  expect_identical(ln_joint_status_record(compact), ln_joint_status_record(large))
  expect_identical(ln_availability_build(list(compact), 1L, "gene"),
    ln_availability_build(list(large), 1L, "gene"))
})

test_that("compact recovery envelopes reproduce full route application without fitting", {
  native <- list(inference_available = FALSE, route = "R2",
    failure_reason = "native_failure", ln_joint = list())
  availability <- ln_availability_build(list(native), 1L, "gene")
  per_gene <- availability
  per_gene$inference_available <- FALSE
  base <- list(
    summary = data.frame(gene_id = 1L, gene = "gene", logFC_a = NA_real_,
      logFC_b = NA_real_, wald_success_a = FALSE, wald_success_b = FALSE),
    test_status = data.frame(gene_id = c(1L, 1L), gene = "gene", test = "wald",
      contrast = c("a", "b"), success = FALSE, failure_reason = "native_failure",
      status = "unavailable", inference = "asymptotic", method = "robust_wald",
      failure_stage = NA_character_),
    overdispersion = data.frame(Subject = NA_real_, Cell = NA_real_),
    convergence = -70, algorithm = 1, availability = availability,
    provenance = list(per_gene = per_gene, mode_diagnostics = list(NULL),
      coefficient_map_internal_to_public = diag(2L)))
  heavy <- matrix(1, 4L, 1000L)
  fixed <- list(route = "fixed_v", result = list(success = TRUE, reason = "certified",
    beta = c(0.1, 0.2), Vbeta = diag(c(0.04, 0.09)), fixed = list(v = 1e-4),
    state = list(full = c(0.1, 0.2, log(1e-4), 0), q = 10, U = heavy),
    trials = list(heavy), influence = heavy), refinement = list(trials = heavy))
  face <- list(route = "identifiable_nuisance_limit", result = list(
    success = TRUE, reason = "certified", beta = c(NA_real_, 0.2),
    Vbeta = matrix(c(NA_real_, NA_real_, NA_real_, 0.09), 2L),
    face = list(estimable = c(FALSE, TRUE), retained = rep(TRUE, 1000L)),
    reduced_full = c(0.2, log(0.3), 0), q = 10,
    reduced = list(result = list(fixed = list(), state = heavy)), U_original = heavy))
  unavailable <- list(route = "withheld", result = list(success = FALSE,
    reason = "reduced_root_worse_than_known_full_candidate"),
    face_audit = list(competitor_audit = heavy), refinement = list(trials = heavy))
  for (answer in list(fixed, face, unavailable)) {
    run <- function(value) suppressWarnings(ln_recovery_apply(base, list(native),
      list(pred = matrix(1, 1L, 2L)), TRUE, FALSE, recovery_results = list(value)))
    compact <- ln_compact_recovery(answer)
    expect_identical(run(compact), run(answer))
    expect_lt(as.numeric(object.size(compact)), as.numeric(object.size(answer)) / 5)
    expect_null(compact$refinement)
    expect_null(compact$face_audit)
    expect_null(compact$result$trials)
    expect_null(compact$result$state$U)
  }
})

test_that("ordinary and disabled LN recovery need no raw recovery inputs", {
  native <- list(success = TRUE, joint = list(candidate_records = "retained"))
  expect_identical(ln_recovery_evaluate(native, NULL, NULL, FALSE),
    list(route = "ordinary_jln", result = native))
  expect_identical(ln_recovery_evaluate(list(success = FALSE), NULL, NULL, FALSE),
    list(route = "withheld", result = list(success = FALSE,
      reason = "recovery_requires_default_box_and_lbfgs")))
})

test_that("per-gene recovery isolates nested fitting RNG state", {
  # The real limiting-face route invokes scrope and doRNG. A deterministic
  # stand-in exercises that side effect without repeating an expensive fit.
  local <- new.env(parent = environment(ln_recovery_evaluate))
  local$ln_recovery_route <- function(...) {
    runif(3L)
    list(route = "withheld", result = list(success = FALSE, reason = "fixture"))
  }
  evaluate <- ln_recovery_evaluate
  environment(evaluate) <- local
  input <- list(pred = matrix(1, 2L, 1L), count = matrix(1, 1L, 2L),
    id = 1:2, offset = c(1, 1))
  set.seed(914L)
  seed <- .Random.seed
  kind <- RNGkind()
  answer <- evaluate(list(success = FALSE), input, 1L, TRUE)
  expect_identical(answer$result$reason, "fixture")
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
  local$ln_recovery_route <- function(...) {
    runif(3L)
    stop("fixture_error")
  }
  answer <- evaluate(list(success = FALSE), input, 1L, TRUE)
  expect_identical(answer$result$reason, "recovery_evaluation_error")
  expect_identical(.Random.seed, seed)
  expect_identical(RNGkind(), kind)
})
