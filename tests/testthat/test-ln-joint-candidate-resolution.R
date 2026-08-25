ln_joint_resolution_test_instance <- list(
  scheme = "synthetic_exact_objective_instance_v1",
  id = "deidentified_fixture"
)

ln_joint_resolution_test_candidate <- function(
    label = "S_direct",
    theta_full_log = c(0, 0, log(1), log(2)),
    objective = 10,
    score = rep(0, 4L),
    score_subject = NULL,
    A = diag(4L),
    required = TRUE,
    core = TRUE,
    certified = NULL,
    objective_instance = ln_joint_resolution_test_instance,
    bounds = list(lower = c(0.1, 0.1), upper = c(10, 10))) {
  theta_full_log <- as.numeric(theta_full_log)
  score <- as.numeric(score)
  A <- as.matrix(A)
  if (is.null(score_subject)) {
    score_subject <- cbind(score / 2, score / 2)
  } else {
    score_subject <- as.matrix(score_subject)
  }
  state <- list(
    full_log = theta_full_log,
    theta_reduced_log = theta_full_log,
    objective_provenance = ln_joint_objective_provenance(),
    U_reduced = score,
    U_subject_reduced = score_subject,
    A_reduced = A,
    raw = list(value = objective)
  )
  donor_additivity <- ln_joint_scaled_max(
    rowSums(score_subject),
    score
  )
  root <- ln_joint_root_metric(score_subject, score)
  newton <- ln_joint_probe_newton_metric(state)
  curvature <- ln_joint_curvature_diagnostics(A)
  dispersion <- exp(theta_full_log[3:4])
  boundary_tolerance <- sqrt(.Machine$double.eps) *
    (1 + pmax(abs(bounds$lower), abs(bounds$upper)))
  interior <- all(
    dispersion > bounds$lower + boundary_tolerance &
      dispersion < bounds$upper - boundary_tolerance
  )
  thresholds <- ln_joint_thresholds()
  checks <- c(
    donor_additivity = isTRUE(
      donor_additivity <= thresholds$donor_additivity
    ),
    root = isTRUE(root <= thresholds$root),
    newton = isTRUE(newton <= thresholds$newton),
    curvature = isTRUE(curvature$positive_definite) &&
      curvature$scaled_min_eigenvalue >= thresholds$curvature_min_eigen,
    condition = isTRUE(
      curvature$scaled_condition_number <= thresholds$condition_number
    ),
    interior = isTRUE(interior),
    constraint = TRUE
  )
  if (is.null(certified)) {
    certified <- all(checks)
  }
  list(
    label = label,
    start_provenance = "synthetic_declared_start",
    required = required,
    core = core,
    n_beta = 2L,
    parameterization = "beta, log(sigma2), log(gamma_nb_size)",
    objective_provenance = ln_joint_objective_provenance(),
    objective_instance = objective_instance,
    bounds = bounds,
    probe_bounds = list(
      lower = c(-100, -100, log(bounds$lower)),
      upper = c(100, 100, log(bounds$upper))
    ),
    probe_map = list(
      schema_version = 1L,
      beta_star = c(0, 0),
      null_basis = diag(2L),
      n_beta = 2L,
      free_dim = 2L
    ),
    finite_endpoint = TRUE,
    certified = isTRUE(certified),
    start = c(0, 0, 1, 2),
    optimizer_endpoint = c(0, 0, 1, 2),
    optimizer_objective = objective,
    polished = FALSE,
    theta_reduced = c(theta_full_log[1:2], exp(theta_full_log[3:4])),
    theta_reduced_log = theta_full_log,
    theta_full_log = theta_full_log,
    objective = objective,
    optimizer = list(convergence = 0L, message = "synthetic"),
    state = state,
    donor_additivity_scaled = donor_additivity,
    root_residual = root,
    newton_displacement = newton,
    curvature = curvature,
    interior = interior,
    constraint_residual = 0,
    checks = checks,
    failure = if (isTRUE(all(checks))) NA_character_ else
      "synthetic_precedence_label_is_not_classification_input"
  )
}

ln_joint_resolution_test_failed <- function(label = "optional_failed",
                                              required = FALSE) {
  bounds <- list(lower = c(0.1, 0.1), upper = c(10, 10))
  list(
    label = label,
    start_provenance = "synthetic_declared_start",
    required = required,
    core = FALSE,
    n_beta = 2L,
    parameterization = "beta, log(sigma2), log(gamma_nb_size)",
    objective_provenance = ln_joint_objective_provenance(),
    objective_instance = ln_joint_resolution_test_instance,
    bounds = bounds,
    probe_bounds = list(
      lower = c(-100, -100, log(bounds$lower)),
      upper = c(100, 100, log(bounds$upper))
    ),
    probe_map = list(
      schema_version = 1L,
      beta_star = c(0, 0),
      null_basis = diag(2L),
      n_beta = 2L,
      free_dim = 2L
    ),
    finite_endpoint = FALSE,
    certified = FALSE,
    objective = NA_real_,
    optimizer = list(convergence = 1L, message = "synthetic failure"),
    failure = "nonfinite_endpoint"
  )
}

ln_joint_resolution_test_core <- function() {
  list(
    ln_joint_resolution_test_candidate("S_direct", objective = 10),
    ln_joint_resolution_test_candidate("S_mid", objective = 10 + 1e-10)
  )
}

ln_joint_resolution_test_review <- function(label = "review") {
  target <- c(0, 0, log(1), log(2))
  delta <- c(2e-5, 0, 0, 0)
  ln_joint_resolution_test_candidate(
    label = label,
    theta_full_log = target + delta,
    objective = 10 + 0.5 * sum(delta^2),
    score = -delta,
    required = FALSE,
    core = FALSE
  )
}

ln_joint_resolution_test_evaluator <- function(selected) {
  target <- selected$theta_full_log
  objective_minimum <- selected$objective
  function(theta_reduced_log) {
    theta_reduced_log <- as.numeric(theta_reduced_log)
    delta <- theta_reduced_log - target
    score <- -delta
    list(
      full_log = theta_reduced_log,
      theta_reduced_log = theta_reduced_log,
      objective_provenance = ln_joint_objective_provenance(),
      U_reduced = score,
      U_subject_reduced = cbind(score / 2, score / 2),
      A_reduced = diag(length(theta_reduced_log)),
      raw = list(value = objective_minimum + 0.5 * sum(delta^2))
    )
  }
}

ln_joint_resolution_test_probe <- function(candidate, theta_reduced_log) {
  step <- as.numeric(solve(
    candidate$state$A_reduced,
    candidate$state$U_reduced
  ))
  score <- rep(0, length(theta_reduced_log))
  list(
    full_log = ln_joint_probe_full_log(theta_reduced_log, candidate),
    theta_reduced_log = as.numeric(theta_reduced_log),
    objective_provenance = ln_joint_objective_provenance(),
    U_reduced = score,
    U_subject_reduced = cbind(score / 2, score / 2),
    A_reduced = candidate$state$A_reduced,
    raw = list(value = as.numeric(
      candidate$objective -
        0.5 * crossprod(step, candidate$state$A_reduced %*% step)
    ))
  )
}

test_that("calibration thresholds are versioned and leave the root certificate unchanged", {
  selected <- ln_joint_thresholds()
  resolution <- ln_joint_resolution_thresholds()

  expect_true(ln_joint_resolution_thresholds_valid(resolution))
  expect_identical(resolution$schema_version, 1L)
  expect_identical(
    resolution$policy_version,
    "joint_ln_candidate_resolution_calibration_v1"
  )
  expect_identical(
    resolution$objective_equivalence,
    selected$hard_stop$objective_agreement
  )
  expect_identical(
    resolution$review_full_state,
    selected$hard_stop$parameter_agreement
  )
  expect_identical(resolution$proposal_root, selected$root)
  expect_identical(resolution$proposal_newton, selected$newton)
  expect_identical(resolution$max_state_evaluations, 1L)
  expect_identical(
    ln_joint_thresholds(),
    selected
  )

  reordered <- resolution[rev(names(resolution))]
  expect_false(ln_joint_resolution_thresholds_valid(reordered))
  loosened_fields <- c(
    "objective_equivalence", "review_full_state", "review_beta_block",
    "review_dispersion_block", "boundary_scaled_distance",
    "curvature_compatibility", "proposal_full_state",
    "proposal_beta_block", "proposal_dispersion_block",
    "proposal_objective_equivalence", "proposal_root", "proposal_newton",
    "contraction_ratio"
  )
  for (field in loosened_fields) {
    changed <- resolution
    changed[[field]] <- resolution[[field]] * (1 + 1e-6)
    expect_false(
      ln_joint_resolution_thresholds_valid(changed),
      info = paste("loosened field", field)
    )
  }
  changed <- resolution
  changed$max_state_evaluations <- 2L
  expect_false(ln_joint_resolution_thresholds_valid(changed))
  changed <- resolution
  changed$policy_version <-
    "joint_ln_candidate_resolution_calibration_v1_loose"
  expect_false(ln_joint_resolution_thresholds_valid(changed))
})

test_that("complete vectors classify every prospective candidate class", {
  thresholds <- ln_joint_resolution_thresholds()
  selected <- ln_joint_resolution_test_core()[[1L]]
  certified <- ln_joint_resolution_test_core()[[2L]]
  dominated <- ln_joint_resolution_test_candidate(
    "dominated",
    objective = 10.1,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  review <- ln_joint_resolution_test_review()
  unresolved <- review
  boundary_theta <- selected$theta_full_log
  boundary_theta[[3L]] <- log(0.1 + 3e-7)
  boundary <- ln_joint_resolution_test_candidate(
    "boundary",
    theta_full_log = boundary_theta,
    objective = 10,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  distinct <- ln_joint_resolution_test_candidate(
    "distinct",
    theta_full_log = review$theta_full_log + c(0, 0.01, 0, 0),
    objective = review$objective,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  hard <- ln_joint_resolution_test_candidate(
    "hard",
    theta_full_log = review$theta_full_log,
    objective = review$objective,
    score = c(-2e-5, 0, 0, 0),
    score_subject = matrix(0, nrow = 4L, ncol = 2L),
    required = FALSE,
    core = FALSE
  )

  expect_identical(
    ln_joint_classify_candidate(
      selected, selected, thresholds, is_selected = TRUE
    )$class,
    "certified_selected"
  )
  expect_identical(
    ln_joint_classify_candidate(certified, selected, thresholds)$class,
    "certified_same_cluster"
  )
  expect_identical(
    ln_joint_classify_candidate(dominated, selected, thresholds)$class,
    "dominated_uncertified"
  )
  expect_identical(
    ln_joint_classify_candidate(review, selected, thresholds)$class,
    "same_basin_review_candidate"
  )
  expect_identical(
    ln_joint_classify_candidate(boundary, selected, thresholds)$class,
    "competitive_boundary"
  )
  expect_identical(
    ln_joint_classify_candidate(distinct, selected, thresholds)$class,
    "competitive_distinct_basin"
  )
  expect_identical(
    ln_joint_classify_candidate(hard, selected, thresholds)$class,
    "competitive_hard_stop"
  )
  expect_identical(
    ln_joint_classify_candidate(hard, selected, thresholds)$reason,
    "failed_checks_not_limited_to_root_and_newton"
  )
  expect_identical(
    ln_joint_classify_candidate(
      ln_joint_resolution_test_failed(), selected, thresholds
    )$class,
    "noncompetitive_failed_attempt"
  )
  expect_identical(
    ln_joint_classify_candidate(
      ln_joint_resolution_test_failed(required = TRUE), selected, thresholds
    )$class,
    "competitive_hard_stop"
  )
})

test_that("bounded review is endpoint-local, immutable, and exactly one evaluation", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  review <- ln_joint_resolution_test_review()
  before <- serialize(review, NULL)
  calls <- 0L
  evaluator <- ln_joint_resolution_test_evaluator(selected)
  counted <- function(theta) {
    calls <<- calls + 1L
    evaluator(theta)
  }

  local <- ln_joint_basin_local_proposal(
    review,
    ln_joint_resolution_thresholds()
  )
  expect_true(local$success)
  expect_equal(local$proposed_theta_reduced_log, selected$theta_full_log,
               tolerance = 1e-15)

  probe <- ln_joint_basin_probe(
    review,
    selected,
    counted,
    ln_joint_resolution_thresholds()
  )
  expect_true(probe$accepted)
  expect_identical(
    probe$termination_reason,
    "same_basin_contraction_established"
  )
  expect_identical(calls, 1L)
  expect_identical(probe$state_evaluations, 1L)
  expect_identical(probe$objective_evaluations, 1L)
  expect_identical(probe$score_evaluations, 1L)
  expect_identical(probe$hessian_evaluations, 1L)
  expect_true(all(unlist(probe$contraction, use.names = FALSE)))
  expect_identical(serialize(review, NULL), before)

  other_selected <- selected
  other_selected$theta_full_log[[2L]] <- 0.05
  other_selected$state$full_log <- other_selected$theta_full_log
  other_selected$state$theta_reduced_log <- other_selected$theta_full_log
  other_probe <- ln_joint_basin_probe(
    review,
    other_selected,
    ln_joint_resolution_test_evaluator(selected),
    ln_joint_resolution_thresholds()
  )
  expect_identical(
    other_probe$proposed$theta_reduced_log,
    probe$proposed$theta_reduced_log
  )
  expect_false(other_probe$accepted)
})

test_that("probe failures do not retry, clamp, or enter candidate selection", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  review <- ln_joint_resolution_test_review()
  outside <- review
  outside$probe_bounds$lower <- outside$state$theta_reduced_log - 1e-8
  outside$probe_bounds$upper <- outside$state$theta_reduced_log + 1e-8
  outside_probe <- ln_joint_basin_probe(
    outside,
    selected,
    function(theta) stop("must not evaluate"),
    ln_joint_resolution_thresholds()
  )
  expect_identical(
    outside_probe$termination_reason,
    "local_proposal_outside_bounds"
  )
  expect_identical(outside_probe$state_evaluations, 0L)

  calls <- 0L
  failed_probe <- ln_joint_basin_probe(
    review,
    selected,
    function(theta) {
      calls <<- calls + 1L
      NULL
    },
    ln_joint_resolution_thresholds()
  )
  expect_false(failed_probe$accepted)
  expect_identical(calls, 1L)
  expect_identical(failed_probe$state_evaluations, 1L)

  candidates <- c(ln_joint_resolution_test_core(), list(review))
  before <- serialize(candidates, NULL)
  resolved <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = ln_joint_resolution_test_probe
  )
  expect_identical(resolved$best_index, 1L)
  expect_identical(serialize(candidates, NULL), before)
  expect_length(resolved$resolution$candidate_resolutions, length(candidates))
})

test_that("default strict route is unchanged while the test seam records nonveto", {
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(ln_joint_resolution_test_review())
  )
  strict <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = ln_joint_resolution_test_probe
  )
  enabled <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = ln_joint_resolution_test_probe,
    enable_prospective = TRUE
  )

  expect_identical(strict$route, "R2")
  expect_identical(strict$failure_reason, "competing_uncertified_candidate")
  expect_identical(strict$resolution$strict_route, "R2")
  expect_identical(
    strict$resolution$strict_failure_reason,
    "competing_uncertified_candidate"
  )
  expect_identical(strict$resolution$prospective_route, "R0")
  expect_true(strict$resolution$would_nonveto)
  expect_identical(strict$resolution$would_route, "R0")
  expect_identical(
    strict$resolution$numerical_resolution,
    "would_same_basin_review_nonveto"
  )
  expect_identical(strict$resolution$mode, "shadow_default_off")
  expect_identical(enabled$route, "R0")
  expect_true(is.na(enabled$failure_reason))
  expect_identical(
    enabled$resolution$mode,
    "calibration_candidate_enabled"
  )
  expect_identical(enabled$best_index, strict$best_index)
  expect_identical(
    enabled$resolution$evaluation_counts$probe_state_evaluations,
    1L
  )
  expect_false(
    enabled$resolution$evaluation_counts$total_engine_counts_available
  )
})

test_that("one relevant veto defeats any resolved same-basin candidate", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  boundary_theta <- selected$theta_full_log
  boundary_theta[[3L]] <- log(0.1 + 3e-7)
  boundary <- ln_joint_resolution_test_candidate(
    "boundary",
    theta_full_log = boundary_theta,
    objective = 10,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(ln_joint_resolution_test_review(), boundary)
  )
  result <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = ln_joint_resolution_test_probe
  )

  expect_false(result$resolution$would_nonveto)
  expect_identical(result$resolution$prospective_route, "R2")
  expect_true(
    "competitive_boundary" %in% result$resolution$global_veto_classes
  )
  expect_identical(
    result$resolution$prospective_failure_reason,
    "competitive_boundary_candidate"
  )

  optional_failed <- c(
    ln_joint_resolution_test_core(),
    list(ln_joint_resolution_test_failed())
  )
  optional_result <- ln_joint_resolve_candidates(optional_failed)
  expect_identical(optional_result$route, "R0")
  expect_false(any(vapply(
    optional_result$resolution$candidate_resolutions,
    `[[`, logical(1), "veto"
  )))

  required_failed <- c(
    ln_joint_resolution_test_core(),
    list(ln_joint_resolution_test_failed(required = TRUE))
  )
  required_result <- ln_joint_resolve_candidates(required_failed)
  expect_identical(required_result$route, "R2")
  expect_true(
    "competitive_hard_stop" %in%
      required_result$resolution$global_veto_classes
  )

  required_dominated <- ln_joint_resolution_test_candidate(
    "required_dominated",
    objective = 10.1,
    score = c(-2e-5, 0, 0, 0),
    required = TRUE,
    core = FALSE
  )
  required_dominated_candidates <- c(
    ln_joint_resolution_test_core(),
    list(required_dominated)
  )
  required_default <- ln_joint_resolve_candidates(
    required_dominated_candidates
  )
  required_enabled <- ln_joint_resolve_candidates(
    required_dominated_candidates,
    enable_prospective = TRUE
  )
  expect_identical(required_default$route, "R0")
  expect_identical(required_default$resolution$prospective_route, "R2")
  expect_identical(
    required_default$resolution$prospective_failure_reason,
    "competitive_hard_stop_candidate"
  )
  expect_identical(required_enabled$route, "R2")
  expect_identical(
    required_enabled$failure_reason,
    "competitive_hard_stop_candidate"
  )
})

test_that("invalid resolution configuration fails closed only prospectively", {
  candidates <- ln_joint_resolution_test_core()
  malformed <- 1
  result <- ln_joint_resolve_candidates(
    candidates,
    resolution_thresholds = malformed
  )
  enabled <- ln_joint_resolve_candidates(
    candidates,
    resolution_thresholds = malformed,
    enable_prospective = TRUE
  )

  expect_identical(result$route, "R0")
  expect_identical(result$resolution$strict_route, "R0")
  expect_identical(result$resolution$prospective_route, "R2")
  expect_identical(
    result$resolution$prospective_failure_reason,
    "competitive_hard_stop_candidate"
  )
  expect_identical(enabled$route, "R2")
  expect_true(is.na(result$resolution$policy_version))
})

test_that("classification uses every check and raw metric, never the scalar label", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  review <- ln_joint_resolution_test_review()
  thresholds <- ln_joint_resolution_thresholds()

  label_changed <- review
  label_changed$failure <- "dispersion_boundary"
  label_changed$optimizer$message <- "desired-significance-target"
  expect_identical(
    ln_joint_classify_candidate(label_changed, selected, thresholds)$class,
    "same_basin_review_candidate"
  )

  simultaneous <- ln_joint_classify_candidate(review, selected, thresholds)
  expect_identical(
    simultaneous$class,
    "same_basin_review_candidate"
  )
  expect_setequal(simultaneous$failed_checks, c("root", "newton"))

  reordered <- review
  reordered$checks <- reordered$checks[rev(names(reordered$checks))]
  expect_identical(
    ln_joint_classify_candidate(reordered, selected, thresholds)$reason,
    "incomplete_check_vector"
  )

  lying_check <- review
  lying_check$checks[["root"]] <- TRUE
  expect_identical(
    ln_joint_classify_candidate(lying_check, selected, thresholds)$reason,
    "check_vector_metric_mismatch"
  )

  missing_metric <- review
  missing_metric$root_residual <- NA_real_
  expect_identical(
    ln_joint_classify_candidate(missing_metric, selected, thresholds)$reason,
    "missing_or_nonfinite_check_metric"
  )

  changed_state <- review
  changed_state$state$U_reduced[[1L]] <- -5e-5
  expect_identical(
    ln_joint_classify_candidate(changed_state, selected, thresholds)$reason,
    "check_metric_record_mismatch"
  )

  changed_objective <- review
  changed_objective$state$raw$value <- 100
  expect_identical(
    ln_joint_classify_candidate(changed_objective, selected, thresholds)$reason,
    "objective_state_mismatch"
  )

  wrong_provenance <- review
  wrong_provenance$objective_instance <- list(id = "other")
  expect_identical(
    ln_joint_classify_candidate(wrong_provenance, selected, thresholds)$reason,
    "fit_provenance_mismatch"
  )
})

test_that("objective, state, boundary, and curvature thresholds are inclusive", {
  thresholds <- ln_joint_resolution_thresholds()
  objective_tolerance <- thresholds$objective_equivalence
  selected_objective <- 10
  positive_objective <- function(gap) {
    selected_objective + gap * (1 + selected_objective) / (1 - gap)
  }
  negative_objective <- function(gap) {
    selected_objective - gap * (1 + selected_objective)
  }
  epsilon <- sqrt(.Machine$double.eps)

  expect_identical(
    ln_joint_objective_comparison(
      positive_objective(objective_tolerance),
      selected_objective,
      objective_tolerance
    )$relation,
    "equivalent"
  )
  expect_identical(
    ln_joint_objective_comparison(
      positive_objective(objective_tolerance * (1 + epsilon)),
      selected_objective,
      objective_tolerance
    )$relation,
    "dominated"
  )
  negative_at <- negative_objective(objective_tolerance)
  representable_negative_gap <- abs(
    (negative_at - selected_objective) /
      (1 + max(abs(negative_at), abs(selected_objective)))
  )
  expect_identical(
    ln_joint_objective_comparison(
      negative_at,
      selected_objective,
      representable_negative_gap
    )$relation,
    "equivalent"
  )
  expect_identical(
    ln_joint_objective_comparison(
      negative_objective(objective_tolerance * (1 + epsilon)),
      selected_objective,
      objective_tolerance
    )$relation,
    "materially_better"
  )

  state_tolerance <- thresholds$review_full_state
  distance_coordinate <- function(distance) distance / (1 - distance)
  at <- ln_joint_resolution_block_metrics(
    c(distance_coordinate(state_tolerance), 0, 0, 0),
    rep(0, 4),
    2L
  )
  above <- ln_joint_resolution_block_metrics(
    c(distance_coordinate(state_tolerance * (1 + epsilon)), 0, 0, 0),
    rep(0, 4),
    2L
  )
  expect_equal(at$full_state, state_tolerance, tolerance = 1e-15)
  expect_gt(above$full_state, state_tolerance)

  boundary_tolerance <- thresholds$boundary_scaled_distance
  lower <- 0.1
  at_boundary <- (lower + boundary_tolerance) /
    (1 - boundary_tolerance)
  above_boundary <- (lower + boundary_tolerance * (1 + epsilon)) /
    (1 - boundary_tolerance * (1 + epsilon))
  at_boundary_metric <- ln_joint_scaled_bound_distance(
    c(0, 0, log(at_boundary), log(2)),
    2L,
    list(lower = c(lower, lower), upper = c(10, 10))
  )
  expect_lte(abs(at_boundary_metric - boundary_tolerance), 2e-17)
  expect_gt(
    ln_joint_scaled_bound_distance(
      c(0, 0, log(above_boundary), log(2)),
      2L,
      list(lower = c(lower, lower), upper = c(10, 10))
    ),
    boundary_tolerance
  )

  selected <- ln_joint_resolution_test_candidate(
    "boundary_selected",
    theta_full_log = c(0, 0, log(at_boundary + 2e-5), log(2)),
    objective = 10
  )
  boundary_candidate <- ln_joint_resolution_test_candidate(
    "boundary_at",
    theta_full_log = c(0, 0, log(at_boundary), log(2)),
    objective = 10,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  at_boundary_threshold <- thresholds
  at_boundary_threshold$boundary_scaled_distance <- at_boundary_metric
  below_boundary_threshold <- at_boundary_threshold
  below_boundary_threshold$boundary_scaled_distance <-
    at_boundary_metric * (1 - epsilon)
  expect_identical(
    ln_joint_classify_candidate(
      boundary_candidate,
      selected,
      at_boundary_threshold
    )$class,
    "competitive_boundary"
  )
  expect_identical(
    ln_joint_classify_candidate(
      boundary_candidate,
      selected,
      below_boundary_threshold
    )$class,
    "same_basin_review_candidate"
  )

  selected <- ln_joint_resolution_test_core()[[1L]]
  curved_A <- diag(4L)
  curved_A[1L, 2L] <- curved_A[2L, 1L] <- 0.2
  curved <- ln_joint_resolution_test_candidate(
    "curved",
    theta_full_log = selected$theta_full_log + c(2e-5, 0, 0, 0),
    objective = 10,
    score = c(-2e-5, 0, 0, 0),
    A = curved_A,
    required = FALSE,
    core = FALSE
  )
  curvature_metric <- ln_joint_curvature_compatibility(
    curved$state$A_reduced,
    selected$state$A_reduced
  )
  at_curvature <- thresholds
  at_curvature$curvature_compatibility <- curvature_metric
  below_curvature <- thresholds
  below_curvature$curvature_compatibility <-
    curvature_metric * (1 - epsilon)
  expect_identical(
    ln_joint_classify_candidate(curved, selected, at_curvature)$class,
    "same_basin_review_candidate"
  )
  expect_identical(
    ln_joint_classify_candidate(curved, selected, below_curvature)$class,
    "competitive_distinct_basin"
  )
})

test_that("full beta and dispersion blocks independently veto distinct basins", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  thresholds <- ln_joint_resolution_thresholds()
  review <- ln_joint_resolution_test_review()
  beta_theta <- review$theta_full_log
  beta_theta[[2L]] <- 0.01
  beta_far <- ln_joint_resolution_test_candidate(
    "beta_far",
    theta_full_log = beta_theta,
    objective = review$objective,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  dispersion_theta <- review$theta_full_log
  dispersion_theta[[4L]] <- dispersion_theta[[4L]] + 0.01
  dispersion_far <- ln_joint_resolution_test_candidate(
    "dispersion_far",
    theta_full_log = dispersion_theta,
    objective = review$objective,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )

  beta_resolution <- ln_joint_classify_candidate(
    beta_far, selected, thresholds
  )
  dispersion_resolution <- ln_joint_classify_candidate(
    dispersion_far, selected, thresholds
  )
  expect_identical(beta_resolution$class, "competitive_distinct_basin")
  expect_gt(
    beta_resolution$metrics$beta_block,
    thresholds$review_beta_block
  )
  expect_identical(
    dispersion_resolution$class,
    "competitive_distinct_basin"
  )
  expect_gt(
    dispersion_resolution$metrics$dispersion_block,
    thresholds$review_dispersion_block
  )
})

test_that("ordering duplicates and optimizer text cannot change global decision", {
  review <- ln_joint_resolution_test_review()
  candidates <- c(ln_joint_resolution_test_core(), list(review))
  permuted <- candidates[c(3L, 2L, 1L)]
  duplicated_candidates <- c(candidates, list(review))

  base <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = ln_joint_resolution_test_probe
  )
  reordered <- ln_joint_resolve_candidates(
    permuted,
    basin_probe = ln_joint_resolution_test_probe
  )
  duplicated <- ln_joint_resolve_candidates(
    duplicated_candidates,
    basin_probe = ln_joint_resolution_test_probe
  )

  expect_identical(base$resolution$prospective_route, "R0")
  expect_identical(
    reordered$resolution$prospective_route,
    base$resolution$prospective_route
  )
  expect_identical(
    duplicated$resolution$prospective_route,
    base$resolution$prospective_route
  )
  expect_identical(
    reordered$resolution$global_veto_classes,
    base$resolution$global_veto_classes
  )
  expect_identical(
    duplicated$resolution$global_veto_classes,
    base$resolution$global_veto_classes
  )

  renamed <- candidates
  renamed[[3L]]$label <- "different_label"
  renamed[[3L]]$optimizer$message <- "different optimizer message"
  renamed_result <- ln_joint_resolve_candidates(
    renamed,
    basin_probe = ln_joint_resolution_test_probe
  )
  expect_identical(
    renamed_result$resolution$prospective_route,
    base$resolution$prospective_route
  )
})

test_that("resolution and additive candidate records serialize losslessly", {
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(ln_joint_resolution_test_review())
  )
  result <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = ln_joint_resolution_test_probe
  )
  restored <- unserialize(serialize(result$resolution, NULL))
  expect_identical(restored, result$resolution)
  expect_identical(
    names(result$resolution),
    c(
      "schema_version", "policy_version", "mode", "thresholds",
      "strict_route", "strict_failure_reason", "prospective_route",
      "prospective_failure_reason", "selected_root_certified",
      "would_nonveto", "would_route", "numerical_resolution",
      "global_veto_classes",
      "candidate_resolutions", "evaluation_counts"
    )
  )
  expect_identical(result$resolution$schema_version, 2L)
  expect_true(all(c(
    "resolution_class", "resolution_veto", "objective_relation",
    "full_state_distance", "beta_block_distance",
    "dispersion_block_distance", "probe_accepted",
    "probe_state_evaluations"
  ) %in% names(result$candidates)))

  records <- ln_joint_candidate_records(
    candidates,
    result$resolution$candidate_resolutions
  )
  expect_length(records, length(candidates))
  expect_identical(
    records[[3L]]$resolution$class,
    "same_basin_review_nonveto"
  )
  expect_identical(records[[3L]]$objective_instance,
                   ln_joint_resolution_test_instance)
  expect_identical(unserialize(serialize(records, NULL)), records)
})

test_that("an evaluator cannot snap the frozen proposal to the selected root", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  target <- selected$theta_full_log
  candidate <- ln_joint_resolution_test_candidate(
    "adversarial_snap",
    theta_full_log = target + c(2e-5, 0, 0, 0),
    objective = 10 + 5e-11,
    score = c(-1e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  local <- ln_joint_basin_local_proposal(
    candidate,
    ln_joint_resolution_thresholds()
  )
  expect_true(local$success)
  expect_false(identical(local$proposed_theta_reduced_log, target))

  evaluations <- 0L
  selected_evaluator <- ln_joint_resolution_test_evaluator(selected)
  adversarial_evaluator <- function(theta_reduced_log) {
    evaluations <<- evaluations + 1L
    selected_evaluator(target)
  }
  probe <- ln_joint_basin_probe(
    candidate,
    selected,
    adversarial_evaluator,
    ln_joint_resolution_thresholds()
  )

  expect_false(probe$accepted)
  expect_identical(
    probe$termination_reason,
    "proposed_state_coordinate_mismatch"
  )
  expect_identical(evaluations, 1L)
  expect_identical(probe$state_evaluations, 1L)
  expect_identical(probe$objective_evaluations, 1L)
  expect_identical(probe$score_evaluations, 1L)
  expect_identical(probe$hessian_evaluations, 1L)
  expect_identical(
    probe$proposed$theta_reduced_log,
    local$proposed_theta_reduced_log
  )
  expect_identical(probe$proposed$evaluated_theta_reduced_log, target)
})

test_that("exact-objective tied roots are order invariant and all are anchors", {
  base <- c(0, 0, log(1), log(2))
  root_a <- ln_joint_resolution_test_candidate(
    "root_a",
    theta_full_log = base,
    objective = 10
  )
  root_b <- ln_joint_resolution_test_candidate(
    "root_b",
    theta_full_log = base + c(9e-7, 0, 0, 0),
    objective = 10
  )
  d <- 9.95e-5
  competitor <- ln_joint_resolution_test_candidate(
    "tied_root_counterexample",
    theta_full_log = base - c(d, 0, 0, 0),
    objective = 10,
    score = c(2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  thresholds <- ln_joint_resolution_thresholds()

  against_a <- ln_joint_classify_candidate(
    competitor,
    root_a,
    thresholds
  )
  against_b <- ln_joint_classify_candidate(
    competitor,
    root_b,
    thresholds
  )
  expect_identical(against_a$class, "same_basin_review_candidate")
  expect_lt(against_a$metrics$full_state, thresholds$review_full_state)
  expect_identical(against_b$class, "competitive_distinct_basin")
  expect_gt(against_b$metrics$full_state, thresholds$review_full_state)

  probe_calls <- 0L
  forbidden_probe <- function(...) {
    probe_calls <<- probe_calls + 1L
    stop("a tied-root incompatibility must veto before probing")
  }
  order_ab <- ln_joint_resolve_candidates(
    list(root_a, root_b, competitor),
    basin_probe = forbidden_probe
  )
  order_ba <- ln_joint_resolve_candidates(
    list(root_b, root_a, competitor),
    basin_probe = forbidden_probe
  )

  expect_identical(probe_calls, 0L)
  expect_identical(
    order_ab$resolution$candidate_resolutions[[3L]]$class,
    "competitive_distinct_basin"
  )
  expect_identical(
    order_ba$resolution$candidate_resolutions[[3L]]$class,
    "competitive_distinct_basin"
  )
  expect_identical(
    order_ab$resolution$global_veto_classes,
    order_ba$resolution$global_veto_classes
  )
  expect_identical(
    order_ab$resolution$prospective_route,
    order_ba$resolution$prospective_route
  )
  expect_identical(order_ab$resolution$prospective_route, "R2")
  expect_identical(
    order_ab$resolution$evaluation_counts$probe_state_evaluations,
    0L
  )
  expect_identical(
    order_ba$resolution$evaluation_counts$probe_state_evaluations,
    0L
  )
})

test_that("an immediate boundary veto skips every bounded review", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  boundary_theta <- selected$theta_full_log
  boundary_theta[[3L]] <- log(0.1 + 3e-7)
  boundary <- ln_joint_resolution_test_candidate(
    "immediate_boundary",
    theta_full_log = boundary_theta,
    objective = 10,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  review <- ln_joint_resolution_test_review("skipped_review")
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(review, boundary)
  )
  probe_calls <- 0L
  forbidden_probe <- function(...) {
    probe_calls <<- probe_calls + 1L
    stop("an immediate substantive veto must skip bounded review")
  }
  result <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = forbidden_probe
  )
  review_record <- result$resolution$candidate_resolutions[[3L]]
  boundary_record <- result$resolution$candidate_resolutions[[4L]]

  expect_identical(probe_calls, 0L)
  expect_identical(boundary_record$class, "competitive_boundary")
  expect_identical(review_record$class, "same_basin_review_unresolved")
  expect_identical(
    review_record$reason,
    "bounded_review_skipped_due_to_substantive_veto"
  )
  expect_identical(
    review_record$probe$termination_reason,
    "skipped_immediate_substantive_veto"
  )
  expect_identical(review_record$probe$state_evaluations, 0L)
  expect_identical(
    result$resolution$evaluation_counts$probe_state_evaluations,
    0L
  )
  expect_identical(
    unserialize(serialize(result$resolution, NULL)),
    result$resolution
  )
})

test_that("fully equivalent review endpoints share one physical probe", {
  review_a <- ln_joint_resolution_test_review("duplicate_review_a")
  review_b <- review_a
  review_b$label <- "duplicate_review_b"
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(review_a, review_b)
  )
  probe_calls <- 0L
  counted_probe <- function(candidate, theta_reduced_log) {
    probe_calls <<- probe_calls + 1L
    ln_joint_resolution_test_probe(
      candidate,
      theta_reduced_log
    )
  }
  result <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = counted_probe
  )
  first <- result$resolution$candidate_resolutions[[3L]]
  second <- result$resolution$candidate_resolutions[[4L]]

  expect_identical(probe_calls, 1L)
  expect_identical(first$class, "same_basin_review_nonveto")
  expect_identical(second$class, "same_basin_review_nonveto")
  expect_false(first$probe_reused)
  expect_true(second$probe_reused)
  expect_identical(first$probe, second$probe)
  expect_identical(
    result$resolution$evaluation_counts$probe_state_evaluations,
    1L
  )
  expect_length(result$resolution$candidate_resolutions, 4L)
  expect_identical(
    result$candidates$label[3:4],
    c("duplicate_review_a", "duplicate_review_b")
  )
  records <- ln_joint_candidate_records(
    candidates,
    result$resolution$candidate_resolutions
  )
  expect_length(records, 4L)
  expect_identical(
    vapply(records[3:4], `[[`, character(1), "label"),
    c("duplicate_review_a", "duplicate_review_b")
  )
})

test_that("curvature compatibility uses the selected root's common scale", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  scaled_A <- 2 * selected$state$A_reduced
  compatibility <- ln_joint_curvature_compatibility(
    scaled_A,
    selected$state$A_reduced
  )
  expect_gt(
    compatibility,
    ln_joint_resolution_thresholds()$curvature_compatibility
  )

  candidate <- ln_joint_resolution_test_candidate(
    "common_scale_curvature",
    theta_full_log = selected$theta_full_log + c(2e-5, 0, 0, 0),
    objective = 10,
    score = c(-2e-5, 0, 0, 0),
    A = scaled_A,
    required = FALSE,
    core = FALSE
  )
  resolution <- ln_joint_classify_candidate(
    candidate,
    selected,
    ln_joint_resolution_thresholds()
  )
  expect_equal(
    resolution$metrics$curvature_compatibility,
    compatibility,
    tolerance = 0
  )
  expect_identical(resolution$class, "competitive_distinct_basin")
})

test_that("dominated derivative failures stay diagnostic but contradictions veto", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  dominated <- ln_joint_resolution_test_candidate(
    "dominated_complete",
    objective = 10.1,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  derivative_failed <- dominated
  derivative_failed$label <- "dominated_derivative_failure"
  derivative_failed$state <- NULL
  derivative_failed$theta_full_log <- NULL
  derivative_failed$checks <- NULL
  derivative_failed$donor_additivity_scaled <- NULL
  derivative_failed$root_residual <- NULL
  derivative_failed$newton_displacement <- NULL
  derivative_failed$curvature <- NULL
  derivative_failed$interior <- NULL
  derivative_failed$constraint_residual <- NULL
  derivative_failed$failure <- "derivative_evaluation_failed"
  diagnostic <- ln_joint_classify_candidate(
    derivative_failed,
    selected,
    ln_joint_resolution_thresholds()
  )
  expect_identical(diagnostic$class, "dominated_uncertified")
  expect_false(diagnostic$veto)
  expect_identical(
    diagnostic$reason,
    "optional_candidate_objectively_dominated"
  )

  contradiction <- dominated
  contradiction$label <- "dominated_raw_check_contradiction"
  contradiction$checks[["root"]] <- TRUE
  hard_stop <- ln_joint_classify_candidate(
    contradiction,
    selected,
    ln_joint_resolution_thresholds()
  )
  expect_identical(hard_stop$class, "competitive_hard_stop")
  expect_true(hard_stop$veto)
  expect_identical(hard_stop$reason, "check_vector_metric_mismatch")
})

test_that("an objectively dominated boundary endpoint remains a nonveto", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  boundary_theta <- selected$theta_full_log
  boundary_theta[[3L]] <- log(0.1 + 3e-7)
  dominated_boundary <- ln_joint_resolution_test_candidate(
    "dominated_boundary",
    theta_full_log = boundary_theta,
    objective = 10.1,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )
  resolution <- ln_joint_classify_candidate(
    dominated_boundary,
    selected,
    ln_joint_resolution_thresholds()
  )

  expect_lte(
    resolution$metrics$boundary_scaled_distance,
    ln_joint_resolution_thresholds()$boundary_scaled_distance
  )
  expect_identical(resolution$objective_relation, "dominated")
  expect_identical(resolution$class, "dominated_uncertified")
  expect_false(resolution$veto)

  global <- ln_joint_resolve_candidates(c(
    ln_joint_resolution_test_core(),
    list(dominated_boundary)
  ))
  expect_identical(global$route, "R0")
  expect_identical(global$resolution$prospective_route, "R0")
  expect_length(global$resolution$global_veto_classes, 0L)
})

test_that("only validated metadata-only optional failures remain nonveto", {
  valid_failed <- ln_joint_resolution_test_failed("validated_optional_failed")
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(valid_failed, list(), 7)
  )
  result <- NULL
  expect_no_error({
    result <- ln_joint_resolve_candidates(candidates)
  })

  expect_identical(nrow(result$candidates), length(candidates))
  expect_length(
    result$resolution$candidate_resolutions,
    length(candidates)
  )
  expect_identical(
    result$resolution$candidate_resolutions[[3L]]$class,
    "noncompetitive_failed_attempt"
  )
  expect_false(result$resolution$candidate_resolutions[[3L]]$veto)
  expect_identical(
    result$resolution$candidate_resolutions[[4L]]$class,
    "competitive_hard_stop"
  )
  expect_identical(
    result$resolution$candidate_resolutions[[4L]]$reason,
    "invalid_candidate_role_flags"
  )
  expect_identical(
    result$resolution$candidate_resolutions[[5L]]$class,
    "competitive_hard_stop"
  )
  expect_identical(
    result$resolution$candidate_resolutions[[5L]]$reason,
    "candidate_not_a_list"
  )
  expect_true(
    "competitive_hard_stop" %in%
      result$resolution$global_veto_classes
  )
})

test_that("nonscalar finite objectives hard-stop without changing row count", {
  empty_objective <- ln_joint_resolution_test_review("empty_objective")
  empty_objective$objective <- numeric(0)
  multiple_objective <- ln_joint_resolution_test_review("multiple_objective")
  multiple_objective$objective <- c(10, 10.1)
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(empty_objective, multiple_objective)
  )
  result <- NULL
  expect_no_error({
    result <- ln_joint_resolve_candidates(candidates)
  })

  expect_identical(nrow(result$candidates), length(candidates))
  expect_length(
    result$resolution$candidate_resolutions,
    length(candidates)
  )
  malformed <- result$resolution$candidate_resolutions[3:4]
  expect_true(all(vapply(
    malformed,
    function(record) identical(record$class, "competitive_hard_stop"),
    logical(1)
  )))
  expect_true(all(vapply(
    malformed,
    function(record) {
      identical(record$reason, "malformed_finite_endpoint_objective")
    },
    logical(1)
  )))
  expect_true(all(is.na(result$candidates$objective[3:4])))
})

test_that("malformed proposed-state callbacks can never claim a nonveto", {
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(ln_joint_resolution_test_review())
  )
  malformed_records <- list(
    list(accepted = TRUE),
    list()
  )

  for (malformed in malformed_records) {
    callback_calls <- 0L
    result <- ln_joint_resolve_candidates(
      candidates,
      basin_probe = function(...) {
        callback_calls <<- callback_calls + 1L
        malformed
      }
    )
    review <- result$resolution$candidate_resolutions[[3L]]

    expect_identical(callback_calls, 1L)
    expect_identical(review$class, "same_basin_review_unresolved")
    expect_true(review$veto)
    expect_identical(
      review$probe$termination_reason,
      "proposed_state_evaluation_failed"
    )
    expect_false(review$probe$accepted)
    expect_identical(review$probe$state_evaluations, 1L)
    expect_identical(
      result$resolution$evaluation_counts$probe_state_evaluations,
      1L
    )
    expect_identical(
      result$resolution$evaluation_counts$probe_objective_evaluations,
      1L
    )
    expect_identical(
      result$resolution$evaluation_counts$probe_score_evaluations,
      1L
    )
    expect_identical(
      result$resolution$evaluation_counts$probe_hessian_evaluations,
      1L
    )
  }
})

test_that("dominated objective-only handling rejects state and stale evidence", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  dominated <- ln_joint_resolution_test_candidate(
    "dominated_integrity",
    objective = 10.1,
    score = c(-2e-5, 0, 0, 0),
    required = FALSE,
    core = FALSE
  )

  atomic_state <- dominated
  atomic_state$state <- 1
  atomic_resolution <- ln_joint_classify_candidate(
    atomic_state,
    selected,
    ln_joint_resolution_thresholds()
  )
  expect_identical(atomic_resolution$class, "competitive_hard_stop")
  expect_identical(atomic_resolution$reason, "invalid_derivative_state")

  stale_evidence <- dominated
  stale_evidence$state <- NULL
  stale_resolution <- ln_joint_classify_candidate(
    stale_evidence,
    selected,
    ln_joint_resolution_thresholds()
  )
  expect_identical(stale_resolution$class, "competitive_hard_stop")
  expect_identical(stale_resolution$reason, "invalid_derivative_state")
})

test_that("candidate and derivative coordinates must be identical", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  contradictory <- ln_joint_resolution_test_review("coordinate_contradiction")
  contradictory$theta_reduced[[1L]] <-
    contradictory$theta_reduced[[1L]] + 1e-3
  contradictory$theta_reduced_log[[1L]] <-
    contradictory$theta_reduced_log[[1L]] + 1e-3

  resolution <- ln_joint_classify_candidate(
    contradictory,
    selected,
    ln_joint_resolution_thresholds()
  )
  expect_identical(resolution$class, "competitive_hard_stop")
  expect_identical(resolution$reason, "invalid_derivative_state")
})

test_that("candidate and selected roots must have the same donor columns", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  mismatch <- ln_joint_resolution_test_review("donor_column_mismatch")
  mismatch$state$U_subject_reduced <- cbind(
    mismatch$state$U_subject_reduced,
    rep(0, nrow(mismatch$state$U_subject_reduced))
  )

  resolution <- ln_joint_classify_candidate(
    mismatch,
    selected,
    ln_joint_resolution_thresholds()
  )
  expect_identical(resolution$class, "competitive_hard_stop")
  expect_identical(resolution$reason, "curvature_dimension_mismatch")
})

test_that("fabricated accepted probe records are only malformed proposed states", {
  thresholds <- ln_joint_resolution_thresholds()
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(ln_joint_resolution_test_review("fabricated_probe_state"))
  )
  empty_evidence <- ln_joint_basin_probe_template(thresholds)
  empty_evidence$accepted <- TRUE
  empty_evidence$termination_reason <-
    "same_basin_contraction_established"
  empty_evidence$state_evaluations <- 1L
  empty_evidence$objective_evaluations <- 1L
  empty_evidence$score_evaluations <- 1L
  empty_evidence$hessian_evaluations <- 1L
  forged_evidence <- empty_evidence
  forged_evidence$initial <- list(forged = TRUE)
  forged_evidence$proposed <- list(forged = TRUE)
  forged_evidence$contraction <- list(
    full_state = TRUE,
    beta_block = TRUE,
    dispersion_block = TRUE,
    root = TRUE,
    newton = TRUE
  )
  forged_evidence$anchor_assessments <- list(list(accepted = TRUE))

  for (fabricated in list(empty_evidence, forged_evidence)) {
    expect_identical(
      names(fabricated),
      names(ln_joint_basin_probe_template(thresholds))
    )
    callback_calls <- 0L
    result <- ln_joint_resolve_candidates(
      candidates,
      basin_probe = function(candidate, theta_reduced_log) {
        callback_calls <<- callback_calls + 1L
        fabricated
      },
      enable_prospective = TRUE
    )
    review <- result$resolution$candidate_resolutions[[3L]]

    expect_identical(callback_calls, 1L)
    expect_identical(review$class, "same_basin_review_unresolved")
    expect_true(review$veto)
    expect_false(review$probe$accepted)
    expect_identical(
      review$probe$termination_reason,
      "proposed_state_evaluation_failed"
    )
    expect_identical(review$probe$state_evaluations, 1L)
    expect_identical(
      result$resolution$evaluation_counts$probe_state_evaluations,
      1L
    )
    expect_false(result$resolution$would_nonveto)
    expect_identical(result$resolution$prospective_route, "R2")
    expect_identical(result$route, "R2")
    expect_identical(
      result$failure_reason,
      "same_basin_review_unresolved"
    )
  }
})

test_that("two distinct review endpoints exhaust the single fixed probe budget", {
  review_a <- ln_joint_resolution_test_review("review_endpoint_a")
  target <- c(0, 0, log(1), log(2))
  delta_b <- c(3e-5, 0, 0, 0)
  review_b <- ln_joint_resolution_test_candidate(
    "review_endpoint_b",
    theta_full_log = target + delta_b,
    objective = 10 + 0.5 * sum(delta_b^2),
    score = -delta_b,
    required = FALSE,
    core = FALSE
  )
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(review_a, review_b)
  )
  callback_calls <- 0L
  result <- ln_joint_resolve_candidates(
    candidates,
    basin_probe = function(candidate, theta_reduced_log) {
      callback_calls <<- callback_calls + 1L
      stop("distinct review endpoints must not receive partial probe effort")
    },
    enable_prospective = TRUE
  )
  reviews <- result$resolution$candidate_resolutions[3:4]

  expect_identical(callback_calls, 0L)
  expect_true(all(vapply(reviews, function(record) {
    identical(record$class, "same_basin_review_unresolved") &&
      isTRUE(record$veto) &&
      identical(
        record$reason,
        "multiple_distinct_review_endpoints_unresolved"
      ) &&
      identical(
        record$probe$termination_reason,
        "multiple_review_endpoints_exceed_fixed_probe_budget"
      ) &&
      identical(record$probe$state_evaluations, 0L)
  }, logical(1))))
  expect_identical(
    result$resolution$evaluation_counts$probe_state_evaluations,
    0L
  )
  expect_identical(
    result$resolution$evaluation_counts$probe_objective_evaluations,
    0L
  )
  expect_identical(result$resolution$prospective_route, "R2")
  expect_identical(result$route, "R2")
})

test_that("every proposal threshold is closed at its deterministic edge", {
  thresholds <- ln_joint_resolution_thresholds()
  edge_multiplier <- c(below = 1 - 1e-4, at = 1, above = 1 + 1e-4)
  target_fields <- c(
    "proposal_full_state", "proposal_beta_block",
    "proposal_dispersion_block", "proposal_objective_equivalence",
    "proposal_root", "proposal_newton"
  )
  for (field in target_fields) {
    target <- thresholds[[field]]
    edge <- target * edge_multiplier
    expect_true(
      ln_joint_probe_contracted_or_target(
        2 * target, edge[["below"]], target,
        thresholds$contraction_ratio
      ),
      info = paste(field, "below")
    )
    expect_true(
      ln_joint_probe_contracted_or_target(
        2 * target, edge[["at"]], target,
        thresholds$contraction_ratio
      ),
      info = paste(field, "at")
    )
    expect_false(
      ln_joint_probe_contracted_or_target(
        2 * target, edge[["above"]], target,
        thresholds$contraction_ratio
      ),
      info = paste(field, "above")
    )
  }

  scaled_coordinate <- function(metric) metric / (1 - metric)
  beta_observed <- vapply(
    thresholds$proposal_beta_block * edge_multiplier,
    function(metric) {
      theta <- c(scaled_coordinate(metric), 0, 0, 0)
      ln_joint_resolution_block_metrics(theta, rep(0, 4L), 2L)$beta_block
    },
    numeric(1)
  )
  full_observed <- vapply(
    thresholds$proposal_full_state * edge_multiplier,
    function(metric) {
      theta <- c(scaled_coordinate(metric), 0, 0, 0)
      ln_joint_resolution_block_metrics(theta, rep(0, 4L), 2L)$full_state
    },
    numeric(1)
  )
  dispersion_observed <- vapply(
    thresholds$proposal_dispersion_block * edge_multiplier,
    function(metric) {
      theta <- c(0, 0, scaled_coordinate(metric), 0)
      ln_joint_resolution_block_metrics(
        theta,
        rep(0, 4L),
        2L
      )$dispersion_block
    },
    numeric(1)
  )
  for (observed in list(
      full_state = full_observed,
      beta_block = beta_observed,
      dispersion_block = dispersion_observed)) {
    expect_lt(observed[["below"]], observed[["at"]])
    expect_lt(observed[["at"]], observed[["above"]])
  }
  expect_equal(
    full_observed[["at"]],
    thresholds$proposal_full_state,
    tolerance = .Machine$double.eps
  )
  expect_equal(
    beta_observed[["at"]],
    thresholds$proposal_beta_block,
    tolerance = .Machine$double.eps
  )
  expect_equal(
    dispersion_observed[["at"]],
    thresholds$proposal_dispersion_block,
    tolerance = .Machine$double.eps
  )

  objective_edges <- lapply(
    thresholds$proposal_objective_equivalence * edge_multiplier,
    function(metric) {
      ln_joint_objective_comparison(
        scaled_coordinate(metric),
        0,
        thresholds$proposal_objective_equivalence
      )
    }
  )
  expect_identical(objective_edges$below$relation, "equivalent")
  expect_identical(objective_edges$at$relation, "equivalent")
  expect_identical(objective_edges$above$relation, "dominated")
  expect_equal(
    objective_edges$at$absolute_gap,
    thresholds$proposal_objective_equivalence,
    tolerance = .Machine$double.eps
  )

  root_observed <- vapply(
    thresholds$proposal_root * edge_multiplier,
    function(metric) {
      ln_joint_root_metric(
        matrix(0, nrow = 4L, ncol = 2L),
        c(metric, 0, 0, 0)
      )
    },
    numeric(1)
  )
  newton_observed <- vapply(
    thresholds$proposal_newton * edge_multiplier,
    function(metric) {
      ln_joint_probe_newton_metric(list(
        theta_reduced_log = rep(0, 4L),
        U_reduced = c(metric, 0, 0, 0),
        A_reduced = diag(4L)
      ))
    },
    numeric(1)
  )
  expect_identical(
    unname(root_observed),
    unname(thresholds$proposal_root * edge_multiplier)
  )
  expect_identical(
    unname(newton_observed),
    unname(thresholds$proposal_newton * edge_multiplier)
  )

  ratio_edge <- thresholds$contraction_ratio * edge_multiplier
  expect_false(ln_joint_probe_contracted_or_target(
    initial = 2, proposed = 1, target = 1, ratio = ratio_edge[["below"]]
  ))
  expect_true(ln_joint_probe_contracted_or_target(
    initial = 2, proposed = 1, target = 1, ratio = ratio_edge[["at"]]
  ))
  expect_true(ln_joint_probe_contracted_or_target(
    initial = 2, proposed = 1, target = 1, ratio = ratio_edge[["above"]]
  ))

  budget_valid <- vapply(0:2, function(budget) {
    changed <- thresholds
    changed$max_state_evaluations <- as.integer(budget)
    ln_joint_resolution_thresholds_valid(changed)
  }, logical(1))
  expect_identical(unname(budget_valid), c(FALSE, TRUE, FALSE))
})

test_that("atomic curvature and absent labels fail safely or serialize as missing", {
  atomic_curvature <- ln_joint_resolution_test_review("atomic_curvature")
  atomic_curvature$curvature <- 1
  atomic_result <- NULL
  expect_no_error({
    atomic_result <- ln_joint_resolve_candidates(c(
      ln_joint_resolution_test_core(),
      list(atomic_curvature)
    ))
  })
  atomic_record <- atomic_result$resolution$candidate_resolutions[[3L]]
  expect_identical(atomic_record$class, "competitive_hard_stop")
  expect_identical(
    atomic_record$reason,
    "missing_or_nonfinite_check_metric"
  )

  missing_label <- ln_joint_resolution_test_review("remove_this_label")
  missing_label$label <- NULL
  missing_result <- NULL
  expect_no_error({
    missing_result <- ln_joint_resolve_candidates(
      c(ln_joint_resolution_test_core(), list(missing_label)),
      basin_probe = ln_joint_resolution_test_probe
    )
  })
  missing_record <- missing_result$resolution$candidate_resolutions[[3L]]
  expect_true(is.na(missing_record$label))
  expect_true(is.na(missing_result$candidates$label[[3L]]))
  expect_identical(missing_record$class, "same_basin_review_nonveto")
  expect_identical(nrow(missing_result$candidates), 3L)
})

test_that("unequal-objective certified anchors are all assessed selected first", {
  selected <- ln_joint_resolution_test_candidate(
    "lower_objective_selected",
    objective = 10
  )
  other <- ln_joint_resolution_test_candidate(
    "higher_objective_certified",
    objective = 10 + 1e-10
  )
  review <- ln_joint_resolution_test_review("all_anchor_review")
  result <- ln_joint_resolve_candidates(
    list(other, selected, review),
    basin_probe = ln_joint_resolution_test_probe
  )
  probe <- result$resolution$candidate_resolutions[[3L]]$probe

  expect_identical(result$best_index, 2L)
  expect_true(probe$accepted)
  expect_length(probe$anchor_assessments, 2L)
  expect_identical(
    vapply(
      probe$anchor_assessments,
      `[[`,
      character(1),
      "anchor_candidate_key"
    ),
    c(
      ln_joint_resolution_candidate_key(selected),
      ln_joint_resolution_candidate_key(other)
    )
  )
  expect_true(all(vapply(
    probe$anchor_assessments,
    `[[`,
    logical(1),
    "accepted"
  )))
  objective_gaps <- vapply(
    probe$anchor_assessments,
    `[[`,
    numeric(1),
    "objective_absolute_gap"
  )
  expect_identical(objective_gaps[[1L]], 0)
  expect_gt(objective_gaps[[2L]], 0)
  expect_lte(
    objective_gaps[[2L]],
    ln_joint_resolution_thresholds()$proposal_objective_equivalence
  )
})

test_that("an echoed reduced proposal cannot hide a snapped full state", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  target <- selected$theta_full_log
  delta <- c(2e-5, 0, 0, 0)
  candidate <- ln_joint_resolution_test_candidate(
    "snapped_full_state",
    theta_full_log = target + delta,
    objective = 10 + 0.5 * sum(delta^2),
    score = -0.75 * delta,
    required = FALSE,
    core = FALSE
  )
  evaluations <- 0L
  evaluator <- function(theta_reduced_log) {
    evaluations <<- evaluations + 1L
    state <- ln_joint_resolution_test_probe(candidate, theta_reduced_log)
    state$full_log <- selected$theta_full_log
    state
  }
  probe <- ln_joint_basin_probe(
    candidate,
    selected,
    evaluator,
    ln_joint_resolution_thresholds(),
    assessment_roots = ln_joint_resolution_test_core()
  )

  expect_false(probe$accepted)
  expect_identical(
    probe$termination_reason,
    "proposed_state_full_map_mismatch"
  )
  expect_identical(evaluations, 1L)
  expect_identical(probe$state_evaluations, 1L)
  expect_identical(
    probe$proposed$theta_reduced_log,
    ln_joint_basin_local_proposal(
      candidate,
      ln_joint_resolution_thresholds()
    )$proposed_theta_reduced_log
  )
})

test_that("a vector label cannot inflate the candidate summary", {
  vector_label <- ln_joint_resolution_test_review("first_label")
  vector_label$label <- c("first_label", "untrusted_extra_label")
  candidates <- c(
    ln_joint_resolution_test_core(),
    list(vector_label)
  )
  result <- NULL
  expect_no_error({
    result <- ln_joint_resolve_candidates(
      candidates,
      basin_probe = ln_joint_resolution_test_probe
    )
  })

  expect_identical(nrow(result$candidates), length(candidates))
  expect_length(
    result$resolution$candidate_resolutions,
    length(candidates)
  )
  expect_identical(result$candidates$label[[3L]], "first_label")
  expect_identical(
    result$resolution$candidate_resolutions[[3L]]$label,
    "first_label"
  )
})

test_that("coherent dominated curvature failures remain optional diagnostics", {
  selected <- ln_joint_resolution_test_core()[[1L]]
  indefinite <- ln_joint_resolution_test_candidate(
    "dominated_indefinite_curvature",
    objective = 10.1,
    score = rep(0, 4L),
    A = diag(c(-1, 1, 1, 1)),
    required = FALSE,
    core = FALSE
  )
  singular <- ln_joint_resolution_test_candidate(
    "dominated_singular_curvature",
    objective = 10.1,
    score = rep(0, 4L),
    A = diag(c(0, 1, 1, 1)),
    required = FALSE,
    core = FALSE
  )

  expect_false(indefinite$curvature$positive_definite)
  expect_identical(indefinite$curvature$scaled_min_eigenvalue, -1)
  expect_identical(indefinite$curvature$scaled_condition_number, Inf)
  expect_false(indefinite$checks[["curvature"]])
  expect_false(indefinite$checks[["condition"]])
  expect_true(all(indefinite$checks[setdiff(
    names(indefinite$checks),
    c("curvature", "condition")
  )]))

  expect_false(singular$curvature$positive_definite)
  expect_identical(singular$curvature$scaled_min_eigenvalue, 0)
  expect_identical(singular$curvature$scaled_condition_number, Inf)
  expect_true(is.na(singular$newton_displacement))
  expect_false(singular$checks[["newton"]])
  expect_false(singular$checks[["curvature"]])
  expect_false(singular$checks[["condition"]])
  expect_true(all(singular$checks[setdiff(
    names(singular$checks),
    c("newton", "curvature", "condition")
  )]))

  indefinite_resolution <- ln_joint_classify_candidate(
    indefinite,
    selected,
    ln_joint_resolution_thresholds()
  )
  singular_resolution <- ln_joint_classify_candidate(
    singular,
    selected,
    ln_joint_resolution_thresholds()
  )
  for (resolution in list(indefinite_resolution, singular_resolution)) {
    expect_identical(resolution$objective_relation, "dominated")
    expect_identical(resolution$class, "dominated_uncertified")
    expect_false(resolution$veto)
    expect_identical(
      resolution$reason,
      "optional_candidate_objectively_dominated"
    )
  }

  global <- ln_joint_resolve_candidates(c(
    ln_joint_resolution_test_core(),
    list(indefinite, singular)
  ))
  expect_identical(global$route, "R0")
  expect_identical(global$resolution$prospective_route, "R0")
  expect_length(global$resolution$global_veto_classes, 0L)
  expect_true(all(!vapply(
    global$resolution$candidate_resolutions[3:4],
    `[[`,
    logical(1),
    "veto"
  )))
})

test_that("vector failure diagnostics cannot inflate a no-root resolution", {
  candidate <- ln_joint_resolution_test_review("only_uncertified_endpoint")
  candidate$failure <- c("a", "b")
  result <- NULL
  expect_no_error({
    result <- ln_joint_resolve_candidates(list(candidate))
  })

  expect_identical(result$best_index, 1L)
  expect_identical(result$route, "R2")
  expect_identical(result$failure_reason, "a")
  expect_identical(result$resolution$strict_route, "R2")
  expect_identical(result$resolution$strict_failure_reason, "a")
  expect_identical(result$resolution$prospective_route, "R2")
  expect_identical(
    result$resolution$candidate_resolutions[[1L]]$class,
    "competitive_hard_stop"
  )
  expect_identical(nrow(result$candidates), 1L)
  expect_identical(result$candidates$failure[[1L]], "a")
  restored <- unserialize(serialize(result, NULL))
  expect_identical(nrow(restored$candidates), 1L)
  expect_identical(restored$candidates, result$candidates)
  expect_identical(restored$resolution, result$resolution)
})

test_that("certified flags cannot override absent finite endpoints", {
  first <- ln_joint_resolution_test_candidate(
    "malformed_certified_without_endpoint_a",
    objective = 10,
    required = FALSE,
    core = TRUE,
    certified = TRUE
  )
  second <- ln_joint_resolution_test_candidate(
    "malformed_certified_without_endpoint_b",
    objective = 10 + 1e-10,
    required = FALSE,
    core = TRUE,
    certified = TRUE
  )
  first$finite_endpoint <- FALSE
  second$finite_endpoint <- FALSE
  candidates <- list(first, second)
  result <- enabled <- NULL
  expect_no_error({
    result <- ln_joint_resolve_candidates(candidates)
    enabled <- ln_joint_resolve_candidates(
      candidates,
      enable_prospective = TRUE
    )
  })

  expect_true(all(vapply(candidates, is.list, logical(1))))
  expect_true(all(vapply(candidates, `[[`, logical(1), "certified")))
  expect_true(all(vapply(candidates, `[[`, logical(1), "core")))
  expect_true(all(vapply(
    candidates,
    ln_joint_finite_scalar_objective,
    logical(1)
  )))
  expect_true(all(!vapply(
    candidates,
    `[[`,
    logical(1),
    "finite_endpoint"
  )))
  expect_true(is.na(result$best_index))
  expect_identical(result$route, "R2")
  expect_identical(result$failure_reason, "no_certified_candidate")
  expect_identical(result$resolution$strict_route, "R2")
  expect_identical(
    result$resolution$strict_failure_reason,
    "no_certified_candidate"
  )
  expect_identical(result$resolution$prospective_route, "R2")
  expect_identical(enabled$route, "R2")
  expect_true(all(vapply(
    result$resolution$candidate_resolutions,
    function(record) {
      identical(record$class, "competitive_hard_stop") && isTRUE(record$veto)
    },
    logical(1)
  )))
  expect_identical(nrow(result$candidates), 2L)
})

test_that("malformed certified field types fail closed without row inflation", {
  base <- ln_joint_resolution_test_candidate(
    "malformed_certified_field",
    required = TRUE,
    core = TRUE,
    certified = TRUE
  )
  mutations <- list(
    n_beta_list = function(x) {
      x$n_beta <- list(2L)
      x
    },
    n_beta_complex = function(x) {
      x$n_beta <- 2 + 0i
      x
    },
    n_beta_matrix = function(x) {
      x$n_beta <- matrix(2L, nrow = 1L, ncol = 1L)
      x
    },
    closure_label = function(x) {
      x$label <- function() "not_a_serializable_scalar_label"
      x
    },
    nested_theta_full_log = function(x) {
      x$theta_full_log <- list(x$theta_full_log)
      x
    },
    nested_theta_reduced_log = function(x) {
      x$theta_reduced_log <- list(x$theta_reduced_log)
      x
    },
    nested_theta_reduced = function(x) {
      x$theta_reduced <- list(x$theta_reduced)
      x
    },
    nested_bounds = function(x) {
      x$bounds <- list(x$bounds)
      x
    },
    nested_probe_bounds = function(x) {
      x$probe_bounds <- list(x$probe_bounds)
      x
    },
    nested_curvature = function(x) {
      x$curvature <- list(x$curvature)
      x
    },
    nested_constraint = function(x) {
      x$constraint_residual <- list(x$constraint_residual)
      x
    },
    inconsistent_natural_theta = function(x) {
      x$theta_reduced[[1L]] <- x$theta_reduced[[1L]] + 0.25
      x
    }
  )

  for (case_name in names(mutations)) {
    candidate <- mutations[[case_name]](base)
    result <- enabled <- NULL
    expect_no_error({
      result <- ln_joint_resolve_candidates(list(candidate))
      enabled <- ln_joint_resolve_candidates(
        list(candidate),
        enable_prospective = TRUE
      )
    })

    expect_identical(result$route, "R2", info = case_name)
    expect_identical(
      result$resolution$strict_route,
      "R2",
      info = case_name
    )
    expect_identical(
      result$resolution$prospective_route,
      "R2",
      info = case_name
    )
    expect_identical(enabled$route, "R2", info = case_name)
    expect_identical(nrow(result$candidates), 1L, info = case_name)
    expect_identical(
      length(result$resolution$candidate_resolutions),
      1L,
      info = case_name
    )
    restored <- unserialize(serialize(result$candidates, NULL))
    expect_identical(nrow(restored), 1L, info = case_name)
  }
})

test_that("callable threshold policy metadata fails closed without evaluation", {
  malformed <- ln_joint_resolution_thresholds()
  malformed$policy_version <- function() "forged_policy"

  template <- empty <- probe <- NULL
  expect_no_error({
    template <- ln_joint_basin_probe_template(malformed)
    empty <- ln_joint_resolve_candidates(
      list(),
      resolution_thresholds = malformed
    )
    probe <- ln_joint_basin_probe(
      ln_joint_resolution_test_review("invalid_policy_probe"),
      ln_joint_resolution_test_core()[[1L]],
      function(theta_reduced_log) {
        stop("an invalid policy must prevent endpoint evaluation")
      },
      resolution_thresholds = malformed
    )
  })

  expect_true(is.na(template$policy_version))
  expect_false(template$accepted)
  expect_identical(template$state_evaluations, 0L)
  expect_identical(empty$route, "R2")
  expect_identical(empty$resolution$strict_route, "R2")
  expect_identical(empty$resolution$prospective_route, "R2")
  expect_true(is.na(empty$resolution$policy_version))
  expect_length(empty$resolution$candidate_resolutions, 0L)
  expect_false(probe$accepted)
  expect_true(is.na(probe$policy_version))
  expect_identical(probe$termination_reason, "invalid_probe_configuration")
  expect_identical(probe$state_evaluations, 0L)
})

test_that("noncanonical top-level candidate inputs become one hard-stop row", {
  environment_input <- new.env(parent = emptyenv())
  function_input <- function() "not_a_candidate_collection"
  data_frame_input <- data.frame(value = 1, stringsAsFactors = FALSE)
  single_candidate_input <- ln_joint_resolution_test_review(
    "unwrapped_single_candidate"
  )
  inputs <- list(
    environment = environment_input,
    function_value = function_input,
    data_frame = data_frame_input,
    unwrapped_candidate = single_candidate_input
  )

  for (case_name in names(inputs)) {
    result <- enabled <- NULL
    expect_no_error({
      result <- ln_joint_resolve_candidates(inputs[[case_name]])
      enabled <- ln_joint_resolve_candidates(
        inputs[[case_name]],
        enable_prospective = TRUE
      )
    })

    expect_identical(result$route, "R2", info = case_name)
    expect_identical(
      result$resolution$strict_route,
      "R2",
      info = case_name
    )
    expect_identical(
      result$resolution$prospective_route,
      "R2",
      info = case_name
    )
    expect_identical(enabled$route, "R2", info = case_name)
    expect_identical(nrow(result$candidates), 1L, info = case_name)
    expect_identical(
      length(result$resolution$candidate_resolutions),
      1L,
      info = case_name
    )
    record <- result$resolution$candidate_resolutions[[1L]]
    expect_identical(record$class, "competitive_hard_stop", info = case_name)
    expect_true(record$veto, info = case_name)
  }
})

test_that("unreviewable strict failures skip bounded review without evaluation", {
  base_theta <- c(0, 0, log(1), log(2))
  parameter_roots <- list(
    ln_joint_resolution_test_candidate(
      "parameter_selected",
      theta_full_log = base_theta,
      objective = 10
    ),
    ln_joint_resolution_test_candidate(
      "parameter_review_zone_root",
      theta_full_log = base_theta + c(2e-5, 0, 0, 0),
      objective = 10
    )
  )
  objective_roots <- list(
    ln_joint_resolution_test_candidate(
      "objective_selected",
      objective = 10
    ),
    ln_joint_resolution_test_candidate(
      "objective_review_zone_root",
      objective = 10 + 2e-7
    )
  )
  insufficient_roots <- list(
    ln_joint_resolution_test_candidate(
      "only_certified_core",
      objective = 10
    )
  )
  scenarios <- list(
    start_parameter_review_zone = c(
      parameter_roots,
      list(ln_joint_resolution_test_review("parameter_skipped_review"))
    ),
    start_objective_review_zone = c(
      objective_roots,
      list(ln_joint_resolution_test_review("objective_skipped_review"))
    ),
    insufficient_certified_core_starts = c(
      insufficient_roots,
      list(ln_joint_resolution_test_review("core_skipped_review"))
    )
  )

  for (failure_name in names(scenarios)) {
    callback_calls <- 0L
    result <- ln_joint_resolve_candidates(
      scenarios[[failure_name]],
      basin_probe = function(candidate, theta_reduced_log) {
        callback_calls <<- callback_calls + 1L
        stop("an unreviewable strict failure must skip bounded review")
      }
    )
    review <- result$resolution$candidate_resolutions[[length(
      scenarios[[failure_name]]
    )]]

    expect_identical(result$route, "R2", info = failure_name)
    expect_identical(result$failure_reason, failure_name)
    expect_identical(
      result$resolution$prospective_route,
      "R2",
      info = failure_name
    )
    expect_identical(callback_calls, 0L, info = failure_name)
    expect_identical(
      review$class,
      "same_basin_review_unresolved",
      info = failure_name
    )
    expect_true(review$veto, info = failure_name)
    expect_identical(
      review$reason,
      "bounded_review_skipped_due_to_strict_failure",
      info = failure_name
    )
    expect_identical(
      review$probe$termination_reason,
      "skipped_strict_failure_not_reviewable",
      info = failure_name
    )
    expect_identical(
      review$probe$state_evaluations,
      0L,
      info = failure_name
    )
    counts <- result$resolution$evaluation_counts
    expect_identical(counts$probe_state_evaluations, 0L, info = failure_name)
    expect_identical(
      counts$probe_objective_evaluations,
      0L,
      info = failure_name
    )
    expect_identical(counts$probe_score_evaluations, 0L, info = failure_name)
    expect_identical(
      counts$probe_hessian_evaluations,
      0L,
      info = failure_name
    )
  }
})

test_that("noncanonical candidate role flags fail closed", {
  flag_names <- c("required", "core", "finite_endpoint", "certified")
  malformed_values <- list(
    numeric = 1,
    missing_logical = NA,
    vector_logical = c(TRUE, FALSE)
  )

  for (flag_name in flag_names) {
    for (value_name in names(malformed_values)) {
      case_name <- paste(flag_name, value_name, sep = "_")
      candidate <- ln_joint_resolution_test_candidate(
        paste0("malformed_role_", case_name),
        required = TRUE,
        core = TRUE,
        certified = TRUE
      )
      candidate[[flag_name]] <- malformed_values[[value_name]]
      result <- enabled <- NULL
      expect_no_error({
        result <- ln_joint_resolve_candidates(list(candidate))
        enabled <- ln_joint_resolve_candidates(
          list(candidate),
          enable_prospective = TRUE
        )
      })

      expect_identical(result$route, "R2", info = case_name)
      expect_identical(
        result$resolution$strict_route,
        "R2",
        info = case_name
      )
      expect_identical(
        result$resolution$prospective_route,
        "R2",
        info = case_name
      )
      expect_identical(enabled$route, "R2", info = case_name)
      expect_identical(nrow(result$candidates), 1L, info = case_name)
      expect_identical(
        length(result$resolution$candidate_resolutions),
        1L,
        info = case_name
      )
      record <- result$resolution$candidate_resolutions[[1L]]
      expect_identical(
        record$class,
        "competitive_hard_stop",
        info = case_name
      )
      expect_true(record$veto, info = case_name)
      restored <- unserialize(serialize(result, NULL))
      expect_identical(nrow(restored$candidates), 1L, info = case_name)
      expect_identical(restored$candidates, result$candidates, info = case_name)
    }
  }
})
