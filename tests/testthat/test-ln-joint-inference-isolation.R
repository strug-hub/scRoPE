ln_isolation_accept_probe <- function(candidate, theta_reduced_log) {
  initial <- candidate$state
  A <- initial$A_reduced
  n_subject <- ncol(initial$U_subject_reduced)
  theta_reduced_log <- as.numeric(theta_reduced_log)
  delta <- theta_reduced_log - initial$theta_reduced_log
  score <- as.numeric(initial$U_reduced - A %*% delta)
  list(
    theta_reduced_log = theta_reduced_log,
    full_log = ln_joint_probe_full_log(theta_reduced_log, candidate),
    objective_provenance = candidate$objective_provenance,
    U_reduced = score,
    U_subject_reduced = matrix(
      rep(score / n_subject, n_subject),
      nrow = length(score)
    ),
    A_reduced = A,
    raw = list(
      value = candidate$objective -
        as.numeric(crossprod(initial$U_reduced, delta)) +
        0.5 * as.numeric(crossprod(delta, A %*% delta))
    )
  )
}

ln_isolation_refresh_candidate <- function(candidate) {
  state <- candidate$state
  thresholds <- ln_joint_thresholds()
  candidate$donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(state$U_subject_reduced),
    state$U_reduced
  )
  candidate$root_residual <- ln_joint_root_metric(
    state$U_subject_reduced,
    state$U_reduced
  )
  candidate$newton_displacement <- ln_joint_probe_newton_metric(state)
  candidate$curvature <- ln_joint_curvature_diagnostics(state$A_reduced)
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

ln_isolation_update_candidate_state <- function(candidate, full_log,
                                                score = NULL) {
  n_beta <- as.integer(candidate$n_beta)
  candidate$theta_full_log <- as.numeric(full_log)
  candidate$state$full_log <- as.numeric(full_log)
  candidate$state$theta_reduced_log <- as.numeric(full_log)
  candidate$theta_reduced_log <- as.numeric(full_log)
  candidate$theta_reduced <- c(
    full_log[seq_len(n_beta)],
    exp(full_log[n_beta + seq_len(2L)])
  )
  candidate$state$beta <- full_log[seq_len(n_beta)]
  candidate$state$sigma2 <- exp(full_log[[n_beta + 1L]])
  candidate$state$gamma <- exp(full_log[[n_beta + 2L]])
  if (!is.null(score)) {
    score <- as.numeric(score)
    n_subject <- ncol(candidate$state$U_subject_reduced)
    candidate$state$U_reduced <- score
    candidate$state$U_subject_reduced <- matrix(
      rep(score / n_subject, n_subject),
      nrow = length(score)
    )
  }
  if (is.list(candidate$state$raw)) {
    candidate$state$raw$value <- as.numeric(candidate$objective)
  }
  ln_isolation_refresh_candidate(candidate)
}

ln_isolation_review_transition <- function(candidate, target = 2e-5) {
  A <- 0.5 * (
    candidate$state$A_reduced + t(candidate$state$A_reduced)
  )
  theta <- as.numeric(candidate$state$theta_reduced_log)
  stopifnot(length(theta) == length(candidate$theta_full_log))
  eig <- eigen(A, symmetric = TRUE)
  direction <- eig$vectors[, which.min(eig$values)]
  direction <- direction /
    max(abs(direction) / (1 + abs(theta)))
  score_base <- as.numeric(-A %*% direction)
  metric_at <- function(scale) {
    score <- scale * score_base
    state <- candidate$state
    state$theta_reduced_log <- theta + scale * direction
    state$U_reduced <- score
    n_subject <- ncol(state$U_subject_reduced)
    state$U_subject_reduced <- matrix(
      rep(score / n_subject, n_subject),
      nrow = length(score)
    )
    max(
      ln_joint_root_metric(state$U_subject_reduced, state$U_reduced),
      ln_joint_probe_newton_metric(state)
    )
  }
  scale <- target / metric_at(1)
  for (iteration in seq_len(4L)) {
    scale <- scale * target / metric_at(scale)
  }
  list(
    full_log = candidate$theta_full_log + scale * direction,
    score = scale * score_base
  )
}

ln_isolation_injected_candidate <- function(selected, kind) {
  candidate <- selected
  candidate$label <- paste0("synthetic_", kind)
  candidate$start_provenance <- "synthetic_test_fixture"
  candidate$required <- FALSE
  candidate$core <- FALSE
  candidate$finite_endpoint <- TRUE
  candidate$certified <- FALSE
  candidate$objective <- selected$objective
  candidate$optimizer$message <- paste0("synthetic_", kind)
  full_log <- as.numeric(selected$theta_full_log)
  n_beta <- as.integer(selected$n_beta)

  if (kind %in% c("same_basin_shadow", "same_basin_unresolved")) {
    transition <- ln_isolation_review_transition(candidate)
    full_log <- transition$full_log
    delta <- full_log - selected$theta_full_log
    candidate$objective <- selected$objective + 0.5 * as.numeric(
      crossprod(delta, selected$state$A_reduced %*% delta)
    )
    candidate <- ln_isolation_update_candidate_state(
      candidate,
      full_log,
      score = transition$score
    )
  } else if (identical(kind, "competitive_boundary")) {
    full_log[[n_beta + 1L]] <- log(candidate$bounds$lower[[1L]])
    candidate <- ln_isolation_update_candidate_state(candidate, full_log)
  } else if (identical(kind, "competitive_distinct_basin")) {
    full_log[[1L]] <- full_log[[1L]] + 0.01
    candidate <- ln_isolation_update_candidate_state(candidate, full_log)
  } else if (identical(kind, "competitive_hard_stop")) {
    full_log[[1L]] <- full_log[[1L]] + 1e-5
    candidate$constraint_residual <- 2 * ln_joint_thresholds()$constraint
    candidate <- ln_isolation_update_candidate_state(candidate, full_log)
  } else {
    stop("Unknown synthetic candidate kind: ", kind)
  }

  if (kind %in% c("same_basin_shadow", "same_basin_unresolved")) {
    failed <- names(candidate$checks)[!candidate$checks]
    stopifnot(
      length(failed) > 0L,
      all(failed %in% c("root", "newton")),
      candidate$root_residual <= ln_joint_thresholds()$hard_stop$root,
      candidate$newton_displacement <= ln_joint_thresholds()$hard_stop$newton
    )
  }
  candidate
}

ln_isolation_fit_with_candidate <- function(
    kind,
    expected_production_enable = TRUE,
    resolver_enable_override = NULL,
    covariance_failure = FALSE,
    additional_tests = character(0)) {
  original_resolver <- getFromNamespace(
    "ln_joint_resolve_candidates",
    "scRoPE"
  )
  original_covariance <- getFromNamespace("ln_joint_covariance", "scRoPE")
  injected_capture <- NULL
  force(kind)
  force(expected_production_enable)
  force(resolver_enable_override)
  force(covariance_failure)
  resolver_wrapper <- function(
      candidates,
      resolution_thresholds = ln_joint_resolution_thresholds(),
      basin_probe = NULL,
      enable_prospective = FALSE) {
    if (!isTRUE(enable_prospective) &&
        isTRUE(expected_production_enable)) {
      return(original_resolver(
        candidates,
        resolution_thresholds = resolution_thresholds,
        basin_probe = basin_probe,
        enable_prospective = FALSE
      ))
    }
    stopifnot(identical(
      isTRUE(enable_prospective),
      isTRUE(expected_production_enable)
    ))
    certified <- which(vapply(candidates, function(candidate) {
      isTRUE(candidate$certified)
    }, logical(1)))
    stopifnot(length(certified) > 0L)
    selected_index <- certified[[which.min(vapply(
      candidates[certified],
      `[[`,
      numeric(1),
      "objective"
    ))]]
    injected <- ln_isolation_injected_candidate(
      candidates[[selected_index]],
      kind
    )
    injected_capture <<- injected
    delegated_probe <- basin_probe
    combined_probe <- function(candidate, theta_reduced_log) {
      if (identical(candidate$label, injected$label)) {
        if (!kind %in% c("same_basin_shadow", "same_basin_unresolved")) {
          stop("Immediate-veto candidate unexpectedly entered basin review")
        }
        if (identical(kind, "same_basin_unresolved")) {
          return(NULL)
        }
        if (!is.function(delegated_probe)) {
          stop("Production basin probe was unavailable")
        }
        return(delegated_probe(
          candidate,
          theta_reduced_log
        ))
      }
      if (is.function(delegated_probe)) {
        return(delegated_probe(
          candidate,
          theta_reduced_log
        ))
      }
      NULL
    }
    original_resolver(
      c(candidates, list(injected)),
      resolution_thresholds = resolution_thresholds,
      basin_probe = combined_probe,
      enable_prospective = if (is.null(resolver_enable_override)) {
        enable_prospective
      } else {
        isTRUE(resolver_enable_override)
      }
    )
  }
  covariance_wrapper <- function(state, nb) {
    if (isTRUE(covariance_failure)) {
      return(list(
        success = FALSE,
        failure_reason = "synthetic_covariance_failure"
      ))
    }
    original_covariance(state, nb)
  }
  local_mocked_bindings(
    ln_joint_resolve_candidates = resolver_wrapper,
    ln_joint_covariance = covariance_wrapper,
    .package = "scRoPE"
  )
  fit <- suppressWarnings(ln_joint_fit_fixture(
    additional_tests = additional_tests
  )$fit)
  attr(fit, "synthetic_injected_candidate") <- injected_capture
  fit
}

ln_isolation_public_inference_columns <- function(summary) {
  grep(
    paste0(
      "^se_|^se_robust_|^p_naive_|^p_robust_|^wald_stat_|",
      "^wald_success_|^wald_failure_|^ci_|_ci_|_lower$|_upper$"
    ),
    names(summary),
    value = TRUE
  )
}

ln_isolation_expect_withheld <- function(fit) {
  columns <- ln_isolation_public_inference_columns(fit$summary)
  numeric_columns <- columns[vapply(
    fit$summary[columns],
    is.numeric,
    logical(1)
  )]
  logical_columns <- columns[vapply(
    fit$summary[columns],
    is.logical,
    logical(1)
  )]
  expect_true(all(is.na(as.matrix(fit$summary[numeric_columns]))))
  expect_true(all(!as.matrix(fit$summary[logical_columns])))
}

ln_isolation_robust_wald_ci <- function(summary, level = 0.95) {
  se_names <- grep("^se_robust_", names(summary), value = TRUE)
  stopifnot(length(se_names) > 0L)
  z <- stats::qnorm(1 - (1 - level) / 2)
  do.call(cbind, lapply(se_names, function(se_name) {
    suffix <- sub("^se_robust_", "", se_name)
    estimate <- summary[[paste0("logFC_", suffix)]]
    se <- summary[[se_name]]
    cbind(lower = estimate - z * se, upper = estimate + z * se)
  }))
}

ln_isolation_synthetic_candidate <- function(
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
    objective_provenance = ln_joint_objective_provenance(),
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
    start_provenance = "synthetic_test_fixture",
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
    objective_provenance = ln_joint_objective_provenance(),
    objective_instance = list(
      scheme = "synthetic_joint_ln_objective_context_v1",
      id = "synthetic-inference-isolation",
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
    donor_additivity_scaled = ln_joint_scaled_max(
      rowSums(score_subject),
      score
    ),
    root_residual = ln_joint_root_metric(score_subject, score),
    newton_displacement = ln_joint_probe_newton_metric(state),
    curvature = ln_joint_curvature_diagnostics(state$A_reduced),
    interior = TRUE,
    constraint_residual = 0,
    checks = setNames(
      rep(TRUE, length(ln_joint_candidate_check_names())),
      ln_joint_candidate_check_names()
    ),
    failure = NA_character_
  )
  ln_isolation_refresh_candidate(candidate)
}

test_that("same-basin shadow is default-off and enabled inference uses only the selected root", {
  baseline <- ln_joint_fit_fixture()$fit
  shadow <- ln_isolation_fit_with_candidate(
    "same_basin_shadow",
    resolver_enable_override = FALSE
  )
  enabled <- ln_isolation_fit_with_candidate("same_basin_shadow")

  baseline_raw <- baseline$diagnostics$raw[[1L]]
  shadow_raw <- shadow$diagnostics$raw[[1L]]
  enabled_raw <- enabled$diagnostics$raw[[1L]]
  shadow_resolution <- shadow_raw$ln_joint$resolution
  enabled_resolution <- enabled_raw$ln_joint$resolution
  shadow_candidate <- Filter(function(record) {
    identical(record$label, "synthetic_same_basin_shadow")
  }, shadow_resolution$candidate_resolutions)

  expect_identical(baseline_raw$route, "R0")
  expect_identical(shadow_raw$route, "R2")
  expect_false(shadow_raw$inference_available)
  expect_identical(shadow_resolution$mode, "shadow_default_off")
  expect_identical(shadow_resolution$strict_route, "R2")
  expect_identical(
    shadow_resolution$strict_failure_reason,
    "competing_uncertified_candidate"
  )
  expect_identical(shadow_resolution$prospective_route, "R0")
  expect_true(shadow_resolution$would_nonveto)
  expect_identical(shadow_resolution$would_route, "R0")
  expect_identical(
    shadow_resolution$numerical_resolution,
    "would_same_basin_review_nonveto"
  )
  expect_length(shadow_candidate, 1L)
  expect_identical(
    shadow_candidate[[1L]]$class,
    "same_basin_review_nonveto"
  )
  expect_true(shadow_candidate[[1L]]$probe$accepted)
  expect_identical(shadow_candidate[[1L]]$probe$state_evaluations, 1L)
  expect_equal(
    shadow_candidate[[1L]]$probe$proposed$theta_reduced_log,
    baseline_raw$ln_joint$state$theta_reduced_log,
    tolerance = 1e-10
  )
  expect_null(shadow_raw$ln_joint$covariance)
  ln_isolation_expect_withheld(shadow)
  expect_true(all(is.na(ln_isolation_robust_wald_ci(shadow$summary))))

  expect_equal(
    shadow_raw$ln_joint$state$full_log,
    baseline_raw$ln_joint$state$full_log,
    tolerance = 1e-10
  )
  expect_equal(
    shadow_raw$ln_joint$state$raw$value,
    baseline_raw$ln_joint$state$raw$value,
    tolerance = 1e-10
  )
  expect_equal(shadow_raw$loglik, baseline_raw$loglik, tolerance = 1e-10)
  coefficient_columns <- grep("^logFC_", names(baseline$summary), value = TRUE)
  expect_equal(
    shadow$summary[coefficient_columns],
    baseline$summary[coefficient_columns],
    tolerance = 1e-10
  )

  expect_identical(enabled_raw$route, "R0")
  expect_true(enabled_raw$inference_available)
  expect_identical(
    enabled_resolution$mode,
    "calibration_candidate_enabled"
  )
  expect_identical(enabled_resolution$strict_route, "R2")
  expect_identical(enabled_resolution$prospective_route, "R0")
  expect_true(enabled_resolution$would_nonveto)
  expect_true(enabled_resolution$selected_root_certified)
  expect_identical(
    enabled_resolution$strict_inference_route,
    "R2"
  )
  expect_identical(
    enabled_resolution$enabled_inference_route,
    "R0"
  )
  expect_true(enabled_resolution$covariance_assessed_from_selected_root)
  expect_true(enabled_resolution$covariance_assessment_success)
  status <- enabled$provenance$per_gene
  expect_identical(enabled$provenance$schema_version, 3L)
  expect_identical(
    enabled$provenance$candidate_resolution_policy,
    "joint_ln_candidate_resolution_calibration_v1"
  )
  expect_identical(
    enabled$provenance$candidate_resolution_mode,
    "calibration_candidate_enabled"
  )
  expect_identical(status$estimator, "joint_ln")
  expect_true(status$certified_root)
  expect_identical(status$inference_status, "available")
  expect_identical(status$strict_inference_status, "withheld")
  expect_identical(status$enabled_inference_status, "available")
  expect_identical(
    status$resolution_status,
    "same_basin_review_resolved"
  )
  expect_identical(status$inference_source, "certified_selected_root")
  expect_true(status$newly_enabled)
  expect_identical(
    status$availability_transition,
    "strict_withheld_to_available"
  )
  expect_true(is.na(status$resolution_veto_classes))
  expect_identical(status$route, status$route_code)
  expect_identical(status$strict_route_code, "R2")
  expect_identical(status$enabled_route_code, "R0")
  expect_identical(
    status$resolution_policy_version,
    "joint_ln_candidate_resolution_calibration_v1"
  )
  expect_identical(
    status$resolution_mode,
    "calibration_candidate_enabled"
  )
  expect_identical(enabled$availability$availability, "inference")
  expect_identical(enabled$availability$inference_status, "available")
  expect_identical(
    enabled$availability$resolution_status,
    "same_basin_review_resolved"
  )
  expect_true(enabled$availability$newly_enabled)
  expect_identical(
    enabled$availability$status_message,
    "Certified inference \u2014 numerical review resolved"
  )
  expect_identical(
    enabled$availability$technical_code,
    "joint_ln_certified_after_same_basin_review"
  )
  injected_record <- attr(enabled, "synthetic_injected_candidate")
  expect_true(is.list(injected_record))
  expect_identical(injected_record$label, "synthetic_same_basin_shadow")
  expect_false(isTRUE(all.equal(
    injected_record$theta_full_log,
    enabled_raw$ln_joint$state$full_log,
    tolerance = 0
  )))
  expect_equal(
    enabled_raw$ln_joint$state$full_log,
    baseline_raw$ln_joint$state$full_log,
    tolerance = 1e-10
  )
  expect_equal(
    enabled_raw$ln_joint$state$raw$value,
    baseline_raw$ln_joint$state$raw$value,
    tolerance = 1e-10
  )
  expect_equal(enabled_raw$loglik, baseline_raw$loglik, tolerance = 1e-10)
  expect_equal(
    enabled_raw$ln_joint$covariance$model_beta,
    baseline_raw$ln_joint$covariance$model_beta,
    tolerance = 1e-10
  )
  expect_equal(
    enabled_raw$ln_joint$covariance$robust_beta,
    baseline_raw$ln_joint$covariance$robust_beta,
    tolerance = 1e-10
  )
  expect_equal(
    enabled_raw$ln_joint$covariance,
    baseline_raw$ln_joint$covariance,
    tolerance = 1e-10
  )
  expect_equal(enabled_raw$naive_cov, baseline_raw$naive_cov, tolerance = 1e-10)
  expect_equal(enabled_raw$robust_cov, baseline_raw$robust_cov, tolerance = 1e-10)
  expect_equal(
    ln_isolation_robust_wald_ci(enabled$summary),
    ln_isolation_robust_wald_ci(baseline$summary),
    tolerance = 1e-10
  )
  inference_columns <- ln_isolation_public_inference_columns(baseline$summary)
  expect_equal(
    enabled$summary[inference_columns],
    baseline$summary[inference_columns],
    tolerance = 1e-10
  )
  restored <- unserialize(serialize(enabled, NULL, version = 3L))
  expect_identical(restored, enabled)
})

test_that("same-basin resolution remains withheld when selected-root covariance fails", {
  failed <- ln_isolation_fit_with_candidate(
    "same_basin_shadow",
    covariance_failure = TRUE
  )
  raw <- failed$diagnostics$raw[[1L]]
  resolution <- raw$ln_joint$resolution
  status <- failed$provenance$per_gene

  expect_identical(resolution$prospective_route, "R0")
  expect_identical(
    resolution$numerical_resolution,
    "would_same_basin_review_nonveto"
  )
  expect_identical(resolution$enabled_inference_route, "R2")
  expect_false(resolution$covariance_assessment_success)
  expect_identical(raw$route, "R2")
  expect_false(raw$inference_available)
  expect_false(raw$ln_joint$covariance$success)
  expect_identical(
    raw$ln_joint$covariance$failure_reason,
    "synthetic_covariance_failure"
  )
  expect_identical(status$resolution_status, "same_basin_review_resolved")
  expect_identical(status$enabled_inference_status, "withheld")
  expect_identical(status$inference_status, "withheld")
  expect_identical(status$inference_source, "none")
  expect_false(status$newly_enabled)
  expect_identical(failed$availability$availability, "unavailable")
  ln_isolation_expect_withheld(failed)
})

test_that("strict-only comparator covariance never leaks through an enabled veto", {
  original_resolver <- getFromNamespace(
    "ln_joint_resolve_candidates",
    "scRoPE"
  )
  resolver_wrapper <- function(
      candidates,
      resolution_thresholds = ln_joint_resolution_thresholds(),
      basin_probe = NULL,
      enable_prospective = FALSE) {
    stopifnot(isTRUE(enable_prospective))
    result <- original_resolver(
      candidates,
      resolution_thresholds = resolution_thresholds,
      basin_probe = basin_probe,
      enable_prospective = enable_prospective
    )
    stopifnot(identical(result$resolution$strict_route, "R0"))
    result$route <- "R2"
    result$failure_reason <- "competitive_hard_stop_candidate"
    result$resolution$prospective_route <- "R2"
    result$resolution$prospective_failure_reason <-
      "competitive_hard_stop_candidate"
    result$resolution$would_nonveto <- FALSE
    result$resolution$would_route <- "R2"
    result$resolution$numerical_resolution <-
      "prospective_veto_competitive_hard_stop"
    result$resolution$global_veto_classes <- "competitive_hard_stop"
    result
  }
  local_mocked_bindings(
    ln_joint_resolve_candidates = resolver_wrapper,
    .package = "scRoPE"
  )

  fit <- suppressWarnings(ln_joint_fit_fixture()$fit)
  raw <- fit$diagnostics$raw[[1L]]
  status <- fit$provenance$per_gene

  expect_identical(raw$route, "R2")
  expect_false(raw$inference_available)
  expect_null(raw$ln_joint$covariance)
  expect_true(raw$ln_joint$strict_comparator_covariance$success)
  expect_identical(
    raw$ln_joint$resolution$strict_inference_route,
    "R0"
  )
  expect_identical(
    raw$ln_joint$resolution$enabled_inference_route,
    "R2"
  )
  expect_identical(status$strict_inference_status, "available")
  expect_identical(status$enabled_inference_status, "withheld")
  expect_identical(status$inference_status, "withheld")
  expect_identical(
    status$availability_transition,
    "strict_available_to_withheld"
  )
  expect_identical(status$inference_source, "none")
  expect_false(status$newly_enabled)
  ln_isolation_expect_withheld(fit)
})

test_that("enabled same-basin Wald leaves Score Pace and LRT pending", {
  fit <- ln_isolation_fit_with_candidate(
    "same_basin_shadow",
    additional_tests = c("score", "lrt")
  )
  status <- fit$provenance$per_gene
  capabilities <- fit$provenance$test_capabilities

  expect_identical(status$inference_status, "available")
  model_and_robust <- grep(
    "^se_|^se_robust_|^p_naive_|^p_robust_|^wald_stat_",
    names(fit$summary),
    value = TRUE
  )
  expect_true(length(model_and_robust) > 0L)
  expect_true(all(is.finite(as.matrix(fit$summary[model_and_robust]))))
  expect_identical(
    capabilities$status[capabilities$test %in% c("score", "lrt")],
    rep("joint_ln_calibration_pending", 2L)
  )
  expect_identical(
    fit$provenance$pace_status,
    "joint_ln_calibration_pending"
  )
  expect_true(all(!fit$lrt_details$success))
  expect_true(all(!fit$lrt_details$score_success))
  expect_true(all(
    fit$lrt_details$failure_reason == "joint_ln_calibration_pending"
  ))
  expect_true(all(
    fit$lrt_details$score_failure == "joint_ln_calibration_pending"
  ))
  expect_true(all(is.na(fit$lrt_details$stat)))
  expect_true(all(is.na(fit$lrt_details$p_value)))
  constrained <- Filter(
    Negate(is.null),
    lapply(fit$diagnostics$raw[[1L]]$lrt, `[[`, "constrained")
  )
  expect_true(length(constrained) > 0L)
  expect_true(all(vapply(constrained, function(item) {
    identical(item$ln_joint$resolution$mode, "shadow_default_off")
  }, logical(1))))
})

test_that("competitive boundary distinct and hard-stop candidates retain withholding", {
  cases <- c(
    competitive_boundary = "competitive_boundary",
    competitive_distinct_basin = "competitive_distinct_basin",
    competitive_hard_stop = "competitive_hard_stop"
  )
  for (kind in names(cases)) {
    fit <- ln_isolation_fit_with_candidate(kind)
    raw <- fit$diagnostics$raw[[1L]]
    resolution <- raw$ln_joint$resolution
    injected <- Filter(function(record) {
      identical(record$label, paste0("synthetic_", kind))
    }, resolution$candidate_resolutions)

    expect_identical(raw$route, "R2", info = kind)
    expect_false(raw$inference_available, info = kind)
    expect_identical(
      resolution$mode,
      "calibration_candidate_enabled",
      info = kind
    )
    expect_identical(resolution$strict_route, "R2", info = kind)
    expect_identical(resolution$prospective_route, "R2", info = kind)
    expect_false(resolution$would_nonveto, info = kind)
    expect_true(cases[[kind]] %in% resolution$global_veto_classes, info = kind)
    expect_length(injected, 1L)
    expect_identical(
      injected[[1L]]$class,
      cases[[kind]],
      info = kind
    )
    expect_true(injected[[1L]]$veto, info = kind)
    expect_null(injected[[1L]]$probe, info = kind)
    expect_null(raw$ln_joint$covariance, info = kind)
    status <- fit$provenance$per_gene
    expect_identical(status$inference_status, "withheld", info = kind)
    expect_identical(
      status$enabled_inference_status,
      "withheld",
      info = kind
    )
    expect_identical(status$inference_source, "none", info = kind)
    expect_false(status$newly_enabled, info = kind)
    expected_status <- c(
      competitive_boundary = "competitive_boundary",
      competitive_distinct_basin = "competitive_distinct_basin",
      competitive_hard_stop = "hard_stop"
    )
    expect_identical(
      status$resolution_status,
      expected_status[[kind]],
      info = kind
    )
    expect_match(
      status$resolution_veto_classes,
      cases[[kind]],
      fixed = TRUE,
      info = kind
    )
    expect_identical(
      resolution$evaluation_counts$probe_state_evaluations,
      0L,
      info = kind
    )
    ln_isolation_expect_withheld(fit)
  }
})

test_that("an unresolved bounded review retains inference withholding", {
  fit <- ln_isolation_fit_with_candidate("same_basin_unresolved")
  raw <- fit$diagnostics$raw[[1L]]
  resolution <- raw$ln_joint$resolution
  injected <- Filter(function(record) {
    identical(record$label, "synthetic_same_basin_unresolved")
  }, resolution$candidate_resolutions)

  expect_identical(raw$route, "R2")
  expect_false(raw$inference_available)
  expect_identical(resolution$prospective_route, "R2")
  expect_false(resolution$would_nonveto)
  expect_true(
    "same_basin_review_unresolved" %in% resolution$global_veto_classes
  )
  expect_length(injected, 1L)
  expect_identical(injected[[1L]]$class, "same_basin_review_unresolved")
  expect_true(injected[[1L]]$veto)
  expect_identical(
    fit$provenance$per_gene$resolution_status,
    "review_unresolved"
  )
  expect_identical(fit$provenance$per_gene$inference_status, "withheld")
  expect_identical(fit$provenance$per_gene$inference_source, "none")
  expect_false(fit$provenance$per_gene$newly_enabled)
  ln_isolation_expect_withheld(fit)
})

test_that("candidate resolution and probe records deep-round-trip and replay", {
  selected <- ln_isolation_synthetic_candidate("certified_a")
  agreeing <- ln_isolation_synthetic_candidate(
    "certified_b",
    objective = 10 + 1e-10,
    theta_full_log = selected$theta_full_log + c(1e-10, 0, 0, 0)
  )
  shadow <- ln_isolation_injected_candidate(selected, "same_basin_shadow")
  candidates <- list(selected, agreeing, shadow)
  candidates_before <- serialize(candidates, connection = NULL, version = 3L)
  resolve <- function() {
    ln_joint_resolve_candidates(
      candidates,
      basin_probe = ln_isolation_accept_probe,
      enable_prospective = FALSE
    )
  }

  first <- resolve()
  second <- resolve()
  expect_identical(
    serialize(candidates, connection = NULL, version = 3L),
    candidates_before
  )
  expect_identical(first, second)
  expect_true(first$resolution$would_nonveto)
  expect_identical(first$resolution$would_route, "R0")
  expect_identical(
    first$resolution$evaluation_counts$probe_state_evaluations,
    1L
  )

  payload <- list(
    resolution = first$resolution,
    immutable_candidates = candidates,
    candidate_records = ln_joint_candidate_records(
      candidates,
      first$resolution$candidate_resolutions
    ),
    candidate_summary = first$candidates
  )
  encoded <- serialize(payload, connection = NULL, version = 3L)
  restored <- unserialize(encoded)
  expect_identical(restored, payload)
  expect_identical(
    restored$candidate_records[[3L]]$resolution$class,
    "same_basin_review_nonveto"
  )
  expect_true(restored$candidate_records[[3L]]$resolution$probe$accepted)
})
