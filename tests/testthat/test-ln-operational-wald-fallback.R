ln_operational_test_ctx <- function(O = diag(2L)) {
  list(
    nb = nrow(O),
    min = c(0.1, 0.1),
    max = c(10, 10),
    coef_internal_to_public = O
  )
}

ln_operational_test_legacy_fit <- function(
    beta = c(0.4, -0.2),
    subject = 1,
    cell = 2,
    convergence = 1L,
    algorithm = 1L,
    rescue = "none") {
  nb <- length(beta)
  vec <- rep(NA_real_, 3L * nb + 5L)
  vec[seq_len(nb)] <- beta
  vec[[nb + 1L]] <- subject
  vec[[nb + 2L]] <- cell
  vec[[3L * nb + 3L]] <- convergence
  vec[[3L * nb + 4L]] <- algorithm
  vec[[3L * nb + 5L]] <- -10
  list(
    vec = vec,
    godambe = NULL,
    debug = list(rescue = rescue)
  )
}

ln_operational_test_gene_result <- function(
    route = "R2",
    failure_reason = "dispersion_boundary",
    checks = c(
      donor_additivity = TRUE,
      root = TRUE,
      newton = TRUE,
      curvature = TRUE,
      condition = TRUE,
      interior = FALSE,
      constraint = TRUE
    )) {
  theta <- c(0.2, -0.1, log(0.5), log(2))
  objective <- 12.5
  list(
    route = route,
    failure_reason = failure_reason,
    loglik = -objective,
    ln_joint = list(
      state = list(full_log = theta),
      candidate_records = list(list(
        label = "selected",
        finite_endpoint = TRUE,
        certified = all(checks),
        objective = objective,
        theta_full_log = theta,
        checks = checks
      ))
    )
  )
}

ln_operational_all_field_names <- function(x) {
  current <- names(x)
  if (is.null(current)) {
    current <- character(0)
  }
  if (!is.list(x)) {
    return(current)
  }
  c(
    current,
    unlist(lapply(x, ln_operational_all_field_names), use.names = FALSE)
  )
}

ln_operational_expect_no_inference_payload <- function(x) {
  field_names <- unique(ln_operational_all_field_names(x))
  forbidden <- paste(
    c(
      "covariance",
      "(^|_)se($|_)",
      "standard_error",
      "wald",
      "(^|_)p($|_)",
      "p_value",
      "pvalue",
      "p_robust",
      "fdr",
      "calibration_eligible",
      "inference_available",
      "inference_success",
      "inferential_success",
      "robust_available"
    ),
    collapse = "|"
  )
  hits <- field_names[grepl(forbidden, field_names, ignore.case = TRUE)]
  expect_identical(hits, character(0))
}

ln_operational_normalise_objective_instance_ids <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  object_names <- names(x)
  for (index in seq_along(x)) {
    value <- x[[index]]
    field <- if (is.null(object_names)) "" else object_names[[index]]
    if (identical(field, "objective_instance") && is.list(value) &&
        identical(
          value$scheme,
          "process_local_exact_joint_ln_objective_context_v1"
        )) {
      value$id <- paste0("<run-local-id-gene-", value$gene_index, ">")
    }
    x[index] <- list(ln_operational_normalise_objective_instance_ids(value))
  }
  x
}

ln_operational_expect_within_fit_objective_identity <- function(fit) {
  raw <- fit$diagnostics$raw
  if (is.null(raw)) {
    return(invisible(NULL))
  }
  tokens <- vector("list", length(raw))
  for (index in seq_along(raw)) {
    gene_fit <- raw[[index]]
    token <- gene_fit$ln_joint$objective_instance
    expect_true(is.list(token))
    expect_identical(
      token$scheme,
      "process_local_exact_joint_ln_objective_context_v1"
    )
    expect_identical(token$gene_index, as.integer(gene_fit$debug$gene_index))
    if (!is.null(gene_fit$ln_joint$state)) {
      expect_identical(
        gene_fit$ln_joint$state$objective_instance,
        token
      )
    }
    if (!is.null(gene_fit$lrt)) {
      for (test_result in gene_fit$lrt) {
        constrained <- test_result$constrained
        if (is.null(constrained) || is.null(constrained$ln_joint)) {
          next
        }
        expect_identical(constrained$ln_joint$objective_instance, token)
        if (!is.null(constrained$ln_joint$state)) {
          expect_identical(
            constrained$ln_joint$state$objective_instance,
            token
          )
        }
        ratio <- constrained$ln_joint$pace_ratio_diagnostic
        if (!is.null(ratio) && !is.null(ratio$objective_instance)) {
          expect_identical(ratio$objective_instance, token)
          expect_identical(
            ratio$endpoints$unrestricted$objective_instance,
            token
          )
          expect_identical(
            ratio$endpoints$constrained$objective_instance,
            token
          )
        }
      }
    }
    tokens[[index]] <- token
  }
  if (length(tokens) > 1L) {
    token_ids <- vapply(tokens, `[[`, character(1), "id")
    expect_identical(length(unique(token_ids)), length(token_ids))
  }
  invisible(NULL)
}

ln_operational_expect_cross_run_identical <- function(left, right) {
  ln_operational_expect_within_fit_objective_identity(left)
  ln_operational_expect_within_fit_objective_identity(right)
  expect_identical(
    ln_operational_normalise_objective_instance_ids(left),
    ln_operational_normalise_objective_instance_ids(right)
  )
}

test_that("operational fallback API is opt-in and none preserves primary", {
  fx <- ln_joint_test_fixture()
  joint_args <- list(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln"
  )

  set.seed(20260721)
  default_none <- expect_no_warning(do.call(scrope, joint_args))
  set.seed(20260721)
  explicit_none <- expect_no_warning(do.call(
    scrope,
    c(joint_args, list(ln_operational_fallback = "none"))
  ))
  ln_operational_expect_cross_run_identical(default_none, explicit_none)
  expect_false("operational_fallback" %in% names(default_none))

  set.seed(20260721)
  enabled <- expect_no_warning(do.call(
    scrope,
    c(joint_args, list(ln_operational_fallback = "legacy_post_pml"))
  ))
  expect_identical(enabled$provenance$per_gene$route, "R0")
  expect_false(enabled$operational_fallback$inference_supported)
  ln_operational_expect_no_inference_payload(enabled$operational_fallback)
  primary_enabled <- enabled
  primary_enabled$operational_fallback <- NULL
  ln_operational_expect_cross_run_identical(primary_enabled, default_none)

  expect_error(
    do.call(scrope, c(
      joint_args,
      list(ln_operational_fallback = "unsupported")
    )),
    "should be one of"
  )

  legacy_args <- joint_args
  legacy_args$ln_estimator <- "legacy_post_pml"
  expect_error(
    do.call(scrope, c(
      legacy_args,
      list(ln_operational_fallback = "legacy_post_pml")
    )),
    "requires ln_estimator='joint_ln'"
  )

  switching_args <- joint_args
  switching_args$allow_per_gene_switch <- TRUE
  expect_error(
    do.call(scrope, c(
      switching_args,
      list(ln_operational_fallback = "legacy_post_pml")
    )),
    "does not support allow_per_gene_switch=TRUE"
  )

  pgmm_args <- joint_args
  pgmm_args$model <- "PGMM"
  pgmm_args$ln_estimator <- NULL
  expect_error(
    do.call(scrope, c(
      pgmm_args,
      list(ln_operational_fallback = "legacy_post_pml")
    )),
    "applies only to model='LN'"
  )
})

test_that("fallback trigger preserves primary R2 boundary selection", {
  eligible <- ln_operational_fallback_eligibility(
    ln_operational_test_gene_result()
  )
  expect_true(eligible$complete_checks)
  expect_true(eligible$selected_state_match)
  expect_true(eligible$selected_objective_match)
  expect_true(eligible$boundary_only)
  expect_true(eligible$eligible)
  expect_identical(
    eligible$eligibility_rule,
    "primary_R2_and_final_failure_reason_dispersion_boundary"
  )
  expect_false(eligible$selected_check_audit_used_for_eligibility)

  boundary_and_root <- ln_operational_test_gene_result()
  boundary_and_root$ln_joint$candidate_records[[1L]]$checks[["root"]] <- FALSE
  boundary_and_root$ln_joint$candidate_records[[1L]]$checks[["newton"]] <- FALSE
  boundary_and_root_eligibility <- ln_operational_fallback_eligibility(
    boundary_and_root
  )
  expect_identical(
    boundary_and_root_eligibility$primary_failure_reason,
    "dispersion_boundary"
  )
  expect_false(boundary_and_root_eligibility$boundary_only)
  expect_true(boundary_and_root_eligibility$eligible)
  expect_true(boundary_and_root_eligibility$joint_root_or_newton_failure)
  expect_identical(
    boundary_and_root_eligibility$joint_failed_checks,
    c("root", "newton", "interior")
  )

  for (reason in c(
    "root_residual_hard_stop",
    "donor_additivity_hard_stop"
  )) {
    expect_false(ln_operational_fallback_eligibility(
      ln_operational_test_gene_result(failure_reason = reason)
    )$eligible)
  }

  incomplete <- ln_operational_test_gene_result()
  incomplete$ln_joint$candidate_records[[1L]]$checks <-
    incomplete$ln_joint$candidate_records[[1L]]$checks[-1L]
  incomplete_eligibility <- ln_operational_fallback_eligibility(incomplete)
  expect_false(incomplete_eligibility$complete_checks)
  expect_true(incomplete_eligibility$eligible)

  state_mismatch <- ln_operational_test_gene_result()
  state_mismatch$ln_joint$state$full_log[[1L]] <- 999
  state_mismatch_eligibility <- ln_operational_fallback_eligibility(
    state_mismatch
  )
  expect_false(state_mismatch_eligibility$selected_state_match)
  expect_true(state_mismatch_eligibility$eligible)

  objective_mismatch <- ln_operational_test_gene_result()
  objective_mismatch$loglik <- -99
  objective_mismatch_eligibility <- ln_operational_fallback_eligibility(
    objective_mismatch
  )
  expect_false(objective_mismatch_eligibility$selected_objective_match)
  expect_true(objective_mismatch_eligibility$eligible)

  missing_candidate <- ln_operational_test_gene_result()
  missing_candidate$ln_joint$candidate_records <- list()
  missing_candidate_eligibility <- ln_operational_fallback_eligibility(
    missing_candidate
  )
  expect_false(missing_candidate_eligibility$complete_checks)
  expect_true(missing_candidate_eligibility$eligible)

  expect_false(ln_operational_fallback_eligibility(
    ln_operational_test_gene_result(route = "R0")
  )$eligible)
})

test_that("legacy assessment returns coefficient and dispersion without covariance", {
  fx <- ln_joint_test_fixture()
  prep <- ln_joint_prepare_fixture(fx)
  prep$ctx$ln_estimator <- "legacy_post_pml"
  prep$ctx$output_re <- FALSE
  prep$ctx$debug_env <- new.env(parent = emptyenv())
  legacy_gene <- fit_gene_unconstrained(1L, prep$posv, prep$ctx)

  legacy_gene$godambe <- NULL
  robust_slots <- 2L * prep$ctx$nb + 2L + seq_len(prep$ctx$nb)
  legacy_gene$vec[robust_slots] <- NA_real_
  assessment <- ln_operational_assess_legacy_fit(legacy_gene, prep$ctx)
  expect_true(assessment$coefficient_available)
  expect_false(assessment$inference_supported)
  expect_identical(
    assessment$inference_status,
    "unsupported_operational_coefficient_only"
  )
  ln_operational_expect_no_inference_payload(assessment)

  public_legacy <- suppressWarnings(scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    keep_diagnostics = TRUE,
    ln_estimator = "legacy_post_pml"
  ))
  terms <- colnames(fx$pred)
  expect_equal(
    assessment$beta_public,
    as.numeric(public_legacy$summary[1L, paste0("logFC_", terms)]),
    tolerance = 1e-10
  )
  expect_equal(
    c(assessment$subject, assessment$cell),
    as.numeric(public_legacy$overdispersion[1L, ]),
    tolerance = 1e-10
  )
})

test_that("legacy assessment transforms coefficients and labels route class", {
  O <- matrix(c(1, 0.5, -0.25, 2), 2L, 2L, byrow = TRUE)
  ctx <- ln_operational_test_ctx(O)
  direct <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(),
    ctx
  )
  rescued <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(rescue = "pql_both"),
    ctx
  )
  boundary <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(subject = ctx$min[[1L]]),
    ctx
  )
  rescued_boundary <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(
      subject = ctx$min[[1L]],
      rescue = "pql_both"
    ),
    ctx
  )

  expect_equal(
    direct$beta_public,
    as.numeric(O %*% c(0.4, -0.2)),
    tolerance = 1e-12
  )
  expect_identical(direct$route_class, "uncertified")
  expect_identical(rescued$route_class, "rescued")
  expect_identical(boundary$route_class, "boundary")
  expect_identical(rescued_boundary$route_class, "rescued_boundary")
  expect_true(all(vapply(
    list(direct, rescued, boundary, rescued_boundary),
    function(x) isTRUE(x$coefficient_available) &&
      identical(x$inference_supported, FALSE),
    logical(1)
  )))
  lapply(
    list(direct, rescued, boundary, rescued_boundary),
    ln_operational_expect_no_inference_payload
  )
})

test_that("legacy operational gates use a coefficient-only failure vocabulary", {
  ctx <- ln_operational_test_ctx()
  invalid_transform_ctx <- ctx
  invalid_transform_ctx$coef_internal_to_public[[1L]] <- Inf
  cases <- list(
    legacy_fit_error = list(
      fit = list(vec = numeric(0)),
      ctx = ctx
    ),
    legacy_nonconvergence = list(
      fit = ln_operational_test_legacy_fit(convergence = -20L),
      ctx = ctx
    ),
    legacy_gradient_convergence_not_accepted = list(
      fit = ln_operational_test_legacy_fit(convergence = -10L),
      ctx = ctx
    ),
    nonfinite_coefficient = list(
      fit = ln_operational_test_legacy_fit(beta = c(Inf, -0.2)),
      ctx = ctx
    ),
    invalid_dispersion = list(
      fit = ln_operational_test_legacy_fit(cell = 0),
      ctx = ctx
    ),
    invalid_coefficient_transform = list(
      fit = ln_operational_test_legacy_fit(),
      ctx = invalid_transform_ctx
    )
  )

  observed <- vapply(cases, function(case) {
    assessment <- ln_operational_assess_legacy_fit(case$fit, case$ctx)
    expect_false(assessment$coefficient_available)
    expect_false(assessment$inference_supported)
    expect_null(assessment$beta_public)
    ln_operational_expect_no_inference_payload(assessment)
    assessment$failure_reason
  }, character(1))
  expect_identical(unname(observed), names(cases))
  expect_setequal(
    unname(unlist(failure_codes$operational_fallback, use.names = FALSE)),
    names(cases)
  )
  expect_error(
    normalise_operational_fallback_failure("other_failure"),
    "Unknown operational fallback failure reason"
  )
})

test_that("sidecar exposes coefficients and diagnostics only for every route", {
  ctx <- ln_operational_test_ctx()
  failed <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(convergence = -20L),
    ctx
  )
  direct <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(),
    ctx
  )
  rescued <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(rescue = "pql_both"),
    ctx
  )
  boundary <- ln_operational_assess_legacy_fit(
    ln_operational_test_legacy_fit(subject = ctx$min[[1L]]),
    ctx
  )

  primary_results <- c(
    list(list(route = "R0"), list(route = "R2")),
    rep(list(list(route = "R2")), 4L)
  )
  eligibility <- list(
    list(
      eligible = FALSE,
      primary_route = "R0",
      primary_failure_reason = NA_character_
    ),
    list(
      eligible = FALSE,
      primary_route = "R2",
      primary_failure_reason = "root_residual_hard_stop"
    ),
    list(
      eligible = TRUE,
      primary_route = "R2",
      primary_failure_reason = "dispersion_boundary"
    ),
    list(
      eligible = TRUE,
      primary_route = "R2",
      primary_failure_reason = "dispersion_boundary"
    ),
    list(
      eligible = TRUE,
      primary_route = "R2",
      primary_failure_reason = "dispersion_boundary"
    ),
    list(
      eligible = TRUE,
      primary_route = "R2",
      primary_failure_reason = "dispersion_boundary"
    )
  )
  fallback_results <- list(
    ln_operational_empty_assessment(),
    ln_operational_empty_assessment(),
    failed,
    direct,
    rescued,
    boundary
  )
  gene_id <- c(8L, 3L, 5L, 1L, 4L, 2L)
  gene <- c("r0", "ineligible", "failed", "direct", "rescued", "boundary")
  terms <- c("second", "first")
  sidecar <- ln_operational_build_sidecar(
    primary_results = primary_results,
    eligibility = eligibility,
    fallback_results = fallback_results,
    gid = gene_id,
    gene_names = gene,
    term_names = terms,
    ctx = ctx
  )

  expect_identical(names(sidecar), c(
    "schema_version", "method", "role", "selection_rule",
    "inference_supported", "inference_status", "per_gene",
    "coefficients", "overdispersion", "diagnostics"
  ))
  expect_identical(sidecar$schema_version, 2L)
  expect_identical(sidecar$method, "legacy_post_pml")
  expect_identical(
    sidecar$role,
    "operational_coefficient_dispersion_diagnostics_only"
  )
  expect_identical(
    sidecar$selection_rule,
    "primary_R2_and_failure_reason_dispersion_boundary"
  )
  expect_false(sidecar$inference_supported)
  expect_identical(
    sidecar$inference_status,
    "unsupported_for_all_fallback_routes"
  )

  expect_identical(names(sidecar$per_gene), c(
    "gene_id", "gene", "primary_route", "primary_failure_reason",
    "fallback_attempted", "fallback_estimator", "fallback_role",
    "fallback_route_class", "fallback_convergence", "fallback_algorithm",
    "fallback_start_provenance", "fallback_boundary",
    "fallback_boundary_component", "fallback_distance_subject_lower",
    "fallback_distance_subject_upper", "fallback_distance_cell_lower",
    "fallback_distance_cell_upper", "coefficient_available",
    "coefficient_source", "fallback_failure_reason", "fallback_status",
    "inference_supported", "inference_status"
  ))
  expect_identical(sidecar$per_gene$gene_id, gene_id)
  expect_identical(sidecar$per_gene$gene, gene)
  expect_identical(
    sidecar$per_gene$fallback_attempted,
    c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
  )
  expect_identical(
    sidecar$per_gene$coefficient_available,
    c(NA, NA, FALSE, TRUE, TRUE, TRUE)
  )
  expect_identical(
    sidecar$per_gene$fallback_status,
    c(
      "not_eligible", "not_eligible", "attempted_failed",
      "coefficient_available", "coefficient_available",
      "coefficient_available"
    )
  )
  expect_identical(
    sidecar$per_gene$fallback_route_class,
    c(NA_character_, NA_character_, "uncertified", "uncertified",
      "rescued", "boundary")
  )
  expect_identical(
    sidecar$per_gene$coefficient_source,
    c(NA_character_, NA_character_, NA_character_,
      rep("legacy_post_pml", 3L))
  )
  expect_identical(
    sidecar$per_gene$fallback_failure_reason,
    c(NA_character_, NA_character_, "legacy_nonconvergence",
      rep(NA_character_, 3L))
  )
  expect_true(all(!sidecar$per_gene$inference_supported))
  expect_true(all(
    sidecar$per_gene$inference_status ==
      "unsupported_operational_coefficient_only"
  ))

  expect_identical(names(sidecar$coefficients), c(
    "gene_id", "gene", "term", "estimate", "coefficient_available",
    "coefficient_source", "fallback_route_class", "fallback_status",
    "fallback_failure_reason", "inference_supported", "inference_status"
  ))
  expect_identical(sidecar$coefficients$gene_id, rep(gene_id, each = 2L))
  expect_identical(sidecar$coefficients$gene, rep(gene, each = 2L))
  expect_identical(sidecar$coefficients$term, rep(terms, times = 6L))
  expect_true(all(is.na(sidecar$coefficients$estimate[1:6])))
  expect_true(all(is.finite(sidecar$coefficients$estimate[7:12])))
  expect_true(all(!sidecar$coefficients$inference_supported))

  expect_identical(names(sidecar$overdispersion), c(
    "gene_id", "gene", "Subject", "Cell", "coefficient_available",
    "fallback_status", "fallback_failure_reason", "inference_supported",
    "inference_status"
  ))
  expect_true(all(is.na(sidecar$overdispersion$Subject[1:3])))
  expect_true(all(is.na(sidecar$overdispersion$Cell[1:3])))
  expect_equal(
    sidecar$overdispersion$Subject[4:6],
    c(direct$subject, rescued$subject, boundary$subject)
  )
  expect_equal(
    sidecar$overdispersion$Cell[4:6],
    c(direct$cell, rescued$cell, boundary$cell)
  )
  expect_true(all(!sidecar$overdispersion$inference_supported))
  expect_identical(
    sidecar$diagnostics$policy,
    "coefficient_dispersion_diagnostics_only"
  )
  ln_operational_expect_no_inference_payload(sidecar)
})

test_that("aggregate warning tolerates missing primary route metadata", {
  sidecar <- list(per_gene = data.frame(
    primary_route = c(NA_character_, "R2"),
    fallback_status = c("not_eligible", "attempted_failed"),
    fallback_boundary = c(NA, NA),
    stringsAsFactors = FALSE
  ))

  expect_warning(
    ln_operational_warn(sidecar),
    "Joint-LN inference was unavailable for 1 gene"
  )
})

ln_operational_force_failed_joint_primary <- function(original_fitter) {
  force(original_fitter)
  function(gene_index, posv, ctx, extra_starts = NULL) {
    result <- original_fitter(
      gene_index = gene_index,
      posv = posv,
      ctx = ctx,
      extra_starts = extra_starts
    )
    records <- result$ln_joint$candidate_records
    selected_index <- which(vapply(records, function(x) {
      !is.null(x$theta_full_log) &&
        identical(
          as.numeric(x$theta_full_log),
          as.numeric(result$ln_joint$state$full_log)
        ) &&
        identical(as.numeric(x$objective), -as.numeric(result$loglik))
    }, logical(1)))[[1L]]
    for (index in seq_along(records)) {
      records[[index]]$certified <- FALSE
      records[[index]]$finite_endpoint <- identical(index, selected_index)
    }
    records[[selected_index]]$checks[] <- TRUE
    records[[selected_index]]$checks[["root"]] <- FALSE
    records[[selected_index]]$checks[["newton"]] <- FALSE
    records[[selected_index]]$checks[["interior"]] <- FALSE
    records[[selected_index]]$failure <- "dispersion_boundary"
    result$ln_joint$candidate_records <- records
    result$route <- "R2"
    result$failure_reason <- "dispersion_boundary"
    result$inference_available <- FALSE
    result$naive_cov[,] <- NA_real_
    result$robust_cov[,] <- NA_real_
    nb <- ctx$nb
    result$vec[nb + 2L + seq_len(2L * nb)] <- NA_real_
    result$vec[[3L * nb + 3L]] <- -70L
    result$debug$conv_code <- -70L
    result$debug$route <- "R2"
    result$debug$failure_reason <- "dispersion_boundary"
    result
  }
}

ln_operational_capture_warnings <- function(expr) {
  messages <- character(0)
  value <- withCallingHandlers(
    force(expr),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = messages)
}

test_that("public fallback is primary-isolated and coefficient-only", {
  fx <- ln_joint_test_fixture()
  original_joint <- getFromNamespace(
    "fit_gene_joint_ln_unconstrained",
    "scRoPE"
  )
  local_mocked_bindings(
    fit_gene_joint_ln_unconstrained =
      ln_operational_force_failed_joint_primary(original_joint),
    .package = "scRoPE"
  )
  args <- list(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    output_re = TRUE,
    additional_tests = c("score", "lrt"),
    lrt_details = TRUE,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln"
  )
  set.seed(20260721)
  primary_capture <- ln_operational_capture_warnings(do.call(
    scrope,
    c(args, list(ln_operational_fallback = "none"))
  ))
  primary <- primary_capture$value
  expect_length(primary_capture$messages, 1L)
  expect_match(
    primary_capture$messages,
    "Joint-LN inference was unavailable for 1 gene"
  )
  expect_match(
    primary_capture$messages,
    "operational coefficients remain unavailable for 1"
  )
  set.seed(20260721)
  enabled_capture <- ln_operational_capture_warnings(do.call(
    scrope,
    c(args, list(ln_operational_fallback = "legacy_post_pml"))
  ))
  enabled <- enabled_capture$value
  expect_length(enabled_capture$messages, 1L)
  expect_match(
    enabled_capture$messages,
    "Joint-LN inference was unavailable for 1 gene"
  )
  expect_match(
    enabled_capture$messages,
    "operational coefficients are available for 1"
  )
  expect_match(
    enabled_capture$messages,
    "operational coefficients remain unavailable for 0"
  )
  expect_match(
    enabled_capture$messages,
    "coefficients, dispersions, and diagnostics only"
  )
  expect_match(enabled_capture$messages, "provides no covariance")

  primary_enabled <- enabled
  primary_enabled$operational_fallback <- NULL
  expect_identical(
    primary_enabled$availability$availability,
    "estimate_only"
  )
  expect_identical(primary$availability$availability, "unavailable")
  primary_enabled$availability <- NULL
  primary_without_availability <- primary
  primary_without_availability$availability <- NULL
  ln_operational_expect_cross_run_identical(
    primary_enabled,
    primary_without_availability
  )
  expect_identical(primary$provenance$per_gene$route, "R2")
  expect_identical(primary$convergence, -70)
  primary_inference_columns <- grep(
    "^se_|^p_naive_|^p_robust_|^wald_stat_",
    names(primary$summary),
    value = TRUE
  )
  expect_true(all(is.na(as.matrix(
    primary$summary[primary_inference_columns]
  ))))
  expect_true(all(!primary$lrt_details$success))
  expect_true(all(!primary$lrt_details$score_success))
  expect_true(all(
    primary$lrt_details$failure_reason == "joint_ln_calibration_pending"
  ))

  sidecar <- enabled$operational_fallback
  per_gene <- sidecar$per_gene
  expect_true(per_gene$fallback_attempted)
  expect_true(per_gene$coefficient_available)
  expect_identical(per_gene$fallback_status, "coefficient_available")
  expect_true(per_gene$fallback_route_class %in%
    c("uncertified", "rescued", "boundary", "rescued_boundary"))
  expect_identical(
    per_gene$fallback_estimator,
    "legacy_post_pml_operational_coefficient"
  )
  expect_identical(
    per_gene$fallback_role,
    "coefficient_dispersion_diagnostics_only"
  )
  expect_false(per_gene$inference_supported)
  expect_false(sidecar$inference_supported)
  expect_false("summary" %in% names(sidecar))
  ln_operational_expect_no_inference_payload(sidecar)

  eligibility <- sidecar$diagnostics$per_gene[[1L]]$eligibility
  expect_true(eligibility$eligible)
  expect_false(eligibility$primary_joint_candidate_checks[["root"]])
  expect_false(eligibility$primary_joint_candidate_checks[["newton"]])
  expect_false(eligibility$primary_joint_candidate_checks[["interior"]])
  expect_true(eligibility$joint_root_or_newton_failure)

  independent <- suppressWarnings(scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    output_re = FALSE,
    keep_diagnostics = TRUE,
    ln_estimator = "legacy_post_pml"
  ))
  terms <- colnames(fx$pred)
  expect_equal(
    sidecar$coefficients$estimate,
    as.numeric(independent$summary[paste0("logFC_", terms)]),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(sidecar$overdispersion[c("Subject", "Cell")]),
    as.numeric(independent$overdispersion[c("Subject", "Cell")]),
    tolerance = 1e-10
  )
})

test_that("legacy gate failure stays attempted without inference-shaped fields", {
  fx <- ln_joint_test_fixture()
  original_joint <- getFromNamespace(
    "fit_gene_joint_ln_unconstrained",
    "scRoPE"
  )
  failed_assessment <- ln_operational_empty_assessment(
    "legacy_gradient_convergence_not_accepted"
  )
  failed_assessment$convergence <- -10L
  failed_assessment$algorithm <- 1L
  failed_assessment$start_provenance <- "direct"
  failed_assessment$route_class <- "uncertified"
  local_mocked_bindings(
    fit_gene_joint_ln_unconstrained =
      ln_operational_force_failed_joint_primary(original_joint),
    ln_operational_attempt_gene = function(...) failed_assessment,
    .package = "scRoPE"
  )
  captured <- ln_operational_capture_warnings(scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln",
    ln_operational_fallback = "legacy_post_pml"
  ))
  expect_length(captured$messages, 1L)
  expect_match(captured$messages, "unavailable for 1 gene")
  expect_match(captured$messages, "coefficients are available for 0")
  expect_match(captured$messages, "coefficients remain unavailable for 1")
  per_gene <- captured$value$operational_fallback$per_gene
  expect_true(per_gene$fallback_attempted)
  expect_false(per_gene$coefficient_available)
  expect_identical(per_gene$fallback_status, "attempted_failed")
  expect_identical(per_gene$fallback_route_class, "uncertified")
  expect_identical(
    per_gene$fallback_failure_reason,
    "legacy_gradient_convergence_not_accepted"
  )
  coefficients <- captured$value$operational_fallback$coefficients
  expect_true(all(!coefficients$coefficient_available))
  expect_true(all(is.na(coefficients$estimate)))
  expect_true(all(!coefficients$inference_supported))
  ln_operational_expect_no_inference_payload(
    captured$value$operational_fallback
  )
})

test_that("aggregate warning reports operational coefficients only", {
  fx <- ln_joint_test_fixture()
  count <- rbind(fx$count, fx$count)
  rownames(count) <- c("gene1", "gene2")
  original_joint <- getFromNamespace(
    "fit_gene_joint_ln_unconstrained",
    "scRoPE"
  )
  local_mocked_bindings(
    fit_gene_joint_ln_unconstrained =
      ln_operational_force_failed_joint_primary(original_joint),
    .package = "scRoPE"
  )
  captured <- ln_operational_capture_warnings(scrope(
    count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln",
    ln_operational_fallback = "legacy_post_pml"
  ))
  fallback_messages <- captured$messages[startsWith(
    captured$messages,
    "Joint-LN inference was unavailable for"
  )]
  expect_length(fallback_messages, 1L)
  expect_match(fallback_messages, "unavailable for 2 genes")
  expect_match(fallback_messages, "coefficients are available for 2")
  expect_match(fallback_messages, "coefficients remain unavailable for 0")
  expect_match(fallback_messages, "no covariance, standard errors")
  expect_true(all(
    captured$value$operational_fallback$per_gene$fallback_status ==
      "coefficient_available"
  ))
  expect_true(all(
    !captured$value$operational_fallback$per_gene$inference_supported
  ))
})

test_that("joint root and Newton failures do not veto coefficient attempt", {
  fit_capture <- ln_operational_capture_warnings(ln_joint_fit_fixture(
    max = c(0.01, 1000),
    ln_operational_fallback = "legacy_post_pml"
  )$fit)
  fit <- fit_capture$value
  expect_length(fit_capture$messages, 1L)
  expect_match(fit_capture$messages, "unavailable for 1 gene")
  expect_match(fit_capture$messages, "coefficients are available for 0")
  expect_match(fit_capture$messages, "coefficients remain unavailable for 1")
  expect_true(fit$operational_fallback$per_gene$fallback_attempted)
  expect_false(fit$operational_fallback$per_gene$coefficient_available)
  expect_identical(
    fit$operational_fallback$per_gene$fallback_status,
    "attempted_failed"
  )
  expect_identical(
    fit$operational_fallback$per_gene$fallback_failure_reason,
    "legacy_nonconvergence"
  )
  eligibility <- fit$operational_fallback$diagnostics$per_gene[[1L]]$eligibility
  expect_false(eligibility$primary_joint_candidate_checks[["root"]])
  expect_false(eligibility$primary_joint_candidate_checks[["newton"]])
  expect_false(eligibility$primary_joint_candidate_checks[["interior"]])
  expect_true(eligibility$joint_root_or_newton_failure)
  ln_operational_expect_no_inference_payload(fit$operational_fallback)
})

test_that("none skips fallback helpers and legacy optional tests are unchanged", {
  fx <- ln_joint_test_fixture()
  local_mocked_bindings(
    ln_operational_fallback_eligibility = function(...) {
      stop("fallback helper should not run")
    },
    .package = "scRoPE"
  )
  expect_no_error(scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    ln_estimator = "joint_ln",
    ln_operational_fallback = "none"
  ))

  legacy_args <- list(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    additional_tests = c("score", "lrt"),
    lrt_details = TRUE,
    keep_diagnostics = TRUE,
    ln_estimator = "legacy_post_pml"
  )
  set.seed(20260721)
  legacy_default <- suppressWarnings(do.call(scrope, legacy_args))
  set.seed(20260721)
  legacy_none <- suppressWarnings(do.call(
    scrope,
    c(legacy_args, list(ln_operational_fallback = "none"))
  ))
  expect_identical(legacy_default, legacy_none)
})

test_that("output modes and design features remain primary-only", {
  fx <- ln_joint_test_fixture()
  donor_factor <- factor(rep(c("a", "b"), length.out = ncol(fx$count)))
  pred <- model.matrix(~ fx$pred[, "x"] + donor_factor)
  pred <- pred[, c(2L, 1L, 3L), drop = FALSE]
  colnames(pred)[[1L]] <- "scaled_x"
  args <- list(
    count = fx$count,
    id = fx$id,
    pred = pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1,
    output_re = TRUE,
    keep_diagnostics = FALSE,
    ln_estimator = "joint_ln"
  )
  set.seed(20260721)
  none <- do.call(scrope, c(args, list(ln_operational_fallback = "none")))
  set.seed(20260721)
  enabled <- do.call(scrope, c(
    args,
    list(ln_operational_fallback = "legacy_post_pml")
  ))
  enabled_primary <- enabled
  enabled_primary$operational_fallback <- NULL
  expect_identical(enabled_primary, none)
  expect_null(enabled$diagnostics)
  expect_true(is.data.frame(enabled$random_effect))
  expect_true(is.list(enabled$operational_fallback$diagnostics))
  expect_false(any(grepl(
    "random_effect|conditional_residual",
    names(enabled$operational_fallback),
    ignore.case = TRUE
  )))
  ln_operational_expect_no_inference_payload(enabled$operational_fallback)
})

test_that("independent runs agree modulo run-local objective IDs", {
  fx <- ln_joint_test_fixture()
  count <- rbind(fx$count, fx$count)
  rownames(count) <- c("gene1", "gene2")
  args <- list(
    count = count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln",
    ln_operational_fallback = "legacy_post_pml"
  )
  set.seed(20260721)
  one <- do.call(scrope, c(args, list(ncore = 1L)))
  set.seed(20260721)
  two <- do.call(scrope, c(args, list(ncore = 1L)))
  ln_operational_expect_cross_run_identical(one, two)
})

test_that("controlled coefficient fallback agrees across independent runs", {
  fx <- ln_joint_test_fixture()
  count <- rbind(fx$count, fx$count)
  rownames(count) <- c("gene1", "gene2")
  original_eligibility <- getFromNamespace(
    "ln_operational_fallback_eligibility",
    "scRoPE"
  )
  local_mocked_bindings(
    ln_operational_fallback_eligibility = function(gene_result) {
      eligibility <- original_eligibility(gene_result)
      eligibility$eligible <- TRUE
      eligibility$primary_route <- "R2"
      eligibility$primary_failure_reason <- "dispersion_boundary"
      eligibility$boundary_only <- FALSE
      eligibility$joint_failed_checks <- c("root", "newton", "interior")
      eligibility$joint_root_or_newton_failure <- TRUE
      eligibility
    },
    .package = "scRoPE"
  )
  args <- list(
    count = count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln",
    ln_operational_fallback = "legacy_post_pml"
  )
  set.seed(20260721)
  one <- suppressWarnings(do.call(scrope, c(args, list(ncore = 1L))))
  set.seed(20260721)
  two <- suppressWarnings(do.call(scrope, c(args, list(ncore = 1L))))
  ln_operational_expect_cross_run_identical(one, two)
  expect_true(all(
    one$operational_fallback$per_gene$fallback_status ==
      "coefficient_available"
  ))
  expect_true(all(
    one$operational_fallback$per_gene$coefficient_available
  ))
  expect_true(all(
    !one$operational_fallback$per_gene$inference_supported
  ))
  ln_operational_expect_no_inference_payload(one$operational_fallback)
})
