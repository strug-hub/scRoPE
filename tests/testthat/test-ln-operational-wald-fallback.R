# These pure helper checks preserve interpretation of historical sidecars.
# Current fitting APIs never create or execute a legacy operational fallback.
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
    paste0(
      "actual_or_strict_comparator_R2_and_",
      "failure_reason_dispersion_boundary"
    )
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
    paste0(
      "actual_or_strict_comparator_R2_and_",
      "failure_reason_dispersion_boundary"
    )
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

test_that("none preserves primary and never invokes legacy fallback helpers", {
  local_mocked_bindings(
    ln_operational_fallback_eligibility = function(...) stop("fallback invoked"),
    ln_operational_attempt_gene = function(...) stop("fallback invoked"),
    .package = "scRoPE"
  )
  fx <- ln_joint_test_fixture()
  args <- list(count = fx$count, id = fx$id, pred = fx$pred, offset = fx$offset,
               verbose = FALSE, cpc = 0, mincp = 1, output_re = TRUE)
  set.seed(20260721)
  default <- do.call(scrope, args)
  set.seed(20260721)
  explicit <- do.call(scrope, c(args, list(ln_operational_fallback = "none")))
  ln_operational_expect_cross_run_identical(default, explicit)
  expect_false("operational_fallback" %in% names(default))
  expect_error(do.call(scrope, c(args, list(
    ln_operational_fallback = "legacy_post_pml"
  ))), "retired")
})

test_that("unavailable joint LN never triggers a hidden legacy fallback", {
  original_joint <- getFromNamespace("fit_gene_joint_ln_unconstrained", "scRoPE")
  local_mocked_bindings(
    fit_gene_joint_ln_unconstrained =
      ln_operational_force_failed_joint_primary(original_joint),
    ln_operational_fallback_eligibility = function(...) stop("fallback invoked"),
    ln_operational_attempt_gene = function(...) stop("fallback invoked"),
    .package = "scRoPE"
  )
  fx <- ln_joint_test_fixture()
  captured <- ln_operational_capture_warnings(scrope(
    fx$count, fx$id, pred = fx$pred, offset = fx$offset,
    verbose = FALSE, cpc = 0, mincp = 1
  ))
  fit <- captured$value
  expect_true(any(grepl("inference was unavailable", captured$messages)))
  expect_identical(fit$availability$availability, "unavailable")
  expect_false(fit$availability$inference_supported)
  expect_false("operational_fallback" %in% names(fit))
  inference_columns <- grep("^(se_|p_)", names(fit$summary), value = TRUE)
  expect_gt(length(inference_columns), 0L)
  expect_true(all(is.na(unlist(fit$summary[inference_columns]))))
})
