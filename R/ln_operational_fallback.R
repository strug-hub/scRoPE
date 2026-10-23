ln_operational_selected_candidate <- function(gene_result) {
  records <- gene_result$ln_joint$candidate_records
  if (is.null(records) || !length(records)) {
    return(NULL)
  }
  certified <- vapply(records, function(x) isTRUE(x$certified), logical(1))
  finite <- vapply(records, function(x) {
    isTRUE(x$finite_endpoint) && length(x$objective) == 1L &&
      is.finite(x$objective)
  }, logical(1))
  eligible_index <- if (any(certified)) which(certified) else which(finite)
  if (!length(eligible_index)) {
    return(NULL)
  }
  objectives <- vapply(records[eligible_index], function(x) {
    as.numeric(x$objective[[1L]])
  }, numeric(1))
  records[[eligible_index[[which.min(objectives)]]]]
}

ln_operational_fallback_eligibility <- function(gene_result) {
  required_checks <- c(
    "donor_additivity", "root", "newton", "curvature", "condition",
    "interior", "constraint"
  )
  selected <- ln_operational_selected_candidate(gene_result)
  checks <- if (is.null(selected)) NULL else selected$checks
  complete_checks <- !is.null(checks) &&
    identical(names(checks), required_checks) &&
    length(checks) == length(required_checks) &&
    all(!is.na(checks))
  selected_state <- if (is.null(selected)) NULL else selected$theta_full_log
  primary_state <- gene_result$ln_joint$state$full_log
  state_match <- !is.null(selected_state) && !is.null(primary_state) &&
    identical(as.numeric(selected_state), as.numeric(primary_state))
  objective_match <- !is.null(selected) &&
    length(selected$objective) == 1L && is.finite(selected$objective) &&
    length(gene_result$loglik) == 1L && is.finite(gene_result$loglik) &&
    identical(as.numeric(selected$objective), -as.numeric(gene_result$loglik))
  boundary_only <- isTRUE(complete_checks) &&
    isTRUE(state_match) && isTRUE(objective_match) &&
    identical(unname(as.logical(checks[["interior"]])), FALSE) &&
    all(as.logical(checks[setdiff(required_checks, "interior")]))
  resolution <- gene_result$ln_joint$resolution
  strict_route <- if (is.list(resolution)) {
    ln_availability_scalar_character(
      resolution$strict_inference_route,
      ln_availability_scalar_character(
        resolution$strict_route,
        gene_result$route
      )
    )
  } else {
    ln_availability_scalar_character(gene_result$route)
  }
  strict_failure <- if (is.list(resolution)) {
    ln_availability_scalar_character(
      resolution$strict_inference_failure_reason,
      ln_availability_scalar_character(
        resolution$strict_failure_reason,
        gene_result$failure_reason
      )
    )
  } else {
    ln_availability_scalar_character(gene_result$failure_reason)
  }
  actual_route <- ln_availability_scalar_character(gene_result$route)
  actual_failure <- ln_availability_scalar_character(
    gene_result$failure_reason
  )
  actual_boundary <- identical(actual_route, "R2") &&
    identical(actual_failure, "dispersion_boundary")
  strict_boundary <- identical(strict_route, "R2") &&
    identical(strict_failure, "dispersion_boundary")
  failed_checks <- if (isTRUE(complete_checks)) {
    names(checks)[!as.logical(checks)]
  } else {
    character(0)
  }
  list(
    eligible = isTRUE(actual_boundary || strict_boundary),
    eligibility_rule =
      paste0(
        "actual_or_strict_comparator_R2_and_",
        "failure_reason_dispersion_boundary"
      ),
    selected_check_audit_used_for_eligibility = FALSE,
    primary_route = if (is.null(gene_result$route)) NA_character_ else
      as.character(gene_result$route[[1L]]),
    primary_failure_reason = if (is.null(gene_result$failure_reason)) {
      NA_character_
    } else {
      as.character(gene_result$failure_reason[[1L]])
    },
    strict_comparator_route = strict_route,
    strict_comparator_failure_reason = strict_failure,
    selected_candidate = selected,
    primary_joint_candidate_checks = checks,
    complete_checks = complete_checks,
    selected_state_match = state_match,
    selected_objective_match = objective_match,
    boundary_only = boundary_only,
    joint_failed_checks = failed_checks,
    joint_root_or_newton_failure = any(
      c("root", "newton") %in% failed_checks
    )
  )
}

ln_operational_boundary_audit <- function(subject, gamma, min, max) {
  subject_bounds <- as.numeric(c(min[[1L]], max[[1L]]))
  gamma_bounds <- as.numeric(c(min[[2L]], max[[2L]]))
  cell_bounds <- c(1 / gamma_bounds[[2L]], 1 / gamma_bounds[[1L]])
  tolerances <- sqrt(.Machine$double.eps) *
    (1 + pmax(abs(as.numeric(min)), abs(as.numeric(max))))

  valid <- all(is.finite(c(subject, gamma))) && subject > 0 && gamma > 0 &&
    subject >= subject_bounds[[1L]] && subject <= subject_bounds[[2L]] &&
    gamma >= gamma_bounds[[1L]] && gamma <= gamma_bounds[[2L]]
  cell <- if (is.finite(gamma) && gamma > 0) 1 / gamma else NA_real_
  components <- character(0)
  if (is.finite(subject)) {
    if (subject <= subject_bounds[[1L]] + tolerances[[1L]]) {
      components <- c(components, "subject_lower")
    }
    if (subject >= subject_bounds[[2L]] - tolerances[[1L]]) {
      components <- c(components, "subject_upper")
    }
  }
  if (is.finite(gamma)) {
    if (gamma >= gamma_bounds[[2L]] - tolerances[[2L]]) {
      components <- c(components, "cell_lower")
    }
    if (gamma <= gamma_bounds[[1L]] + tolerances[[2L]]) {
      components <- c(components, "cell_upper")
    }
  }
  subject_range <- diff(subject_bounds)
  cell_range <- diff(cell_bounds)
  list(
    valid = valid,
    interior = isTRUE(valid) && !length(components),
    boundary = isTRUE(valid) && length(components) > 0L,
    component = if (length(components)) paste(components, collapse = ";") else
      NA_character_,
    cell = cell,
    distance_subject_lower = (subject - subject_bounds[[1L]]) / subject_range,
    distance_subject_upper = (subject_bounds[[2L]] - subject) / subject_range,
    distance_cell_lower = (cell - cell_bounds[[1L]]) / cell_range,
    distance_cell_upper = (cell_bounds[[2L]] - cell) / cell_range,
    natural_tolerance_subject = tolerances[[1L]],
    natural_tolerance_gamma = tolerances[[2L]],
    subject_bounds = subject_bounds,
    gamma_bounds = gamma_bounds,
    cell_bounds = cell_bounds
  )
}

ln_operational_legacy_start <- function(legacy_fit) {
  rescue <- legacy_fit$debug$rescue
  if (!is.null(rescue) && length(rescue) && !is.na(rescue[[1L]]) &&
      !identical(as.character(rescue[[1L]]), "none")) {
    "pql_start"
  } else {
    "direct"
  }
}

ln_operational_route_class <- function(start_provenance, boundary) {
  rescued <- identical(start_provenance, "pql_start")
  at_boundary <- isTRUE(boundary$boundary)
  if (rescued && at_boundary) {
    return("rescued_boundary")
  }
  if (rescued) {
    return("rescued")
  }
  if (at_boundary) {
    return("boundary")
  }
  "uncertified"
}

ln_operational_empty_assessment <- function(failure_reason = NA_character_) {
  list(
    coefficient_available = FALSE,
    failure_reason = normalise_operational_fallback_failure(failure_reason),
    convergence = NA_integer_,
    algorithm = NA_integer_,
    start_provenance = NA_character_,
    route_class = NA_character_,
    boundary = NA,
    boundary_component = NA_character_,
    distance_subject_lower = NA_real_,
    distance_subject_upper = NA_real_,
    distance_cell_lower = NA_real_,
    distance_cell_upper = NA_real_,
    beta_public = NULL,
    subject = NA_real_,
    cell = NA_real_,
    inference_supported = FALSE,
    inference_status = "unsupported_operational_coefficient_only",
    diagnostics = list()
  )
}

ln_operational_assess_legacy_fit <- function(legacy_fit, ctx) {
  nb <- ctx$nb
  expected_length <- 3L * nb + 5L
  if (is.null(legacy_fit$vec) || length(legacy_fit$vec) < expected_length) {
    return(ln_operational_empty_assessment("legacy_fit_error"))
  }
  vec <- as.numeric(legacy_fit$vec)
  beta_internal <- vec[seq_len(nb)]
  subject <- vec[[nb + 1L]]
  cell <- vec[[nb + 2L]]
  gamma <- if (is.finite(cell) && cell > 0) 1 / cell else NA_real_
  convergence <- as.integer(vec[[3L * nb + 3L]])
  algorithm <- as.integer(vec[[3L * nb + 4L]])
  start_provenance <- ln_operational_legacy_start(legacy_fit)
  boundary <- ln_operational_boundary_audit(
    subject, gamma, ctx$min, ctx$max
  )
  result <- ln_operational_empty_assessment()
  result$convergence <- convergence
  result$algorithm <- algorithm
  result$start_provenance <- start_provenance
  result$route_class <- ln_operational_route_class(start_provenance, boundary)
  result$boundary <- boundary$boundary
  result$boundary_component <- boundary$component
  result$distance_subject_lower <- boundary$distance_subject_lower
  result$distance_subject_upper <- boundary$distance_subject_upper
  result$distance_cell_lower <- boundary$distance_cell_lower
  result$distance_cell_upper <- boundary$distance_cell_upper

  fail <- function(reason, diagnostics = list()) {
    result$failure_reason <- normalise_operational_fallback_failure(reason)
    result$diagnostics <- c(result$diagnostics, diagnostics)
    result
  }
  if (is.na(convergence) || convergence < 0L) {
    reason <- if (identical(convergence, -10L)) {
      "legacy_gradient_convergence_not_accepted"
    } else {
      "legacy_nonconvergence"
    }
    return(fail(reason, list(boundary = boundary)))
  }
  if (any(!is.finite(beta_internal))) {
    return(fail("nonfinite_coefficient", list(boundary = boundary)))
  }
  if (!isTRUE(boundary$valid)) {
    return(fail("invalid_dispersion", list(boundary = boundary)))
  }
  O <- as.matrix(ctx$coef_internal_to_public)
  if (!identical(dim(O), c(nb, nb)) || any(!is.finite(O))) {
    return(fail("invalid_coefficient_transform", list(boundary = boundary)))
  }
  beta_public <- as.numeric(O %*% beta_internal)
  if (any(!is.finite(beta_public))) {
    return(fail("nonfinite_coefficient", list(boundary = boundary)))
  }

  result$coefficient_available <- TRUE
  result$failure_reason <- NA_character_
  result$beta_public <- beta_public
  result$subject <- subject
  result$cell <- cell
  result$diagnostics <- list(
    boundary = boundary,
    route_class = result$route_class,
    convergence = convergence,
    algorithm = algorithm,
    start_provenance = start_provenance,
    rescue_label = if (is.null(legacy_fit$debug$rescue) ||
        !length(legacy_fit$debug$rescue)) {
      NA_character_
    } else {
      as.character(legacy_fit$debug$rescue[[1L]])
    },
    coefficient_internal = beta_internal,
    coefficient_public = beta_public,
    dispersion = list(subject = subject, cell = cell, gamma_nb_size = gamma),
    inference_supported = FALSE,
    inference_status = result$inference_status
  )
  result
}

ln_operational_attempt_gene <- function(gene_index, posv, ctx) {
  stop(paste0(
    "The legacy-LN operational fallback is retired. Use ",
    "ln_operational_fallback='none'. For historical reproduction, use the ",
    "exact preserved scRoPE version recorded for the analysis ",
    "(0.2.5.9005 or earlier)."
  ), call. = FALSE)
}

ln_operational_build_sidecar <- function(primary_results, eligibility,
                                         fallback_results, gid, gene_names,
                                         term_names, ctx) {
  n_gene <- length(primary_results)
  n_term <- length(term_names)
  per_gene <- data.frame(
    gene_id = as.integer(gid),
    gene = as.character(gene_names),
    primary_route = vapply(eligibility, `[[`, character(1), "primary_route"),
    primary_failure_reason = vapply(
      eligibility, `[[`, character(1), "primary_failure_reason"
    ),
    fallback_attempted = rep(FALSE, n_gene),
    fallback_estimator = rep(NA_character_, n_gene),
    fallback_role = rep(NA_character_, n_gene),
    fallback_route_class = rep(NA_character_, n_gene),
    fallback_convergence = rep(NA_integer_, n_gene),
    fallback_algorithm = rep(NA_integer_, n_gene),
    fallback_start_provenance = rep(NA_character_, n_gene),
    fallback_boundary = rep(NA, n_gene),
    fallback_boundary_component = rep(NA_character_, n_gene),
    fallback_distance_subject_lower = rep(NA_real_, n_gene),
    fallback_distance_subject_upper = rep(NA_real_, n_gene),
    fallback_distance_cell_lower = rep(NA_real_, n_gene),
    fallback_distance_cell_upper = rep(NA_real_, n_gene),
    coefficient_available = rep(NA, n_gene),
    coefficient_source = rep(NA_character_, n_gene),
    fallback_failure_reason = rep(NA_character_, n_gene),
    fallback_status = rep("not_eligible", n_gene),
    inference_supported = rep(FALSE, n_gene),
    inference_status = rep(
      "unsupported_operational_coefficient_only",
      n_gene
    ),
    stringsAsFactors = FALSE
  )

  coefficient_rows <- vector("list", n_gene)
  overdispersion <- data.frame(
    gene_id = as.integer(gid),
    gene = as.character(gene_names),
    Subject = rep(NA_real_, n_gene),
    Cell = rep(NA_real_, n_gene),
    coefficient_available = rep(NA, n_gene),
    fallback_status = rep("not_eligible", n_gene),
    fallback_failure_reason = rep(NA_character_, n_gene),
    inference_supported = rep(FALSE, n_gene),
    inference_status = rep(
      "unsupported_operational_coefficient_only",
      n_gene
    ),
    stringsAsFactors = FALSE
  )
  diagnostic_records <- vector("list", n_gene)
  names(diagnostic_records) <- as.character(gene_names)

  for (index in seq_len(n_gene)) {
    attempted <- isTRUE(eligibility[[index]]$eligible)
    assessment <- fallback_results[[index]]
    if (attempted) {
      per_gene$fallback_attempted[[index]] <- TRUE
      per_gene$fallback_estimator[[index]] <-
        "legacy_post_pml_operational_coefficient"
      per_gene$fallback_role[[index]] <-
        "coefficient_dispersion_diagnostics_only"
      per_gene$fallback_route_class[[index]] <- assessment$route_class
      per_gene$fallback_convergence[[index]] <- assessment$convergence
      per_gene$fallback_algorithm[[index]] <- assessment$algorithm
      per_gene$fallback_start_provenance[[index]] <- assessment$start_provenance
      per_gene$fallback_boundary[[index]] <- assessment$boundary
      per_gene$fallback_boundary_component[[index]] <-
        assessment$boundary_component
      per_gene$fallback_distance_subject_lower[[index]] <-
        assessment$distance_subject_lower
      per_gene$fallback_distance_subject_upper[[index]] <-
        assessment$distance_subject_upper
      per_gene$fallback_distance_cell_lower[[index]] <-
        assessment$distance_cell_lower
      per_gene$fallback_distance_cell_upper[[index]] <-
        assessment$distance_cell_upper
      per_gene$coefficient_available[[index]] <-
        isTRUE(assessment$coefficient_available)
      per_gene$fallback_failure_reason[[index]] <- assessment$failure_reason
      per_gene$fallback_status[[index]] <- if (
          isTRUE(assessment$coefficient_available)) {
        "coefficient_available"
      } else {
        "attempted_failed"
      }
      overdispersion$coefficient_available[[index]] <-
        per_gene$coefficient_available[[index]]
      overdispersion$fallback_status[[index]] <-
        per_gene$fallback_status[[index]]
      overdispersion$fallback_failure_reason[[index]] <-
        assessment$failure_reason
      if (isTRUE(assessment$coefficient_available)) {
        per_gene$coefficient_source[[index]] <- "legacy_post_pml"
        overdispersion$Subject[[index]] <- assessment$subject
        overdispersion$Cell[[index]] <- assessment$cell
      }
    }

    available <- if (attempted) {
      isTRUE(assessment$coefficient_available)
    } else {
      NA
    }
    coefficient_rows[[index]] <- data.frame(
      gene_id = rep(as.integer(gid[[index]]), n_term),
      gene = rep(as.character(gene_names[[index]]), n_term),
      term = as.character(term_names),
      estimate = if (isTRUE(available)) assessment$beta_public else
        rep(NA_real_, n_term),
      coefficient_available = rep(available, n_term),
      coefficient_source = rep(
        if (isTRUE(available)) "legacy_post_pml" else NA_character_,
        n_term
      ),
      fallback_route_class = rep(
        if (attempted) assessment$route_class else NA_character_,
        n_term
      ),
      fallback_status = rep(per_gene$fallback_status[[index]], n_term),
      fallback_failure_reason = rep(
        if (attempted && !isTRUE(available)) assessment$failure_reason else
          NA_character_,
        n_term
      ),
      inference_supported = rep(FALSE, n_term),
      inference_status = rep(
        "unsupported_operational_coefficient_only",
        n_term
      ),
      stringsAsFactors = FALSE
    )
    diagnostic_records[[index]] <- list(
      eligibility = eligibility[[index]],
      attempted = attempted,
      assessment = if (attempted) assessment$diagnostics else NULL
    )
  }

  list(
    schema_version = 2L,
    method = "legacy_post_pml",
    role = "operational_coefficient_dispersion_diagnostics_only",
    selection_rule =
      paste0(
        "actual_or_strict_comparator_R2_and_",
        "failure_reason_dispersion_boundary"
      ),
    inference_supported = FALSE,
    inference_status = "unsupported_for_all_fallback_routes",
    per_gene = per_gene,
    coefficients = do.call(rbind, coefficient_rows),
    overdispersion = overdispersion,
    diagnostics = list(
      policy = "coefficient_dispersion_diagnostics_only",
      boundary_tolerance = "sqrt(.Machine$double.eps)*(1+pmax(abs(min),abs(max)))",
      per_gene = diagnostic_records
    )
  )
}

ln_operational_warn <- function(sidecar) {
  per_gene <- sidecar$per_gene
  r2 <- !is.na(per_gene$primary_route) & per_gene$primary_route == "R2"
  n_r2 <- sum(r2)
  if (!n_r2) {
    return(invisible(NULL))
  }
  available <- r2 & per_gene$fallback_status == "coefficient_available"
  boundary_available <- available & !is.na(per_gene$fallback_boundary) &
    per_gene$fallback_boundary
  warning(sprintf(
    paste0(
      "Joint-LN inference was unavailable for %d gene%s. Independent legacy ",
      "post-PML operational coefficients are available for %d; operational ",
      "coefficients remain unavailable for %d. Available coefficients on or ",
      "near a dispersion boundary: %d. The fallback sidecar contains ",
      "coefficients, dispersions, and diagnostics only for rescued, boundary, ",
      "or otherwise uncertified legacy routes. It provides no covariance, ",
      "standard errors, Wald statistics, p-values, inferential success, or ",
      "FDR-eligible output."
    ),
    n_r2,
    if (n_r2 == 1L) "" else "s",
    sum(available),
    n_r2 - sum(available),
    sum(boundary_available)
  ), call. = FALSE)
  invisible(NULL)
}
