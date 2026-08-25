ln_availability_scalar_character <- function(value, default = NA_character_) {
  if (is.null(value) || !length(value)) {
    return(default)
  }
  value <- as.character(value[[1L]])
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  value
}

ln_joint_resolution_status <- function(gene_result, resolution) {
  if (!is.list(resolution) && isTRUE(gene_result$inference_available)) {
    return("not_needed")
  }
  classes <- if (is.list(resolution) &&
      !is.null(resolution$global_veto_classes)) {
    values <- as.character(resolution$global_veto_classes)
    unique(values[!is.na(values) & nzchar(values)])
  } else {
    character(0)
  }
  if ("competitive_boundary" %in% classes) {
    return("competitive_boundary")
  }
  if ("competitive_distinct_basin" %in% classes) {
    return("competitive_distinct_basin")
  }
  if ("competitive_hard_stop" %in% classes) {
    return("hard_stop")
  }
  if ("same_basin_review_unresolved" %in% classes) {
    return("review_unresolved")
  }
  if (is.list(resolution) && isTRUE(resolution$would_nonveto) &&
      identical(
        ln_availability_scalar_character(resolution$prospective_route),
        "R0"
      )) {
    return("same_basin_review_resolved")
  }
  if (is.list(resolution) && identical(
      ln_availability_scalar_character(resolution$strict_route),
      "R0"
    )) {
    return("not_needed")
  }

  failure <- ln_availability_scalar_character(gene_result$failure_reason)
  if (is.na(failure) && is.list(resolution)) {
    failure <- ln_availability_scalar_character(
      resolution$strict_failure_reason
    )
  }
  failure_lower <- tolower(failure)
  if (!is.na(failure_lower) && grepl("boundary", failure_lower, fixed = TRUE)) {
    return("competitive_boundary")
  }
  if (!is.na(failure_lower) && grepl("distinct", failure_lower, fixed = TRUE)) {
    return("competitive_distinct_basin")
  }
  if (!is.na(failure_lower) && grepl("unresolved", failure_lower, fixed = TRUE)) {
    return("review_unresolved")
  }
  hard_stop_failures <- c(
    "required_start_failed", "singular_joint_curvature",
    "covariance_identity_failed", "derivative_evaluation_failed",
    "derivative_audit_failed"
  )
  if (!is.na(failure_lower) && (
      grepl("hard_stop", failure_lower, fixed = TRUE) ||
      failure_lower %in% hard_stop_failures)) {
    return("hard_stop")
  }
  "certification_failure"
}

ln_joint_status_record <- function(gene_result) {
  resolution <- gene_result$ln_joint$resolution
  inference_available <- isTRUE(gene_result$inference_available)
  actual_route <- ln_availability_scalar_character(
    gene_result$route,
    if (inference_available) "R0" else "R2"
  )
  strict_route <- if (is.list(resolution) &&
      "strict_inference_route" %in% names(resolution)) {
    ln_availability_scalar_character(resolution$strict_inference_route)
  } else if (is.list(resolution)) {
    ln_availability_scalar_character(resolution$strict_route, actual_route)
  } else {
    actual_route
  }
  enabled_route <- if (is.list(resolution) &&
      "enabled_inference_route" %in% names(resolution)) {
    ln_availability_scalar_character(resolution$enabled_inference_route)
  } else if (is.list(resolution)) {
    ln_availability_scalar_character(resolution$prospective_route, actual_route)
  } else {
    actual_route
  }
  strict_failure <- if (is.list(resolution) &&
      "strict_inference_failure_reason" %in% names(resolution)) {
    ln_availability_scalar_character(
      resolution$strict_inference_failure_reason
    )
  } else if (is.list(resolution)) {
    ln_availability_scalar_character(resolution$strict_failure_reason)
  } else {
    ln_availability_scalar_character(gene_result$failure_reason)
  }
  enabled_failure <- if (is.list(resolution) &&
      "enabled_inference_failure_reason" %in% names(resolution)) {
    ln_availability_scalar_character(
      resolution$enabled_inference_failure_reason
    )
  } else if (is.list(resolution)) {
    ln_availability_scalar_character(resolution$prospective_failure_reason)
  } else {
    ln_availability_scalar_character(gene_result$failure_reason)
  }
  selected_root_certified <- if (is.list(resolution) &&
      is.logical(resolution$selected_root_certified) &&
      length(resolution$selected_root_certified) == 1L &&
      !is.na(resolution$selected_root_certified)) {
    isTRUE(resolution$selected_root_certified)
  } else {
    inference_available
  }
  newly_enabled <- inference_available && identical(strict_route, "R2") &&
    identical(enabled_route, "R0") && is.list(resolution) &&
    isTRUE(resolution$would_nonveto)
  transition <- if (identical(strict_route, "R0") &&
      identical(enabled_route, "R0")) {
    "unchanged_available"
  } else if (identical(strict_route, "R2") &&
      identical(enabled_route, "R0")) {
    "strict_withheld_to_available"
  } else if (identical(strict_route, "R0") &&
      identical(enabled_route, "R2")) {
    "strict_available_to_withheld"
  } else {
    "unchanged_withheld"
  }
  veto_classes <- if (is.list(resolution) &&
      !is.null(resolution$global_veto_classes)) {
    values <- sort(unique(as.character(resolution$global_veto_classes)))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values)) paste(values, collapse = ";") else NA_character_
  } else {
    NA_character_
  }

  list(
    estimator = "joint_ln",
    certified_root = selected_root_certified,
    inference_status = if (inference_available) "available" else "withheld",
    strict_inference_status = if (is.na(strict_route)) {
      "not_assessed"
    } else if (identical(strict_route, "R0")) {
      "available"
    } else {
      "withheld"
    },
    enabled_inference_status = if (is.na(enabled_route)) {
      "not_assessed"
    } else if (identical(enabled_route, "R0")) {
      "available"
    } else {
      "withheld"
    },
    resolution_status = ln_joint_resolution_status(gene_result, resolution),
    inference_source = if (inference_available) {
      "certified_selected_root"
    } else {
      "none"
    },
    newly_enabled = newly_enabled,
    availability_transition = transition,
    resolution_veto_classes = veto_classes,
    resolution_policy_version = if (is.list(resolution)) {
      ln_availability_scalar_character(resolution$policy_version)
    } else {
      NA_character_
    },
    resolution_mode = if (is.list(resolution)) {
      ln_availability_scalar_character(resolution$mode)
    } else {
      NA_character_
    },
    route_code = actual_route,
    strict_route_code = strict_route,
    enabled_route_code = enabled_route,
    strict_failure_reason = strict_failure,
    enabled_failure_reason = enabled_failure
  )
}

ln_availability_strict_failure <- function(gene_result) {
  failure <- ln_availability_scalar_character(gene_result$failure_reason)
  if (!is.na(failure)) {
    return(failure)
  }
  resolution <- gene_result$ln_joint$resolution
  failure <- if (is.list(resolution)) {
    ln_availability_scalar_character(resolution$strict_failure_reason)
  } else {
    NA_character_
  }
  if (is.na(failure)) "no_certified_candidate" else failure
}

ln_availability_resolution_classes <- function(gene_result) {
  resolution <- gene_result$ln_joint$resolution
  if (!is.list(resolution) || is.null(resolution$global_veto_classes)) {
    return(character(0))
  }
  classes <- as.character(resolution$global_veto_classes)
  unique(classes[!is.na(classes) & nzchar(classes)])
}

ln_availability_unavailable_reason <- function(gene_result, technical_code) {
  classes <- ln_availability_resolution_classes(gene_result)
  code <- tolower(technical_code)

  boundary <- any(grepl("boundary", classes, fixed = TRUE)) ||
    grepl("boundary", code, fixed = TRUE)
  if (isTRUE(boundary)) {
    return(list(
      reason_group = "competitive_boundary",
      status_message = paste0(
        "Inference withheld \u2014 a boundary solution could not be ruled out"
      )
    ))
  }

  distinct <- any(grepl("distinct", classes, fixed = TRUE)) ||
    grepl("distinct", code, fixed = TRUE)
  if (isTRUE(distinct)) {
    return(list(
      reason_group = "competitive_solutions",
      status_message = paste0(
        "Inference withheld \u2014 materially different solutions remain competitive"
      )
    ))
  }

  hard_stop <- any(grepl("hard_stop", classes, fixed = TRUE)) ||
    grepl("hard_stop", code, fixed = TRUE) ||
    code %in% c(
      "required_start_failed",
      "singular_joint_curvature",
      "covariance_identity_failed",
      "derivative_evaluation_failed",
      "derivative_audit_failed"
    )
  if (isTRUE(hard_stop)) {
    return(list(
      reason_group = "hard_stop",
      status_message = "A reliable solution could not be certified"
    ))
  }

  list(
    reason_group = "certification_failure",
    status_message = "A reliable solution could not be certified"
  )
}

ln_availability_fallback_available <- function(sidecar, gid, gene_names) {
  n_gene <- length(gid)
  available <- rep(FALSE, n_gene)
  if (is.null(sidecar)) {
    return(available)
  }
  if (!is.list(sidecar) ||
      !identical(sidecar$schema_version, 2L) ||
      !identical(sidecar$method, "legacy_post_pml") ||
      !identical(
        sidecar$role,
        "operational_coefficient_dispersion_diagnostics_only"
      ) ||
      !identical(sidecar$inference_supported, FALSE)) {
    return(available)
  }

  per_gene <- sidecar$per_gene
  required <- c(
    "gene_id", "gene", "fallback_attempted", "fallback_estimator",
    "coefficient_available", "coefficient_source", "fallback_status",
    "inference_supported"
  )
  if (!is.data.frame(per_gene) ||
      !all(required %in% names(per_gene)) ||
      nrow(per_gene) != n_gene ||
      !identical(as.character(per_gene$gene_id), as.character(gid)) ||
      !identical(as.character(per_gene$gene), as.character(gene_names))) {
    return(available)
  }

  coefficients <- sidecar$coefficients
  coefficient_required <- c(
    "gene_id", "gene", "term", "estimate", "coefficient_available",
    "coefficient_source", "fallback_status", "inference_supported"
  )
  if (!is.data.frame(coefficients) ||
      !all(coefficient_required %in% names(coefficients))) {
    return(available)
  }
  coefficient_terms <- as.character(coefficients$term)
  if (!nrow(coefficients) || any(is.na(coefficient_terms)) ||
      any(!nzchar(coefficient_terms))) {
    return(available)
  }
  common_terms <- sort(unique(coefficient_terms))
  rectangular_payload <- vapply(seq_len(n_gene), function(index) {
    rows <- which(
      as.character(coefficients$gene_id) == as.character(gid[[index]]) &
        as.character(coefficients$gene) == as.character(gene_names[[index]])
    )
    length(rows) == length(common_terms) &&
      !anyDuplicated(coefficient_terms[rows]) &&
      identical(sort(coefficient_terms[rows]), common_terms)
  }, logical(1))
  known_payload_row <- vapply(seq_len(nrow(coefficients)), function(index) {
    isTRUE(any(
      as.character(gid) == as.character(coefficients$gene_id[[index]]) &
        as.character(gene_names) == as.character(coefficients$gene[[index]])
    ))
  }, logical(1))
  if (!all(rectangular_payload) || !all(known_payload_row)) {
    return(available)
  }

  vapply(seq_len(n_gene), function(index) {
    payload_rows <- which(
      as.character(coefficients$gene_id) == as.character(gid[[index]]) &
        as.character(coefficients$gene) == as.character(gene_names[[index]])
    )
    payload_valid <- length(payload_rows) > 0L &&
      all(!is.na(coefficients$term[payload_rows])) &&
      all(nzchar(as.character(coefficients$term[payload_rows]))) &&
      is.numeric(coefficients$estimate) &&
      all(is.finite(coefficients$estimate[payload_rows])) &&
      all(coefficients$coefficient_available[payload_rows] %in% TRUE) &&
      all(coefficients$coefficient_source[payload_rows] == "legacy_post_pml") &&
      all(coefficients$fallback_status[payload_rows] ==
        "coefficient_available") &&
      all(coefficients$inference_supported[payload_rows] %in% FALSE)

    isTRUE(payload_valid) &&
      isTRUE(per_gene$fallback_attempted[[index]]) &&
      identical(
        per_gene$fallback_estimator[[index]],
        "legacy_post_pml_operational_coefficient"
      ) &&
      isTRUE(per_gene$coefficient_available[[index]]) &&
      identical(per_gene$coefficient_source[[index]], "legacy_post_pml") &&
      identical(per_gene$fallback_status[[index]], "coefficient_available") &&
      identical(per_gene$inference_supported[[index]], FALSE)
  }, logical(1))
}

ln_availability_build <- function(primary_results, gid, gene_names,
                                  operational_fallback = NULL) {
  n_gene <- length(primary_results)
  if (length(gid) != n_gene || length(gene_names) != n_gene) {
    stop("Internal joint-LN availability inputs are misaligned.")
  }

  availability <- rep("unavailable", n_gene)
  reason_group <- rep("certification_failure", n_gene)
  status_message <- rep("A reliable solution could not be certified", n_gene)
  inference_supported <- rep(FALSE, n_gene)
  estimate_source <- rep(NA_character_, n_gene)
  attention_flag <- rep(TRUE, n_gene)
  technical_code <- rep("no_certified_candidate", n_gene)
  status_records <- lapply(primary_results, ln_joint_status_record)

  fallback_available <- ln_availability_fallback_available(
    operational_fallback,
    gid,
    gene_names
  )

  for (index in seq_len(n_gene)) {
    primary_available <- isTRUE(primary_results[[index]]$inference_available)
    if (primary_available) {
      availability[[index]] <- "inference"
      reason_group[[index]] <- "certified_inference"
      status_message[[index]] <- if (isTRUE(
          status_records[[index]]$newly_enabled)) {
        "Certified inference \u2014 numerical review resolved"
      } else {
        "Certified inference"
      }
      inference_supported[[index]] <- TRUE
      estimate_source[[index]] <- "joint_ln"
      attention_flag[[index]] <- FALSE
      technical_code[[index]] <- if (isTRUE(
          status_records[[index]]$newly_enabled)) {
        "joint_ln_certified_after_same_basin_review"
      } else {
        "joint_ln_certified"
      }
      next
    }

    if (isTRUE(fallback_available[[index]])) {
      availability[[index]] <- "estimate_only"
      reason_group[[index]] <- "operational_estimate_only"
      status_message[[index]] <-
        "Estimate only \u2014 uncertainty is not available"
      estimate_source[[index]] <- "legacy_post_pml"
      technical_code[[index]] <-
        "legacy_post_pml_coefficient_available"
      next
    }

    technical_code[[index]] <-
      ln_availability_strict_failure(primary_results[[index]])
    reason <- ln_availability_unavailable_reason(
      primary_results[[index]],
      technical_code[[index]]
    )
    reason_group[[index]] <- reason$reason_group
    status_message[[index]] <- reason$status_message
  }

  data.frame(
    gene_id = as.integer(gid),
    gene = as.character(gene_names),
    estimator = vapply(status_records, `[[`, character(1), "estimator"),
    certified_root = vapply(
      status_records, `[[`, logical(1), "certified_root"
    ),
    inference_status = vapply(
      status_records, `[[`, character(1), "inference_status"
    ),
    strict_inference_status = vapply(
      status_records, `[[`, character(1), "strict_inference_status"
    ),
    enabled_inference_status = vapply(
      status_records, `[[`, character(1), "enabled_inference_status"
    ),
    resolution_status = vapply(
      status_records, `[[`, character(1), "resolution_status"
    ),
    inference_source = vapply(
      status_records, `[[`, character(1), "inference_source"
    ),
    newly_enabled = vapply(
      status_records, `[[`, logical(1), "newly_enabled"
    ),
    availability_transition = vapply(
      status_records, `[[`, character(1), "availability_transition"
    ),
    resolution_policy_version = vapply(
      status_records, `[[`, character(1), "resolution_policy_version"
    ),
    availability = availability,
    reason_group = reason_group,
    status_message = status_message,
    inference_supported = inference_supported,
    estimate_source = estimate_source,
    attention_flag = attention_flag,
    technical_code = technical_code,
    stringsAsFactors = FALSE
  )
}

ln_availability_warn <- function(availability,
                                 operational_fallback_requested = FALSE) {
  if (!is.data.frame(availability) || !nrow(availability)) {
    return(invisible(NULL))
  }
  required <- c("availability", "attention_flag")
  if (!all(required %in% names(availability))) {
    return(invisible(NULL))
  }

  attention <- availability$attention_flag %in% TRUE
  n_unavailable_inference <- sum(attention)
  if (!n_unavailable_inference) {
    return(invisible(NULL))
  }
  n_gene <- nrow(availability)
  n_inference <- sum(availability$availability == "inference", na.rm = TRUE)
  n_estimate_only <- sum(
    availability$availability == "estimate_only",
    na.rm = TRUE
  )
  n_unavailable <- sum(
    availability$availability == "unavailable",
    na.rm = TRUE
  )

  message <- sprintf(
    paste0(
      "Joint-LN inference was unavailable for %d gene%s. Certified primary ",
      "inference is available for %d of %d."
    ),
    n_unavailable_inference,
    if (n_unavailable_inference == 1L) "" else "s",
    n_inference,
    n_gene
  )
  if (isTRUE(operational_fallback_requested)) {
    message <- paste0(
      message,
      sprintf(
        paste0(
          " Independent legacy post-PML operational coefficients are ",
          "available for %d; operational coefficients remain unavailable ",
          "for %d."
        ),
        n_estimate_only,
        n_unavailable
      ),
      paste0(
        " The fallback sidecar contains coefficients, dispersions, and ",
        "diagnostics only and provides no covariance, standard errors, Wald ",
        "statistics, p-values, inferential success, or FDR-eligible output."
      )
    )
  } else {
    message <- paste0(
      message,
      sprintf(
        paste0(
          " No operational fallback was requested; operational coefficients ",
          "remain unavailable for %d."
        ),
        n_unavailable
      )
    )
  }
  message <- paste0(
    message,
    " See fit$availability for one row per gene/request."
  )
  warning(message, call. = FALSE)
  invisible(NULL)
}
