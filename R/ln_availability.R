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
      status_message[[index]] <- "Certified inference"
      inference_supported[[index]] <- TRUE
      estimate_source[[index]] <- "joint_ln"
      attention_flag[[index]] <- FALSE
      technical_code[[index]] <- "joint_ln_certified"
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
