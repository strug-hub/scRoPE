ln_availability_test_primary <- function(
    inference_available = FALSE,
    failure_reason = "no_certified_candidate",
    resolution = NULL) {
  list(
    inference_available = inference_available,
    failure_reason = failure_reason,
    vec = c(0.4, 0.2, 0.3),
    ln_joint = list(
      state = list(full_log = c(0.4, log(0.2), log(0.3))),
      resolution = resolution
    )
  )
}

ln_availability_test_sidecar <- function(gid, gene_names, available) {
  n_gene <- length(gid)
  n_term <- 2L
  data <- data.frame(
    gene_id = as.integer(gid),
    gene = as.character(gene_names),
    fallback_attempted = rep(TRUE, n_gene),
    fallback_estimator = rep(
      "legacy_post_pml_operational_coefficient",
      n_gene
    ),
    coefficient_available = as.logical(available),
    coefficient_source = ifelse(available, "legacy_post_pml", NA_character_),
    fallback_status = ifelse(
      available,
      "coefficient_available",
      "attempted_failed"
    ),
    inference_supported = rep(FALSE, n_gene),
    stringsAsFactors = FALSE
  )
  coefficients <- data.frame(
    gene_id = rep(as.integer(gid), each = n_term),
    gene = rep(as.character(gene_names), each = n_term),
    term = rep(c("term_a", "term_b"), times = n_gene),
    estimate = unlist(lapply(seq_len(n_gene), function(index) {
      if (isTRUE(available[[index]])) {
        c(index / 10, -index / 20)
      } else {
        rep(NA_real_, n_term)
      }
    }), use.names = FALSE),
    coefficient_available = rep(as.logical(available), each = n_term),
    coefficient_source = rep(
      ifelse(available, "legacy_post_pml", NA_character_),
      each = n_term
    ),
    fallback_status = rep(
      ifelse(available, "coefficient_available", "attempted_failed"),
      each = n_term
    ),
    inference_supported = rep(FALSE, n_gene * n_term),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = 2L,
    method = "legacy_post_pml",
    role = "operational_coefficient_dispersion_diagnostics_only",
    inference_supported = FALSE,
    per_gene = data,
    coefficients = coefficients,
    diagnostics = list(preserved = TRUE)
  )
}

ln_availability_capture_warnings <- function(expr) {
  messages <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = messages)
}

test_that("availability has one typed row per gene rather than contrast", {
  primary <- list(
    ln_availability_test_primary(inference_available = TRUE),
    ln_availability_test_primary(failure_reason = "dispersion_boundary")
  )
  availability <- ln_availability_build(
    primary_results = primary,
    gid = c(7L, 9L),
    gene_names = c("gene_a", "gene_b")
  )

  expect_identical(names(availability), c(
    "gene_id", "gene", "estimator", "certified_root",
    "inference_status", "strict_inference_status",
    "enabled_inference_status", "resolution_status", "inference_source",
    "newly_enabled", "availability_transition",
    "resolution_policy_version", "availability", "reason_group",
    "status_message", "inference_supported", "estimate_source",
    "attention_flag", "technical_code"
  ))
  expect_identical(nrow(availability), 2L)
  expect_type(availability$gene_id, "integer")
  expect_type(availability$gene, "character")
  expect_type(availability$inference_supported, "logical")
  expect_type(availability$attention_flag, "logical")
  expect_identical(availability$estimator, rep("joint_ln", 2L))
  expect_identical(availability$certified_root, c(TRUE, FALSE))
  expect_identical(
    availability$inference_status,
    c("available", "withheld")
  )
  expect_identical(
    availability$strict_inference_status,
    c("available", "withheld")
  )
  expect_identical(
    availability$enabled_inference_status,
    c("available", "withheld")
  )
  expect_identical(
    availability$resolution_status,
    c("not_needed", "competitive_boundary")
  )
  expect_identical(
    availability$inference_source,
    c("certified_selected_root", "none")
  )
  expect_identical(availability$newly_enabled, c(FALSE, FALSE))
  expect_identical(
    availability$availability_transition,
    c("unchanged_available", "unchanged_withheld")
  )
  expect_identical(
    availability$availability,
    c("inference", "unavailable")
  )
  expect_identical(
    availability$reason_group,
    c("certified_inference", "competitive_boundary")
  )
  expect_identical(availability$estimate_source, c("joint_ln", NA_character_))
  expect_identical(availability$inference_supported, c(TRUE, FALSE))
  expect_identical(availability$attention_flag, c(FALSE, TRUE))
  expect_identical(
    availability$technical_code,
    c("joint_ln_certified", "dispersion_boundary")
  )
  expect_identical(unserialize(serialize(availability, NULL)), availability)
})

test_that("primary then explicit fallback precedence is fail closed", {
  gid <- 1:3
  gene_names <- paste0("gene_", gid)
  primary <- list(
    ln_availability_test_primary(inference_available = TRUE),
    ln_availability_test_primary(
      failure_reason = "competing_uncertified_candidate"
    ),
    ln_availability_test_primary(
      failure_reason = "competing_uncertified_candidate"
    )
  )
  sidecar <- ln_availability_test_sidecar(
    gid,
    gene_names,
    available = c(TRUE, TRUE, FALSE)
  )
  availability <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = sidecar
  )

  expect_identical(
    availability$availability,
    c("inference", "estimate_only", "unavailable")
  )
  expect_identical(
    availability$reason_group,
    c(
      "certified_inference",
      "operational_estimate_only",
      "certification_failure"
    )
  )
  expect_identical(
    availability$estimate_source,
    c("joint_ln", "legacy_post_pml", NA_character_)
  )
  expect_identical(
    availability$technical_code,
    c(
      "joint_ln_certified",
      "legacy_post_pml_coefficient_available",
      "competing_uncertified_candidate"
    )
  )

  invalid_sidecar <- sidecar
  invalid_sidecar$inference_supported <- TRUE
  fail_closed <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = invalid_sidecar
  )
  expect_identical(
    fail_closed$availability,
    c("inference", "unavailable", "unavailable")
  )
  expect_true(all(is.na(fail_closed$estimate_source[-1L])))

  missing_payload <- sidecar
  missing_payload$coefficients <- NULL
  missing_payload_result <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = missing_payload
  )
  expect_identical(
    missing_payload_result$availability,
    c("inference", "unavailable", "unavailable")
  )

  incomplete_rectangle <- sidecar
  missing_term_row <- which(
    incomplete_rectangle$coefficients$gene_id == gid[[2L]] &
      incomplete_rectangle$coefficients$term == "term_b"
  )[[1L]]
  incomplete_rectangle$coefficients <-
    incomplete_rectangle$coefficients[-missing_term_row, , drop = FALSE]
  incomplete_rectangle_result <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = incomplete_rectangle
  )
  expect_false(any(
    incomplete_rectangle_result$availability == "estimate_only"
  ))
  expect_identical(
    incomplete_rectangle_result$availability,
    c("inference", "unavailable", "unavailable")
  )

  for (identity_field in c("gene_id", "gene")) {
    missing_identity <- sidecar
    missing_identity$coefficients[[identity_field]][[1L]] <- if (
        identical(identity_field, "gene_id")) {
      NA_integer_
    } else {
      NA_character_
    }
    missing_identity_result <- NULL
    expect_no_error({
      missing_identity_result <- ln_availability_build(
        primary,
        gid,
        gene_names,
        operational_fallback = missing_identity
      )
    })
    expect_false(any(
      missing_identity_result$availability == "estimate_only"
    ))
    expect_identical(
      missing_identity_result$availability,
      c("inference", "unavailable", "unavailable")
    )
  }

  nonfinite_payload <- sidecar
  nonfinite_row <- which(
    nonfinite_payload$coefficients$gene_id == gid[[2L]]
  )[[1L]]
  nonfinite_payload$coefficients$estimate[[nonfinite_row]] <- NA_real_
  nonfinite_payload_result <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = nonfinite_payload
  )
  expect_identical(
    nonfinite_payload_result$availability,
    c("inference", "unavailable", "unavailable")
  )

  wrong_source <- sidecar
  wrong_source$coefficients$coefficient_source[
    wrong_source$coefficients$gene_id == gid[[2L]]
  ] <- "untrusted_source"
  wrong_source_result <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = wrong_source
  )
  expect_identical(
    wrong_source_result$availability,
    c("inference", "unavailable", "unavailable")
  )

  wrong_status <- sidecar
  wrong_status$coefficients$fallback_status[
    wrong_status$coefficients$gene_id == gid[[2L]]
  ] <- "attempted_failed"
  wrong_status_result <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = wrong_status
  )
  expect_identical(
    wrong_status_result$availability,
    c("inference", "unavailable", "unavailable")
  )

  inference_claim <- sidecar
  inference_claim$coefficients$inference_supported[
    inference_claim$coefficients$gene_id == gid[[2L]]
  ] <- TRUE
  inference_claim_result <- ln_availability_build(
    primary,
    gid,
    gene_names,
    operational_fallback = inference_claim
  )
  expect_identical(
    inference_claim_result$availability,
    c("inference", "unavailable", "unavailable")
  )
})

test_that("strict reason groups use shadow diagnostics without enabling them", {
  primary <- list(
    ln_availability_test_primary(
      failure_reason = "competing_uncertified_candidate",
      resolution = list(
        strict_failure_reason = "competing_uncertified_candidate",
        would_nonveto = TRUE,
        would_route = "R0",
        global_veto_classes = character(0)
      )
    ),
    ln_availability_test_primary(
      failure_reason = "competing_uncertified_candidate",
      resolution = list(global_veto_classes = "competitive_boundary")
    ),
    ln_availability_test_primary(
      failure_reason = "competing_uncertified_candidate",
      resolution = list(global_veto_classes = "competitive_distinct_basin")
    ),
    ln_availability_test_primary(
      failure_reason = "derivative_audit_failed"
    )
  )
  availability <- ln_availability_build(
    primary,
    gid = 1:4,
    gene_names = paste0("gene_", 1:4)
  )

  expect_true(all(availability$availability == "unavailable"))
  expect_true(all(!availability$inference_supported))
  expect_true(all(is.na(availability$estimate_source)))
  expect_identical(availability$reason_group, c(
    "certification_failure",
    "competitive_boundary",
    "competitive_solutions",
    "hard_stop"
  ))
  expect_match(
    availability$status_message[[2L]],
    "boundary solution"
  )
  expect_match(
    availability$status_message[[3L]],
    "materially different solutions"
  )
})

test_that("availability warning is consolidated and preserves other warnings", {
  availability <- data.frame(
    availability = c("inference", "estimate_only", "unavailable"),
    attention_flag = c(FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  captured <- ln_availability_capture_warnings({
    warning("Unrelated input warning.", call. = FALSE)
    ln_availability_warn(
      availability,
      operational_fallback_requested = TRUE
    )
  })

  expect_length(captured$messages, 2L)
  availability_messages <- captured$messages[startsWith(
    captured$messages,
    "Joint-LN inference was unavailable for"
  )]
  expect_length(availability_messages, 1L)
  expect_match(availability_messages, "unavailable for 2 genes")
  expect_match(availability_messages, "inference is available for 1 of 3")
  expect_match(availability_messages, "coefficients are available for 1")
  expect_match(availability_messages, "coefficients remain unavailable for 1")
  expect_match(availability_messages, "provides no covariance")
  expect_true("Unrelated input warning." %in% captured$messages)

  primary_only <- ln_availability_capture_warnings(
    ln_availability_warn(availability)
  )
  expect_length(primary_only$messages, 1L)
  expect_false(grepl(
    "legacy post-PML",
    primary_only$messages,
    fixed = TRUE
  ))

  expect_no_warning(ln_availability_warn(data.frame(
    availability = rep("inference", 2L),
    attention_flag = rep(FALSE, 2L),
    stringsAsFactors = FALSE
  )))
})

test_that("availability stays joint-LN scoped", {
  joint <- expect_no_warning(ln_joint_fit_fixture()$fit)
  expect_identical(joint$availability$availability, "inference")
  expect_true(joint$availability$inference_supported)

  fx <- ln_joint_test_fixture()
  legacy <- scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    keep_diagnostics = FALSE,
    ln_estimator = "legacy_post_pml"
  )
  expect_false("availability" %in% names(legacy))
})

test_that("compact descriptive provenance does not require raw diagnostics", {
  fx <- ln_joint_test_fixture()
  fit <- scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    keep_diagnostics = FALSE,
    ln_estimator = "joint_ln"
  )
  compact <- fit$provenance$per_gene
  required <- c(
    "estimator", "certified_root", "inference_status",
    "strict_inference_status", "enabled_inference_status",
    "resolution_status", "inference_source", "newly_enabled",
    "availability_transition", "resolution_veto_classes",
    "resolution_policy_version", "resolution_mode", "route_code",
    "strict_route_code", "enabled_route_code", "strict_failure_reason",
    "enabled_failure_reason"
  )

  expect_true(all(required %in% names(compact)))
  expect_identical(
    fit$provenance$candidate_resolution_policy,
    unique(compact$resolution_policy_version)
  )
  expect_identical(
    fit$provenance$candidate_resolution_mode,
    unique(compact$resolution_mode)
  )
  expect_identical(unserialize(serialize(compact, NULL)), compact)
  expect_identical(
    fit$availability$inference_status,
    compact$inference_status
  )
  expect_identical(
    fit$availability$resolution_status,
    compact$resolution_status
  )
})

test_that("joint-LN output retains typed missing inference and expert detail", {
  captured <- ln_availability_capture_warnings(
    ln_joint_fit_fixture(max = c(0.01, 1000))$fit
  )
  fit <- captured$value

  expect_s3_class(fit$availability, "data.frame")
  expect_identical(nrow(fit$availability), nrow(fit$overdispersion))
  expect_identical(fit$availability$availability, "unavailable")
  expect_identical(fit$availability$estimate_source, NA_character_)
  expect_length(grep("^p", names(fit$availability)), 0L)

  inference_columns <- grep(
    "^se_|^p_naive_|^p_robust_|^wald_stat_",
    names(fit$summary),
    value = TRUE
  )
  expect_true(length(inference_columns) > 0L)
  inference_values <- as.matrix(fit$summary[inference_columns])
  expect_true(all(is.na(inference_values)))
  expect_false(any(inference_values == 1, na.rm = TRUE))
  expect_true(is.list(fit$diagnostics))
  expect_true(is.list(fit$diagnostics$raw))
  expect_s3_class(fit$provenance$per_gene, "data.frame")

  availability_messages <- captured$messages[startsWith(
    captured$messages,
    "Joint-LN inference was unavailable for"
  )]
  expect_length(availability_messages, 1L)
})
