# Estimator-aligned joint-LN fitting ---------------------------------------

ln_joint_thresholds <- function() {
  list(
    donor_additivity = 1e-10,
    objective_gradient = 1e-5,
    sensitivity_score = 1e-4,
    finite_difference_stability = 1e-4,
    root = 1e-6,
    newton = 1e-6,
    parameter_agreement = 1e-6,
    objective_agreement = 1e-8,
    constraint = 1e-10,
    curvature_min_eigen = 1e-8,
    condition_number = 1e10,
    covariance_identity = 1e-8,
    natural_log_coordinate = 1e-7,
    internal_public_coordinate = 1e-7,
    mode_score = 1e-10,
    mode_newton = 1e-10,
    hard_stop = list(
      donor_additivity = 1e-8,
      objective_gradient = 1e-4,
      sensitivity_score = 1e-3,
      finite_difference_stability = 1e-3,
      root = 1e-4,
      newton = 1e-4,
      parameter_agreement = 1e-4,
      objective_agreement = 1e-6,
      constraint = 1e-8,
      condition_number = 1e12,
      covariance_identity = 1e-7,
      natural_log_coordinate = 1e-6,
      internal_public_coordinate = 1e-6,
      mode_score = 1e-8,
      mode_newton = 1e-8
    )
  )
}

ln_joint_candidate_check_names <- function() {
  c(
    "donor_additivity", "root", "newton", "curvature", "condition",
    "interior", "constraint"
  )
}

# Shadow-only candidate-resolution thresholds.  The existing selected-root
# certificate above is deliberately not changed.  Review-zone thresholds reuse
# its pre-existing hard-stop scales; a proposed basin map must return to the
# stricter selected-root scales before it can be recorded as a prospective
# nonveto.
ln_joint_resolution_thresholds <- function() {
  selected <- ln_joint_thresholds()
  list(
    schema_version = 1L,
    policy_version = "joint_ln_candidate_resolution_shadow_v1",
    objective_equivalence = selected$hard_stop$objective_agreement,
    review_full_state = selected$hard_stop$parameter_agreement,
    review_beta_block = selected$hard_stop$parameter_agreement,
    review_dispersion_block = selected$hard_stop$parameter_agreement,
    boundary_scaled_distance = selected$root,
    curvature_compatibility =
      selected$hard_stop$finite_difference_stability,
    proposal_full_state = selected$parameter_agreement,
    proposal_beta_block = selected$parameter_agreement,
    proposal_dispersion_block = selected$parameter_agreement,
    proposal_objective_equivalence = selected$objective_agreement,
    proposal_root = selected$root,
    proposal_newton = selected$newton,
    contraction_ratio = 0.5,
    max_state_evaluations = 1L
  )
}

ln_joint_resolution_thresholds_valid <- function(x) {
  expected <- ln_joint_resolution_thresholds()
  required <- c(
    "schema_version", "policy_version", "objective_equivalence",
    "review_full_state", "review_beta_block", "review_dispersion_block",
    "boundary_scaled_distance", "curvature_compatibility",
    "proposal_full_state", "proposal_beta_block",
    "proposal_dispersion_block", "proposal_objective_equivalence",
    "proposal_root", "proposal_newton", "contraction_ratio",
    "max_state_evaluations"
  )
  if (!is.list(x) || !identical(names(x), required) ||
      !identical(x, expected) ||
      !identical(x$schema_version, 1L) ||
      length(x$policy_version) != 1L || is.na(x$policy_version) ||
      !nzchar(x$policy_version)) {
    return(FALSE)
  }
  numeric_names <- setdiff(
    required,
    c("schema_version", "policy_version", "max_state_evaluations")
  )
  numeric_fields_valid <- vapply(x[numeric_names], function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value)
  }, logical(1))
  if (!all(numeric_fields_valid)) {
    return(FALSE)
  }
  numeric_values <- unlist(x[numeric_names], use.names = FALSE)
  all(is.finite(numeric_values)) && all(numeric_values > 0) &&
    x$contraction_ratio < 1 &&
    identical(x$max_state_evaluations, 1L)
}

ln_joint_objective_provenance <- function() {
  list(
    id = "joint_ln_ptmg_working_objective_v1",
    estimator = "joint_ln",
    criterion = "ptmg_ll",
    derivative_criterion = "ptmg_ll_der_hes4",
    objective_convention = "minimized_q_equals_negative_working_log_objective",
    score_convention = "maximizing_U_equals_negative_gradient_of_q",
    objective_coordinates = "beta, sigma2, gamma_nb_size",
    derivative_coordinates = "beta, log(sigma2), log(gamma_nb_size)",
    donor_unit = "independent_subject",
    subject_additive = TRUE,
    corrected_zero_based_positive_cell_index = TRUE
  )
}

ln_joint_objective_instance_registry <- local({
  registry <- new.env(parent = emptyenv())
  counter <- 0
  max_entries <- 16L

  matching_keys <- function(gene_index) {
    if (!ln_joint_positive_integer_scalar(gene_index)) {
      return(character(0))
    }
    keys <- ls(registry, all.names = TRUE)
    if (!length(keys)) {
      return(character(0))
    }
    matches <- vapply(keys, function(key) {
      entry <- get(key, envir = registry, inherits = FALSE)
      identical(entry$token$gene_index, as.integer(gene_index))
    }, logical(1))
    keys[matches]
  }

  matching_key <- function(gene_index) {
    keys <- matching_keys(gene_index)
    if (!length(keys)) {
      return(NULL)
    }
    sequences <- vapply(keys, function(key) {
      get(key, envir = registry, inherits = FALSE)$sequence
    }, numeric(1))
    keys[[which.max(sequences)]]
  }

  function(action = c(
    "register",
    "lookup",
    "store_unrestricted",
    "get_unrestricted",
    "remove_matching",
    "size",
    "clear"
  ), gene_index = NULL, posv = NULL, ctx = NULL, fit = NULL,
  instance = NULL) {
    action <- match.arg(action)
    if (identical(action, "size")) {
      return(length(ls(registry, all.names = TRUE)))
    }
    if (identical(action, "clear")) {
      keys <- ls(registry, all.names = TRUE)
      if (length(keys)) {
        rm(list = keys, envir = registry)
      }
      return(invisible(TRUE))
    }
    if (!ln_joint_positive_integer_scalar(gene_index)) {
      return(NULL)
    }
    if (identical(action, "remove_matching")) {
      keys <- matching_keys(gene_index)
      if (length(keys)) {
        rm(list = keys, envir = registry)
      }
      return(invisible(TRUE))
    }
    if (identical(action, "register")) {
      counter <<- counter + 1
      key <- paste0("objective_instance_", Sys.getpid(), "_", counter)
      token <- list(
        scheme = "process_local_exact_joint_ln_objective_context_v1",
        id = paste(
          Sys.getpid(),
          format(counter, scientific = FALSE),
          sep = "|"
        ),
        gene_index = as.integer(gene_index)
      )
      assign(
        key,
        list(
          token = token,
          sequence = counter,
          unrestricted_fit = NULL
        ),
        envir = registry
      )
      existing <- ls(registry, all.names = TRUE)
      if (length(existing) > max_entries) {
        sequences <- vapply(existing, function(existing_key) {
          get(existing_key, envir = registry, inherits = FALSE)$sequence
        }, numeric(1))
        remove_keys <- existing[order(sequences)][seq_len(
          length(existing) - max_entries
        )]
        rm(list = remove_keys, envir = registry)
      }
      return(token)
    }
    key <- matching_key(gene_index)
    if (is.null(key)) {
      return(NULL)
    }
    entry <- get(key, envir = registry, inherits = FALSE)
    if (identical(action, "lookup")) {
      return(entry$token)
    }
    if (identical(action, "store_unrestricted")) {
      if (!is.list(fit)) {
        return(NULL)
      }
      entry$unrestricted_fit <- fit
      assign(key, entry, envir = registry)
      return(entry$token)
    }
    entry$unrestricted_fit
  }
})

ln_joint_scaled_max <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (!length(a) || length(a) != length(b) || any(!is.finite(c(a, b)))) {
    return(NA_real_)
  }
  max(abs(a - b) / (1 + pmax(abs(a), abs(b))))
}

ln_joint_relative_frobenius <- function(a, b) {
  if (is.null(a) || is.null(b) || !all(dim(a) == dim(b)) ||
      any(!is.finite(c(a, b)))) {
    return(NA_real_)
  }
  as.numeric(norm(a - b, "F") / max(norm(b, "F"), .Machine$double.eps))
}

ln_joint_real_vector <- function(value, length = NULL, finite = FALSE) {
  valid <- is.atomic(value) && typeof(value) %in% c("integer", "double") &&
    is.null(dim(value))
  if (!valid || (!is.null(length) && base::length(value) != length)) {
    return(FALSE)
  }
  !isTRUE(finite) || all(is.finite(value))
}

ln_joint_real_matrix <- function(value, finite = FALSE) {
  valid <- is.matrix(value) &&
    typeof(value) %in% c("integer", "double")
  if (!valid) {
    return(FALSE)
  }
  !isTRUE(finite) || all(is.finite(value))
}

ln_joint_positive_integer_scalar <- function(value, allow_zero = FALSE) {
  ln_joint_real_vector(value, length = 1L, finite = TRUE) &&
    isTRUE(value == as.integer(value)) &&
    if (isTRUE(allow_zero)) value >= 0 else value >= 1
}

ln_joint_scalar_logical <- function(value) {
  is.logical(value) && length(value) == 1L && !is.na(value)
}

ln_joint_resolution_candidate_key <- function(candidate) {
  if (!is.list(candidate) || is.data.frame(candidate)) {
    return("joint_ln_candidate_v1|malformed_nonlist_candidate")
  }
  numeric_key <- function(x) {
    if (is.null(x) || !length(x)) {
      return("<missing>")
    }
    if (!ln_joint_real_vector(x)) {
      return("<malformed>")
    }
    paste(format(as.numeric(x), digits = 17, scientific = TRUE), collapse = ";")
  }
  paste(
    "joint_ln_candidate_v1",
    isTRUE(candidate$required),
    isTRUE(candidate$core),
    isTRUE(candidate$finite_endpoint),
    isTRUE(candidate$certified),
    numeric_key(candidate$objective),
    numeric_key(candidate$theta_full_log),
    sep = "|"
  )
}

ln_joint_finite_scalar_objective <- function(candidate) {
  is.list(candidate) && !is.data.frame(candidate) &&
    ln_joint_real_vector(candidate$objective, length = 1L, finite = TRUE)
}

ln_joint_candidate_label <- function(candidate, fallback = NA_character_) {
  if (!is.list(candidate) || is.data.frame(candidate) ||
      is.null(candidate$label) ||
      !is.atomic(candidate$label) ||
      !length(candidate$label)) {
    return(fallback)
  }
  label <- as.character(candidate$label[[1L]])
  if (length(label) != 1L || is.na(label) || !nzchar(label)) {
    fallback
  } else {
    label
  }
}

ln_joint_diagnostic_character <- function(value,
                                          fallback = NA_character_) {
  if (is.null(value) || !is.atomic(value) || !length(value)) {
    return(fallback)
  }
  result <- as.character(value[[1L]])
  if (length(result) != 1L || is.na(result)) fallback else result
}

ln_joint_diagnostic_numeric <- function(value, fallback = NA_real_) {
  if (!ln_joint_real_vector(value) || !length(value)) {
    return(fallback)
  }
  as.numeric(value[[1L]])
}

ln_joint_diagnostic_integer <- function(value, fallback = NA_integer_) {
  if (!ln_joint_real_vector(value) || !length(value) ||
      !is.finite(value[[1L]])) {
    return(fallback)
  }
  as.integer(value[[1L]])
}

ln_joint_diagnostic_vector_string <- function(value,
                                              fallback = NA_character_) {
  if (!ln_joint_real_vector(value)) {
    return(fallback)
  }
  paste(format(value, digits = 17), collapse = ";")
}

ln_joint_probe_map_record <- function(map, n_beta) {
  list(
    schema_version = 1L,
    beta_star = as.numeric(map$beta_star),
    null_basis = unname(as.matrix(map$K)),
    n_beta = as.integer(n_beta),
    free_dim = as.integer(map$free_dim)
  )
}

ln_joint_probe_map_valid <- function(probe_map, n_beta,
                                     reduced_dimension = NULL) {
  if (!is.list(probe_map) ||
      !identical(
        names(probe_map),
        c("schema_version", "beta_star", "null_basis", "n_beta", "free_dim")
      ) ||
      !identical(probe_map$schema_version, 1L) ||
      !ln_joint_positive_integer_scalar(n_beta) ||
      !identical(probe_map$n_beta, as.integer(n_beta)) ||
      !ln_joint_positive_integer_scalar(probe_map$free_dim, allow_zero = TRUE) ||
      !is.integer(probe_map$free_dim) ||
      !ln_joint_real_vector(
        probe_map$beta_star,
        length = as.integer(n_beta),
        finite = TRUE
      ) ||
      !ln_joint_real_matrix(probe_map$null_basis, finite = TRUE) ||
      !identical(
        dim(probe_map$null_basis),
        c(as.integer(n_beta), probe_map$free_dim)
      )) {
    return(FALSE)
  }
  if (!is.null(reduced_dimension) &&
      (!ln_joint_positive_integer_scalar(reduced_dimension) ||
        as.integer(reduced_dimension) != probe_map$free_dim + 2L)) {
    return(FALSE)
  }
  TRUE
}

ln_joint_probe_full_log <- function(theta_reduced_log, candidate) {
  if (!ln_joint_real_vector(theta_reduced_log, finite = TRUE) ||
      !is.list(candidate) || is.data.frame(candidate) ||
      !ln_joint_probe_map_valid(
        candidate$probe_map,
        candidate$n_beta,
        length(theta_reduced_log)
      )) {
    return(NULL)
  }
  theta_reduced_log <- as.numeric(theta_reduced_log)
  probe_map <- candidate$probe_map
  free_dim <- probe_map$free_dim
  zeta <- if (free_dim) {
    theta_reduced_log[seq_len(free_dim)]
  } else {
    numeric(0)
  }
  beta <- if (free_dim) {
    as.numeric(probe_map$beta_star + probe_map$null_basis %*% zeta)
  } else {
    as.numeric(probe_map$beta_star)
  }
  c(beta, theta_reduced_log[free_dim + seq_len(2L)])
}

ln_joint_objective_only_endpoint_record <- function(candidate) {
  if (!is.list(candidate) || is.data.frame(candidate) ||
      !is.null(candidate$state)) {
    return(FALSE)
  }
  derivative_fields <- c(
    "theta_full_log", "checks", "donor_additivity_scaled",
    "root_residual", "newton_displacement", "curvature", "interior",
    "constraint_residual"
  )
  all(vapply(derivative_fields, function(field) {
    is.null(candidate[[field]])
  }, logical(1)))
}

ln_joint_optional_failed_attempt_integrity <- function(candidate,
                                                       selected = NULL) {
  if (!is.list(candidate) || is.data.frame(candidate) ||
      !identical(candidate$required, FALSE) ||
      !identical(candidate$core, FALSE) ||
      !identical(candidate$finite_endpoint, FALSE) ||
      !identical(candidate$certified, FALSE) ||
      !is.null(candidate$state)) {
    return(FALSE)
  }
  n_beta <- candidate$n_beta
  if (!ln_joint_positive_integer_scalar(n_beta) ||
      !identical(
        candidate$parameterization,
        "beta, log(sigma2), log(gamma_nb_size)"
      ) ||
      !identical(
        candidate$objective_provenance,
        ln_joint_objective_provenance()
      ) ||
      is.null(candidate$objective_instance)) {
    return(FALSE)
  }
  bounds_ok <- is.list(candidate$bounds) &&
    identical(names(candidate$bounds), c("lower", "upper")) &&
    ln_joint_real_vector(candidate$bounds$lower, length = 2L, finite = TRUE) &&
    ln_joint_real_vector(candidate$bounds$upper, length = 2L, finite = TRUE) &&
    all(candidate$bounds$lower > 0) &&
    all(candidate$bounds$upper > candidate$bounds$lower)
  probe_bounds_ok <- is.list(candidate$probe_bounds) &&
    identical(names(candidate$probe_bounds), c("lower", "upper")) &&
    ln_joint_real_vector(candidate$probe_bounds$lower, finite = TRUE) &&
    length(candidate$probe_bounds$lower) >= 2L &&
    ln_joint_real_vector(candidate$probe_bounds$upper, finite = TRUE) &&
    length(candidate$probe_bounds$upper) ==
      length(candidate$probe_bounds$lower) &&
    all(candidate$probe_bounds$upper > candidate$probe_bounds$lower)
  probe_dimension <- if (is.list(candidate$probe_bounds)) {
    length(candidate$probe_bounds$lower)
  } else {
    NA_integer_
  }
  probe_map_ok <- ln_joint_probe_map_valid(
    candidate$probe_map,
    candidate$n_beta,
    probe_dimension
  )
  if (!bounds_ok || !probe_bounds_ok || !probe_map_ok) {
    return(FALSE)
  }
  if (!is.null(selected)) {
    if (!is.list(selected) ||
        !ln_joint_positive_integer_scalar(selected$n_beta) ||
        !identical(as.integer(candidate$n_beta),
                   as.integer(selected$n_beta)) ||
        !identical(candidate$parameterization, selected$parameterization) ||
        !identical(candidate$objective_provenance,
                   selected$objective_provenance) ||
        !identical(candidate$objective_instance,
                   selected$objective_instance) ||
        !identical(candidate$bounds, selected$bounds) ||
        !identical(candidate$probe_bounds, selected$probe_bounds) ||
        !identical(candidate$probe_map, selected$probe_map)) {
      return(FALSE)
    }
  }
  TRUE
}

ln_joint_objective_comparison <- function(candidate, selected, tolerance) {
  candidate <- as.numeric(candidate)
  selected <- as.numeric(selected)
  if (length(candidate) != 1L || length(selected) != 1L ||
      !is.finite(candidate) || !is.finite(selected) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    return(list(
      valid = FALSE,
      relation = "not_comparable",
      signed_gap = NA_real_,
      absolute_gap = NA_real_
    ))
  }
  signed_gap <- (candidate - selected) /
    (1 + max(abs(candidate), abs(selected)))
  relation <- if (signed_gap > tolerance) {
    "dominated"
  } else if (signed_gap < -tolerance) {
    "materially_better"
  } else {
    "equivalent"
  }
  list(
    valid = TRUE,
    relation = relation,
    signed_gap = signed_gap,
    absolute_gap = abs(signed_gap)
  )
}

ln_joint_resolution_block_metrics <- function(candidate, selected, n_beta) {
  candidate <- as.numeric(candidate)
  selected <- as.numeric(selected)
  if (length(n_beta) != 1L || !is.finite(n_beta) || n_beta < 1L) {
    return(list(
      full_state = NA_real_, beta_block = NA_real_,
      dispersion_block = NA_real_
    ))
  }
  n_beta <- as.integer(n_beta)
  expected <- n_beta + 2L
  if (length(candidate) != expected || length(selected) != expected ||
      any(!is.finite(c(candidate, selected)))) {
    return(list(
      full_state = NA_real_, beta_block = NA_real_,
      dispersion_block = NA_real_
    ))
  }
  list(
    full_state = ln_joint_scaled_max(candidate, selected),
    beta_block = ln_joint_scaled_max(
      candidate[seq_len(n_beta)],
      selected[seq_len(n_beta)]
    ),
    dispersion_block = ln_joint_scaled_max(
      candidate[n_beta + seq_len(2L)],
      selected[n_beta + seq_len(2L)]
    )
  )
}

ln_joint_scaled_bound_distance <- function(full_log, n_beta, bounds) {
  full_log <- as.numeric(full_log)
  if (!is.list(bounds) || !identical(names(bounds), c("lower", "upper")) ||
      length(bounds$lower) != 2L || length(bounds$upper) != 2L ||
      any(!is.finite(c(bounds$lower, bounds$upper))) ||
      any(bounds$lower <= 0) || any(bounds$upper <= bounds$lower) ||
      length(n_beta) != 1L || !is.finite(n_beta)) {
    return(NA_real_)
  }
  n_beta <- as.integer(n_beta)
  if (length(full_log) != n_beta + 2L || any(!is.finite(full_log))) {
    return(NA_real_)
  }
  dispersion <- exp(full_log[n_beta + seq_len(2L)])
  if (any(!is.finite(dispersion))) {
    return(NA_real_)
  }
  lower <- as.numeric(bounds$lower)
  upper <- as.numeric(bounds$upper)
  lower_distance <- (dispersion - lower) /
    (1 + pmax(abs(dispersion), abs(lower)))
  upper_distance <- (upper - dispersion) /
    (1 + pmax(abs(dispersion), abs(upper)))
  min(c(lower_distance, upper_distance))
}

ln_joint_scaled_curvature <- function(A, scale = NULL) {
  A <- as.matrix(A)
  if (!length(A) || nrow(A) != ncol(A) || any(!is.finite(A))) {
    return(NULL)
  }
  A <- 0.5 * (A + t(A))
  if (is.null(scale)) {
    scale <- sqrt(pmax(abs(diag(A)), .Machine$double.eps))
  } else {
    scale <- as.numeric(scale)
    if (length(scale) != nrow(A) || any(!is.finite(scale)) ||
        any(scale <= 0)) {
      return(NULL)
    }
  }
  A / outer(scale, scale)
}

ln_joint_curvature_compatibility <- function(candidate, selected) {
  selected_matrix <- as.matrix(selected)
  if (!length(selected_matrix) ||
      nrow(selected_matrix) != ncol(selected_matrix) ||
      any(!is.finite(selected_matrix))) {
    return(NA_real_)
  }
  selected_matrix <- 0.5 * (selected_matrix + t(selected_matrix))
  common_scale <- sqrt(pmax(
    abs(diag(selected_matrix)),
    .Machine$double.eps
  ))
  candidate <- ln_joint_scaled_curvature(candidate, common_scale)
  selected <- ln_joint_scaled_curvature(selected_matrix, common_scale)
  if (is.null(candidate) || is.null(selected) ||
      !identical(dim(candidate), dim(selected))) {
    return(NA_real_)
  }
  as.numeric(
    norm(candidate - selected, "F") /
      (1 + max(norm(candidate, "F"), norm(selected, "F")))
  )
}

ln_joint_complete_candidate_checks <- function(checks) {
  expected <- ln_joint_candidate_check_names()
  if (is.null(checks) || !is.atomic(checks) ||
      !identical(names(checks), expected) || length(checks) != length(expected)) {
    return(list(valid = FALSE, values = NULL, failed = character(0)))
  }
  values <- as.logical(checks)
  if (any(is.na(values))) {
    return(list(valid = FALSE, values = NULL, failed = character(0)))
  }
  names(values) <- expected
  list(valid = TRUE, values = values, failed = expected[!values])
}

ln_joint_candidate_objective_integrity <- function(candidate, selected) {
  fail <- function(reason) list(valid = FALSE, reason = reason)
  if (!is.list(candidate) || is.data.frame(candidate) ||
      !is.list(selected) || is.data.frame(selected)) {
    return(fail("candidate_not_a_list"))
  }
  role_fields <- c("required", "core", "finite_endpoint", "certified")
  if (!all(vapply(candidate[role_fields], ln_joint_scalar_logical,
                  logical(1)))) {
    return(fail("invalid_candidate_role_flags"))
  }
  n_beta <- candidate$n_beta
  selected_n_beta <- selected$n_beta
  if (!ln_joint_positive_integer_scalar(n_beta) ||
      !ln_joint_positive_integer_scalar(selected_n_beta) ||
      as.integer(selected_n_beta) != as.integer(n_beta)) {
    return(fail("beta_dimension_mismatch"))
  }
  if (!identical(
      candidate$parameterization,
      "beta, log(sigma2), log(gamma_nb_size)")) {
    return(fail("parameterization_mismatch"))
  }
  expected_provenance <- ln_joint_objective_provenance()
  if (!identical(candidate$objective_provenance, expected_provenance)) {
    return(fail("objective_provenance_mismatch"))
  }
  if (is.null(candidate$objective_instance) ||
      !identical(candidate$objective_provenance,
                 selected$objective_provenance) ||
      !identical(candidate$objective_instance,
                 selected$objective_instance)) {
    return(fail("fit_provenance_mismatch"))
  }
  bounds_ok <- is.list(candidate$bounds) &&
    identical(names(candidate$bounds), c("lower", "upper")) &&
    ln_joint_real_vector(candidate$bounds$lower, length = 2L, finite = TRUE) &&
    ln_joint_real_vector(candidate$bounds$upper, length = 2L, finite = TRUE) &&
    all(candidate$bounds$lower > 0) &&
    all(candidate$bounds$upper > candidate$bounds$lower) &&
    identical(candidate$bounds, selected$bounds)
  if (!bounds_ok) {
    return(fail("fit_coordinate_or_bounds_mismatch"))
  }
  if (!ln_joint_real_vector(candidate$theta_reduced, finite = TRUE) ||
      !ln_joint_real_vector(candidate$theta_reduced_log, finite = TRUE) ||
      !ln_joint_real_vector(selected$theta_reduced, finite = TRUE)) {
    return(fail("invalid_objective_endpoint"))
  }
  theta_natural <- as.numeric(candidate$theta_reduced)
  theta_log <- as.numeric(candidate$theta_reduced_log)
  selected_reduced_dimension <- length(selected$theta_reduced)
  expected_theta_log <- c(
    head(theta_natural, -2L),
    log(tail(theta_natural, 2L))
  )
  if (selected_reduced_dimension < 2L ||
      length(theta_natural) != selected_reduced_dimension ||
      length(theta_log) != length(theta_natural) ||
      any(tail(theta_natural, 2L) <= 0) ||
      !is.finite(ln_joint_scaled_max(theta_log, expected_theta_log)) ||
      ln_joint_scaled_max(theta_log, expected_theta_log) >
        16 * .Machine$double.eps) {
    return(fail("invalid_objective_endpoint"))
  }
  if (!ln_joint_probe_map_valid(
      candidate$probe_map,
      candidate$n_beta,
      length(theta_log)) ||
      !identical(candidate$probe_map, selected$probe_map)) {
    return(fail("fit_reduced_to_full_map_mismatch"))
  }
  if (!ln_joint_finite_scalar_objective(candidate)) {
    return(fail("objective_not_comparable"))
  }
  list(valid = TRUE, reason = NA_character_)
}

ln_joint_candidate_integrity <- function(
    candidate, selected = NULL,
    allow_nonfinite_failure_metrics = FALSE) {
  fail <- function(reason) list(valid = FALSE, reason = reason)
  if (!is.list(candidate) || is.data.frame(candidate)) {
    return(fail("candidate_not_a_list"))
  }
  role_fields <- c("required", "core", "finite_endpoint", "certified")
  if (!all(vapply(candidate[role_fields], ln_joint_scalar_logical,
                  logical(1)))) {
    return(fail("invalid_candidate_role_flags"))
  }
  n_beta <- candidate$n_beta
  if (!ln_joint_positive_integer_scalar(n_beta)) {
    return(fail("invalid_beta_dimension"))
  }
  n_beta <- as.integer(n_beta)
  if (!ln_joint_real_vector(
      candidate$theta_full_log,
      length = n_beta + 2L,
      finite = TRUE)) {
    return(fail("invalid_full_state"))
  }
  if (!identical(
      candidate$parameterization,
      "beta, log(sigma2), log(gamma_nb_size)")) {
    return(fail("parameterization_mismatch"))
  }
  expected_provenance <- ln_joint_objective_provenance()
  if (!identical(candidate$objective_provenance, expected_provenance)) {
    return(fail("objective_provenance_mismatch"))
  }
  if (is.null(candidate$objective_instance)) {
    return(fail("objective_instance_missing"))
  }
  bounds_ok <- is.list(candidate$bounds) &&
    identical(names(candidate$bounds), c("lower", "upper")) &&
    ln_joint_real_vector(candidate$bounds$lower, length = 2L, finite = TRUE) &&
    ln_joint_real_vector(candidate$bounds$upper, length = 2L, finite = TRUE) &&
    all(candidate$bounds$lower > 0) &&
    all(candidate$bounds$upper > candidate$bounds$lower)
  if (!bounds_ok) {
    return(fail("invalid_bounds"))
  }
  if (!ln_joint_real_vector(candidate$theta_reduced, finite = TRUE) ||
      !ln_joint_real_vector(candidate$theta_reduced_log, finite = TRUE) ||
      length(candidate$theta_reduced) < 2L ||
      length(candidate$theta_reduced_log) !=
        length(candidate$theta_reduced) ||
      any(tail(candidate$theta_reduced, 2L) <= 0)) {
    return(fail("invalid_objective_endpoint"))
  }
  expected_theta_log <- c(
    head(candidate$theta_reduced, -2L),
    log(tail(candidate$theta_reduced, 2L))
  )
  endpoint_coordinate_error <- ln_joint_scaled_max(
    candidate$theta_reduced_log,
    expected_theta_log
  )
  if (!is.finite(endpoint_coordinate_error) ||
      endpoint_coordinate_error > 16 * .Machine$double.eps) {
    return(fail("invalid_objective_endpoint"))
  }
  checks <- ln_joint_complete_candidate_checks(candidate$checks)
  if (!isTRUE(checks$valid)) {
    return(fail("incomplete_check_vector"))
  }
  state <- candidate$state
  if (!is.list(state) ||
      !ln_joint_real_matrix(state$A_reduced, finite = TRUE) ||
      nrow(state$A_reduced) != ncol(state$A_reduced) ||
      !ln_joint_real_vector(
        state$theta_reduced_log,
        length = nrow(state$A_reduced),
        finite = TRUE
      ) ||
      !ln_joint_real_vector(
        candidate$theta_reduced_log,
        length = nrow(state$A_reduced),
        finite = TRUE
      ) ||
      !identical(as.numeric(state$theta_reduced_log),
                 as.numeric(candidate$theta_reduced_log)) ||
      !ln_joint_real_vector(
        state$U_reduced,
        length = nrow(state$A_reduced),
        finite = TRUE
      ) ||
      !ln_joint_real_matrix(state$U_subject_reduced, finite = TRUE) ||
      nrow(state$U_subject_reduced) != nrow(state$A_reduced) ||
      ncol(state$U_subject_reduced) < 1L ||
      !ln_joint_real_vector(
        state$full_log,
        length = n_beta + 2L,
        finite = TRUE
      ) ||
      !identical(as.numeric(state$full_log),
                 as.numeric(candidate$theta_full_log)) ||
      !identical(state$objective_provenance,
                 candidate$objective_provenance) ||
      !is.list(state$raw) ||
      !ln_joint_real_vector(state$raw$value, length = 1L, finite = TRUE)) {
    return(fail("invalid_derivative_state"))
  }
  probe_bounds <- candidate$probe_bounds
  if (!is.list(probe_bounds) ||
      !identical(names(probe_bounds), c("lower", "upper")) ||
      !ln_joint_real_vector(
        probe_bounds$lower,
        length = nrow(state$A_reduced),
        finite = TRUE
      ) ||
      !ln_joint_real_vector(
        probe_bounds$upper,
        length = nrow(state$A_reduced),
        finite = TRUE
      ) ||
      any(probe_bounds$upper <= probe_bounds$lower)) {
    return(fail("invalid_probe_bounds"))
  }
  expected_full_log <- ln_joint_probe_full_log(
    candidate$theta_reduced_log,
    candidate
  )
  if (is.null(expected_full_log) ||
      !identical(as.numeric(candidate$theta_full_log),
                 as.numeric(expected_full_log))) {
    return(fail("reduced_to_full_map_mismatch"))
  }

  # Fail closed when a serialized check vector contradicts the immutable raw
  # endpoint record.  The prospective classifier must consume the complete
  # vector, but it must not trust hand-edited or partially serialized flags.
  selected_thresholds <- ln_joint_thresholds()
  raw_donor_additivity <- ln_joint_scaled_max(
    rowSums(state$U_subject_reduced),
    state$U_reduced
  )
  raw_root <- ln_joint_root_metric(
    state$U_subject_reduced,
    state$U_reduced
  )
  raw_newton <- ln_joint_probe_newton_metric(state)
  raw_curvature <- ln_joint_curvature_diagnostics(state$A_reduced)
  raw_dispersion <- exp(state$full_log[n_beta + seq_len(2L)])
  boundary_tolerance <- sqrt(.Machine$double.eps) *
    (1 + pmax(abs(candidate$bounds$lower), abs(candidate$bounds$upper)))
  raw_interior <- all(
    raw_dispersion > candidate$bounds$lower + boundary_tolerance &
      raw_dispersion < candidate$bounds$upper - boundary_tolerance
  )
  if (!is.list(candidate$curvature)) {
    return(fail("missing_or_nonfinite_check_metric"))
  }
  metric_scalars <- list(
    donor_additivity = candidate$donor_additivity_scaled,
    root = candidate$root_residual,
    newton = candidate$newton_displacement,
    curvature_min_eigen = candidate$curvature$scaled_min_eigenvalue,
    condition_number = candidate$curvature$scaled_condition_number,
    constraint = candidate$constraint_residual
  )
  finite_required_metrics <- list(
    donor_additivity = candidate$donor_additivity_scaled,
    root = candidate$root_residual,
    curvature_min_eigen = candidate$curvature$scaled_min_eigenvalue,
    constraint = candidate$constraint_residual
  )
  scalar_finite <- function(value) {
    ln_joint_real_vector(value, length = 1L, finite = TRUE)
  }
  scalar_missing <- function(value) {
    ln_joint_real_vector(value, length = 1L) && is.na(value)
  }
  scalar_positive_infinite <- function(value) {
    ln_joint_real_vector(value, length = 1L) && is.infinite(value) &&
      value > 0
  }
  newton_metric_valid <- scalar_finite(candidate$newton_displacement) ||
    (isTRUE(allow_nonfinite_failure_metrics) &&
      scalar_missing(candidate$newton_displacement))
  condition_metric_valid <- scalar_finite(
    candidate$curvature$scaled_condition_number
  ) || (isTRUE(allow_nonfinite_failure_metrics) &&
    scalar_positive_infinite(candidate$curvature$scaled_condition_number))
  raw_newton_valid <- scalar_finite(raw_newton) ||
    (isTRUE(allow_nonfinite_failure_metrics) && scalar_missing(raw_newton))
  raw_condition_valid <- scalar_finite(
    raw_curvature$scaled_condition_number
  ) || (isTRUE(allow_nonfinite_failure_metrics) &&
    scalar_positive_infinite(raw_curvature$scaled_condition_number))
  if (length(metric_scalars) != 6L ||
      !all(vapply(finite_required_metrics, scalar_finite, logical(1))) ||
      !isTRUE(newton_metric_valid) || !isTRUE(condition_metric_valid) ||
      !is.logical(candidate$interior) || length(candidate$interior) != 1L ||
      is.na(candidate$interior) ||
      !is.logical(candidate$curvature$positive_definite) ||
      length(candidate$curvature$positive_definite) != 1L ||
      is.na(candidate$curvature$positive_definite) ||
      !is.finite(raw_donor_additivity) || !is.finite(raw_root) ||
      !isTRUE(raw_newton_valid) ||
      !is.finite(raw_curvature$scaled_min_eigenvalue) ||
      !isTRUE(raw_condition_valid) ||
      length(raw_dispersion) != 2L || any(!is.finite(raw_dispersion))) {
    return(fail("missing_or_nonfinite_check_metric"))
  }
  raw_metrics_match <-
    identical(as.numeric(candidate$donor_additivity_scaled),
              as.numeric(raw_donor_additivity)) &&
    identical(as.numeric(candidate$root_residual), as.numeric(raw_root)) &&
    identical(as.numeric(candidate$newton_displacement),
              as.numeric(raw_newton)) &&
    identical(
      as.numeric(candidate$curvature$scaled_min_eigenvalue),
      as.numeric(raw_curvature$scaled_min_eigenvalue)
    ) &&
    identical(
      as.numeric(candidate$curvature$scaled_condition_number),
      as.numeric(raw_curvature$scaled_condition_number)
    ) &&
    identical(isTRUE(candidate$curvature$positive_definite),
              isTRUE(raw_curvature$positive_definite)) &&
    identical(isTRUE(candidate$interior), isTRUE(raw_interior))
  if (!isTRUE(raw_metrics_match)) {
    return(fail("check_metric_record_mismatch"))
  }
  expected_checks <- c(
    donor_additivity = is.finite(raw_donor_additivity) &&
      raw_donor_additivity <=
      selected_thresholds$donor_additivity,
    root = is.finite(raw_root) && raw_root <= selected_thresholds$root,
    newton = is.finite(raw_newton) &&
      raw_newton <= selected_thresholds$newton,
    curvature = isTRUE(raw_curvature$positive_definite) &&
      raw_curvature$scaled_min_eigenvalue >=
        selected_thresholds$curvature_min_eigen,
    condition = is.finite(raw_curvature$scaled_condition_number) &&
      raw_curvature$scaled_condition_number <=
      selected_thresholds$condition_number,
    interior = raw_interior,
    constraint = is.finite(candidate$constraint_residual) &&
      candidate$constraint_residual <=
      selected_thresholds$constraint
  )
  if (!identical(unname(checks$values), unname(expected_checks))) {
    return(fail("check_vector_metric_mismatch"))
  }
  objective_alignment <- ln_joint_objective_comparison(
    candidate$objective,
    state$raw$value,
    selected_thresholds$hard_stop$objective_agreement
  )
  if (!isTRUE(objective_alignment$valid) ||
      !identical(objective_alignment$relation, "equivalent")) {
    return(fail("objective_state_mismatch"))
  }
  if (!is.null(selected)) {
    selected_n_beta <- selected$n_beta
    if (!ln_joint_positive_integer_scalar(selected_n_beta) ||
        as.integer(selected_n_beta) != n_beta) {
      return(fail("beta_dimension_mismatch"))
    }
    if (!identical(candidate$objective_provenance,
                   selected$objective_provenance) ||
        !identical(candidate$objective_instance,
                   selected$objective_instance)) {
      return(fail("fit_provenance_mismatch"))
    }
    if (!identical(candidate$parameterization, selected$parameterization) ||
        !identical(candidate$bounds, selected$bounds) ||
        !identical(candidate$probe_map, selected$probe_map)) {
      return(fail("fit_coordinate_or_bounds_mismatch"))
    }
    selected_state <- selected$state
    if (!is.list(selected_state)) {
      return(fail("selected_derivative_state_invalid"))
    }
    selected_A <- selected_state$A_reduced
    selected_U_subject <- selected_state$U_subject_reduced
    if (!is.matrix(selected_A) ||
        !identical(dim(state$A_reduced), dim(selected_A)) ||
        !is.matrix(selected_U_subject) ||
        ncol(state$U_subject_reduced) != ncol(selected_U_subject)) {
      return(fail("curvature_dimension_mismatch"))
    }
  }
  list(
    valid = TRUE,
    reason = NA_character_,
    n_beta = n_beta,
    checks = checks$values,
    failed_checks = checks$failed
  )
}

ln_joint_candidate_resolution_record <- function(candidate) {
  candidate_list <- if (is.list(candidate) && !is.data.frame(candidate)) {
    candidate
  } else {
    list()
  }
  label <- ln_joint_candidate_label(candidate_list)
  list(
    candidate_key = ln_joint_resolution_candidate_key(candidate),
    label = label,
    required = isTRUE(candidate_list$required),
    core = isTRUE(candidate_list$core),
    class = "competitive_hard_stop",
    veto = TRUE,
    reason = "unclassified_candidate",
    failed_checks = character(0),
    objective_relation = "not_comparable",
    objective_signed_gap = NA_real_,
    objective_absolute_gap = NA_real_,
    metrics = list(
      full_state = NA_real_,
      beta_block = NA_real_,
      dispersion_block = NA_real_,
      boundary_scaled_distance = NA_real_,
      curvature_compatibility = NA_real_
    ),
    probe = NULL,
    probe_reused = FALSE
  )
}

ln_joint_finish_candidate_resolution <- function(record, class, veto, reason) {
  record$class <- class
  record$veto <- isTRUE(veto)
  record$reason <- reason
  record
}

ln_joint_classify_candidate <- function(candidate, selected,
                                        resolution_thresholds,
                                        is_selected = FALSE) {
  record <- ln_joint_candidate_resolution_record(candidate)
  if (!is.list(candidate) || is.data.frame(candidate)) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE,
      "candidate_not_a_list"
    ))
  }
  role_fields <- c("required", "core", "finite_endpoint", "certified")
  if (!all(vapply(candidate[role_fields], ln_joint_scalar_logical,
                  logical(1)))) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE,
      "invalid_candidate_role_flags"
    ))
  }
  finite_comparable <- isTRUE(candidate$finite_endpoint) &&
    ln_joint_finite_scalar_objective(candidate)
  if (!finite_comparable) {
    if (ln_joint_optional_failed_attempt_integrity(candidate, selected)) {
      return(ln_joint_finish_candidate_resolution(
        record, "noncompetitive_failed_attempt", FALSE,
        "optional_candidate_has_no_comparable_endpoint"
      ))
    }
    if (isTRUE(candidate$required) || isTRUE(candidate$certified) ||
        isTRUE(candidate$finite_endpoint)) {
      return(ln_joint_finish_candidate_resolution(
        record, "competitive_hard_stop", TRUE,
        if (isTRUE(candidate$required)) {
          "required_candidate_has_no_comparable_endpoint"
        } else {
          "malformed_finite_endpoint_objective"
        }
      ))
    }
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE,
      "malformed_optional_failed_attempt"
    ))
  }

  objective_integrity <- ln_joint_candidate_objective_integrity(
    candidate,
    selected
  )
  if (!isTRUE(objective_integrity$valid)) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE, objective_integrity$reason
    ))
  }
  objective <- ln_joint_objective_comparison(
    candidate$objective,
    selected$objective,
    resolution_thresholds$objective_equivalence
  )
  record$objective_relation <- objective$relation
  record$objective_signed_gap <- objective$signed_gap
  record$objective_absolute_gap <- objective$absolute_gap
  dominated_optional <- !isTRUE(is_selected) &&
    !isTRUE(candidate$certified) && !isTRUE(candidate$required) &&
    identical(objective$relation, "dominated")
  if (isTRUE(dominated_optional) &&
      ln_joint_objective_only_endpoint_record(candidate) &&
      identical(objective$relation, "dominated")) {
    return(ln_joint_finish_candidate_resolution(
      record, "dominated_uncertified", FALSE,
      "optional_candidate_objectively_dominated"
    ))
  }

  integrity <- ln_joint_candidate_integrity(
    candidate,
    selected,
    allow_nonfinite_failure_metrics = isTRUE(dominated_optional)
  )
  if (!isTRUE(integrity$valid)) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE, integrity$reason
    ))
  }
  record$failed_checks <- integrity$failed_checks
  blocks <- ln_joint_resolution_block_metrics(
    candidate$theta_full_log,
    selected$theta_full_log,
    integrity$n_beta
  )
  record$metrics$full_state <- blocks$full_state
  record$metrics$beta_block <- blocks$beta_block
  record$metrics$dispersion_block <- blocks$dispersion_block
  record$metrics$boundary_scaled_distance <-
    ln_joint_scaled_bound_distance(
      candidate$theta_full_log,
      integrity$n_beta,
      candidate$bounds
    )
  record$metrics$curvature_compatibility <-
    ln_joint_curvature_compatibility(
      candidate$state$A_reduced,
      selected$state$A_reduced
    )

  if (isTRUE(is_selected)) {
    if (!isTRUE(candidate$certified) || length(integrity$failed_checks)) {
      return(ln_joint_finish_candidate_resolution(
        record, "competitive_hard_stop", TRUE,
        "selected_candidate_certificate_inconsistent"
      ))
    }
    return(ln_joint_finish_candidate_resolution(
      record, "certified_selected", FALSE, "selected_certified_root"
    ))
  }
  if (isTRUE(candidate$certified)) {
    if (length(integrity$failed_checks)) {
      return(ln_joint_finish_candidate_resolution(
        record, "competitive_hard_stop", TRUE,
        "certified_candidate_check_vector_inconsistent"
      ))
    }
    return(ln_joint_finish_candidate_resolution(
      record, "certified_same_cluster", FALSE,
      "certified_candidate_handled_by_existing_consensus"
    ))
  }
  if (!isTRUE(objective$valid) ||
      any(!is.finite(unlist(record$metrics, use.names = FALSE)))) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE,
      "competitive_metrics_invalid"
    ))
  }
  if (isTRUE(candidate$required)) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE,
      "required_candidate_not_certified"
    ))
  }
  if (identical(objective$relation, "dominated")) {
    return(ln_joint_finish_candidate_resolution(
      record, "dominated_uncertified", FALSE,
      "optional_candidate_objectively_dominated"
    ))
  }
  if (record$metrics$boundary_scaled_distance <=
      resolution_thresholds$boundary_scaled_distance) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_boundary", TRUE,
      "objective_competitive_candidate_near_boundary"
    ))
  }
  if (identical(objective$relation, "materially_better")) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_distinct_basin", TRUE,
      "candidate_objective_materially_better_than_selected"
    ))
  }
  distinct <- record$metrics$full_state >
      resolution_thresholds$review_full_state ||
    record$metrics$beta_block > resolution_thresholds$review_beta_block ||
    record$metrics$dispersion_block >
      resolution_thresholds$review_dispersion_block ||
    record$metrics$curvature_compatibility >
      resolution_thresholds$curvature_compatibility
  if (distinct) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_distinct_basin", TRUE,
      "objective_competitive_candidate_has_distinct_state_or_curvature"
    ))
  }

  allowed_failures <- c("root", "newton")
  if (!length(integrity$failed_checks) ||
      any(!integrity$failed_checks %in% allowed_failures)) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE,
      "failed_checks_not_limited_to_root_and_newton"
    ))
  }
  selected_thresholds <- ln_joint_thresholds()
  root_review <- isTRUE(integrity$checks[["root"]]) ||
    (is.finite(candidate$root_residual) &&
      candidate$root_residual > selected_thresholds$root &&
      candidate$root_residual <= selected_thresholds$hard_stop$root)
  newton_review <- isTRUE(integrity$checks[["newton"]]) ||
    (is.finite(candidate$newton_displacement) &&
      candidate$newton_displacement > selected_thresholds$newton &&
      candidate$newton_displacement <= selected_thresholds$hard_stop$newton)
  if (!root_review || !newton_review) {
    return(ln_joint_finish_candidate_resolution(
      record, "competitive_hard_stop", TRUE,
      "root_or_newton_metric_outside_review_zone"
    ))
  }
  ln_joint_finish_candidate_resolution(
    record, "same_basin_review_candidate", TRUE,
    "eligible_for_bounded_same_basin_review"
  )
}

ln_joint_classify_candidate_against_anchors <- function(
    candidate, selected, assessment_roots, resolution_thresholds,
    is_selected = FALSE) {
  if (!is.list(candidate) || is.data.frame(candidate) ||
      isTRUE(is_selected) ||
      isTRUE(candidate$certified) ||
      !is.list(assessment_roots) || length(assessment_roots) <= 1L) {
    return(ln_joint_classify_candidate(
      candidate,
      selected,
      resolution_thresholds,
      is_selected = is_selected
    ))
  }
  records <- lapply(assessment_roots, function(anchor) {
    ln_joint_classify_candidate(
      candidate,
      anchor,
      resolution_thresholds,
      is_selected = FALSE
    )
  })
  priority <- c(
    certified_selected = 0L,
    certified_same_cluster = 0L,
    noncompetitive_failed_attempt = 1L,
    dominated_uncertified = 2L,
    same_basin_review_candidate = 3L,
    competitive_distinct_basin = 4L,
    competitive_boundary = 5L,
    competitive_hard_stop = 6L
  )
  scores <- vapply(records, function(record) {
    score <- unname(priority[record$class])
    if (length(score) != 1L || is.na(score)) 6L else as.integer(score)
  }, integer(1))
  records[[which.max(scores)]]
}

ln_joint_probe_equivalence_record <- function(candidate) {
  state <- if (is.list(candidate$state)) candidate$state else list()
  list(
    required = isTRUE(candidate$required),
    core = isTRUE(candidate$core),
    finite_endpoint = isTRUE(candidate$finite_endpoint),
    certified = isTRUE(candidate$certified),
    n_beta = candidate$n_beta,
    parameterization = candidate$parameterization,
    objective_provenance = candidate$objective_provenance,
    objective_instance = candidate$objective_instance,
    bounds = candidate$bounds,
    probe_bounds = candidate$probe_bounds,
    probe_map = candidate$probe_map,
    theta_reduced = candidate$theta_reduced,
    theta_reduced_log = candidate$theta_reduced_log,
    theta_full_log = candidate$theta_full_log,
    objective = candidate$objective,
    donor_additivity_scaled = candidate$donor_additivity_scaled,
    root_residual = candidate$root_residual,
    newton_displacement = candidate$newton_displacement,
    curvature = candidate$curvature,
    interior = candidate$interior,
    constraint_residual = candidate$constraint_residual,
    checks = candidate$checks,
    state = list(
      full_log = state$full_log,
      theta_reduced_log = state$theta_reduced_log,
      objective_provenance = state$objective_provenance,
      U_reduced = state$U_reduced,
      U_subject_reduced = state$U_subject_reduced,
      A_reduced = state$A_reduced,
      raw_value = if (is.list(state$raw)) state$raw$value else NULL
    )
  )
}

ln_joint_probe_equivalent <- function(candidate, reference) {
  identical(
    ln_joint_probe_equivalence_record(candidate),
    ln_joint_probe_equivalence_record(reference)
  )
}

ln_joint_basin_probe_template <- function(resolution_thresholds) {
  policy_version <- if (is.list(resolution_thresholds) &&
      is.character(resolution_thresholds$policy_version) &&
      length(resolution_thresholds$policy_version) == 1L &&
      !is.na(resolution_thresholds$policy_version)) {
    resolution_thresholds$policy_version
  } else {
    NA_character_
  }
  max_state_evaluations <- if (is.list(resolution_thresholds) &&
      ln_joint_positive_integer_scalar(
        resolution_thresholds$max_state_evaluations,
        allow_zero = TRUE
      )) {
    as.integer(resolution_thresholds$max_state_evaluations)
  } else {
    NA_integer_
  }
  list(
    schema_version = 1L,
    policy_version = policy_version,
    accepted = FALSE,
    termination_reason = "not_started",
    evaluation_budget = list(
      max_state_evaluations = max_state_evaluations
    ),
    state_evaluations = 0L,
    objective_evaluations = 0L,
    score_evaluations = 0L,
    hessian_evaluations = 0L,
    initial = NULL,
    proposed = NULL,
    contraction = NULL,
    anchor_assessments = list()
  )
}

# This proposal is deliberately endpoint-local.  It has no selected-root
# argument and therefore cannot steer, interpolate, snap, or replace the
# uncertified endpoint.
ln_joint_basin_local_proposal <- function(candidate, resolution_thresholds) {
  fail <- function(reason) list(success = FALSE, reason = reason)
  if (!ln_joint_resolution_thresholds_valid(resolution_thresholds) ||
      !is.list(candidate) || is.data.frame(candidate) ||
      !is.list(candidate$state)) {
    return(fail("invalid_probe_input"))
  }
  state <- candidate$state
  theta <- as.numeric(state$theta_reduced_log)
  U <- as.numeric(state$U_reduced)
  A <- as.matrix(state$A_reduced)
  bounds <- candidate$probe_bounds
  dimension_ok <- length(theta) > 0L && length(U) == length(theta) &&
    all(dim(A) == c(length(theta), length(theta))) &&
    is.list(bounds) && length(bounds$lower) == length(theta) &&
    length(bounds$upper) == length(theta)
  if (!dimension_ok || any(!is.finite(c(
      theta, U, A, bounds$lower, bounds$upper)))) {
    return(fail("invalid_local_derivative_state"))
  }
  step <- tryCatch(solve(A, U), error = function(e) NULL)
  if (is.null(step) || length(step) != length(theta) ||
      any(!is.finite(step))) {
    return(fail("local_newton_solve_failed"))
  }
  proposed <- theta + as.numeric(step)
  if (any(!is.finite(proposed))) {
    return(fail("nonfinite_local_proposal"))
  }
  if (any(proposed <= bounds$lower) || any(proposed >= bounds$upper)) {
    return(list(
      success = FALSE,
      reason = "local_proposal_outside_bounds",
      initial_theta_reduced_log = theta,
      local_step = as.numeric(step),
      proposed_theta_reduced_log = proposed
    ))
  }
  list(
    success = TRUE,
    reason = NA_character_,
    initial_theta_reduced_log = theta,
    local_step = as.numeric(step),
    proposed_theta_reduced_log = proposed
  )
}

ln_joint_probe_newton_metric <- function(state) {
  U <- as.numeric(state$U_reduced)
  A <- as.matrix(state$A_reduced)
  theta <- as.numeric(state$theta_reduced_log)
  if (!length(U) || length(theta) != length(U) ||
      !all(dim(A) == c(length(U), length(U))) ||
      any(!is.finite(c(U, A, theta)))) {
    return(NA_real_)
  }
  step <- tryCatch(solve(A, U), error = function(e) NULL)
  if (is.null(step) || any(!is.finite(step))) {
    return(NA_real_)
  }
  max(abs(step) / (1 + abs(theta)))
}

ln_joint_probe_contracted_or_target <- function(initial, proposed, target,
                                                ratio) {
  if (any(!is.finite(c(initial, proposed, target, ratio))) ||
      target <= 0 || ratio <= 0 || ratio >= 1 || initial < 0 || proposed < 0) {
    return(FALSE)
  }
  if (initial <= target) {
    proposed <= target
  } else {
    proposed <= target && proposed <= ratio * initial
  }
}

ln_joint_basin_probe <- function(candidate, selected, evaluate_state,
                                 resolution_thresholds =
                                   ln_joint_resolution_thresholds(),
                                 assessment_roots = NULL) {
  result <- ln_joint_basin_probe_template(resolution_thresholds)
  if (!ln_joint_resolution_thresholds_valid(resolution_thresholds) ||
      !is.function(evaluate_state)) {
    result$termination_reason <- "invalid_probe_configuration"
    return(result)
  }
  if (is.null(assessment_roots)) {
    assessment_roots <- list(selected)
  }
  if (!is.list(assessment_roots) || !length(assessment_roots) ||
      any(!vapply(assessment_roots, is.list, logical(1)))) {
    result$termination_reason <- "invalid_probe_assessment_roots"
    return(result)
  }
  proposal <- ln_joint_basin_local_proposal(candidate, resolution_thresholds)
  if (!isTRUE(proposal$success)) {
    result$termination_reason <- proposal$reason
    result$initial <- list(
      theta_reduced_log = proposal$initial_theta_reduced_log
    )
    result$proposed <- list(
      theta_reduced_log = proposal$proposed_theta_reduced_log,
      local_step = proposal$local_step
    )
    return(result)
  }

  # The certified state is first consulted only after the local proposal is
  # frozen.  It is used solely to assess the already proposed map.
  proposed_state <- tryCatch(
    evaluate_state(proposal$proposed_theta_reduced_log),
    error = function(e) NULL
  )
  result$state_evaluations <- 1L
  result$objective_evaluations <- 1L
  result$score_evaluations <- 1L
  result$hessian_evaluations <- 1L
  if (!is.list(proposed_state) || !is.list(proposed_state$raw) ||
      length(proposed_state$raw$value) != 1L ||
      !is.finite(proposed_state$raw$value) ||
      !identical(proposed_state$objective_provenance,
                 candidate$objective_provenance)) {
    result$termination_reason <- "proposed_state_evaluation_failed"
    result$initial <- list(
      theta_reduced_log = proposal$initial_theta_reduced_log
    )
    result$proposed <- list(
      theta_reduced_log = proposal$proposed_theta_reduced_log,
      local_step = proposal$local_step
    )
    return(result)
  }
  if (!identical(
      as.numeric(proposed_state$theta_reduced_log),
      as.numeric(proposal$proposed_theta_reduced_log))) {
    result$termination_reason <- "proposed_state_coordinate_mismatch"
    result$initial <- list(
      theta_reduced_log = proposal$initial_theta_reduced_log
    )
    result$proposed <- list(
      theta_reduced_log = proposal$proposed_theta_reduced_log,
      local_step = proposal$local_step,
      evaluated_theta_reduced_log =
        as.numeric(proposed_state$theta_reduced_log)
    )
    return(result)
  }

  expected_proposed_full_log <- ln_joint_probe_full_log(
    proposal$proposed_theta_reduced_log,
    candidate
  )
  if (is.null(expected_proposed_full_log) ||
      !identical(
        as.numeric(proposed_state$full_log),
        as.numeric(expected_proposed_full_log)
      )) {
    result$termination_reason <- "proposed_state_full_map_mismatch"
    result$initial <- list(
      theta_reduced_log = proposal$initial_theta_reduced_log
    )
    result$proposed <- list(
      theta_reduced_log = proposal$proposed_theta_reduced_log,
      local_step = proposal$local_step,
      evaluated_full_log = as.numeric(proposed_state$full_log),
      expected_full_log = as.numeric(expected_proposed_full_log)
    )
    return(result)
  }

  n_beta <- as.integer(candidate$n_beta)
  proposed_full_log <- as.numeric(proposed_state$full_log)
  proposed_blocks <- ln_joint_resolution_block_metrics(
    proposed_full_log,
    selected$theta_full_log,
    n_beta
  )
  initial_blocks <- ln_joint_resolution_block_metrics(
    candidate$theta_full_log,
    selected$theta_full_log,
    n_beta
  )
  proposed_objective <- ln_joint_objective_comparison(
    proposed_state$raw$value,
    selected$objective,
    resolution_thresholds$proposal_objective_equivalence
  )
  initial_root <- as.numeric(candidate$root_residual)
  initial_newton <- as.numeric(candidate$newton_displacement)
  proposed_root <- ln_joint_root_metric(
    proposed_state$U_subject_reduced,
    proposed_state$U_reduced
  )
  proposed_newton <- ln_joint_probe_newton_metric(proposed_state)
  proposed_boundary <- ln_joint_scaled_bound_distance(
    proposed_full_log,
    n_beta,
    candidate$bounds
  )
  proposed_curvature_compatibility <- ln_joint_curvature_compatibility(
    proposed_state$A_reduced,
    selected$state$A_reduced
  )
  proposed_curvature <- ln_joint_curvature_diagnostics(
    proposed_state$A_reduced
  )
  donor_additivity <- ln_joint_scaled_max(
    rowSums(as.matrix(proposed_state$U_subject_reduced)),
    as.numeric(proposed_state$U_reduced)
  )

  contraction <- list(
    full_state = ln_joint_probe_contracted_or_target(
      initial_blocks$full_state,
      proposed_blocks$full_state,
      resolution_thresholds$proposal_full_state,
      resolution_thresholds$contraction_ratio
    ),
    beta_block = ln_joint_probe_contracted_or_target(
      initial_blocks$beta_block,
      proposed_blocks$beta_block,
      resolution_thresholds$proposal_beta_block,
      resolution_thresholds$contraction_ratio
    ),
    dispersion_block = ln_joint_probe_contracted_or_target(
      initial_blocks$dispersion_block,
      proposed_blocks$dispersion_block,
      resolution_thresholds$proposal_dispersion_block,
      resolution_thresholds$contraction_ratio
    ),
    root = ln_joint_probe_contracted_or_target(
      initial_root,
      proposed_root,
      resolution_thresholds$proposal_root,
      resolution_thresholds$contraction_ratio
    ),
    newton = ln_joint_probe_contracted_or_target(
      initial_newton,
      proposed_newton,
      resolution_thresholds$proposal_newton,
      resolution_thresholds$contraction_ratio
    )
  )
  selected_thresholds <- ln_joint_thresholds()
  state_valid <- length(proposed_full_log) == n_beta + 2L &&
    all(is.finite(c(
      proposed_full_log,
      unlist(proposed_blocks, use.names = FALSE),
      proposed_root,
      proposed_newton,
      proposed_boundary,
      proposed_curvature_compatibility,
      donor_additivity
    ))) &&
    isTRUE(proposed_objective$valid) &&
    identical(proposed_objective$relation, "equivalent") &&
    proposed_boundary > resolution_thresholds$boundary_scaled_distance &&
    proposed_curvature_compatibility <=
      resolution_thresholds$curvature_compatibility &&
    isTRUE(proposed_curvature$positive_definite) &&
    proposed_curvature$scaled_min_eigenvalue >=
      selected_thresholds$curvature_min_eigen &&
    proposed_curvature$scaled_condition_number <=
      selected_thresholds$condition_number &&
    donor_additivity <= selected_thresholds$donor_additivity &&
    all(unlist(contraction, use.names = FALSE))

  result$anchor_assessments <- lapply(assessment_roots, function(anchor) {
    anchor_proposed_blocks <- ln_joint_resolution_block_metrics(
      proposed_full_log,
      anchor$theta_full_log,
      n_beta
    )
    anchor_initial_blocks <- ln_joint_resolution_block_metrics(
      candidate$theta_full_log,
      anchor$theta_full_log,
      n_beta
    )
    anchor_objective <- ln_joint_objective_comparison(
      proposed_state$raw$value,
      anchor$objective,
      resolution_thresholds$proposal_objective_equivalence
    )
    anchor_curvature <- ln_joint_curvature_compatibility(
      proposed_state$A_reduced,
      anchor$state$A_reduced
    )
    anchor_contraction <- list(
      full_state = ln_joint_probe_contracted_or_target(
        anchor_initial_blocks$full_state,
        anchor_proposed_blocks$full_state,
        resolution_thresholds$proposal_full_state,
        resolution_thresholds$contraction_ratio
      ),
      beta_block = ln_joint_probe_contracted_or_target(
        anchor_initial_blocks$beta_block,
        anchor_proposed_blocks$beta_block,
        resolution_thresholds$proposal_beta_block,
        resolution_thresholds$contraction_ratio
      ),
      dispersion_block = ln_joint_probe_contracted_or_target(
        anchor_initial_blocks$dispersion_block,
        anchor_proposed_blocks$dispersion_block,
        resolution_thresholds$proposal_dispersion_block,
        resolution_thresholds$contraction_ratio
      )
    )
    accepted <- all(is.finite(c(
      unlist(anchor_initial_blocks, use.names = FALSE),
      unlist(anchor_proposed_blocks, use.names = FALSE),
      anchor_curvature
    ))) &&
      isTRUE(anchor_objective$valid) &&
      identical(anchor_objective$relation, "equivalent") &&
      anchor_curvature <= resolution_thresholds$curvature_compatibility &&
      all(unlist(anchor_contraction, use.names = FALSE))
    list(
      anchor_candidate_key = ln_joint_resolution_candidate_key(anchor),
      accepted = isTRUE(accepted),
      objective_relation = anchor_objective$relation,
      objective_absolute_gap = anchor_objective$absolute_gap,
      initial = anchor_initial_blocks,
      proposed = anchor_proposed_blocks,
      curvature_compatibility = anchor_curvature,
      contraction = anchor_contraction
    )
  })
  state_valid <- isTRUE(state_valid) &&
    all(vapply(result$anchor_assessments, function(assessment) {
      isTRUE(assessment$accepted)
    }, logical(1)))

  result$initial <- list(
    theta_reduced_log = proposal$initial_theta_reduced_log,
    objective = as.numeric(candidate$objective),
    full_state_distance = initial_blocks$full_state,
    beta_block_distance = initial_blocks$beta_block,
    dispersion_block_distance = initial_blocks$dispersion_block,
    root_residual = initial_root,
    newton_displacement = initial_newton
  )
  result$proposed <- list(
    theta_reduced_log = proposal$proposed_theta_reduced_log,
    local_step = proposal$local_step,
    full_log = proposed_full_log,
    objective = as.numeric(proposed_state$raw$value),
    objective_relation = proposed_objective$relation,
    objective_absolute_gap = proposed_objective$absolute_gap,
    full_state_distance = proposed_blocks$full_state,
    beta_block_distance = proposed_blocks$beta_block,
    dispersion_block_distance = proposed_blocks$dispersion_block,
    boundary_scaled_distance = proposed_boundary,
    curvature_compatibility = proposed_curvature_compatibility,
    curvature_positive_definite =
      isTRUE(proposed_curvature$positive_definite),
    curvature_scaled_min_eigenvalue =
      as.numeric(proposed_curvature$scaled_min_eigenvalue),
    curvature_scaled_condition_number =
      as.numeric(proposed_curvature$scaled_condition_number),
    root_residual = proposed_root,
    newton_displacement = proposed_newton,
    donor_additivity_scaled = donor_additivity
  )
  result$contraction <- contraction
  result$accepted <- isTRUE(state_valid)
  result$termination_reason <- if (isTRUE(state_valid)) {
    "same_basin_contraction_established"
  } else {
    "same_basin_contraction_not_established"
  }
  result
}

ln_joint_basin_probe_record_valid <- function(
    probe, candidate, assessment_roots, resolution_thresholds) {
  template <- ln_joint_basin_probe_template(resolution_thresholds)
  expected_anchor_count <- length(assessment_roots)
  if (!is.list(probe) || !identical(names(probe), names(template)) ||
      !identical(probe$schema_version, 1L) ||
      !identical(probe$policy_version,
                 resolution_thresholds$policy_version) ||
      !is.logical(probe$accepted) || length(probe$accepted) != 1L ||
      is.na(probe$accepted) ||
      !is.character(probe$termination_reason) ||
      length(probe$termination_reason) != 1L ||
      is.na(probe$termination_reason) ||
      !nzchar(probe$termination_reason) ||
      !identical(
        probe$evaluation_budget,
        list(max_state_evaluations = 1L)
      ) ||
      !is.numeric(expected_anchor_count) ||
      length(expected_anchor_count) != 1L ||
      !is.finite(expected_anchor_count) || expected_anchor_count < 1L) {
    return(FALSE)
  }
  count_names <- c(
    "state_evaluations", "objective_evaluations",
    "score_evaluations", "hessian_evaluations"
  )
  counts_valid <- vapply(count_names, function(name) {
    value <- probe[[name]]
    is.integer(value) && length(value) == 1L && !is.na(value) &&
      value >= 0L && value <= 1L
  }, logical(1))
  if (!all(counts_valid)) {
    return(FALSE)
  }
  counts <- vapply(count_names, function(name) probe[[name]], integer(1))
  if (length(unique(counts)) != 1L) {
    return(FALSE)
  }
  if (!isTRUE(probe$accepted)) {
    return(TRUE)
  }
  finite_scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }
  initial_names <- c(
    "theta_reduced_log", "objective", "full_state_distance",
    "beta_block_distance", "dispersion_block_distance", "root_residual",
    "newton_displacement"
  )
  proposed_names <- c(
    "theta_reduced_log", "local_step", "full_log", "objective",
    "objective_relation", "objective_absolute_gap", "full_state_distance",
    "beta_block_distance", "dispersion_block_distance",
    "boundary_scaled_distance", "curvature_compatibility",
    "curvature_positive_definite", "curvature_scaled_min_eigenvalue",
    "curvature_scaled_condition_number", "root_residual",
    "newton_displacement", "donor_additivity_scaled"
  )
  if (!identical(probe$termination_reason,
                 "same_basin_contraction_established") ||
      !all(counts == 1L) || !is.list(probe$initial) ||
      !is.list(probe$proposed) || !is.list(probe$contraction) ||
      !all(initial_names %in% names(probe$initial)) ||
      !all(proposed_names %in% names(probe$proposed)) ||
      !identical(
        names(probe$contraction),
        c("full_state", "beta_block", "dispersion_block", "root", "newton")
      ) ||
      !all(vapply(probe$contraction, identical, logical(1), TRUE)) ||
      !is.list(probe$anchor_assessments) ||
      length(probe$anchor_assessments) != as.integer(expected_anchor_count)) {
    return(FALSE)
  }
  proposal <- ln_joint_basin_local_proposal(
    candidate,
    resolution_thresholds
  )
  if (!isTRUE(proposal$success) ||
      !identical(as.numeric(probe$initial$theta_reduced_log),
                 as.numeric(proposal$initial_theta_reduced_log)) ||
      !identical(as.numeric(probe$proposed$theta_reduced_log),
                 as.numeric(proposal$proposed_theta_reduced_log)) ||
      !identical(as.numeric(probe$proposed$local_step),
                 as.numeric(proposal$local_step)) ||
      !identical(as.numeric(probe$initial$objective),
                 as.numeric(candidate$objective)) ||
      !identical(as.numeric(probe$initial$root_residual),
                 as.numeric(candidate$root_residual)) ||
      !identical(as.numeric(probe$initial$newton_displacement),
                 as.numeric(candidate$newton_displacement)) ||
      !identical(
        as.numeric(probe$proposed$full_log),
        as.numeric(ln_joint_probe_full_log(
          proposal$proposed_theta_reduced_log,
          candidate
        ))
      ) ||
      length(probe$proposed$full_log) != as.integer(candidate$n_beta) + 2L ||
      any(!is.finite(c(
        probe$proposed$theta_reduced_log,
        probe$proposed$local_step,
        probe$proposed$full_log
      )))) {
    return(FALSE)
  }
  proposed_scalar_names <- setdiff(
    proposed_names,
    c(
      "theta_reduced_log", "local_step", "full_log", "objective_relation",
      "curvature_positive_definite"
    )
  )
  if (!all(vapply(probe$proposed[proposed_scalar_names],
                  finite_scalar, logical(1))) ||
      !identical(probe$proposed$objective_relation, "equivalent") ||
      !identical(probe$proposed$curvature_positive_definite, TRUE) ||
      probe$proposed$objective_absolute_gap >
        resolution_thresholds$proposal_objective_equivalence ||
      probe$proposed$boundary_scaled_distance <=
        resolution_thresholds$boundary_scaled_distance ||
      probe$proposed$curvature_compatibility >
        resolution_thresholds$curvature_compatibility ||
      probe$proposed$curvature_scaled_min_eigenvalue <
        ln_joint_thresholds()$curvature_min_eigen ||
      probe$proposed$curvature_scaled_condition_number >
        ln_joint_thresholds()$condition_number ||
      probe$proposed$root_residual > resolution_thresholds$proposal_root ||
      probe$proposed$newton_displacement >
        resolution_thresholds$proposal_newton ||
      probe$proposed$donor_additivity_scaled >
        ln_joint_thresholds()$donor_additivity) {
    return(FALSE)
  }
  expected_root_contraction <- ln_joint_probe_contracted_or_target(
    probe$initial$root_residual,
    probe$proposed$root_residual,
    resolution_thresholds$proposal_root,
    resolution_thresholds$contraction_ratio
  )
  expected_newton_contraction <- ln_joint_probe_contracted_or_target(
    probe$initial$newton_displacement,
    probe$proposed$newton_displacement,
    resolution_thresholds$proposal_newton,
    resolution_thresholds$contraction_ratio
  )
  if (!isTRUE(expected_root_contraction) ||
      !isTRUE(expected_newton_contraction) ||
      !identical(probe$contraction$root, expected_root_contraction) ||
      !identical(probe$contraction$newton, expected_newton_contraction)) {
    return(FALSE)
  }
  anchors_valid <- vapply(seq_along(assessment_roots), function(index) {
    anchor <- assessment_roots[[index]]
    assessment <- probe$anchor_assessments[[index]]
    expected_initial <- ln_joint_resolution_block_metrics(
      candidate$theta_full_log,
      anchor$theta_full_log,
      candidate$n_beta
    )
    if (!is.list(assessment) || !isTRUE(assessment$accepted) ||
        !identical(
          assessment$anchor_candidate_key,
          ln_joint_resolution_candidate_key(anchor)
        ) ||
        !identical(assessment$objective_relation, "equivalent") ||
        !is.list(assessment$initial) ||
        !identical(assessment$initial, expected_initial) ||
        !is.list(assessment$proposed) ||
        !all(vapply(assessment$proposed, finite_scalar, logical(1))) ||
        !finite_scalar(assessment$objective_absolute_gap) ||
        assessment$objective_absolute_gap >
          resolution_thresholds$proposal_objective_equivalence ||
        !finite_scalar(assessment$curvature_compatibility) ||
        assessment$curvature_compatibility >
          resolution_thresholds$curvature_compatibility ||
        !is.list(assessment$contraction) ||
        !identical(
          names(assessment$contraction),
          c("full_state", "beta_block", "dispersion_block")
        )) {
      return(FALSE)
    }
    targets <- c(
      full_state = resolution_thresholds$proposal_full_state,
      beta_block = resolution_thresholds$proposal_beta_block,
      dispersion_block = resolution_thresholds$proposal_dispersion_block
    )
    expected_contraction <- mapply(
      ln_joint_probe_contracted_or_target,
      initial = assessment$initial[names(targets)],
      proposed = assessment$proposed[names(targets)],
      target = targets,
      MoreArgs = list(ratio = resolution_thresholds$contraction_ratio),
      SIMPLIFY = TRUE,
      USE.NAMES = TRUE
    )
    all(unlist(
      assessment$proposed[names(targets)],
      use.names = TRUE
    ) <= targets) &&
      all(expected_contraction) &&
      identical(
        unname(unlist(assessment$contraction, use.names = FALSE)),
        unname(as.logical(expected_contraction))
      )
  }, logical(1))
  if (!all(anchors_valid)) {
    return(FALSE)
  }
  primary_assessment <- probe$anchor_assessments[[1L]]
  identical(
    c(
      full_state = probe$initial$full_state_distance,
      beta_block = probe$initial$beta_block_distance,
      dispersion_block = probe$initial$dispersion_block_distance
    ),
    unlist(primary_assessment$initial, use.names = TRUE)
  ) &&
    identical(
      c(
        full_state = probe$proposed$full_state_distance,
        beta_block = probe$proposed$beta_block_distance,
        dispersion_block = probe$proposed$dispersion_block_distance
      ),
      unlist(primary_assessment$proposed, use.names = TRUE)
    ) &&
    identical(
      as.numeric(probe$proposed$objective_absolute_gap),
      as.numeric(primary_assessment$objective_absolute_gap)
    ) &&
    identical(
      as.numeric(probe$proposed$curvature_compatibility),
      as.numeric(primary_assessment$curvature_compatibility)
    )
}

ln_joint_map <- function(nb, L = NULL, b = NULL) {
  if (is.null(L)) {
    beta_star <- rep(0, nb)
    K <- diag(nb)
    constrained <- FALSE
  } else {
    constr <- constraint_reparameterization(L, b)
    beta_star <- as.numeric(constr$beta_star)
    K <- constr$null_basis
    if (is.null(K)) {
      K <- matrix(0, nrow = nb, ncol = 0L)
    }
    constrained <- TRUE
  }
  free_dim <- ncol(K)
  compose <- function(zeta) {
    if (!free_dim) beta_star else as.numeric(beta_star + K %*% zeta)
  }
  project <- function(beta) {
    if (!free_dim) {
      numeric(0)
    } else {
      as.numeric(qr.solve(K, as.numeric(beta) - beta_star, tol = 1e-10))
    }
  }
  derivative_map <- matrix(0, nrow = nb + 2L, ncol = free_dim + 2L)
  if (free_dim) {
    derivative_map[seq_len(nb), seq_len(free_dim)] <- K
  }
  derivative_map[nb + 1L, free_dim + 1L] <- 1
  derivative_map[nb + 2L, free_dim + 2L] <- 1
  list(
    constrained = constrained,
    beta_star = beta_star,
    K = K,
    free_dim = free_dim,
    compose = compose,
    project = project,
    derivative_map = derivative_map,
    L = L,
    b = b
  )
}

ln_joint_full_natural <- function(theta_reduced, map) {
  free_dim <- map$free_dim
  zeta <- if (free_dim) theta_reduced[seq_len(free_dim)] else numeric(0)
  dispersion <- theta_reduced[free_dim + seq_len(2L)]
  c(map$compose(zeta), dispersion)
}

ln_joint_objective <- function(theta_reduced, gene_index, posv, ctx, map) {
  full <- ln_joint_full_natural(theta_reduced, map)
  ptmg_ll(
    full,
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    n_one = posv$n_onetwo,
    ytwo = posv$ytwo,
    fam = ctx$id,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[gene_index, ],
    posind = ctx$posind[[gene_index]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k
  )
}

ln_joint_gradient <- function(theta_reduced, gene_index, posv, ctx, map) {
  full <- ln_joint_full_natural(theta_reduced, map)
  gradient_full <- ptmg_der(
    full,
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    n_one = posv$n_onetwo,
    ytwo = posv$ytwo,
    fam = ctx$id,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[gene_index, ],
    posind = ctx$posind[[gene_index]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k
  )
  as.numeric(crossprod(map$derivative_map, gradient_full))
}

ln_joint_derivative_state <- function(theta_reduced_log, gene_index, posv, ctx, map) {
  free_dim <- map$free_dim
  zeta <- if (free_dim) theta_reduced_log[seq_len(free_dim)] else numeric(0)
  beta <- map$compose(zeta)
  full_log <- c(beta, theta_reduced_log[free_dim + seq_len(2L)])
  raw <- ptmg_ll_der_hes4(
    full_log,
    ctx$pred,
    ctx$offset,
    posv$Y,
    posv$n_onetwo,
    posv$ytwo,
    ctx$fid,
    ctx$cumsumy[gene_index, ],
    ctx$posind[[gene_index]],
    posv$posindy,
    ctx$nb,
    ctx$nind,
    ctx$k,
    correct_posindy = TRUE
  )
  T <- map$derivative_map
  A_full <- 0.5 * (raw$hessian + t(raw$hessian))
  U_full <- -as.numeric(raw$gradient)
  U_subject_full <- -as.matrix(raw$per_subject_gradients)
  A_reduced <- crossprod(T, A_full %*% T)
  A_reduced <- 0.5 * (A_reduced + t(A_reduced))
  U_reduced <- as.numeric(crossprod(T, U_full))
  U_subject_reduced <- crossprod(T, U_subject_full)
  list(
    full_log = full_log,
    theta_reduced_log = theta_reduced_log,
    objective_provenance = ln_joint_objective_provenance(),
    beta = beta,
    sigma2 = exp(full_log[ctx$nb + 1L]),
    gamma = exp(full_log[ctx$nb + 2L]),
    q_gradient_full = as.numeric(raw$gradient),
    q_gradient_reduced = -U_reduced,
    U_full = U_full,
    U_reduced = U_reduced,
    U_subject_full = U_subject_full,
    U_subject_reduced = U_subject_reduced,
    A_full = A_full,
    A_reduced = A_reduced,
    raw = raw
  )
}

ln_joint_curvature_diagnostics <- function(A) {
  A <- 0.5 * (A + t(A))
  if (!length(A) || any(!is.finite(A))) {
    return(list(
      scaled_min_eigenvalue = NA_real_,
      scaled_condition_number = Inf,
      positive_definite = FALSE
    ))
  }
  d <- pmax(abs(diag(A)), sqrt(.Machine$double.eps))
  scaled <- A / sqrt(outer(d, d))
  eig <- tryCatch(
    eigen(scaled, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eig) || any(!is.finite(eig))) {
    return(list(
      scaled_min_eigenvalue = NA_real_,
      scaled_condition_number = Inf,
      positive_definite = FALSE
    ))
  }
  min_eig <- min(eig)
  cond <- if (min_eig > 0) max(eig) / min_eig else Inf
  list(
    scaled_min_eigenvalue = min_eig,
    scaled_condition_number = cond,
    positive_definite = min_eig > 0,
    scaled_eigenvalues = eig
  )
}

ln_joint_score_matrix_diagnostics <- function(x,
                                              condition_limit =
                                                ln_joint_thresholds()$condition_number,
                                              symmetry_tolerance =
                                                ln_joint_thresholds()$hard_stop$
                                                  covariance_identity) {
  unavailable <- function(nrow = NA_integer_, ncol = NA_integer_) {
    list(
      nrow = as.integer(nrow),
      ncol = as.integer(ncol),
      finite = FALSE,
      square = FALSE,
      symmetry_error = NA_real_,
      symmetry_acceptable = FALSE,
      rank = 0L,
      rank_tolerance = NA_real_,
      condition_number = Inf,
      scaled_condition_number = Inf,
      eigenvalues = NULL,
      scaled_eigenvalues = NULL,
      positive_definite = FALSE,
      full_rank = FALSE,
      condition_acceptable = FALSE,
      usable = FALSE
    )
  }
  if (is.null(x)) {
    return(unavailable())
  }
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  if (!length(x) || nr != nc || any(!is.finite(x))) {
    out <- unavailable(nr, nc)
    out$finite <- length(x) > 0L && all(is.finite(x))
    out$square <- nr == nc
    return(out)
  }
  symmetry_error <- as.numeric(
    norm(x - t(x), "F") / (1 + norm(x, "F"))
  )
  symmetric_x <- 0.5 * (x + t(x))
  singular_values <- tryCatch(
    svd(symmetric_x, nu = 0L, nv = 0L)$d,
    error = function(e) NULL
  )
  eigenvalues <- tryCatch(
    eigen(symmetric_x, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  diagonal_scale <- sqrt(pmax(
    abs(diag(symmetric_x)),
    .Machine$double.eps
  ))
  scaled_x <- symmetric_x / outer(diagonal_scale, diagonal_scale)
  scaled_singular_values <- tryCatch(
    svd(scaled_x, nu = 0L, nv = 0L)$d,
    error = function(e) NULL
  )
  scaled_eigenvalues <- tryCatch(
    eigen(scaled_x, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(singular_values) || is.null(eigenvalues) ||
      is.null(scaled_singular_values) || is.null(scaled_eigenvalues) ||
      any(!is.finite(c(
        singular_values,
        eigenvalues,
        scaled_singular_values,
        scaled_eigenvalues
      )))) {
    out <- unavailable(nr, nc)
    out$finite <- TRUE
    out$square <- TRUE
    out$symmetry_error <- symmetry_error
    return(out)
  }
  max_singular <- max(singular_values)
  rank_tolerance <- max(nr, nc) * .Machine$double.eps * max_singular
  rank <- if (max_singular > 0) {
    sum(singular_values > rank_tolerance)
  } else {
    0L
  }
  scaled_max_singular <- max(scaled_singular_values)
  scaled_rank_tolerance <- max(nr, nc) * .Machine$double.eps *
    scaled_max_singular
  scaled_rank <- if (scaled_max_singular > 0) {
    sum(scaled_singular_values > scaled_rank_tolerance)
  } else {
    0L
  }
  condition_number <- if (rank == nr && min(singular_values) > 0) {
    max_singular / min(singular_values)
  } else {
    Inf
  }
  scaled_condition_number <- if (scaled_rank == nr &&
      min(scaled_singular_values) > 0) {
    scaled_max_singular / min(scaled_singular_values)
  } else {
    Inf
  }
  eigen_tolerance <- nr * .Machine$double.eps *
    max(abs(eigenvalues))
  positive_definite <- min(eigenvalues) > eigen_tolerance
  full_rank <- rank == nr
  symmetry_acceptable <- is.finite(symmetry_error) &&
    symmetry_error <= symmetry_tolerance
  condition_acceptable <- is.finite(scaled_condition_number) &&
    scaled_condition_number <= condition_limit
  list(
    nrow = as.integer(nr),
    ncol = as.integer(nc),
    finite = TRUE,
    square = TRUE,
    symmetry_error = symmetry_error,
    symmetry_acceptable = symmetry_acceptable,
    rank = as.integer(rank),
    scaled_rank = as.integer(scaled_rank),
    rank_tolerance = rank_tolerance,
    scaled_rank_tolerance = scaled_rank_tolerance,
    condition_number = condition_number,
    scaled_condition_number = scaled_condition_number,
    eigenvalues = eigenvalues,
    scaled_eigenvalues = scaled_eigenvalues,
    positive_definite = positive_definite,
    full_rank = full_rank,
    condition_acceptable = condition_acceptable,
    usable = symmetry_acceptable && full_rank && positive_definite &&
      condition_acceptable
  )
}

ln_joint_score_coordinate_map <- function(L, nb, null_basis = NULL,
                                          interest_basis = NULL) {
  thresholds <- ln_joint_thresholds()
  L <- as.matrix(L)
  r <- nrow(L)
  observed_rank <- if (length(L) && all(is.finite(L))) {
    qr(L)$rank
  } else {
    0L
  }
  unavailable <- function(reason) {
    list(
      success = FALSE,
      failure_reason = reason,
      requested_rows = as.integer(r),
      rank = as.integer(observed_rank)
    )
  }
  if (!length(L) || ncol(L) != nb || r < 1L || r > nb ||
      any(!is.finite(L))) {
    return(unavailable("invalid_hypothesis_matrix"))
  }
  if (observed_rank != r) {
    return(unavailable("hypothesis_rank_deficient"))
  }
  qrt <- qr(t(L))
  if (qrt$rank != r) {
    return(unavailable("hypothesis_rank_deficient"))
  }
  complete_basis <- qr.Q(qrt, complete = TRUE)
  if (is.null(interest_basis)) {
    interest_basis <- complete_basis[, seq_len(r), drop = FALSE]
  } else {
    interest_basis <- as.matrix(interest_basis)
  }
  nuisance_beta_dim <- nb - r
  if (is.null(null_basis)) {
    null_basis <- if (nuisance_beta_dim) {
      complete_basis[, r + seq_len(nuisance_beta_dim), drop = FALSE]
    } else {
      matrix(0, nrow = nb, ncol = 0L)
    }
  } else {
    null_basis <- as.matrix(null_basis)
  }
  if (!all(dim(interest_basis) == c(nb, r)) ||
      any(!is.finite(interest_basis))) {
    return(unavailable("invalid_interest_basis"))
  }
  if (!all(dim(null_basis) == c(nb, nuisance_beta_dim)) ||
      any(!is.finite(null_basis))) {
    return(unavailable("invalid_null_basis"))
  }
  tested_map <- L %*% interest_basis
  if (qr(tested_map)$rank != r) {
    return(unavailable("interest_basis_rank_deficient"))
  }
  null_residual_scaled <- if (nuisance_beta_dim) {
    max(abs(L %*% null_basis)) /
      (1 + norm(L, "F") * norm(null_basis, "F"))
  } else {
    0
  }
  if (!is.finite(null_residual_scaled) ||
      null_residual_scaled > thresholds$hard_stop$constraint) {
    return(unavailable("null_basis_not_in_hypothesis_kernel"))
  }
  beta_map <- cbind(interest_basis, null_basis)
  if (qr(beta_map)$rank != nb) {
    return(unavailable("interest_nuisance_map_rank_deficient"))
  }
  beta_map_condition <- tryCatch(
    kappa(beta_map, exact = TRUE),
    error = function(e) Inf
  )
  if (!is.finite(beta_map_condition) ||
      beta_map_condition > thresholds$hard_stop$condition_number) {
    return(unavailable("interest_nuisance_map_ill_conditioned"))
  }
  full_dim <- nb + 2L
  full_map <- matrix(0, nrow = full_dim, ncol = full_dim)
  full_map[seq_len(nb), seq_len(nb)] <- beta_map
  full_map[nb + 1L, nb + 1L] <- 1
  full_map[nb + 2L, nb + 2L] <- 1
  list(
    success = TRUE,
    failure_reason = NA_character_,
    rank = as.integer(r),
    interest_basis = interest_basis,
    null_basis = null_basis,
    beta = beta_map,
    full = full_map,
    tested_hypothesis_map = tested_map,
    null_residual_scaled = null_residual_scaled,
    beta_map_condition_number = beta_map_condition,
    interest_crossprod = crossprod(interest_basis),
    null_crossprod = crossprod(null_basis)
  )
}

ln_joint_robust_score_diagnostic <- function(constrained_fit,
                                             null_basis = NULL,
                                             interest_basis = NULL) {
  thresholds <- ln_joint_thresholds()
  empty_matrix_diagnostics <- ln_joint_score_matrix_diagnostics(NULL)
  result <- list(
    diagnostic_available = FALSE,
    failure_reason = NA_character_,
    public_available = FALSE,
    public_status = "joint_ln_calibration_pending",
    statistic = NA_real_,
    df = NA_integer_,
    theoretical_df = NA_integer_,
    theoretical_reference = paste(
      "chi-square first-order reference only under smooth interior",
      "independent-donor fixed-dimension regularity"
    ),
    score_sign = "maximizing_U_equals_negative_gradient_of_q",
    no_numerical_repair = TRUE,
    hypothesis = NULL,
    coordinate_map = NULL,
    U_coordinate = NULL,
    U_subject_coordinate = NULL,
    U_eff = NULL,
    U_eff_subject = NULL,
    H_psi_chi = NULL,
    H_chi_chi = NULL,
    J_eff = NULL,
    projection = NULL,
    diagnostics = list(
      constraint_residual_scaled = NA_real_,
      nuisance_root_scaled = NA_real_,
      full_donor_additivity_scaled = NA_real_,
      coordinate_donor_additivity_scaled = NA_real_,
      efficient_donor_additivity_scaled = NA_real_,
      efficient_vs_tested_score_scaled = NA_real_,
      J_eff_projection_identity_scaled = NA_real_,
      score_sign_aggregate_scaled = NA_real_,
      score_sign_donor_scaled = NA_real_,
      nuisance_root_scale = paste(
        "max_j abs(sum_i U_ji) /",
        "(1 + sqrt(sum_i U_ji^2))"
      ),
      H_coordinate_symmetry_scaled = NA_real_,
      H_chi_chi = empty_matrix_diagnostics,
      J_eff = empty_matrix_diagnostics
    )
  )
  fail <- function(reason) {
    result$failure_reason <- reason
    result
  }

  constraint <- constrained_fit$constraint
  L <- constraint$L
  b <- constraint$b
  if (is.null(L) || is.null(b)) {
    return(fail("missing_hypothesis"))
  }
  L <- as.matrix(L)
  b <- as.numeric(b)
  if (!nrow(L) || length(b) != nrow(L) ||
      any(!is.finite(c(L, b)))) {
    return(fail("invalid_hypothesis"))
  }
  hypothesis_rank <- qr(L)$rank
  result$hypothesis <- list(
    L = L,
    b = b,
    requested_rows = as.integer(nrow(L)),
    rank = as.integer(hypothesis_rank)
  )
  if (hypothesis_rank != nrow(L)) {
    return(fail("hypothesis_rank_deficient"))
  }
  result$df <- result$theoretical_df <- as.integer(nrow(L))

  if (!identical(constrained_fit$route, "R0")) {
    reason <- if (identical(constrained_fit$failure_reason,
                            "dispersion_boundary")) {
      "constrained_boundary"
    } else {
      "constrained_route_not_R0"
    }
    return(fail(reason))
  }
  if (!isTRUE(constrained_fit$inference_available) ||
      !isTRUE(constrained_fit$lrt_gate$converged)) {
    return(fail("constrained_endpoint_uncertified"))
  }
  if (!isTRUE(constrained_fit$lrt_gate$interior) ||
      !isTRUE(constrained_fit$lrt_gate$nuisance_interior)) {
    return(fail("constrained_boundary"))
  }
  state <- constrained_fit$ln_joint$state
  if (is.null(state) || is.null(state$beta)) {
    return(fail("missing_constrained_score_state"))
  }
  nb <- length(state$beta)
  full_dim <- nb + 2L
  U_full <- as.numeric(state$U_full)
  U_subject_full <- as.matrix(state$U_subject_full)
  H_full <- as.matrix(state$A_full)
  state_dimensions_ok <- length(U_full) == full_dim &&
    nrow(U_subject_full) == full_dim && ncol(U_subject_full) > 0L &&
    all(dim(H_full) == c(full_dim, full_dim))
  if (!state_dimensions_ok || any(!is.finite(c(
    state$beta,
    U_full,
    U_subject_full,
    H_full
  )))) {
    return(fail("invalid_constrained_score_state"))
  }
  state_full_log <- as.numeric(state$full_log)
  state_coordinates_match <- length(state_full_log) == full_dim &&
    identical(
      state_full_log[seq_len(nb)],
      as.numeric(state$beta)
    ) &&
    length(state$sigma2) == 1L && length(state$gamma) == 1L &&
    identical(
      as.numeric(state$sigma2),
      exp(state_full_log[nb + 1L])
    ) &&
    identical(
      as.numeric(state$gamma),
      exp(state_full_log[nb + 2L])
    )
  if (!isTRUE(state_coordinates_match)) {
    return(fail("constrained_state_coordinate_mismatch"))
  }
  returned_theta <- constrained_fit$theta
  returned_state_matches <- is.list(returned_theta) &&
    identical(as.numeric(returned_theta$beta), as.numeric(state$beta)) &&
    identical(as.numeric(returned_theta$subVar), as.numeric(state$sigma2)) &&
    identical(as.numeric(returned_theta$cellVar), as.numeric(state$gamma)) &&
    identical(
      as.numeric(constrained_fit$beta_constrained),
      as.numeric(state$beta)
    )
  if (!isTRUE(returned_state_matches)) {
    return(fail("constrained_return_state_mismatch"))
  }
  candidate_records <- constrained_fit$ln_joint$candidate_records
  certified_state_match <- length(candidate_records) > 0L && any(vapply(
    candidate_records,
    function(candidate) {
      isTRUE(candidate$certified) &&
        !is.null(candidate$theta_full_log) &&
        identical(
          as.numeric(candidate$theta_full_log),
          state_full_log
        )
    },
    logical(1)
  ))
  if (!certified_state_match) {
    return(fail("constrained_endpoint_certification_missing"))
  }

  map <- ln_joint_score_coordinate_map(
    L,
    nb,
    null_basis = if (is.null(null_basis)) constraint$null_basis else null_basis,
    interest_basis = interest_basis
  )
  if (!isTRUE(map$success)) {
    return(fail(map$failure_reason))
  }
  result$coordinate_map <- map

  constraint_residual <- as.numeric(L %*% state$beta - b)
  constraint_denom <- 1 + max(abs(b)) +
    max(rowSums(abs(L))) * max(abs(state$beta))
  constraint_residual_scaled <- max(abs(constraint_residual)) /
    constraint_denom
  result$diagnostics$constraint_residual_scaled <-
    constraint_residual_scaled
  if (!is.finite(constraint_residual_scaled) ||
      constraint_residual_scaled > thresholds$constraint) {
    return(fail("constraint_not_satisfied"))
  }

  score_sign_aggregate_scaled <- ln_joint_scaled_max(
    U_full,
    -state$q_gradient_full
  )
  score_sign_donor_scaled <- ln_joint_scaled_max(
    U_subject_full,
    -state$raw$per_subject_gradients
  )
  full_donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_subject_full),
    U_full
  )
  result$diagnostics$score_sign_aggregate_scaled <-
    score_sign_aggregate_scaled
  result$diagnostics$score_sign_donor_scaled <- score_sign_donor_scaled
  result$diagnostics$full_donor_additivity_scaled <-
    full_donor_additivity_scaled
  if (!all(is.finite(c(
    score_sign_aggregate_scaled,
    score_sign_donor_scaled
  ))) || max(score_sign_aggregate_scaled, score_sign_donor_scaled) >
      thresholds$donor_additivity) {
    return(fail("maximizing_score_sign_mismatch"))
  }
  if (!is.finite(full_donor_additivity_scaled) ||
      full_donor_additivity_scaled > thresholds$donor_additivity) {
    return(fail("full_score_donor_additivity_failed"))
  }

  T <- map$full
  U_coordinate <- as.numeric(crossprod(T, U_full))
  U_subject_coordinate <- crossprod(T, U_subject_full)
  H_coordinate_raw <- crossprod(T, H_full %*% T)
  H_coordinate_symmetry_scaled <- as.numeric(
    norm(H_coordinate_raw - t(H_coordinate_raw), "F") /
      (1 + norm(H_coordinate_raw, "F"))
  )
  result$diagnostics$H_coordinate_symmetry_scaled <-
    H_coordinate_symmetry_scaled
  if (!is.finite(H_coordinate_symmetry_scaled) ||
      H_coordinate_symmetry_scaled >
        thresholds$hard_stop$covariance_identity) {
    return(fail("transformed_sensitivity_asymmetric"))
  }
  H_coordinate <- H_coordinate_raw
  H_coordinate <- 0.5 * (H_coordinate + t(H_coordinate))
  r <- result$df
  interest <- seq_len(r)
  nuisance <- r + seq_len(full_dim - r)
  H_psi_chi <- H_coordinate[interest, nuisance, drop = FALSE]
  H_chi_chi <- H_coordinate[nuisance, nuisance, drop = FALSE]
  coordinate_donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_subject_coordinate),
    U_coordinate
  )
  nuisance_root_scaled <- ln_joint_root_metric(
    U_subject_coordinate[nuisance, , drop = FALSE],
    U_coordinate[nuisance]
  )
  result$U_coordinate <- U_coordinate
  result$U_subject_coordinate <- U_subject_coordinate
  result$H_psi_chi <- H_psi_chi
  result$H_chi_chi <- H_chi_chi
  result$diagnostics$coordinate_donor_additivity_scaled <-
    coordinate_donor_additivity_scaled
  result$diagnostics$nuisance_root_scaled <- nuisance_root_scaled
  if (!is.finite(coordinate_donor_additivity_scaled) ||
      coordinate_donor_additivity_scaled > thresholds$donor_additivity) {
    return(fail("transformed_score_donor_additivity_failed"))
  }
  if (!is.finite(nuisance_root_scaled) ||
      nuisance_root_scaled > thresholds$root) {
    return(fail("constrained_nuisance_root_not_certified"))
  }

  H_chi_diagnostics <- ln_joint_score_matrix_diagnostics(H_chi_chi)
  result$diagnostics$H_chi_chi <- H_chi_diagnostics
  if (!isTRUE(H_chi_diagnostics$symmetry_acceptable)) {
    return(fail("H_chi_chi_asymmetric"))
  }
  if (!isTRUE(H_chi_diagnostics$full_rank)) {
    return(fail("H_chi_chi_rank_deficient"))
  }
  if (!isTRUE(H_chi_diagnostics$positive_definite)) {
    return(fail("H_chi_chi_not_positive_definite"))
  }
  if (!isTRUE(H_chi_diagnostics$condition_acceptable)) {
    return(fail("H_chi_chi_ill_conditioned"))
  }
  nuisance_adjustment <- tryCatch(
    t(solve(t(H_chi_chi), t(H_psi_chi))),
    error = function(e) NULL
  )
  if (is.null(nuisance_adjustment) ||
      any(!is.finite(nuisance_adjustment))) {
    return(fail("H_chi_chi_solve_failed"))
  }
  projection <- cbind(diag(r), -nuisance_adjustment)
  U_eff <- as.numeric(projection %*% U_coordinate)
  U_eff_subject <- projection %*% U_subject_coordinate
  J_eff <- tcrossprod(U_eff_subject)
  J_coordinate <- tcrossprod(U_subject_coordinate)
  J_eff_projected <- projection %*% J_coordinate %*% t(projection)
  efficient_donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_eff_subject),
    U_eff
  )
  efficient_vs_tested_score_scaled <- ln_joint_scaled_max(
    U_eff,
    U_coordinate[interest]
  )
  J_eff_projection_identity_scaled <- ln_joint_relative_frobenius(
    J_eff,
    J_eff_projected
  )
  result$projection <- projection
  result$U_eff <- U_eff
  result$U_eff_subject <- U_eff_subject
  result$J_eff <- J_eff
  result$diagnostics$efficient_donor_additivity_scaled <-
    efficient_donor_additivity_scaled
  result$diagnostics$efficient_vs_tested_score_scaled <-
    efficient_vs_tested_score_scaled
  result$diagnostics$J_eff_projection_identity_scaled <-
    J_eff_projection_identity_scaled
  if (!is.finite(efficient_donor_additivity_scaled) ||
      efficient_donor_additivity_scaled > thresholds$donor_additivity) {
    return(fail("efficient_score_donor_additivity_failed"))
  }
  if (!is.finite(efficient_vs_tested_score_scaled) ||
      efficient_vs_tested_score_scaled > thresholds$root) {
    return(fail("efficient_score_not_equal_tested_score_at_null"))
  }
  if (!is.finite(J_eff_projection_identity_scaled) ||
      J_eff_projection_identity_scaled > thresholds$covariance_identity) {
    return(fail("J_eff_projection_identity_failed"))
  }

  J_eff_diagnostics <- ln_joint_score_matrix_diagnostics(J_eff)
  result$diagnostics$J_eff <- J_eff_diagnostics
  if (!isTRUE(J_eff_diagnostics$symmetry_acceptable)) {
    return(fail("J_eff_asymmetric"))
  }
  if (!isTRUE(J_eff_diagnostics$full_rank)) {
    return(fail("J_eff_rank_deficient"))
  }
  if (!isTRUE(J_eff_diagnostics$positive_definite)) {
    return(fail("J_eff_not_positive_definite"))
  }
  if (!isTRUE(J_eff_diagnostics$condition_acceptable)) {
    return(fail("J_eff_ill_conditioned"))
  }
  J_eff_solution <- tryCatch(
    solve(J_eff, U_eff),
    error = function(e) NULL
  )
  if (is.null(J_eff_solution) || any(!is.finite(J_eff_solution))) {
    return(fail("J_eff_solve_failed"))
  }
  statistic <- as.numeric(crossprod(U_eff, J_eff_solution))
  if (length(statistic) != 1L || !is.finite(statistic) || statistic < 0) {
    return(fail("nonfinite_or_negative_score_statistic"))
  }
  result$statistic <- statistic
  result$diagnostic_available <- TRUE
  result$failure_reason <- NA_character_
  result
}

ln_joint_ratio_endpoint_audit <- function(fit,
                                          endpoint = c(
                                            "unrestricted",
                                            "constrained"
                                          ),
                                          objective_context = NULL) {
  endpoint <- match.arg(endpoint)
  thresholds <- ln_joint_thresholds()
  fit_is_list <- is.list(fit)
  fit_route <- if (fit_is_list) fit$route else NULL
  fit_debug <- if (fit_is_list && is.list(fit$debug)) fit$debug else list()
  fit_gene_index <- tryCatch(
    suppressWarnings(as.integer(fit_debug$gene_index)),
    error = function(e) NA_integer_
  )
  if (length(fit_gene_index) != 1L || is.na(fit_gene_index)) {
    fit_gene_index <- NA_integer_
  }
  result <- list(
    success = FALSE,
    failure_reason = NA_character_,
    endpoint = endpoint,
    route = fit_route,
    inference_available = fit_is_list && isTRUE(fit$inference_available),
    gene_index = fit_gene_index,
    n_beta = NA_integer_,
    n_donors = NA_integer_,
    full_log = NULL,
    state_objective = NA_real_,
    stored_state_objective = NA_real_,
    reevaluated_objective = NA_real_,
    candidate_objective = NA_real_,
    loglik = NA_real_,
    loglik_ln = NA_real_,
    certified_candidate_index = NA_integer_,
    objective_provenance = NULL,
    objective_instance = NULL,
    diagnostics = list(
      candidate_objective_state_scaled = NA_real_,
      loglik_objective_scaled = NA_real_,
      loglik_ln_objective_scaled = NA_real_,
      score_sign_scaled = NA_real_,
      donor_additivity_scaled = NA_real_,
      objective_recheck_scaled = NA_real_,
      full_log_recheck_scaled = NA_real_,
      U_full_recheck_scaled = NA_real_,
      U_subject_full_recheck_scaled = NA_real_,
      A_full_recheck_scaled = NA_real_,
      U_reduced_recheck_scaled = NA_real_,
      U_subject_reduced_recheck_scaled = NA_real_,
      A_reduced_recheck_scaled = NA_real_,
      root_residual = NA_real_,
      newton_displacement = NA_real_,
      curvature = NULL
    )
  )
  fail <- function(reason) {
    result$failure_reason <- paste0(endpoint, "_", reason)
    result
  }

  if (!fit_is_list) {
    return(fail("fit_invalid"))
  }
  if (!identical(fit$route, "R0") || !isTRUE(fit$inference_available)) {
    return(fail("route_not_R0"))
  }
  expected_method <- if (identical(endpoint, "unrestricted")) {
    "joint_ln"
  } else {
    "joint_ln_constrained"
  }
  if (!identical(fit_debug$method, expected_method)) {
    return(fail("estimator_mismatch"))
  }
  if (!is.list(fit$ln_joint)) {
    return(fail("joint_state_invalid"))
  }
  state <- fit$ln_joint$state
  if (!is.list(state) || !is.list(state$raw) || is.null(state$beta) ||
      is.null(state$raw$value)) {
    return(fail("state_missing"))
  }
  nb <- length(state$beta)
  full_dim <- nb + 2L
  U_full <- as.numeric(state$U_full)
  U_subject_full <- as.matrix(state$U_subject_full)
  A_full <- as.matrix(state$A_full)
  result$n_beta <- as.integer(nb)
  result$n_donors <- as.integer(ncol(U_subject_full))
  state_dimensions_ok <- length(U_full) == full_dim &&
    nrow(U_subject_full) == full_dim && ncol(U_subject_full) > 0L &&
    all(dim(A_full) == c(full_dim, full_dim))
  if (!state_dimensions_ok || any(!is.finite(c(
    state$beta,
    U_full,
    U_subject_full,
    A_full,
    state$raw$value
  )))) {
    return(fail("state_invalid"))
  }
  U_reduced <- as.numeric(state$U_reduced)
  U_subject_reduced <- as.matrix(state$U_subject_reduced)
  A_reduced <- as.matrix(state$A_reduced)
  theta_reduced_log <- as.numeric(state$theta_reduced_log)
  reduced_dim <- length(U_reduced)
  reduced_dimensions_ok <- reduced_dim > 0L &&
    nrow(U_subject_reduced) == reduced_dim &&
    ncol(U_subject_reduced) == ncol(U_subject_full) &&
    all(dim(A_reduced) == c(reduced_dim, reduced_dim)) &&
    length(theta_reduced_log) == reduced_dim
  if (!reduced_dimensions_ok || any(!is.finite(c(
    U_reduced,
    U_subject_reduced,
    A_reduced,
    theta_reduced_log
  )))) {
    return(fail("reduced_state_invalid"))
  }
  result$diagnostics$root_residual <- ln_joint_root_metric(
    U_subject_reduced,
    U_reduced
  )
  result$diagnostics$curvature <- ln_joint_curvature_diagnostics(A_reduced)
  newton_step <- tryCatch(
    solve(A_reduced, U_reduced),
    error = function(e) NULL
  )
  result$diagnostics$newton_displacement <- if (!is.null(newton_step) &&
      all(is.finite(newton_step))) {
    max(abs(newton_step) / (1 + abs(theta_reduced_log)))
  } else {
    NA_real_
  }
  if (!is.finite(result$diagnostics$root_residual) ||
      result$diagnostics$root_residual > thresholds$root) {
    return(fail("selected_root_not_certified"))
  }
  if (!is.finite(result$diagnostics$newton_displacement) ||
      result$diagnostics$newton_displacement > thresholds$newton) {
    return(fail("selected_newton_displacement_not_certified"))
  }
  curvature <- result$diagnostics$curvature
  if (!isTRUE(curvature$positive_definite) ||
      !is.finite(curvature$scaled_min_eigenvalue) ||
      curvature$scaled_min_eigenvalue < thresholds$curvature_min_eigen) {
    return(fail("selected_curvature_not_certified"))
  }
  if (!is.finite(curvature$scaled_condition_number) ||
      curvature$scaled_condition_number > thresholds$condition_number) {
    return(fail("selected_curvature_condition_not_certified"))
  }

  expected_provenance <- ln_joint_objective_provenance()
  fit_provenance <- fit$ln_joint$objective_provenance
  state_provenance <- state$objective_provenance
  if (!identical(fit_provenance, expected_provenance) ||
      !identical(state_provenance, expected_provenance)) {
    return(fail("objective_provenance_mismatch"))
  }
  fit_instance <- fit$ln_joint$objective_instance
  state_instance <- state$objective_instance
  if (is.null(fit_instance) || is.null(state_instance)) {
    return(fail("objective_instance_missing"))
  }
  if (!identical(fit_instance, state_instance)) {
    return(fail("objective_instance_mismatch"))
  }
  if (!identical(
    fit$ln_joint$parameterization,
    expected_provenance$derivative_coordinates
  ) || !identical(fit$ln_joint$score_sign, "maximizing")) {
    return(fail("coordinate_or_score_convention_mismatch"))
  }
  result$objective_provenance <- expected_provenance
  result$objective_instance <- fit_instance

  state_full_log <- as.numeric(state$full_log)
  state_coordinates_match <- length(state_full_log) == full_dim &&
    identical(
      state_full_log[seq_len(nb)],
      as.numeric(state$beta)
    ) &&
    length(state$sigma2) == 1L && length(state$gamma) == 1L &&
    identical(as.numeric(state$sigma2), exp(state_full_log[nb + 1L])) &&
    identical(as.numeric(state$gamma), exp(state_full_log[nb + 2L]))
  if (!isTRUE(state_coordinates_match)) {
    return(fail("state_coordinate_mismatch"))
  }
  result$full_log <- state_full_log
  returned_theta <- fit$theta
  returned_state_matches <- is.list(returned_theta) &&
    identical(as.numeric(returned_theta$beta), as.numeric(state$beta)) &&
    identical(as.numeric(returned_theta$subVar), as.numeric(state$sigma2)) &&
    identical(as.numeric(returned_theta$cellVar), as.numeric(state$gamma))
  if (identical(endpoint, "constrained")) {
    returned_state_matches <- returned_state_matches &&
      identical(as.numeric(fit$beta_constrained), as.numeric(state$beta))
  }
  if (!isTRUE(returned_state_matches)) {
    return(fail("return_state_mismatch"))
  }

  if (!is.list(objective_context) ||
      !is.list(objective_context$ctx) ||
      !is.list(objective_context$posv) ||
      length(objective_context$gene_index) != 1L ||
      !is.finite(objective_context$gene_index) ||
      !identical(
        as.integer(objective_context$gene_index),
        result$gene_index
      )) {
    return(fail("objective_context_missing_or_mismatched"))
  }
  ctx <- objective_context$ctx
  posv <- objective_context$posv
  if (length(ctx$nb) != 1L || !is.finite(ctx$nb) ||
      as.integer(ctx$nb) != nb) {
    return(fail("objective_context_dimension_mismatch"))
  }
  endpoint_map <- tryCatch(
    if (identical(endpoint, "unrestricted")) {
      ln_joint_map(nb)
    } else {
      if (!is.list(fit$constraint) || is.null(fit$constraint$L) ||
          is.null(fit$constraint$b)) {
        stop("constraint map missing")
      }
      ln_joint_map(nb, fit$constraint$L, fit$constraint$b)
    },
    error = function(e) NULL
  )
  if (is.null(endpoint_map) ||
      length(theta_reduced_log) != endpoint_map$free_dim + 2L) {
    return(fail("objective_map_invalid"))
  }
  reevaluated_state <- tryCatch(
    ln_joint_derivative_state(
      theta_reduced_log,
      result$gene_index,
      posv,
      ctx,
      endpoint_map
    ),
    error = function(e) NULL
  )
  if (!is.list(reevaluated_state) || !is.list(reevaluated_state$raw) ||
      !is.finite(reevaluated_state$raw$value) ||
      !identical(
        reevaluated_state$objective_provenance,
        expected_provenance
      )) {
    return(fail("objective_reevaluation_failed"))
  }
  result$reevaluated_objective <- as.numeric(reevaluated_state$raw$value)
  result$stored_state_objective <- as.numeric(state$raw$value)
  result$diagnostics$objective_recheck_scaled <- ln_joint_scaled_max(
    result$stored_state_objective,
    result$reevaluated_objective
  )
  result$diagnostics$full_log_recheck_scaled <- ln_joint_scaled_max(
    state_full_log,
    reevaluated_state$full_log
  )
  result$diagnostics$U_full_recheck_scaled <- ln_joint_scaled_max(
    U_full,
    reevaluated_state$U_full
  )
  result$diagnostics$U_subject_full_recheck_scaled <- ln_joint_scaled_max(
    U_subject_full,
    reevaluated_state$U_subject_full
  )
  result$diagnostics$A_full_recheck_scaled <- ln_joint_scaled_max(
    A_full,
    reevaluated_state$A_full
  )
  result$diagnostics$U_reduced_recheck_scaled <- ln_joint_scaled_max(
    U_reduced,
    reevaluated_state$U_reduced
  )
  result$diagnostics$U_subject_reduced_recheck_scaled <-
    ln_joint_scaled_max(
      U_subject_reduced,
      reevaluated_state$U_subject_reduced
    )
  result$diagnostics$A_reduced_recheck_scaled <- ln_joint_scaled_max(
    A_reduced,
    reevaluated_state$A_reduced
  )
  recheck_diagnostics <- unlist(
    result$diagnostics[c(
      "full_log_recheck_scaled",
      "U_full_recheck_scaled",
      "U_subject_full_recheck_scaled",
      "A_full_recheck_scaled",
      "U_reduced_recheck_scaled",
      "U_subject_reduced_recheck_scaled",
      "A_reduced_recheck_scaled"
    )],
    use.names = FALSE
  )
  if (!is.finite(result$diagnostics$objective_recheck_scaled) ||
      result$diagnostics$objective_recheck_scaled >
        thresholds$objective_agreement) {
    return(fail("objective_reevaluation_mismatch"))
  }
  if (length(recheck_diagnostics) != 7L ||
      any(!is.finite(recheck_diagnostics)) ||
      result$diagnostics$full_log_recheck_scaled >
        thresholds$natural_log_coordinate ||
      any(recheck_diagnostics[-1L] > thresholds$covariance_identity)) {
    return(fail("derivative_state_reevaluation_mismatch"))
  }

  candidate_records <- fit$ln_joint$candidate_records
  if (!is.list(candidate_records) || !length(candidate_records) ||
      any(!vapply(candidate_records, is.list, logical(1)))) {
    return(fail("candidate_records_invalid"))
  }
  state_matches <- which(vapply(
    candidate_records,
    function(candidate) {
      isTRUE(candidate$certified) &&
        isTRUE(candidate$interior) &&
        !is.null(candidate$theta_full_log) &&
        identical(as.numeric(candidate$theta_full_log), state_full_log)
    },
    logical(1)
  ))
  if (!length(state_matches)) {
    return(fail("endpoint_certification_missing"))
  }
  provenance_matches <- state_matches[vapply(
    candidate_records[state_matches],
    function(candidate) {
      identical(candidate$objective_provenance, expected_provenance)
    },
    logical(1)
  )]
  if (!length(provenance_matches)) {
    return(fail("candidate_objective_provenance_mismatch"))
  }

  state_objective <- result$reevaluated_objective
  candidate_errors <- vapply(
    candidate_records[provenance_matches],
    function(candidate) {
      ln_joint_scaled_max(candidate$objective, state_objective)
    },
    numeric(1)
  )
  objective_matches <- provenance_matches[
    is.finite(candidate_errors) &
      candidate_errors <= thresholds$objective_agreement
  ]
  if (!length(objective_matches)) {
    result$diagnostics$candidate_objective_state_scaled <-
      if (any(is.finite(candidate_errors))) {
        min(candidate_errors[is.finite(candidate_errors)])
      } else {
        NA_real_
      }
    return(fail("objective_state_mismatch"))
  }
  candidate_index <- objective_matches[[1L]]
  candidate_objective <- as.numeric(
    candidate_records[[candidate_index]]$objective
  )
  result$state_objective <- state_objective
  result$candidate_objective <- candidate_objective
  result$certified_candidate_index <- as.integer(candidate_index)
  result$diagnostics$candidate_objective_state_scaled <-
    ln_joint_scaled_max(candidate_objective, state_objective)

  loglik <- as.numeric(fit$loglik)
  loglik_ln <- as.numeric(fit$loglik_ln)
  result$loglik <- loglik
  result$loglik_ln <- loglik_ln
  result$diagnostics$loglik_objective_scaled <-
    ln_joint_scaled_max(loglik, -state_objective)
  result$diagnostics$loglik_ln_objective_scaled <-
    ln_joint_scaled_max(loglik_ln, -state_objective)
  if (length(loglik) != 1L || length(loglik_ln) != 1L ||
      any(!is.finite(c(loglik, loglik_ln))) ||
      !is.finite(result$diagnostics$loglik_objective_scaled) ||
      !is.finite(result$diagnostics$loglik_ln_objective_scaled) ||
      result$diagnostics$loglik_objective_scaled >
        thresholds$objective_agreement ||
      result$diagnostics$loglik_ln_objective_scaled >
        thresholds$objective_agreement) {
    return(fail("loglik_objective_mismatch"))
  }

  result$diagnostics$score_sign_scaled <- ln_joint_scaled_max(
    U_full,
    -state$q_gradient_full
  )
  result$diagnostics$donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_subject_full),
    U_full
  )
  if (!all(is.finite(c(
    result$diagnostics$score_sign_scaled,
    result$diagnostics$donor_additivity_scaled
  ))) || result$diagnostics$score_sign_scaled >
      thresholds$donor_additivity) {
    return(fail("maximizing_score_sign_mismatch"))
  }
  if (result$diagnostics$donor_additivity_scaled >
      thresholds$donor_additivity) {
    return(fail("donor_additivity_failed"))
  }

  result$success <- TRUE
  result$failure_reason <- NA_character_
  result
}

ln_joint_pace_ratio_diagnostic <- function(unrestricted_fit,
                                           constrained_fit,
                                           null_basis = NULL,
                                           interest_basis = NULL,
                                           objective_context = NULL) {
  empty_matrix_diagnostics <- ln_joint_score_matrix_diagnostics(NULL)
  result <- list(
    diagnostic_available = FALSE,
    failure_reason = NA_character_,
    public_available = FALSE,
    public_status = "joint_ln_calibration_pending",
    statistic = NA_real_,
    adjusted_ratio = NA_real_,
    df = NA_integer_,
    theoretical_df = NA_integer_,
    theoretical_reference = paste(
      "chi-square first-order reference only under smooth interior",
      "independent-donor fixed-dimension regularity"
    ),
    method = "pace_invariant_joint_ln_working_objective_ratio_audit",
    no_numerical_repair = TRUE,
    wP = NA_real_,
    wU = NA_real_,
    qP = NA_real_,
    multiplier = NA_real_,
    scale = NA_real_,
    U_eff = NULL,
    U_eff_subject = NULL,
    H_eff = NULL,
    J_eff = NULL,
    coordinate_map = NULL,
    projection = NULL,
    objective_provenance = NULL,
    objective_instance = NULL,
    endpoints = list(unrestricted = NULL, constrained = NULL),
    score_diagnostic = NULL,
    diagnostics = list(
      same_criterion = FALSE,
      same_objective_instance = FALSE,
      same_objective = FALSE,
      same_gene = FALSE,
      same_dimension = FALSE,
      same_donor_count = FALSE,
      wP_loglik_identity_scaled = NA_real_,
      H_eff_projection_identity_scaled = NA_real_,
      wU_score_identity_scaled = NA_real_,
      wP_minus_qP = NA_real_,
      wP_over_qP = NA_real_,
      adjusted_minus_wU = NA_real_,
      pace_decomposition_identity_scaled = NA_real_,
      H_eff = empty_matrix_diagnostics,
      J_eff = empty_matrix_diagnostics
    )
  )
  fail <- function(reason) {
    result$failure_reason <- reason
    result
  }

  if (!is.list(unrestricted_fit)) {
    return(fail("unrestricted_fit_invalid"))
  }
  if (!length(unrestricted_fit)) {
    return(fail("unrestricted_pair_not_registered"))
  }
  if (!is.list(constrained_fit) || !length(constrained_fit)) {
    return(fail("constrained_fit_invalid"))
  }
  fit_instance <- function(fit) {
    if (!is.list(fit) || !is.list(fit$ln_joint)) {
      return(NULL)
    }
    fit$ln_joint$objective_instance
  }
  unrestricted_instance <- fit_instance(unrestricted_fit)
  constrained_instance <- fit_instance(constrained_fit)
  if (!is.null(unrestricted_instance) && !is.null(constrained_instance) &&
      !identical(unrestricted_instance, constrained_instance)) {
    return(fail("endpoint_objective_instance_mismatch"))
  }
  unrestricted_endpoint <- ln_joint_ratio_endpoint_audit(
    unrestricted_fit,
    "unrestricted",
    objective_context = objective_context
  )
  constrained_endpoint <- ln_joint_ratio_endpoint_audit(
    constrained_fit,
    "constrained",
    objective_context = objective_context
  )
  result$endpoints$unrestricted <- unrestricted_endpoint
  result$endpoints$constrained <- constrained_endpoint
  if (!isTRUE(unrestricted_endpoint$success)) {
    return(fail(unrestricted_endpoint$failure_reason))
  }
  if (!isTRUE(constrained_endpoint$success)) {
    return(fail(constrained_endpoint$failure_reason))
  }

  result$diagnostics$same_criterion <- identical(
    unrestricted_endpoint$objective_provenance,
    constrained_endpoint$objective_provenance
  )
  result$diagnostics$same_objective_instance <- identical(
    unrestricted_endpoint$objective_instance,
    constrained_endpoint$objective_instance
  )
  result$diagnostics$same_objective <-
    isTRUE(result$diagnostics$same_criterion) &&
    isTRUE(result$diagnostics$same_objective_instance)
  result$diagnostics$same_gene <- identical(
    unrestricted_endpoint$gene_index,
    constrained_endpoint$gene_index
  ) && !is.na(unrestricted_endpoint$gene_index)
  result$diagnostics$same_dimension <- identical(
    unrestricted_endpoint$n_beta,
    constrained_endpoint$n_beta
  )
  result$diagnostics$same_donor_count <- identical(
    unrestricted_endpoint$n_donors,
    constrained_endpoint$n_donors
  )
  if (!isTRUE(result$diagnostics$same_criterion)) {
    return(fail("endpoint_criterion_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_objective_instance)) {
    return(fail("endpoint_objective_instance_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_gene)) {
    return(fail("endpoint_gene_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_dimension)) {
    return(fail("endpoint_dimension_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_donor_count)) {
    return(fail("endpoint_donor_count_mismatch"))
  }
  result$objective_provenance <-
    unrestricted_endpoint$objective_provenance
  result$objective_instance <- unrestricted_endpoint$objective_instance

  wP <- 2 * (
    constrained_endpoint$state_objective -
      unrestricted_endpoint$state_objective
  )
  wP_loglik <- 2 * (
    unrestricted_endpoint$loglik_ln - constrained_endpoint$loglik_ln
  )
  result$wP <- wP
  result$diagnostics$wP_loglik_identity_scaled <- ln_joint_scaled_max(
    wP,
    wP_loglik
  )
  if (!is.finite(wP) ||
      !is.finite(result$diagnostics$wP_loglik_identity_scaled) ||
      result$diagnostics$wP_loglik_identity_scaled >
        ln_joint_thresholds()$objective_agreement) {
    return(fail("working_objective_ratio_identity_failed"))
  }
  if (wP < 0) {
    return(fail("negative_working_objective_ratio"))
  }

  score <- ln_joint_robust_score_diagnostic(
    constrained_fit,
    null_basis = null_basis,
    interest_basis = interest_basis
  )
  result$score_diagnostic <- score
  result$df <- result$theoretical_df <- score$df
  if (!isTRUE(score$diagnostic_available)) {
    return(fail(paste0(
      "constrained_score_",
      if (is.null(score$failure_reason) || is.na(score$failure_reason)) {
        "unavailable"
      } else {
        score$failure_reason
      }
    )))
  }
  result$coordinate_map <- score$coordinate_map
  result$projection <- score$projection
  result$U_eff <- score$U_eff
  result$U_eff_subject <- score$U_eff_subject
  result$J_eff <- score$J_eff
  result$diagnostics$J_eff <- score$diagnostics$J_eff

  state <- constrained_fit$ln_joint$state
  T <- score$coordinate_map$full
  H_coordinate <- crossprod(T, state$A_full %*% T)
  H_coordinate <- 0.5 * (H_coordinate + t(H_coordinate))
  r <- score$df
  interest <- seq_len(r)
  nuisance <- r + seq_len(nrow(H_coordinate) - r)
  H_psi_psi <- H_coordinate[interest, interest, drop = FALSE]
  H_chi_psi <- H_coordinate[nuisance, interest, drop = FALSE]
  H_chi_solution <- tryCatch(
    solve(score$H_chi_chi, H_chi_psi),
    error = function(e) NULL
  )
  if (is.null(H_chi_solution) || any(!is.finite(H_chi_solution))) {
    return(fail("H_chi_chi_solve_failed_for_H_eff"))
  }
  H_eff_schur <- H_psi_psi - score$H_psi_chi %*% H_chi_solution
  H_eff_projected <- score$projection %*% H_coordinate %*%
    t(score$projection)
  result$diagnostics$H_eff_projection_identity_scaled <-
    ln_joint_relative_frobenius(H_eff_schur, H_eff_projected)
  if (!is.finite(
    result$diagnostics$H_eff_projection_identity_scaled
  ) || result$diagnostics$H_eff_projection_identity_scaled >
      ln_joint_thresholds()$covariance_identity) {
    return(fail("H_eff_projection_identity_failed"))
  }
  H_eff <- 0.5 * (H_eff_schur + t(H_eff_schur))
  H_eff_diagnostics <- ln_joint_score_matrix_diagnostics(H_eff)
  result$H_eff <- H_eff
  result$diagnostics$H_eff <- H_eff_diagnostics
  if (!isTRUE(H_eff_diagnostics$symmetry_acceptable)) {
    return(fail("H_eff_asymmetric"))
  }
  if (!isTRUE(H_eff_diagnostics$full_rank)) {
    return(fail("H_eff_rank_deficient"))
  }
  if (!isTRUE(H_eff_diagnostics$positive_definite)) {
    return(fail("H_eff_not_positive_definite"))
  }
  if (!isTRUE(H_eff_diagnostics$condition_acceptable)) {
    return(fail("H_eff_ill_conditioned"))
  }

  qP_solution <- tryCatch(
    solve(H_eff, score$U_eff),
    error = function(e) NULL
  )
  wU_solution <- tryCatch(
    solve(score$J_eff, score$U_eff),
    error = function(e) NULL
  )
  if (is.null(qP_solution) || any(!is.finite(qP_solution))) {
    return(fail("H_eff_solve_failed"))
  }
  if (is.null(wU_solution) || any(!is.finite(wU_solution))) {
    return(fail("J_eff_solve_failed"))
  }
  qP <- as.numeric(crossprod(score$U_eff, qP_solution))
  wU_reconstructed <- as.numeric(crossprod(
    score$U_eff,
    wU_solution
  ))
  wU <- as.numeric(score$statistic)
  result$qP <- qP
  result$wU <- wU
  result$diagnostics$wU_score_identity_scaled <- ln_joint_scaled_max(
    wU,
    wU_reconstructed
  )
  if (!is.finite(result$diagnostics$wU_score_identity_scaled) ||
      result$diagnostics$wU_score_identity_scaled >
        ln_joint_thresholds()$covariance_identity) {
    return(fail("wU_score_identity_failed"))
  }
  if (!is.finite(qP) || qP <= .Machine$double.eps) {
    return(fail("qP_nonpositive_or_numerically_zero"))
  }
  if (!is.finite(wU) || wU <= .Machine$double.eps) {
    return(fail("wU_nonpositive_or_numerically_zero"))
  }

  multiplier <- wU / qP
  adjusted_ratio <- multiplier * wP
  result$diagnostics$wP_minus_qP <- wP - qP
  result$diagnostics$wP_over_qP <- wP / qP
  result$diagnostics$adjusted_minus_wU <- adjusted_ratio - wU
  pace_decomposition <- wU + multiplier * (wP - qP)
  result$diagnostics$pace_decomposition_identity_scaled <-
    ln_joint_scaled_max(adjusted_ratio, pace_decomposition)
  if (!is.finite(multiplier) || multiplier <= 0) {
    return(fail("nonpositive_pace_multiplier"))
  }
  if (!is.finite(adjusted_ratio) || adjusted_ratio < 0) {
    return(fail("invalid_adjusted_working_objective_ratio"))
  }
  if (!is.finite(
    result$diagnostics$pace_decomposition_identity_scaled
  ) || result$diagnostics$pace_decomposition_identity_scaled >
      ln_joint_thresholds()$covariance_identity) {
    return(fail("pace_decomposition_identity_failed"))
  }

  result$multiplier <- result$scale <- multiplier
  result$adjusted_ratio <- result$statistic <- adjusted_ratio
  result$diagnostic_available <- TRUE
  result$failure_reason <- NA_character_
  result
}

ln_joint_root_metric <- function(U_subject, U) {
  if (is.null(U_subject) || !is.matrix(U_subject) || any(!is.finite(U_subject)) ||
      any(!is.finite(U))) {
    return(NA_real_)
  }
  denom <- 1 + sqrt(rowSums(U_subject^2))
  max(abs(U) / denom)
}

ln_joint_fd_gradient <- function(fn, theta, step_multiplier = 1) {
  theta <- as.numeric(theta)
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  vapply(seq_along(theta), function(jj) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    (fn(plus) - fn(minus)) / (2 * step[jj])
  }, numeric(1))
}

ln_joint_fd_jacobian <- function(fn, theta, step_multiplier = 1) {
  theta <- as.numeric(theta)
  base <- as.numeric(fn(theta))
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  out <- matrix(NA_real_, nrow = length(base), ncol = length(theta))
  for (jj in seq_along(theta)) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    out[, jj] <- (as.numeric(fn(plus)) - as.numeric(fn(minus))) /
      (2 * step[jj])
  }
  out
}

ln_joint_donor_objective_components <- function(theta_reduced_log, gene_index,
                                                 posv, ctx, map) {
  free_dim <- map$free_dim
  theta_reduced <- c(
    if (free_dim) theta_reduced_log[seq_len(free_dim)] else numeric(0),
    exp(theta_reduced_log[free_dim + seq_len(2L)])
  )
  full <- ln_joint_full_natural(theta_reduced, map)
  beta <- full[seq_len(ctx$nb)]
  sigma2 <- full[ctx$nb + 1L]
  gamma <- full[ctx$nb + 2L]
  exp_sigma2 <- exp(sigma2)
  alpha <- 1 / (exp_sigma2 - 1)
  lambda <- 1 / (sqrt(exp_sigma2) * (exp_sigma2 - 1))
  components <- lapply(seq_len(ctx$k), function(subject) {
    first_zero <- ctx$fid[subject] - 1L
    last_zero <- ctx$fid[subject + 1L] - 2L
    cell_index <- seq.int(first_zero + 1L, last_zero + 1L)
    positive <- which(
      posv$posindy >= first_zero & posv$posindy <= last_zero
    )
    donor_y <- posv$Y[positive]
    donor_posindy <- posv$posindy[positive] - first_zero
    donor_sum <- sum(donor_y)
    donor_n <- length(cell_index)
    eta <- as.numeric(
      ctx$offset[cell_index] +
        ctx$pred[cell_index, , drop = FALSE] %*% beta
    )
    mu <- exp(eta)
    sum_mu <- sum(mu)
    total <- sum_mu + lambda
    donor_weight <- (donor_sum + alpha) / total
    log_cell_total <- log(donor_weight * mu + gamma)
    fn_components <- c(
      positive_linear = if (length(positive))
        sum(eta[donor_posindy + 1L] * donor_y) else 0,
      donor_log_total = -(donor_sum + alpha) * log(total),
      alpha_log_lambda = alpha * log(lambda),
      cell_gamma_normalizer = donor_n * gamma * log(gamma),
      donor_mean_term = donor_weight * sum_mu,
      cell_log_normalizer = -gamma * sum(log_cell_total),
      positive_cell_log = if (length(positive))
        -sum(donor_y * log_cell_total[donor_posindy + 1L]) else 0,
      donor_lgamma = if (donor_sum > 0)
        lgamma(donor_sum + alpha) - lgamma(alpha) else 0,
      cell_lgamma = if (length(positive))
        sum(lgamma(donor_y + gamma) - lgamma(gamma)) else 0
    )
    -fn_components
  })
  out <- do.call(rbind, components)
  rownames(out) <- paste0("donor_", seq_len(ctx$k))
  out
}

ln_joint_donor_objective_values <- function(theta_reduced_log, gene_index,
                                             posv, ctx, map) {
  rowSums(ln_joint_donor_objective_components(
    theta_reduced_log,
    gene_index,
    posv,
    ctx,
    map
  ))
}

ln_joint_fd_gradient_donor_matrix <- function(fn, theta,
                                              step_multiplier = 1) {
  theta <- as.numeric(theta)
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  donor_count <- nrow(as.matrix(fn(theta)))
  out <- vapply(seq_along(theta), function(jj) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    rowSums(
      (as.matrix(fn(plus)) - as.matrix(fn(minus))) / (2 * step[jj])
    )
  }, numeric(donor_count))
  t(out)
}

ln_joint_fd_gradient_donorwise <- function(fn, theta, step_multiplier = 1) {
  rowSums(ln_joint_fd_gradient_donor_matrix(fn, theta, step_multiplier))
}

ln_joint_fd_sensitivity_donorwise <- function(fn, theta,
                                               step_multiplier = 1) {
  theta <- as.numeric(theta)
  base <- as.matrix(fn(theta))
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  out <- matrix(NA_real_, nrow = nrow(base), ncol = length(theta))
  for (jj in seq_along(theta)) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    donor_jacobian <- (as.matrix(fn(plus)) - as.matrix(fn(minus))) /
      (2 * step[jj])
    out[, jj] <- -rowSums(donor_jacobian)
  }
  out
}

ln_joint_fd_audit <- function(theta_reduced_log, gene_index, posv, ctx, map) {
  thresholds <- ln_joint_thresholds()
  free_dim <- map$free_dim
  to_natural_reduced <- function(z) {
    c(
      if (free_dim) z[seq_len(free_dim)] else numeric(0),
      exp(z[free_dim + seq_len(2L)])
    )
  }
  objective_log <- function(z) {
    ln_joint_objective(to_natural_reduced(z), gene_index, posv, ctx, map)
  }
  q_gradient_log <- function(z) {
    ln_joint_derivative_state(z, gene_index, posv, ctx, map)$q_gradient_reduced
  }
  maximizing_score_log <- function(z) {
    ln_joint_derivative_state(z, gene_index, posv, ctx, map)$U_reduced
  }
  donor_objective_log <- function(z) {
    ln_joint_donor_objective_components(z, gene_index, posv, ctx, map)
  }
  donor_score_log <- function(z) {
    ln_joint_derivative_state(z, gene_index, posv, ctx, map)$U_subject_reduced
  }

  state <- ln_joint_derivative_state(
    theta_reduced_log,
    gene_index,
    posv,
    ctx,
    map
  )
  fd_gradient_global_h <- ln_joint_fd_gradient(
    objective_log,
    theta_reduced_log,
    1
  )
  fd_gradient_global_h2 <- ln_joint_fd_gradient(
    objective_log,
    theta_reduced_log,
    0.5
  )
  fd_donor_gradient_h <- ln_joint_fd_gradient_donor_matrix(
    donor_objective_log,
    theta_reduced_log,
    1
  )
  fd_donor_gradient_h2 <- ln_joint_fd_gradient_donor_matrix(
    donor_objective_log,
    theta_reduced_log,
    0.5
  )
  fd_gradient_h <- rowSums(fd_donor_gradient_h)
  fd_gradient_h2 <- rowSums(fd_donor_gradient_h2)
  fd_sensitivity_global_h <- -ln_joint_fd_jacobian(
    maximizing_score_log,
    theta_reduced_log,
    1
  )
  fd_sensitivity_global_h2 <- -ln_joint_fd_jacobian(
    maximizing_score_log,
    theta_reduced_log,
    0.5
  )
  fd_sensitivity_h <- ln_joint_fd_sensitivity_donorwise(
    donor_score_log,
    theta_reduced_log,
    1
  )
  fd_sensitivity_h2 <- ln_joint_fd_sensitivity_donorwise(
    donor_score_log,
    theta_reduced_log,
    0.5
  )

  full_natural <- c(state$beta, state$sigma2, state$gamma)
  natural_gradient <- function(theta) {
    ptmg_der(
      theta,
      X = ctx$pred,
      offset = ctx$offset,
      Y = posv$Y,
      n_one = posv$n_onetwo,
      ytwo = posv$ytwo,
      fam = ctx$id,
      fid = ctx$fid,
      cumsumy = ctx$cumsumy[gene_index, ],
      posind = ctx$posind[[gene_index]],
      posindy = posv$posindy,
      nb = ctx$nb,
      nind = ctx$nind,
      k = ctx$k
    )
  }
  q_gradient_natural <- natural_gradient(full_natural)
  natural_sensitivity_fd <- ln_joint_fd_jacobian(
    natural_gradient,
    full_natural,
    1
  )
  natural_sensitivity_fd_h2 <- ln_joint_fd_jacobian(
    natural_gradient,
    full_natural,
    0.5
  )
  natural_sensitivity <- state$raw$hessian_natural
  q_gradient_natural_analytic <- state$raw$gradient_natural
  J <- diag(c(rep(1, ctx$nb), state$sigma2, state$gamma))
  q_gradient_log_from_natural <- as.numeric(crossprod(
    J,
    q_gradient_natural
  ))
  q_hessian_log_from_natural <- crossprod(J, natural_sensitivity %*% J)
  q_hessian_log_from_natural[ctx$nb + 1L, ctx$nb + 1L] <-
    q_hessian_log_from_natural[ctx$nb + 1L, ctx$nb + 1L] +
    q_gradient_natural[ctx$nb + 1L] * state$sigma2
  q_hessian_log_from_natural[ctx$nb + 2L, ctx$nb + 2L] <-
    q_hessian_log_from_natural[ctx$nb + 2L, ctx$nb + 2L] +
    q_gradient_natural[ctx$nb + 2L] * state$gamma

  metrics <- c(
    objective_value = abs(objective_log(theta_reduced_log) - state$raw$value) /
      (1 + max(abs(objective_log(theta_reduced_log)), abs(state$raw$value))),
    objective_additivity = abs(
      sum(donor_objective_log(theta_reduced_log)) -
        objective_log(theta_reduced_log)
    ) / (1 + abs(objective_log(theta_reduced_log))),
    objective_gradient = ln_joint_scaled_max(
      state$q_gradient_reduced,
      fd_gradient_h
    ),
    donor_objective_gradient = ln_joint_scaled_max(
      -state$U_subject_reduced,
      fd_donor_gradient_h
    ),
    sensitivity_score = ln_joint_scaled_max(
      state$A_reduced,
      fd_sensitivity_h
    ),
    fd_gradient_stability = ln_joint_scaled_max(fd_gradient_h, fd_gradient_h2),
    fd_donor_gradient_stability = ln_joint_scaled_max(
      fd_donor_gradient_h,
      fd_donor_gradient_h2
    ),
    fd_sensitivity_stability = ln_joint_scaled_max(
      fd_sensitivity_h,
      fd_sensitivity_h2
    ),
    donor_additivity = ln_joint_scaled_max(
      rowSums(state$U_subject_reduced),
      state$U_reduced
    ),
    aggregate_formula_vs_donor_sum = ln_joint_scaled_max(
      state$raw$gradient_aggregate_formula,
      state$q_gradient_full
    ),
    natural_log_gradient = ln_joint_scaled_max(
      state$q_gradient_full,
      q_gradient_log_from_natural
    ),
    natural_log_sensitivity = ln_joint_scaled_max(
      state$A_full,
      q_hessian_log_from_natural
    ),
    natural_gradient_implementation = ln_joint_scaled_max(
      q_gradient_natural_analytic,
      q_gradient_natural
    ),
    natural_sensitivity_score = ln_joint_scaled_max(
      natural_sensitivity,
      natural_sensitivity_fd
    ),
    fd_natural_sensitivity_stability = ln_joint_scaled_max(
      natural_sensitivity_fd,
      natural_sensitivity_fd_h2
    )
  )
  target <- c(
    objective_value = thresholds$objective_agreement,
    objective_additivity = thresholds$objective_agreement,
    objective_gradient = thresholds$objective_gradient,
    donor_objective_gradient = thresholds$objective_gradient,
    sensitivity_score = thresholds$sensitivity_score,
    fd_gradient_stability = thresholds$finite_difference_stability,
    fd_donor_gradient_stability = thresholds$finite_difference_stability,
    fd_sensitivity_stability = thresholds$finite_difference_stability,
    donor_additivity = thresholds$donor_additivity,
    aggregate_formula_vs_donor_sum = thresholds$donor_additivity,
    natural_log_gradient = thresholds$natural_log_coordinate,
    natural_log_sensitivity = thresholds$natural_log_coordinate,
    natural_gradient_implementation = thresholds$objective_gradient,
    natural_sensitivity_score = thresholds$sensitivity_score,
    fd_natural_sensitivity_stability =
      thresholds$finite_difference_stability
  )
  hard <- c(
    objective_value = thresholds$hard_stop$objective_agreement,
    objective_additivity = thresholds$hard_stop$objective_agreement,
    objective_gradient = thresholds$hard_stop$objective_gradient,
    donor_objective_gradient = thresholds$hard_stop$objective_gradient,
    sensitivity_score = thresholds$hard_stop$sensitivity_score,
    fd_gradient_stability = thresholds$hard_stop$finite_difference_stability,
    fd_donor_gradient_stability =
      thresholds$hard_stop$finite_difference_stability,
    fd_sensitivity_stability = thresholds$hard_stop$finite_difference_stability,
    donor_additivity = thresholds$hard_stop$donor_additivity,
    aggregate_formula_vs_donor_sum = thresholds$hard_stop$donor_additivity,
    natural_log_gradient = thresholds$hard_stop$natural_log_coordinate,
    natural_log_sensitivity = thresholds$hard_stop$natural_log_coordinate,
    natural_gradient_implementation = thresholds$hard_stop$objective_gradient,
    natural_sensitivity_score = thresholds$hard_stop$sensitivity_score,
    fd_natural_sensitivity_stability =
      thresholds$hard_stop$finite_difference_stability
  )
  list(
    metrics = metrics,
    roundoff_metrics = c(
      global_objective_gradient = ln_joint_scaled_max(
        state$q_gradient_reduced,
        fd_gradient_global_h
      ),
      global_objective_gradient_stability = ln_joint_scaled_max(
        fd_gradient_global_h,
        fd_gradient_global_h2
      ),
      global_sensitivity_score = ln_joint_scaled_max(
        state$A_reduced,
        fd_sensitivity_global_h
      ),
      global_sensitivity_stability = ln_joint_scaled_max(
        fd_sensitivity_global_h,
        fd_sensitivity_global_h2
      )
    ),
    target = target,
    hard_stop = hard,
    target_pass = is.finite(metrics) & metrics <= target,
    hard_stop_pass = is.finite(metrics) & metrics <= hard,
    fd_gradient_h = fd_gradient_h,
    fd_gradient_h2 = fd_gradient_h2,
    fd_donor_gradient_h = fd_donor_gradient_h,
    fd_donor_gradient_h2 = fd_donor_gradient_h2,
    fd_gradient_global_h = fd_gradient_global_h,
    fd_gradient_global_h2 = fd_gradient_global_h2,
    fd_sensitivity_h = fd_sensitivity_h,
    fd_sensitivity_h2 = fd_sensitivity_h2,
    fd_sensitivity_global_h = fd_sensitivity_global_h,
    fd_sensitivity_global_h2 = fd_sensitivity_global_h2,
    analytic_q_gradient = state$q_gradient_reduced,
    analytic_sensitivity = state$A_reduced,
    natural_q_gradient = q_gradient_natural,
    natural_sensitivity = natural_sensitivity,
    natural_sensitivity_fd = natural_sensitivity_fd,
    natural_sensitivity_fd_h2 = natural_sensitivity_fd_h2,
    q_gradient_log_from_natural = q_gradient_log_from_natural,
    q_hessian_log_from_natural = q_hessian_log_from_natural
  )
}

ln_joint_constraint_metric <- function(beta, map) {
  if (!isTRUE(map$constrained)) {
    return(0)
  }
  L <- as.matrix(map$L)
  b <- as.numeric(map$b)
  residual <- as.numeric(L %*% beta - b)
  denom <- 1 + max(abs(b)) + max(rowSums(abs(L))) * max(abs(beta))
  max(abs(residual)) / denom
}

ln_joint_candidate_state <- function(theta_reduced, objective, optimizer, label,
                                     start_provenance, required, core,
                                     gene_index, posv, ctx, map,
                                     objective_instance = NULL) {
  thresholds <- ln_joint_thresholds()
  free_dim <- map$free_dim
  metadata <- list(
    n_beta = as.integer(ctx$nb),
    parameterization = "beta, log(sigma2), log(gamma_nb_size)",
    objective_provenance = ln_joint_objective_provenance(),
    objective_instance = objective_instance,
    bounds = list(
      lower = as.numeric(ctx$min),
      upper = as.numeric(ctx$max)
    ),
    probe_bounds = list(
      lower = c(if (free_dim) rep(-100, free_dim) else numeric(0),
                log(ctx$min)),
      upper = c(if (free_dim) rep(100, free_dim) else numeric(0),
                log(ctx$max))
    ),
    probe_map = ln_joint_probe_map_record(map, ctx$nb)
  )
  finite_endpoint <- length(theta_reduced) == free_dim + 2L &&
    all(is.finite(theta_reduced)) && is.finite(objective) &&
    all(theta_reduced[free_dim + seq_len(2L)] > 0)
  if (!finite_endpoint) {
    return(list(
      label = label,
      start_provenance = start_provenance,
      required = required,
      core = core,
      n_beta = metadata$n_beta,
      parameterization = metadata$parameterization,
      objective_provenance = metadata$objective_provenance,
      objective_instance = metadata$objective_instance,
      bounds = metadata$bounds,
      probe_bounds = metadata$probe_bounds,
      probe_map = metadata$probe_map,
      finite_endpoint = FALSE,
      certified = FALSE,
      objective = objective,
      optimizer = optimizer,
      failure = "nonfinite_endpoint"
    ))
  }
  theta_log <- c(
    if (free_dim) theta_reduced[seq_len(free_dim)] else numeric(0),
    log(theta_reduced[free_dim + seq_len(2L)])
  )
  state <- tryCatch(
    ln_joint_derivative_state(theta_log, gene_index, posv, ctx, map),
    error = function(e) NULL
  )
  if (is.null(state)) {
    return(list(
      label = label,
      start_provenance = start_provenance,
      required = required,
      core = core,
      n_beta = metadata$n_beta,
      parameterization = metadata$parameterization,
      objective_provenance = metadata$objective_provenance,
      objective_instance = metadata$objective_instance,
      bounds = metadata$bounds,
      probe_bounds = metadata$probe_bounds,
      probe_map = metadata$probe_map,
      finite_endpoint = TRUE,
      certified = FALSE,
      theta_reduced = theta_reduced,
      theta_reduced_log = theta_log,
      objective = objective,
      optimizer = optimizer,
      failure = "derivative_evaluation_failed"
    ))
  }
  donor_sum <- rowSums(state$U_subject_reduced)
  donor_additivity <- ln_joint_scaled_max(donor_sum, state$U_reduced)
  root <- ln_joint_root_metric(state$U_subject_reduced, state$U_reduced)
  curvature <- ln_joint_curvature_diagnostics(state$A_reduced)
  newton_step <- tryCatch(
    solve(state$A_reduced, state$U_reduced),
    error = function(e) rep(NA_real_, length(state$U_reduced))
  )
  newton <- if (all(is.finite(newton_step))) {
    max(abs(newton_step) / (1 + abs(theta_log)))
  } else {
    NA_real_
  }
  dispersion <- theta_reduced[free_dim + seq_len(2L)]
  boundary_tolerance <- sqrt(.Machine$double.eps) *
    (1 + pmax(abs(ctx$min), abs(ctx$max)))
  interior <- all(
    dispersion > ctx$min + boundary_tolerance &
      dispersion < ctx$max - boundary_tolerance
  )
  constraint <- ln_joint_constraint_metric(state$beta, map)
  checks <- c(
    donor_additivity = is.finite(donor_additivity) &&
      donor_additivity <= thresholds$donor_additivity,
    root = is.finite(root) && root <= thresholds$root,
    newton = is.finite(newton) && newton <= thresholds$newton,
    curvature = isTRUE(curvature$positive_definite) &&
      is.finite(curvature$scaled_min_eigenvalue) &&
      curvature$scaled_min_eigenvalue >= thresholds$curvature_min_eigen,
    condition = is.finite(curvature$scaled_condition_number) &&
      curvature$scaled_condition_number <= thresholds$condition_number,
    interior = interior,
    constraint = is.finite(constraint) && constraint <= thresholds$constraint
  )
  failure <- NA_character_
  if (!checks[["donor_additivity"]]) {
    failure <- if (is.finite(donor_additivity) &&
        donor_additivity <= thresholds$hard_stop$donor_additivity) {
      "donor_additivity_target_failed"
    } else {
      "donor_additivity_hard_stop"
    }
  } else if (!checks[["interior"]]) {
    failure <- "dispersion_boundary"
  } else if (!checks[["root"]]) {
    failure <- if (is.finite(root) && root <= thresholds$hard_stop$root) {
      "root_residual_review_zone"
    } else {
      "root_residual_hard_stop"
    }
  } else if (!checks[["newton"]]) {
    failure <- if (is.finite(newton) && newton <= thresholds$hard_stop$newton) {
      "newton_displacement_review_zone"
    } else {
      "newton_displacement_hard_stop"
    }
  } else if (!isTRUE(curvature$positive_definite) ||
      !is.finite(curvature$scaled_min_eigenvalue) ||
      curvature$scaled_min_eigenvalue <= 0) {
    failure <- "nonpositive_joint_curvature"
  } else if (!checks[["curvature"]]) {
    failure <- "curvature_review_zone"
  } else if (!checks[["condition"]]) {
    failure <- if (is.finite(curvature$scaled_condition_number) &&
        curvature$scaled_condition_number <= thresholds$hard_stop$condition_number) {
      "curvature_condition_review_zone"
    } else {
      "curvature_condition_unusable"
    }
  } else if (!checks[["constraint"]]) {
    failure <- if (is.finite(constraint) &&
        constraint <= thresholds$hard_stop$constraint) {
      "constraint_review_zone"
    } else {
      "constraint_hard_stop"
    }
  }
  list(
    label = label,
    start_provenance = start_provenance,
    required = required,
    core = core,
    n_beta = metadata$n_beta,
    parameterization = metadata$parameterization,
    objective_instance = metadata$objective_instance,
    bounds = metadata$bounds,
    probe_bounds = metadata$probe_bounds,
    probe_map = metadata$probe_map,
    finite_endpoint = TRUE,
    certified = all(checks),
    theta_reduced = theta_reduced,
    theta_reduced_log = theta_log,
    theta_full_log = state$full_log,
    objective_provenance = state$objective_provenance,
    objective = objective,
    optimizer = optimizer,
    state = state,
    donor_additivity_scaled = donor_additivity,
    root_residual = root,
    newton_displacement = newton,
    curvature = curvature,
    interior = interior,
    constraint_residual = constraint,
    checks = checks,
    failure = failure
  )
}

ln_joint_polish <- function(theta_reduced, gene_index, posv, ctx, map,
                            lower, upper, max_iter = 25L) {
  free_dim <- map$free_dim
  z <- c(
    if (free_dim) theta_reduced[seq_len(free_dim)] else numeric(0),
    log(theta_reduced[free_dim + seq_len(2L)])
  )
  lower_log <- c(if (free_dim) lower[seq_len(free_dim)] else numeric(0), log(ctx$min))
  upper_log <- c(if (free_dim) upper[seq_len(free_dim)] else numeric(0), log(ctx$max))
  objective_at <- function(zz) {
    natural <- c(
      if (free_dim) zz[seq_len(free_dim)] else numeric(0),
      exp(zz[free_dim + seq_len(2L)])
    )
    ln_joint_objective(natural, gene_index, posv, ctx, map)
  }
  current_obj <- objective_at(z)
  for (iter in seq_len(max_iter)) {
    state <- tryCatch(
      ln_joint_derivative_state(z, gene_index, posv, ctx, map),
      error = function(e) NULL
    )
    if (is.null(state)) break
    metric <- ln_joint_root_metric(state$U_subject_reduced, state$U_reduced)
    if (is.finite(metric) && metric <= ln_joint_thresholds()$root / 10) break
    step <- tryCatch(
      solve(state$A_reduced, state$q_gradient_reduced),
      error = function(e) NULL
    )
    if (is.null(step) || any(!is.finite(step))) break
    accepted <- FALSE
    for (half in 0:20) {
      candidate <- z - step / (2^half)
      if (any(candidate <= lower_log) || any(candidate >= upper_log)) next
      candidate_obj <- tryCatch(objective_at(candidate), error = function(e) NA_real_)
      candidate_state <- tryCatch(
        ln_joint_derivative_state(candidate, gene_index, posv, ctx, map),
        error = function(e) NULL
      )
      if (!is.finite(candidate_obj) || is.null(candidate_state)) next
      candidate_metric <- ln_joint_root_metric(
        candidate_state$U_subject_reduced,
        candidate_state$U_reduced
      )
      objective_ok <- candidate_obj <= current_obj +
        1e-12 * (1 + abs(current_obj))
      if (objective_ok && is.finite(candidate_metric) &&
          (!is.finite(metric) || candidate_metric < metric)) {
        z <- candidate
        current_obj <- candidate_obj
        accepted <- TRUE
        break
      }
    }
    if (!accepted) break
  }
  list(
    theta_reduced = c(
      if (free_dim) z[seq_len(free_dim)] else numeric(0),
      exp(z[free_dim + seq_len(2L)])
    ),
    objective = current_obj
  )
}

ln_joint_clamp <- function(x, lower, upper) {
  gap <- 1e-6 * (upper - lower)
  pmin(pmax(x, lower + gap), upper - gap)
}

ln_joint_pql_start <- function(best, gene_index, posv, ctx, map) {
  full <- ln_joint_full_natural(best$theta_reduced, map)
  beta <- full[seq_len(ctx$nb)]
  sigma2 <- full[ctx$nb + 1L]
  gamma <- full[ctx$nb + 2L]
  gni <- (ctx$nind / ctx$k) * gamma
  cv2p <- if (ctx$ncell > 0L) {
    get_cv(ctx$offset, ctx$pred, beta, ctx$cell_ind, ctx$ncell, ctx$nind)
  } else {
    ctx$cv2
  }
  kappa_obs <- gni / (1 + cv2p)
  pair_trigger <- gni < ctx$cutoff_cell
  sigma_trigger <- !pair_trigger &&
    ((kappa_obs < 20) || ((kappa_obs < ctx$kappa) &&
      (sigma2 < (8 / max(kappa_obs, 1e-8)))))
  if (!pair_trigger && !sigma_trigger) {
    return(list(triggered = FALSE, status = "not_triggered"))
  }
  ord <- if ((posv$mct * (ctx$nind / ctx$k)) < 3) 3 else 1
  common <- list(
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    n_one = posv$n_onetwo,
    ytwo = posv$ytwo,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[gene_index, ],
    posind = ctx$posind[[gene_index]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k,
    betas = beta,
    reml = 0,
    eps = ctx$eps,
    ord = ord,
    intcol = ctx$intcol
  )
  dispersion <- tryCatch({
    if (pair_trigger) {
      fit <- do.call(
        bobyqa,
        c(list(x0 = c(sigma2, gamma), fn = pql_ll,
               lower = ctx$min, upper = ctx$max), common)
      )
      as.numeric(fit$par[seq_len(2L)])
    } else {
      fit <- do.call(
        nlminb,
        c(list(start = sigma2, objective = pql_gamma_ll,
               lower = ctx$min[1L], upper = ctx$max[1L], gamma = gamma), common)
      )
      c(as.numeric(fit$par[[1L]]), gamma)
    }
  }, error = function(e) NULL)
  if (is.null(dispersion) || length(dispersion) != 2L ||
      any(!is.finite(dispersion))) {
    return(list(
      triggered = TRUE,
      status = "generation_failed",
      method = if (pair_trigger) "pql_pair" else "pql_sigma"
    ))
  }
  list(
    triggered = TRUE,
    status = "generated",
    method = if (pair_trigger) "pql_pair" else "pql_sigma",
    start = c(map$project(beta), dispersion),
    gni = gni,
    kappa_obs = kappa_obs,
    cv2p = cv2p
  )
}

ln_joint_candidate_summary <- function(candidates, candidate_resolutions = NULL) {
  do.call(rbind, lapply(seq_along(candidates), function(index) {
    x <- candidates[[index]]
    if (!is.list(x) || is.data.frame(x)) {
      x <- list(
        label = NA_character_,
        start_provenance = NA_character_,
        optimizer = list()
      )
    }
    if (!is.list(x$optimizer)) x$optimizer <- list()
    curvature <- if (is.list(x$curvature)) x$curvature else list()
    resolution <- if (is.list(candidate_resolutions) &&
        length(candidate_resolutions) >= index) {
      candidate_resolutions[[index]]
    } else {
      NULL
    }
    probe <- if (is.list(resolution)) resolution$probe else NULL
    data.frame(
      label = ln_joint_candidate_label(x),
      start_provenance = ln_joint_diagnostic_character(x$start_provenance),
      required = isTRUE(x$required),
      core = isTRUE(x$core),
      finite_endpoint = isTRUE(x$finite_endpoint),
      certified = isTRUE(x$certified),
      optimizer_convergence = ln_joint_diagnostic_integer(
        x$optimizer$convergence
      ),
      optimizer_message = ln_joint_diagnostic_character(
        x$optimizer$message
      ),
      start = ln_joint_diagnostic_vector_string(x$start),
      optimizer_endpoint = ln_joint_diagnostic_vector_string(
        x$optimizer_endpoint
      ),
      optimizer_objective = ln_joint_diagnostic_numeric(
        x$optimizer_objective
      ),
      polished = isTRUE(x$polished),
      theta_reduced = ln_joint_diagnostic_vector_string(x$theta_reduced),
      theta_full_log = ln_joint_diagnostic_vector_string(x$theta_full_log),
      objective = if (ln_joint_finite_scalar_objective(x)) {
        as.numeric(x$objective)
      } else {
        NA_real_
      },
      root_residual = ln_joint_diagnostic_numeric(x$root_residual),
      newton_displacement = ln_joint_diagnostic_numeric(
        x$newton_displacement
      ),
      donor_additivity_scaled = ln_joint_diagnostic_numeric(
        x$donor_additivity_scaled
      ),
      scaled_min_eigenvalue = ln_joint_diagnostic_numeric(
        curvature$scaled_min_eigenvalue
      ),
      scaled_condition_number = ln_joint_diagnostic_numeric(
        curvature$scaled_condition_number
      ),
      failure = ln_joint_diagnostic_character(x$failure),
      resolution_class = if (is.list(resolution)) {
        ln_joint_diagnostic_character(resolution$class)
      } else {
        NA_character_
      },
      resolution_veto = if (is.list(resolution)) isTRUE(resolution$veto) else NA,
      objective_relation = if (is.list(resolution)) {
        ln_joint_diagnostic_character(resolution$objective_relation)
      } else {
        NA_character_
      },
      objective_signed_gap = if (is.list(resolution)) {
        ln_joint_diagnostic_numeric(resolution$objective_signed_gap)
      } else {
        NA_real_
      },
      full_state_distance = if (is.list(resolution)) {
        ln_joint_diagnostic_numeric(resolution$metrics$full_state)
      } else {
        NA_real_
      },
      beta_block_distance = if (is.list(resolution)) {
        ln_joint_diagnostic_numeric(resolution$metrics$beta_block)
      } else {
        NA_real_
      },
      dispersion_block_distance = if (is.list(resolution)) {
        ln_joint_diagnostic_numeric(resolution$metrics$dispersion_block)
      } else {
        NA_real_
      },
      boundary_scaled_distance = if (is.list(resolution)) {
        ln_joint_diagnostic_numeric(
          resolution$metrics$boundary_scaled_distance
        )
      } else {
        NA_real_
      },
      curvature_compatibility = if (is.list(resolution)) {
        ln_joint_diagnostic_numeric(
          resolution$metrics$curvature_compatibility
        )
      } else {
        NA_real_
      },
      probe_accepted = if (is.list(probe)) isTRUE(probe$accepted) else NA,
      probe_termination_reason = if (is.list(probe)) {
        ln_joint_diagnostic_character(probe$termination_reason)
      } else {
        NA_character_
      },
      probe_state_evaluations = if (is.list(probe)) {
        ln_joint_diagnostic_integer(probe$state_evaluations)
      } else {
        NA_integer_
      },
      probe_reused = if (is.list(resolution)) {
        isTRUE(resolution$probe_reused)
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  }))
}

ln_joint_candidate_records <- function(candidates, candidate_resolutions = NULL) {
  lapply(seq_along(candidates), function(index) {
    x <- candidates[[index]]
    if (!is.list(x) || is.data.frame(x)) {
      x <- list(
        label = NA_character_,
        start_provenance = NA_character_,
        optimizer = list()
      )
    }
    if (!is.list(x$optimizer)) x$optimizer <- list()
    resolution <- if (is.list(candidate_resolutions) &&
        length(candidate_resolutions) >= index) {
      candidate_resolutions[[index]]
    } else {
      NULL
    }
    list(
      label = ln_joint_candidate_label(x),
      start_provenance = ln_joint_diagnostic_character(x$start_provenance),
      required = isTRUE(x$required),
      core = isTRUE(x$core),
      start = x$start,
      optimizer_endpoint = x$optimizer_endpoint,
      optimizer_objective = x$optimizer_objective,
      polished = isTRUE(x$polished),
      theta_reduced = x$theta_reduced,
      theta_reduced_log = x$theta_reduced_log,
      theta_full_log = x$theta_full_log,
      n_beta = x$n_beta,
      parameterization = x$parameterization,
      objective_provenance = x$objective_provenance,
      objective_instance = x$objective_instance,
      bounds = x$bounds,
      probe_bounds = x$probe_bounds,
      probe_map = x$probe_map,
      finite_endpoint = isTRUE(x$finite_endpoint),
      certified = isTRUE(x$certified),
      objective = x$objective,
      optimizer_convergence = x$optimizer$convergence,
      optimizer_message = x$optimizer$message,
      donor_additivity_scaled = x$donor_additivity_scaled,
      root_residual = x$root_residual,
      newton_displacement = x$newton_displacement,
      curvature = x$curvature,
      interior = x$interior,
      constraint_residual = x$constraint_residual,
      checks = x$checks,
      failure = x$failure,
      resolution = resolution
    )
  })
}

ln_joint_resolve_candidates <- function(
    candidates,
    resolution_thresholds = ln_joint_resolution_thresholds(),
    basin_probe = NULL,
    enable_prospective = FALSE) {
  looks_like_one_candidate <- is.list(candidates) &&
    any(c("finite_endpoint", "certified", "objective") %in% names(candidates))
  if (!is.list(candidates) || is.data.frame(candidates) ||
      isTRUE(looks_like_one_candidate)) {
    candidates <- list(candidates)
  }
  thresholds <- ln_joint_thresholds()
  certified <- which(vapply(candidates, function(x) {
    is.list(x) && isTRUE(x$finite_endpoint) && isTRUE(x$certified) &&
      ln_joint_finite_scalar_objective(x)
  }, logical(1)))
  finite <- which(vapply(candidates, function(x) {
    is.list(x) && isTRUE(x$finite_endpoint) &&
      ln_joint_finite_scalar_objective(x)
  }, logical(1)))
  best_idx <- if (length(certified)) {
    certified[which.min(vapply(candidates[certified], `[[`, numeric(1), "objective"))]
  } else if (length(finite)) {
    finite[which.min(vapply(candidates[finite], `[[`, numeric(1), "objective"))]
  } else {
    NA_integer_
  }
  assessment_indices <- if (is.na(best_idx)) {
    integer(0)
  } else if (best_idx %in% certified) {
    c(best_idx, setdiff(certified, best_idx))
  } else {
    best_idx
  }
  assessment_roots <- candidates[assessment_indices]
  failure <- NA_character_
  parameter_agreement <- NA_real_
  objective_agreement <- NA_real_
  competing <- character(0)
  dominated <- character(0)
  certified_cluster_valid <- FALSE
  required_failed <- vapply(candidates, function(x) {
    is.list(x) && isTRUE(x$required) &&
      (!isTRUE(x$finite_endpoint) ||
        !ln_joint_finite_scalar_objective(x))
  }, logical(1))
  if (any(required_failed)) {
    failure <- "required_start_failed"
  } else if (!length(certified)) {
    failure <- if (!is.na(best_idx)) {
      ln_joint_diagnostic_character(
        candidates[[best_idx]]$failure,
        fallback = "no_certified_candidate"
      )
    } else {
      "no_certified_candidate"
    }
  } else {
    best <- candidates[[best_idx]]
    certified_integrity <- lapply(candidates[certified], function(x) {
      tryCatch({
        objective_integrity <- ln_joint_candidate_objective_integrity(x, best)
        state_integrity <- ln_joint_candidate_integrity(x, best)
        list(
          valid = isTRUE(objective_integrity$valid) &&
            isTRUE(state_integrity$valid),
          failed_checks = if (isTRUE(state_integrity$valid)) {
            state_integrity$failed_checks
          } else {
            character(0)
          }
        )
      }, error = function(e) {
        list(valid = FALSE, failed_checks = character(0))
      })
    })
    certified_valid <- vapply(certified_integrity, function(integrity) {
      isTRUE(integrity$valid) && !length(integrity$failed_checks)
    }, logical(1))
    certified_cluster_valid <- length(certified_valid) > 0L &&
      all(certified_valid)
    if (!isTRUE(certified_cluster_valid)) {
      best_position <- match(best_idx, certified)
      failure <- if (length(best_position) == 1L &&
          !is.na(best_position) && !isTRUE(certified_valid[[best_position]])) {
        "selected_candidate_certificate_inconsistent"
      } else {
        "certified_candidate_check_vector_inconsistent"
      }
    }
    if (isTRUE(certified_cluster_valid)) {
      core_certified <- sum(vapply(candidates, function(x) {
      is.list(x) && isTRUE(x$core) && isTRUE(x$finite_endpoint) &&
        isTRUE(x$certified) &&
        ln_joint_finite_scalar_objective(x)
      }, logical(1)))
      if (core_certified < 2L) {
        failure <- "insufficient_certified_core_starts"
      }
      parameter_agreement <- max(vapply(candidates[certified], function(x) {
        ln_joint_scaled_max(x$theta_full_log, best$theta_full_log)
      }, numeric(1)))
      objective_agreement <- max(vapply(candidates[certified], function(x) {
        abs(x$objective - best$objective) /
          (1 + max(abs(x$objective), abs(best$objective)))
      }, numeric(1)))
      if (is.na(failure) && (!is.finite(parameter_agreement) ||
          parameter_agreement > thresholds$parameter_agreement)) {
        failure <- if (is.finite(parameter_agreement) &&
            parameter_agreement <= thresholds$hard_stop$parameter_agreement) {
          "start_parameter_review_zone"
        } else {
          "start_parameter_hard_stop"
        }
      }
      if (is.na(failure) && (!is.finite(objective_agreement) ||
          objective_agreement > thresholds$objective_agreement)) {
        failure <- if (is.finite(objective_agreement) &&
            objective_agreement <= thresholds$hard_stop$objective_agreement) {
          "start_objective_review_zone"
        } else {
          "start_objective_hard_stop"
        }
      }
      uncertified_finite <- which(vapply(candidates, function(x) {
        is.list(x) && isTRUE(x$finite_endpoint) &&
          ln_joint_finite_scalar_objective(x) && !isTRUE(x$certified)
      }, logical(1)))
      if (length(uncertified_finite)) {
        classifications <- lapply(candidates[uncertified_finite], function(x) {
          diagnostic_label <- ln_joint_candidate_label(
            x,
            fallback = ln_joint_resolution_candidate_key(x)
          )
          relative_gap <- (x$objective - best$objective) /
            (1 + max(abs(x$objective), abs(best$objective)))
          if (!is.finite(relative_gap) ||
              relative_gap <= thresholds$hard_stop$objective_agreement) {
            list(competing = diagnostic_label, dominated = NA_character_)
          } else {
            list(competing = NA_character_, dominated = diagnostic_label)
          }
        })
        competing <- vapply(classifications, `[[`, character(1), "competing")
        competing <- na.omit(competing)
        dominated <- vapply(classifications, `[[`, character(1), "dominated")
        dominated <- na.omit(dominated)
        if (is.na(failure) && length(competing)) {
          failure <- "competing_uncertified_candidate"
        }
      }
    }
  }
  strict_route <- if (is.na(failure)) "R0" else "R2"
  strict_failure <- failure

  threshold_configuration_valid <-
    ln_joint_resolution_thresholds_valid(resolution_thresholds)
  candidate_resolutions <- vector("list", length(candidates))
  physical_probe_records <- list()
  if (length(candidates)) {
    if (!threshold_configuration_valid || is.na(best_idx) ||
        !isTRUE(certified_cluster_valid)) {
      candidate_resolutions <- lapply(candidates, function(candidate) {
        record <- ln_joint_candidate_resolution_record(candidate)
        if (ln_joint_optional_failed_attempt_integrity(candidate)) {
          ln_joint_finish_candidate_resolution(
            record, "noncompetitive_failed_attempt", FALSE,
            "optional_candidate_has_no_comparable_endpoint"
          )
        } else {
          ln_joint_finish_candidate_resolution(
            record, "competitive_hard_stop", TRUE,
            if (threshold_configuration_valid) {
              "selected_candidate_unavailable"
            } else {
              "resolution_threshold_configuration_invalid"
            }
          )
        }
      })
    } else {
      selected <- candidates[[best_idx]]
      candidate_resolutions <- lapply(seq_along(candidates), function(index) {
        tryCatch(
          ln_joint_classify_candidate_against_anchors(
            candidates[[index]],
            selected,
            assessment_roots,
            resolution_thresholds,
            is_selected = identical(index, best_idx)
          ),
          error = function(e) {
            ln_joint_finish_candidate_resolution(
              ln_joint_candidate_resolution_record(candidates[[index]]),
              "competitive_hard_stop",
              TRUE,
              "malformed_candidate_record"
            )
          }
        )
      })
      review_indices <- which(vapply(candidate_resolutions, function(x) {
        identical(x$class, "same_basin_review_candidate")
      }, logical(1)))
      if (length(review_indices)) {
        strict_failure_reviewable <- identical(
          strict_failure,
          "competing_uncertified_candidate"
        )
        substantive_veto <- any(vapply(
          candidate_resolutions,
          function(resolution) {
            isTRUE(resolution$veto) &&
              !identical(
                resolution$class,
                "same_basin_review_candidate"
              )
          },
          logical(1)
        ))
        if (!isTRUE(strict_failure_reviewable)) {
          for (index in review_indices) {
            probe <- ln_joint_basin_probe_template(resolution_thresholds)
            probe$termination_reason <-
              "skipped_strict_failure_not_reviewable"
            candidate_resolutions[[index]]$probe <- probe
            candidate_resolutions[[index]] <-
              ln_joint_finish_candidate_resolution(
                candidate_resolutions[[index]],
                "same_basin_review_unresolved",
                TRUE,
                "bounded_review_skipped_due_to_strict_failure"
              )
          }
        } else if (isTRUE(substantive_veto)) {
          for (index in review_indices) {
            probe <- ln_joint_basin_probe_template(resolution_thresholds)
            probe$termination_reason <-
              "skipped_immediate_substantive_veto"
            candidate_resolutions[[index]]$probe <- probe
            candidate_resolutions[[index]] <-
              ln_joint_finish_candidate_resolution(
                candidate_resolutions[[index]],
                "same_basin_review_unresolved",
                TRUE,
                "bounded_review_skipped_due_to_substantive_veto"
              )
          }
        } else {
          review_representatives <- integer(0)
          for (index in review_indices) {
            already_represented <- if (length(review_representatives)) {
              any(vapply(review_representatives, function(reference_index) {
                ln_joint_probe_equivalent(
                  candidates[[index]],
                  candidates[[reference_index]]
                )
              }, logical(1)))
            } else {
              FALSE
            }
            if (!already_represented) {
              review_representatives <- c(review_representatives, index)
            }
          }
          if (length(review_representatives) > 1L) {
            for (index in review_indices) {
              probe <- ln_joint_basin_probe_template(resolution_thresholds)
              probe$termination_reason <-
                "multiple_review_endpoints_exceed_fixed_probe_budget"
              candidate_resolutions[[index]]$probe <- probe
              candidate_resolutions[[index]] <-
                ln_joint_finish_candidate_resolution(
                  candidate_resolutions[[index]],
                  "same_basin_review_unresolved",
                  TRUE,
                  "multiple_distinct_review_endpoints_unresolved"
                )
            }
          } else {
            cached_candidates <- list()
            cached_probes <- list()
            for (index in review_indices) {
            reused_index <- if (length(cached_candidates)) {
              match(TRUE, vapply(cached_candidates, function(reference) {
                ln_joint_probe_equivalent(candidates[[index]], reference)
              }, logical(1)), nomatch = 0L)
            } else {
              0L
            }
            if (reused_index > 0L) {
              probe <- cached_probes[[reused_index]]
              candidate_resolutions[[index]]$probe_reused <- TRUE
            } else {
              probe_callback_available <- is.function(basin_probe)
              probe <- if (probe_callback_available) {
                tryCatch(
                  ln_joint_basin_probe(
                    candidate = candidates[[index]],
                    selected = selected,
                    evaluate_state = function(theta_reduced_log) {
                      basin_probe(candidates[[index]], theta_reduced_log)
                    },
                    resolution_thresholds = resolution_thresholds,
                    assessment_roots = assessment_roots
                  ),
                  error = function(e) NULL
                )
              } else {
                NULL
              }
              if (!ln_joint_basin_probe_record_valid(
                  probe,
                  candidates[[index]],
                  assessment_roots,
                  resolution_thresholds)) {
                probe <- ln_joint_basin_probe_template(
                  resolution_thresholds
                )
                probe$termination_reason <- if (probe_callback_available) {
                  "malformed_or_failed_bounded_review_record"
                } else {
                  "bounded_review_not_run"
                }
              }
              cached_candidates[[length(cached_candidates) + 1L]] <-
                candidates[[index]]
              cached_probes[[length(cached_probes) + 1L]] <- probe
              physical_probe_records[[length(physical_probe_records) + 1L]] <-
                probe
            }
            candidate_resolutions[[index]]$probe <- probe
            if (isTRUE(probe$accepted)) {
              candidate_resolutions[[index]] <-
                ln_joint_finish_candidate_resolution(
                  candidate_resolutions[[index]],
                  "same_basin_review_nonveto",
                  FALSE,
                  "bounded_same_basin_review_resolved"
                )
            } else {
              candidate_resolutions[[index]] <-
                ln_joint_finish_candidate_resolution(
                  candidate_resolutions[[index]],
                  "same_basin_review_unresolved",
                  TRUE,
                  "bounded_same_basin_review_unresolved"
                )
            }
            }
          }
        }
      }
    }
  }

  veto_classes <- if (length(candidate_resolutions)) {
    vapply(candidate_resolutions, function(x) {
      if (isTRUE(x$veto)) x$class else NA_character_
    }, character(1))
  } else {
    character(0)
  }
  global_veto_classes <- sort(unique(na.omit(veto_classes)))
  nonveto_count <- sum(vapply(candidate_resolutions, function(x) {
    identical(x$class, "same_basin_review_nonveto")
  }, logical(1)))
  would_nonveto <- identical(strict_route, "R2") &&
    identical(strict_failure, "competing_uncertified_candidate") &&
    nonveto_count > 0L && !length(global_veto_classes)
  priority <- c(
    competitive_boundary = "competitive_boundary_candidate",
    competitive_distinct_basin = "competitive_distinct_basin_candidate",
    competitive_hard_stop = "competitive_hard_stop_candidate",
    same_basin_review_unresolved = "same_basin_review_unresolved"
  )
  matched_veto <- names(priority)[names(priority) %in% global_veto_classes]
  prospective_veto_failure <- if (length(matched_veto)) {
    unname(priority[[matched_veto[[1L]]]])
  } else {
    NA_character_
  }
  prospective_route <- if (length(global_veto_classes)) {
    "R2"
  } else if (isTRUE(would_nonveto)) {
    "R0"
  } else {
    strict_route
  }
  prospective_failure <- if (length(global_veto_classes)) {
    if (!is.na(prospective_veto_failure)) {
      prospective_veto_failure
    } else {
      "competitive_hard_stop_candidate"
    }
  } else if (isTRUE(would_nonveto)) {
    NA_character_
  } else {
    strict_failure
  }
  numerical_resolution <- if (isTRUE(would_nonveto)) {
    "would_same_basin_review_nonveto"
  } else if ("competitive_boundary" %in% global_veto_classes) {
    "prospective_veto_competitive_boundary"
  } else if ("competitive_distinct_basin" %in% global_veto_classes) {
    "prospective_veto_competitive_distinct_basin"
  } else if ("competitive_hard_stop" %in% global_veto_classes) {
    "prospective_veto_competitive_hard_stop"
  } else if ("same_basin_review_unresolved" %in% global_veto_classes) {
    "prospective_veto_same_basin_review_unresolved"
  } else if (identical(strict_route, "R0")) {
    "strict_certified"
  } else {
    "strict_failure_not_candidate_reviewable"
  }
  probe_records <- Filter(is.list, physical_probe_records)
  evaluation_counts <- list(
    total_engine_counts_available = FALSE,
    total_engine_counts_reason =
      "existing_optimizer_and_polish_calls_are_not_instrumented",
    probe_state_evaluations = sum(vapply(probe_records, function(x) {
      as.integer(x$state_evaluations)
    }, integer(1))),
    probe_objective_evaluations = sum(vapply(probe_records, function(x) {
      as.integer(x$objective_evaluations)
    }, integer(1))),
    probe_score_evaluations = sum(vapply(probe_records, function(x) {
      as.integer(x$score_evaluations)
    }, integer(1))),
    probe_hessian_evaluations = sum(vapply(probe_records, function(x) {
      as.integer(x$hessian_evaluations)
    }, integer(1)))
  )
  shadow_resolution <- list(
    schema_version = 1L,
    policy_version = if (isTRUE(threshold_configuration_valid)) {
      resolution_thresholds$policy_version
    } else {
      NA_character_
    },
    mode = if (isTRUE(enable_prospective)) {
      "test_enablement"
    } else {
      "shadow_default_off"
    },
    thresholds = resolution_thresholds,
    strict_route = strict_route,
    strict_failure_reason = strict_failure,
    prospective_route = prospective_route,
    prospective_failure_reason = prospective_failure,
    would_nonveto = isTRUE(would_nonveto),
    would_route = prospective_route,
    numerical_resolution = numerical_resolution,
    global_veto_classes = as.character(global_veto_classes),
    candidate_resolutions = candidate_resolutions,
    evaluation_counts = evaluation_counts
  )
  route <- if (isTRUE(enable_prospective)) prospective_route else strict_route
  reported_failure <- if (isTRUE(enable_prospective)) {
    prospective_failure
  } else {
    strict_failure
  }
  list(
    route = route,
    failure_reason = reported_failure,
    best_index = best_idx,
    parameter_agreement = parameter_agreement,
    objective_agreement = objective_agreement,
    competing_uncertified = as.character(competing),
    dominated_uncertified = as.character(dominated),
    resolution = shadow_resolution,
    candidates = ln_joint_candidate_summary(candidates, candidate_resolutions)
  )
}

ln_joint_fit_engine <- function(
    gene_index, posv, ctx, map, extra_starts = NULL,
    objective_instance = NULL,
    resolution_thresholds = ln_joint_resolution_thresholds(),
    enable_prospective = FALSE) {
  free_dim <- map$free_dim
  lower <- c(if (free_dim) rep(-100, free_dim) else numeric(0), ctx$min)
  upper <- c(if (free_dim) rep(100, free_dim) else numeric(0), ctx$max)
  betae <- rep(0, ctx$nb)
  betae[ctx$intcol] <- log(posv$mct) - ctx$moffset
  beta_starts <- list(
    S_direct = betae,
    S_mid = betae,
    S_up = replace(betae, ctx$intcol, betae[ctx$intcol] + log(2)),
    S_down = replace(betae, ctx$intcol, betae[ctx$intcol] - log(2))
  )
  dispersion_starts <- list(
    S_direct = c(1, ctx$cell_init),
    S_mid = sqrt(ctx$min * ctx$max),
    S_up = c(1, ctx$cell_init),
    S_down = c(1, ctx$cell_init)
  )
  core_starts <- Map(function(beta, dispersion) {
    ln_joint_clamp(c(map$project(beta), dispersion), lower, upper)
  }, beta_starts, dispersion_starts)
  keys <- vapply(core_starts, function(x) paste(format(x, digits = 17), collapse = "|"), character(1))
  core_starts <- core_starts[!duplicated(keys)]

  objective <- function(theta) ln_joint_objective(theta, gene_index, posv, ctx, map)
  gradient <- function(theta) ln_joint_gradient(theta, gene_index, posv, ctx, map)
  controls <- list(eval.max = 500L, iter.max = 300L, rel.tol = 1e-12,
                   x.tol = 1e-10, xf.tol = 1e-12)

  run_start <- function(start, label, start_provenance, required, core,
                        optimizer = "nlminb") {
    fit <- tryCatch({
      if (identical(optimizer, "lbfgs")) {
        ans <- lbfgs(
          start,
          function(theta) list(objective = objective(theta), gradient = gradient(theta)),
          lower = lower,
          upper = upper,
          control = list(ftol_abs = ctx$eps)
        )
        list(
          par = ans$par,
          objective = objective(ans$par),
          convergence = ans$convergence,
          message = ans$message
        )
      } else if (identical(optimizer, "trust")) {
        log_start <- c(
          if (free_dim) start[seq_len(free_dim)] else numeric(0),
          log(start[free_dim + seq_len(2L)])
        )
        ans <- trust(
          objfun = function(z) {
            state <- ln_joint_derivative_state(z, gene_index, posv, ctx, map)
            natural <- c(
              if (free_dim) z[seq_len(free_dim)] else numeric(0),
              exp(z[free_dim + seq_len(2L)])
            )
            list(
              value = objective(natural),
              gradient = state$q_gradient_reduced,
              hessian = state$A_reduced
            )
          },
          parinit = log_start,
          rinit = 2,
          rmax = 100,
          iterlim = 100
        )
        par <- c(
          if (free_dim) ans$argument[seq_len(free_dim)] else numeric(0),
          exp(ans$argument[free_dim + seq_len(2L)])
        )
        list(
          par = par,
          objective = objective(par),
          convergence = if (isTRUE(ans$converged)) 0L else 1L,
          message = "trust"
        )
      } else {
        nlminb(
          start = start,
          objective = objective,
          gradient = gradient,
          lower = lower,
          upper = upper,
          control = controls
        )
      }
    }, error = function(e) {
      list(par = rep(NA_real_, length(start)), objective = NA_real_,
           convergence = NA_integer_, message = conditionMessage(e))
    })
    endpoint <- fit$par
    endpoint_objective <- fit$objective
    optimizer_endpoint <- endpoint
    optimizer_objective <- endpoint_objective
    polished <- FALSE
    if (length(endpoint) == length(start) && all(is.finite(endpoint)) &&
        all(endpoint[free_dim + seq_len(2L)] > 0)) {
      polish <- tryCatch(
        ln_joint_polish(endpoint, gene_index, posv, ctx, map, lower, upper),
        error = function(e) NULL
      )
      if (!is.null(polish) && is.finite(polish$objective)) {
        endpoint <- polish$theta_reduced
        endpoint_objective <- polish$objective
        polished <- !isTRUE(all.equal(
          endpoint,
          optimizer_endpoint,
          tolerance = 0,
          check.attributes = FALSE
        ))
      }
    }
    candidate <- ln_joint_candidate_state(
      endpoint,
      endpoint_objective,
      optimizer = fit,
      label = label,
      start_provenance = start_provenance,
      required = required,
      core = core,
      gene_index = gene_index,
      posv = posv,
      ctx = ctx,
      map = map,
      objective_instance = objective_instance
    )
    candidate$start <- as.numeric(start)
    candidate$optimizer_endpoint <- optimizer_endpoint
    candidate$optimizer_objective <- optimizer_objective
    candidate$polished <- polished
    candidate
  }

  candidates <- Map(function(start, label) {
    run_start(start, label, "direct", required = TRUE, core = TRUE)
  }, core_starts, names(core_starts))

  primary <- run_start(
    core_starts[[1L]],
    "primary_optimizer",
    "direct",
    required = FALSE,
    core = FALSE,
    optimizer = ctx$opt
  )
  candidates <- c(candidates, list(primary))

  finite_for_pql <- which(vapply(candidates, function(x) {
    isTRUE(x$finite_endpoint) && is.finite(x$objective)
  }, logical(1)))
  pql <- list(triggered = FALSE, status = "not_evaluated")
  if (length(finite_for_pql)) {
    prelim <- candidates[[finite_for_pql[
      which.min(vapply(candidates[finite_for_pql], `[[`, numeric(1), "objective"))
    ]]]
    pql <- ln_joint_pql_start(prelim, gene_index, posv, ctx, map)
    if (isTRUE(pql$triggered) && identical(pql$status, "generated")) {
      candidates <- c(candidates, list(run_start(
        ln_joint_clamp(pql$start, lower, upper),
        "pql_completed_start",
        "pql_start",
        required = TRUE,
        core = FALSE
      )))
    } else if (isTRUE(pql$triggered)) {
      candidates <- c(candidates, list(list(
        label = "pql_completed_start",
        start_provenance = "pql_start",
        required = TRUE,
        core = FALSE,
        n_beta = as.integer(ctx$nb),
        parameterization = "beta, log(sigma2), log(gamma_nb_size)",
        objective_provenance = ln_joint_objective_provenance(),
        objective_instance = objective_instance,
        bounds = list(
          lower = as.numeric(ctx$min),
          upper = as.numeric(ctx$max)
        ),
        probe_bounds = list(
          lower = c(if (free_dim) rep(-100, free_dim) else numeric(0),
                    log(ctx$min)),
          upper = c(if (free_dim) rep(100, free_dim) else numeric(0),
                    log(ctx$max))
        ),
        probe_map = ln_joint_probe_map_record(map, ctx$nb),
        finite_endpoint = FALSE,
        certified = FALSE,
        objective = NA_real_,
        optimizer = list(convergence = NA_integer_),
        failure = pql$status
      )))
    }
  }

  if (!is.null(extra_starts)) {
    for (ii in seq_along(extra_starts)) {
      start_full <- as.numeric(extra_starts[[ii]])
      if (length(start_full) != ctx$nb + 2L || any(!is.finite(start_full))) next
      start_reduced <- c(map$project(start_full[seq_len(ctx$nb)]),
                         start_full[ctx$nb + seq_len(2L)])
      candidates <- c(candidates, list(run_start(
        ln_joint_clamp(start_reduced, lower, upper),
        paste0("declared_extra_", ii),
        "other_declared_start",
        required = FALSE,
        core = FALSE
      )))
    }
  }

  basin_probe <- function(candidate, theta_reduced_log) {
    ln_joint_derivative_state(
      theta_reduced_log,
      gene_index,
      posv,
      ctx,
      map
    )
  }
  resolution <- ln_joint_resolve_candidates(
    candidates,
    resolution_thresholds = resolution_thresholds,
    basin_probe = basin_probe,
    enable_prospective = enable_prospective
  )
  best <- if (!is.na(resolution$best_index)) candidates[[resolution$best_index]] else NULL
  list(
    route = resolution$route,
    failure_reason = resolution$failure_reason,
    best = best,
    start_provenance = if (!is.null(best)) best$start_provenance else NA_character_,
    candidates = candidates,
    diagnostics = resolution,
    derivative_audit = NULL,
    pql = pql,
    map = map
  )
}

ln_joint_covariance <- function(state, nb) {
  A <- state$A_full
  U_subject <- state$U_subject_full
  B <- U_subject %*% t(U_subject)
  A_inv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(A_inv) || any(!is.finite(A_inv))) {
    return(list(success = FALSE, failure_reason = "singular_joint_curvature"))
  }
  model_full <- 0.5 * (A_inv + t(A_inv))
  robust_full <- A_inv %*% B %*% t(A_inv)
  robust_full <- 0.5 * (robust_full + t(robust_full))
  influence <- A_inv %*% U_subject
  direct_influence <- influence %*% t(influence)
  direct_error <- ln_joint_relative_frobenius(robust_full, direct_influence)
  profile_available <- FALSE
  profile_error <- NA_real_
  profile_beta <- NULL
  B_inv <- tryCatch(solve(B), error = function(e) NULL)
  if (!is.null(B_inv) && all(is.finite(B_inv))) {
    G <- t(A) %*% B_inv %*% A
    idx_beta <- seq_len(nb)
    idx_rho <- nb + seq_len(nrow(A) - nb)
    G_rr_inv <- tryCatch(solve(G[idx_rho, idx_rho, drop = FALSE]),
                         error = function(e) NULL)
    if (!is.null(G_rr_inv)) {
      G_beta <- G[idx_beta, idx_beta, drop = FALSE] -
        G[idx_beta, idx_rho, drop = FALSE] %*% G_rr_inv %*%
        G[idx_rho, idx_beta, drop = FALSE]
      profile_beta <- tryCatch(solve(G_beta), error = function(e) NULL)
      if (!is.null(profile_beta) && all(is.finite(profile_beta))) {
        profile_available <- TRUE
        profile_error <- ln_joint_relative_frobenius(
          robust_full[idx_beta, idx_beta, drop = FALSE],
          profile_beta
        )
      }
    }
  }
  thresholds <- ln_joint_thresholds()
  beta_model <- model_full[seq_len(nb), seq_len(nb), drop = FALSE]
  beta_robust <- robust_full[seq_len(nb), seq_len(nb), drop = FALSE]
  success <- is.finite(direct_error) &&
    direct_error <= thresholds$covariance_identity &&
    (!profile_available || (is.finite(profile_error) &&
      profile_error <= thresholds$covariance_identity)) &&
    all(is.finite(c(beta_model, beta_robust))) &&
    all(diag(beta_model) > 0) && all(diag(beta_robust) > 0)
  list(
    success = success,
    failure_reason = if (success) NA_character_ else "covariance_identity_failed",
    A = A,
    B = B,
    model_full = model_full,
    robust_full = robust_full,
    model_beta = beta_model,
    robust_beta = beta_robust,
    influence = influence,
    direct_influence_error = direct_error,
    profile_available = profile_available,
    profile_beta = profile_beta,
    profile_error = profile_error
  )
}

ln_joint_fixed_modes <- function(state, gene_index, posv, ctx) {
  unavailable <- function(status, method = "fixed_beta_pml") {
    list(
      values = rep(NA_real_, ctx$k),
      method = method,
      status = status,
      success = FALSE,
      max_score_residual = NA_real_,
      max_newton_displacement = NA_real_,
      donor_status = rep(status, ctx$k),
      alpha = NA_real_,
      lambda = NA_real_
    )
  }
  if (!isTRUE(ctx$output_re)) {
    return(unavailable("not_requested", method = "none"))
  }
  if (is.null(state)) {
    return(unavailable("fit_unavailable"))
  }
  result <- tryCatch(
    opt_pml_fixed_beta_modes(
      ctx$pred,
      ctx$offset,
      posv$Y,
      ctx$fid - 1L,
      ctx$cumsumy[gene_index, ],
      posv$posindy,
      ctx$nind,
      ctx$k,
      state$beta,
      c(state$sigma2, state$gamma),
      ln_joint_thresholds()$mode_score,
      200L
    ),
    error = function(e) NULL
  )
  if (is.null(result)) {
    return(unavailable("solver_error"))
  }
  max_score <- as.numeric(result$max_score_residual)
  max_newton <- as.numeric(result$max_newton_displacement)
  donor_status <- as.character(result$status)
  success <- length(result$logw) == ctx$k &&
    all(is.finite(result$logw)) &&
    all(donor_status == "ok") &&
    is.finite(max_score) && max_score <= ln_joint_thresholds()$mode_score &&
    is.finite(max_newton) && max_newton <= ln_joint_thresholds()$mode_newton
  list(
    values = if (success) as.numeric(result$logw) else rep(NA_real_, ctx$k),
    method = "fixed_beta_pml",
    status = if (success) "ok" else "gate_failed",
    success = success,
    max_score_residual = max_score,
    max_newton_displacement = max_newton,
    donor_status = donor_status,
    alpha = as.numeric(result$alpha),
    lambda = as.numeric(result$lambda),
    parameterization = result$parameterization
  )
}

fit_gene_joint_ln_unconstrained <- function(gene_index, posv, ctx,
                                             extra_starts = NULL,
                                             resolution_thresholds =
                                               ln_joint_resolution_thresholds(),
                                             enable_prospective = FALSE) {
  objective_instance <- ln_joint_objective_instance_registry(
    "register",
    gene_index,
    posv,
    ctx
  )
  map <- ln_joint_map(ctx$nb)
  engine <- ln_joint_fit_engine(
    gene_index,
    posv,
    ctx,
    map,
    extra_starts,
    objective_instance = objective_instance,
    resolution_thresholds = resolution_thresholds,
    enable_prospective = enable_prospective
  )
  best <- engine$best
  state <- if (!is.null(best)) best$state else NULL
  if (!is.null(state)) {
    state$objective_instance <- objective_instance
  }
  route <- engine$route
  failure <- engine$failure_reason
  covariance <- NULL
  if (identical(route, "R0")) {
    covariance <- ln_joint_covariance(state, ctx$nb)
    if (!isTRUE(covariance$success)) {
      route <- "R2"
      failure <- covariance$failure_reason
    }
  }
  inference_available <- identical(route, "R0")
  mode <- ln_joint_fixed_modes(if (inference_available) state else NULL,
                               gene_index, posv, ctx)
  beta <- if (!is.null(state)) state$beta else rep(NA_real_, ctx$nb)
  sigma2 <- if (!is.null(state)) state$sigma2 else NA_real_
  gamma <- if (!is.null(state)) state$gamma else NA_real_
  model_cov <- if (inference_available) covariance$model_beta else
    matrix(NA_real_, ctx$nb, ctx$nb)
  robust_cov <- if (inference_available) covariance$robust_beta else
    matrix(NA_real_, ctx$nb, ctx$nb)
  loglik <- if (!is.null(best) && is.finite(best$objective)) -best$objective else NA_real_
  conv_code <- if (inference_available) 1L else -70L
  restemp <- c(
    beta,
    sigma2,
    if (is.finite(gamma) && gamma > 0) 1 / gamma else NA_real_,
    diag(model_cov),
    diag(robust_cov),
    conv_code,
    1L,
    loglik
  )
  if (ctx$output_re) {
    restemp <- c(restemp, mode$values)
  }
  derivative_object <- if (!is.null(state)) {
    list(
      gradient = state$U_full,
      hessian = state$A_full,
      observed_info = state$A_full,
      per_subject_gradients = state$U_subject_full
    )
  } else NULL
  godambe <- if (!is.null(derivative_object))
    extract_godambe_components(derivative_object, ctx$nb) else NULL
  theta <- list(beta = beta, subVar = sigma2, cellVar = gamma)
  result <- list(
    vec = restemp,
    godambe = godambe,
    loglik = loglik,
    loglik_ln = loglik,
    theta = theta,
    theta_ln = theta,
    theta_ln_rescued = theta,
    theta_unrescued = theta,
    score = if (!is.null(state)) state$U_full else rep(NA_real_, ctx$nb + 2L),
    naive_cov = model_cov,
    robust_cov = robust_cov,
    route = route,
    start_provenance = engine$start_provenance,
    inference_available = inference_available,
    failure_reason = failure,
    mode_method = mode$method,
    mode_status = mode$status,
    mode = mode,
    ln_joint = list(
      parameterization = "beta, log(sigma2), log(gamma_nb_size)",
      score_sign = "maximizing",
      objective_provenance = if (!is.null(state)) {
        state$objective_provenance
      } else {
        ln_joint_objective_provenance()
      },
      objective_instance = objective_instance,
      state = state,
      covariance = covariance,
      candidates = engine$diagnostics$candidates,
      candidate_records = ln_joint_candidate_records(
        engine$candidates,
        engine$diagnostics$resolution$candidate_resolutions
      ),
      resolution = engine$diagnostics$resolution,
      parameter_agreement = engine$diagnostics$parameter_agreement,
      objective_agreement = engine$diagnostics$objective_agreement,
      derivative_audit = engine$derivative_audit,
      pql_start = engine$pql
    ),
    ln_refinement = list(
      refined = TRUE,
      grad_norm_start = NA_real_,
      grad_norm_final = if (!is.null(best)) best$root_residual else NA_real_,
      extra_attempts = if (is.null(extra_starts)) 0L else length(extra_starts),
      extra_success = 0L
    ),
    debug = list(
      method = "joint_ln",
      gene_index = gene_index,
      fit_code = 1L,
      conv_LN = conv_code,
      conv_code = conv_code,
      route = route,
      failure_reason = failure
    )
  )
  ln_joint_objective_instance_registry(
    "store_unrestricted",
    gene_index,
    posv,
    ctx,
    fit = result
  )
  result
}

fit_gene_joint_ln_constrained <- function(
    gene_index, posv, ctx, L, b,
    resolution_thresholds = ln_joint_resolution_thresholds(),
    enable_prospective = FALSE) {
  objective_instance <- ln_joint_objective_instance_registry(
    "lookup",
    gene_index,
    posv,
    ctx
  )
  map <- ln_joint_map(ctx$nb, L, b)
  engine <- ln_joint_fit_engine(
    gene_index,
    posv,
    ctx,
    map,
    objective_instance = objective_instance,
    resolution_thresholds = resolution_thresholds,
    enable_prospective = enable_prospective
  )
  best <- engine$best
  state <- if (!is.null(best)) best$state else NULL
  if (!is.null(state)) {
    state$objective_instance <- objective_instance
  }
  route <- engine$route
  failure <- engine$failure_reason
  beta <- if (!is.null(state)) state$beta else rep(NA_real_, ctx$nb)
  sigma2 <- if (!is.null(state)) state$sigma2 else NA_real_
  gamma <- if (!is.null(state)) state$gamma else NA_real_
  loglik <- if (!is.null(best) && is.finite(best$objective)) -best$objective else NA_real_
  inference_available <- identical(route, "R0")
  derivative_object <- if (!is.null(state)) {
    list(
      gradient = state$U_full,
      hessian = state$A_full,
      observed_info = state$A_full,
      per_subject_gradients = state$U_subject_full
    )
  } else NULL
  godambe <- if (!is.null(derivative_object))
    extract_godambe_components(derivative_object, ctx$nb) else NULL
  theta <- list(beta = beta, subVar = sigma2, cellVar = gamma)
  conv_code <- if (inference_available) 1L else -70L
  restemp <- c(
    beta,
    sigma2,
    if (is.finite(gamma) && gamma > 0) 1 / gamma else NA_real_,
    rep(NA_real_, 2L * ctx$nb),
    conv_code,
    1L,
    loglik
  )
  if (ctx$output_re) restemp <- c(restemp, rep(NA_real_, ctx$k))
  reduced <- if (!is.null(state)) list(
    parameterization = "zeta, log(sigma^2), log(gamma_nb_size)",
    map = map$derivative_map,
    theta = best$theta_reduced_log,
    theta_red = best$theta_reduced_log,
    gamma = if (map$free_dim) best$theta_reduced_log[seq_len(map$free_dim)] else numeric(0),
    score = state$U_reduced,
    g_red = state$U_reduced,
    observed_info = state$A_reduced,
    H_red = state$A_reduced,
    per_subject_gradients = state$U_subject_reduced
  ) else NULL
  result <- list(
    vec = restemp,
    godambe = godambe,
    beta_constrained = beta,
    loglik = loglik,
    loglik_ln = loglik,
    loglik_ln_pre = loglik,
    loglik_ln_rescued = loglik,
    theta = theta,
    theta_ln = theta,
    theta_ln_rescued = theta,
    theta_unrescued = theta,
    constraint = list(
      L = L,
      b = b,
      beta_star = map$beta_star,
      null_basis = map$K,
      free_dim = map$free_dim
    ),
    reduced = reduced,
    diagnostics = if (!is.null(best)) list(
      reduced_gradient_max_abs = max(abs(state$U_reduced)),
      H_red_positive_definite = isTRUE(best$curvature$positive_definite),
      H_red_rcond = if (is.finite(best$curvature$scaled_condition_number))
        1 / best$curvature$scaled_condition_number else 0,
      max_abs_newton_step_scaled = best$newton_displacement,
      constraint_residual_scaled = best$constraint_residual
    ) else list(),
    score = if (!is.null(state)) state$U_full else rep(NA_real_, ctx$nb + 2L),
    robust_cov = matrix(NA_real_, ctx$nb, ctx$nb),
    route = route,
    start_provenance = engine$start_provenance,
    inference_available = inference_available,
    failure_reason = failure,
    lrt_gate = list(
      converged = inference_available,
      interior = if (!is.null(best)) isTRUE(best$interior) else FALSE,
      nuisance_interior = if (!is.null(best)) isTRUE(best$interior) else FALSE,
      constraint_satisfied = if (!is.null(best))
        best$constraint_residual <= ln_joint_thresholds()$constraint else FALSE,
      calibration_available = FALSE,
      calibration_status = "joint_ln_calibration_pending"
    ),
    ln_joint = list(
      parameterization = "beta, log(sigma2), log(gamma_nb_size)",
      score_sign = "maximizing",
      objective_provenance = if (!is.null(state)) {
        state$objective_provenance
      } else {
        ln_joint_objective_provenance()
      },
      objective_instance = objective_instance,
      state = state,
      candidates = engine$diagnostics$candidates,
      candidate_records = ln_joint_candidate_records(
        engine$candidates,
        engine$diagnostics$resolution$candidate_resolutions
      ),
      resolution = engine$diagnostics$resolution,
      parameter_agreement = engine$diagnostics$parameter_agreement,
      objective_agreement = engine$diagnostics$objective_agreement,
      derivative_audit = engine$derivative_audit
    ),
    debug = list(
      method = "joint_ln_constrained",
      gene_index = gene_index,
      fit_code = 1L,
      conv_LN = conv_code,
      conv_code = conv_code,
      free_dim = map$free_dim,
      route = route,
      failure_reason = failure
    )
  )
  result$ln_joint$score_diagnostic <-
    ln_joint_robust_score_diagnostic(result)
  unrestricted_fit <- ln_joint_objective_instance_registry(
    "get_unrestricted",
    gene_index,
    posv,
    ctx
  )
  result$ln_joint$pace_ratio_diagnostic <-
    ln_joint_pace_ratio_diagnostic(
      if (is.list(unrestricted_fit)) unrestricted_fit else list(),
      result,
      objective_context = list(
        gene_index = as.integer(gene_index),
        posv = posv,
        ctx = ctx
      )
    )
  result
}
