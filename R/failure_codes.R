# Canonical labels for adjusted-test failure reasons
failure_codes <- list(
  operational_fallback = list(
    legacy_fit_error = "legacy_fit_error",
    legacy_nonconvergence = "legacy_nonconvergence",
    legacy_gradient_convergence_not_accepted =
      "legacy_gradient_convergence_not_accepted",
    nonfinite_coefficient = "nonfinite_coefficient",
    invalid_dispersion = "invalid_dispersion",
    invalid_coefficient_transform = "invalid_coefficient_transform"
  ),
  wald = list(
    missing_godambe = "missing_godambe",
    missing_blocks = "missing_blocks",
    singular_S = "singular_S",
    singular_G_ll = "singular_G_ll",
    singular_G_beta = "singular_G_beta",
    singular_cov = "singular_cov",
    undefined_stat = "undefined_stat",
    singular_godambe = "singular_godambe",
    required_start_failed = "required_start_failed",
    no_certified_candidate = "no_certified_candidate",
    insufficient_certified_core_starts = "insufficient_certified_core_starts",
    start_parameter_review_zone = "start_parameter_review_zone",
    start_parameter_hard_stop = "start_parameter_hard_stop",
    start_objective_review_zone = "start_objective_review_zone",
    start_objective_hard_stop = "start_objective_hard_stop",
    competing_uncertified_candidate = "competing_uncertified_candidate",
    singular_joint_curvature = "singular_joint_curvature",
    covariance_identity_failed = "covariance_identity_failed",
    donor_additivity_target_failed = "donor_additivity_target_failed",
    donor_additivity_hard_stop = "donor_additivity_hard_stop",
    root_residual_review_zone = "root_residual_review_zone",
    root_residual_hard_stop = "root_residual_hard_stop",
    newton_displacement_review_zone = "newton_displacement_review_zone",
    newton_displacement_hard_stop = "newton_displacement_hard_stop",
    nonpositive_joint_curvature = "nonpositive_joint_curvature",
    curvature_review_zone = "curvature_review_zone",
    curvature_condition_review_zone = "curvature_condition_review_zone",
    curvature_condition_unusable = "curvature_condition_unusable",
    dispersion_boundary = "dispersion_boundary",
    constraint_review_zone = "constraint_review_zone",
    constraint_hard_stop = "constraint_hard_stop",
    derivative_evaluation_failed = "derivative_evaluation_failed",
    derivative_audit_failed = "derivative_audit_failed",
    objective_value_review_zone = "objective_value_review_zone",
    objective_value_hard_stop = "objective_value_hard_stop",
    objective_additivity_review_zone = "objective_additivity_review_zone",
    objective_additivity_hard_stop = "objective_additivity_hard_stop",
    objective_gradient_review_zone = "objective_gradient_review_zone",
    objective_gradient_hard_stop = "objective_gradient_hard_stop",
    sensitivity_score_review_zone = "sensitivity_score_review_zone",
    sensitivity_score_hard_stop = "sensitivity_score_hard_stop",
    fd_gradient_stability_review_zone = "fd_gradient_stability_review_zone",
    fd_gradient_stability_hard_stop = "fd_gradient_stability_hard_stop",
    fd_sensitivity_stability_review_zone = "fd_sensitivity_stability_review_zone",
    fd_sensitivity_stability_hard_stop = "fd_sensitivity_stability_hard_stop",
    natural_log_gradient_review_zone = "natural_log_gradient_review_zone",
    natural_log_gradient_hard_stop = "natural_log_gradient_hard_stop",
    natural_log_sensitivity_review_zone = "natural_log_sensitivity_review_zone",
    natural_log_sensitivity_hard_stop = "natural_log_sensitivity_hard_stop",
    natural_gradient_implementation_review_zone = "natural_gradient_implementation_review_zone",
    natural_gradient_implementation_hard_stop = "natural_gradient_implementation_hard_stop",
    natural_sensitivity_score_review_zone = "natural_sensitivity_score_review_zone",
    natural_sensitivity_score_hard_stop = "natural_sensitivity_score_hard_stop"
  ),
  lrt = list(
    negative_raw_lr = "negative_raw_lr",
    singular_Gpsi = "singular_Gpsi",
    singular_Hpsi = "singular_Hpsi",
    singular_Hlambda = "singular_Hlambda",
    nonconverged_constrained = "nonconverged_constrained",
    constraint_not_satisfied = "constraint_not_satisfied",
    boundary_constrained = "boundary_constrained",
    nonpositive_curvature = "nonpositive_curvature",
    nonpositive_scale = "nonpositive_scale",
    undefined_scale = "undefined_scale",
    wald_fallback_failed = "wald_fallback_failed",
    missing_profile_lr = "missing_profile_lr",
    undefined_adjustment = "undefined_adjustment",
    joint_ln_calibration_pending = "joint_ln_calibration_pending"
  ),
  score = list(
    singular_Spsi = "singular_Spsi",
    undefined_stat = "undefined_stat",
    joint_ln_calibration_pending = "joint_ln_calibration_pending"
  )
)

# Helper to normalise failure reasons to the known vocabulary
normalise_failure <- function(reason, type = c("wald", "lrt", "score")) {
  type <- match.arg(type)
  if (is.null(reason)) {
    return(NA_character_)
  }

  reason_vec <- as.character(reason)
  original_na <- is.na(reason)
  reason_vec[original_na] <- NA_character_
  out <- rep(NA_character_, length(reason_vec))
  missing_idx <- is.na(reason_vec) | (!original_na & !nzchar(reason_vec))
  if (all(missing_idx)) {
    return(if (length(out) == 1) out[[1]] else out)
  }

  known_idx <- reason_vec %in% failure_codes[[type]]
  out[known_idx] <- reason_vec[known_idx]
  other_idx <- !(known_idx | missing_idx)
  out[other_idx] <- paste0("other_", reason_vec[other_idx])

  if (length(out) == 1) {
    return(out[[1]])
  }
  out
}

normalise_operational_fallback_failure <- function(reason) {
  if (is.null(reason) || !length(reason) || is.na(reason[[1L]]) ||
      !nzchar(as.character(reason[[1L]]))) {
    return(NA_character_)
  }
  reason <- as.character(reason[[1L]])
  known <- unname(unlist(failure_codes$operational_fallback, use.names = FALSE))
  if (!reason %in% known) {
    stop("Unknown operational fallback failure reason: ", reason)
  }
  reason
}
