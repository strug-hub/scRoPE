# Public assembly only: fitting, native statistics and numerical gates remain
# in ln_joint.R. The initial opt-in scope matches the retained scalar-null
# study, including its accepted unrestricted-primary requirement for Score.
ln_joint_public_optional_tests <- function(unrestricted_fit, constrained_fit,
                                          L, b, intercept_indices,
                                          coef_internal_to_public) {
  L <- as.matrix(L)
  b <- as.numeric(b)
  df <- as.integer(nrow(L))
  score <- list(
    method = "joint_ln_robust_score", df = df,
    stat = NA_real_, pval = NA_real_, ridge = NA_real_, success = FALSE,
    status = "unavailable", inference = "experimental_asymptotic",
    failure_reason = NA_character_, failure_stage = NA_character_
  )
  lrt <- list(
    method = "joint_ln_pace_adjusted_ratio", df = df,
    stat = NA_real_, pval = NA_real_, wP_adjusted = NA_real_,
    success = FALSE, adjusted_available = FALSE,
    status = "unavailable", inference = "experimental_asymptotic",
    adjusted_failure_reason = NA_character_, failure_reason = NA_character_,
    failure_stage = NA_character_, fallback = FALSE
  )
  reason_or <- function(value, fallback) {
    if (is.character(value) && length(value) == 1L && !is.na(value) &&
        nzchar(value)) value else fallback
  }
  fail_both <- function(reason, stage) {
    score$failure_reason <- reason
    score$failure_stage <- stage
    lrt$failure_reason <- lrt$adjusted_failure_reason <- reason
    lrt$failure_stage <- stage
    list(score = score, lrt = lrt)
  }
  selected <- which(L != 0)
  scope_ok <- nrow(L) == 1L && ncol(L) > 0L &&
    all(is.finite(L)) && length(b) == 1L && is.finite(b) && b == 0 &&
    length(selected) == 1L && !selected %in% intercept_indices
  if (!isTRUE(scope_ok)) {
    return(fail_both("joint_ln_optional_scope_unsupported", "scope"))
  }
  primary_ok <- identical(unrestricted_fit$route, "R0") &&
    isTRUE(unrestricted_fit$inference_available) &&
    isTRUE(unrestricted_fit$ln_joint$resolution$selected_root_certified)
  if (!primary_ok) {
    return(fail_both(reason_or(unrestricted_fit$failure_reason,
                               "joint_ln_primary_unavailable"), "primary"))
  }
  gate <- constrained_fit$lrt_gate
  constrained_ok <- identical(constrained_fit$route, "R0") &&
    isTRUE(constrained_fit$inference_available) && isTRUE(gate$converged) &&
    isTRUE(gate$interior) && isTRUE(gate$nuisance_interior) &&
    isTRUE(gate$constraint_satisfied)
  if (!constrained_ok) {
    return(fail_both(reason_or(constrained_fit$failure_reason,
                               "joint_ln_constrained_unavailable"),
                     "constrained"))
  }
  expected <- ln_joint_objective_provenance()
  primary_joint <- unrestricted_fit$ln_joint
  null_joint <- constrained_fit$ln_joint
  if (!identical(primary_joint$objective_provenance, expected) ||
      !identical(null_joint$objective_provenance, expected)) {
    return(fail_both("endpoint_criterion_mismatch", "objective"))
  }
  if (is.null(primary_joint$objective_instance) ||
      !identical(primary_joint$objective_instance,
                 null_joint$objective_instance)) {
    return(fail_both("endpoint_objective_instance_mismatch", "objective"))
  }
  map <- as.matrix(coef_internal_to_public)
  map_ok <- identical(dim(map), c(ncol(L), ncol(L))) &&
    all(is.finite(map))
  if (!map_ok || !identical(as.numeric(constrained_fit$constraint$L),
                            as.numeric(L %*% map)) ||
      !identical(as.numeric(constrained_fit$constraint$b), b)) {
    return(fail_both("joint_ln_public_constraint_mismatch", "constraint"))
  }

  native_score <- null_joint$score_diagnostic
  native_pace <- null_joint$pace_ratio_diagnostic
  matrix_positive <- function(matrix, diagnostics) {
    if (is.null(matrix)) return(NA)
    value <- diagnostics$positive_definite
    if (is.logical(value) && length(value) == 1L) value else NA
  }
  # Detail flags describe the native effective matrices, not the constrained
  # nuisance Hessian. An uncomputed native quantity remains unknown.
  lrt$gate <- list(
    profiled_curvature_positive = matrix_positive(
      native_pace$H_eff, native_pace$diagnostics$H_eff
    ),
    profiled_variability_positive = matrix_positive(
      native_pace$J_eff, native_pace$diagnostics$J_eff
    ),
    scale_positive = if (length(native_pace$multiplier) == 1L &&
        is.finite(native_pace$multiplier)) {
      native_pace$multiplier > 0
    } else NA
  )
  usable <- function(x) {
    isTRUE(x$diagnostic_available) &&
      isTRUE(x$no_numerical_repair) &&
      is.numeric(x$statistic) && length(x$statistic) == 1L &&
      is.finite(x$statistic) && x$statistic >= 0 &&
      identical(x$df, 1L) && identical(x$theoretical_df, 1L) &&
      identical(x$theoretical_reference, paste(
        "chi-square first-order reference only under smooth interior",
        "independent-donor fixed-dimension regularity"
      ))
  }
  if (usable(native_score)) {
    score$stat <- native_score$statistic
    score$pval <- stats::pchisq(score$stat, df = 1L, lower.tail = FALSE)
    score$success <- TRUE
    score$status <- "experimental_asymptotic"
  } else {
    score$failure_reason <- reason_or(native_score$failure_reason,
                                      "joint_ln_score_diagnostic_unavailable")
    score$failure_stage <- "diagnostic"
  }
  if (usable(native_pace) && isTRUE(score$success) &&
      identical(native_pace$score_diagnostic, native_score)) {
    lrt$wP <- native_pace$wP
    lrt$wU <- native_pace$wU
    lrt$qP <- native_pace$qP
    lrt$wP_minus_qP <- native_pace$diagnostics$wP_minus_qP
    lrt$wP_over_qP <- native_pace$diagnostics$wP_over_qP
    lrt$scale <- native_pace$multiplier
    lrt$stat <- lrt$wP_adjusted <- native_pace$statistic
    lrt$pval <- stats::pchisq(lrt$stat, df = 1L, lower.tail = FALSE)
    lrt$success <- lrt$adjusted_available <- TRUE
    lrt$status <- "experimental_asymptotic"
  } else {
    lrt$failure_reason <- lrt$adjusted_failure_reason <-
      reason_or(native_pace$failure_reason,
                if (isTRUE(native_pace$diagnostic_available)) {
                  "joint_ln_pace_score_identity_mismatch"
                } else {
                  "joint_ln_pace_diagnostic_unavailable"
                })
    lrt$failure_stage <- "diagnostic"
  }
  list(score = score, lrt = lrt)
}
