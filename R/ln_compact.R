# Per-gene retention boundary. Run optional tests and recovery while native
# candidate/donor evidence is still available, then return only report inputs.
ln_compact_fields <- function(x, fields) {
  if (is.null(x)) return(NULL)
  x[intersect(fields, names(x))]
}

ln_compact_optional <- function(entry, keep_details = FALSE) {
  out <- ln_compact_fields(entry, c(
    "df", "stat", "pval", "success", "failure_reason", "scale", "wP", "wU",
    "qP", "method", "status", "failure_stage", "ridge", "fallback",
    "score_stat", "score_pval", "score_ridge", "score_success", "score_failure",
    "score_method", "score_status", "score_failure_stage"
  ))
  if (isTRUE(keep_details)) {
    out$raw <- ln_compact_fields(entry$raw, c(
      "wP_minus_qP", "wP_over_qP", "adjusted_available",
      "adjusted_failure_reason", "gate"
    ))
    # The detail table consumes only these scalars, never the null fit or
    # derivative matrices attached to optional-test diagnostic objects.
    out$raw$gate <- ln_compact_fields(entry$raw$gate, c(
      "converged", "convergence_basis", "unrestricted_route",
      "unrestricted_failure_reason", "constrained_route",
      "constrained_failure_reason", "legacy_converged", "legacy_conv_code",
      "information_scaled_eligible", "information_scaled_convergence_ok",
      "information_scaled_accepted", "H_red_positive_definite",
      "diagnostics_finite", "H_red_rcond", "reduced_gradient_max_abs",
      "max_abs_newton_step_scaled", "gradient_quadratic_lr_bound",
      "newton_step_scaled_ok", "gradient_quadratic_lr_bound_ok",
      "reduced_parameterization", "constraint_satisfied",
      "constraint_residual_max_abs", "nuisance_hessian_invertible", "interior",
      "nuisance_interior", "profiled_curvature_positive",
      "profiled_variability_positive", "scale_positive"
    ))
  }
  out
}

ln_compact_recovery <- function(answer) {
  f <- answer$result
  result <- ln_compact_fields(f, c("success", "reason"))
  if (isTRUE(f$success) && !identical(answer$route, "ordinary_jln")) {
    result$beta <- f$beta
    result$Vbeta <- f$Vbeta
    if (identical(answer$route, "identifiable_nuisance_limit")) {
      result$face <- ln_compact_fields(f$face, "estimable")
      result$reduced_full <- f$reduced_full
      result$q <- f$q
      result$reduced <- list(result = ln_compact_fields(f$reduced$result, "fixed"))
    } else {
      result$state <- ln_compact_fields(f$state, c("full", "q"))
      result$fixed <- f$fixed
    }
  }
  list(route = answer$route, result = result)
}

ln_finalize_gene <- function(fit, keep_diagnostics, keep_lrt_details = FALSE,
                             recovery_input = NULL, gene_index = NULL,
                             recovery_enabled = FALSE) {
  recovery <- NULL
  if (!is.null(recovery_input)) {
    recovery <- ln_recovery_evaluate(
      list(success = isTRUE(fit$inference_available), joint = fit$ln_joint),
      recovery_input, gene_index, recovery_enabled
    )
    if (!isTRUE(keep_diagnostics)) recovery <- ln_compact_recovery(recovery)
  }
  if (!isTRUE(keep_diagnostics)) {
    out <- ln_compact_fields(fit, c(
      "vec", "route", "start_provenance", "inference_available",
      "failure_reason", "mode_method", "mode_status"
    ))
    # Preserve the exact native status inputs. Candidate rows and donor-level
    # mode placeholders are diagnostics, not inputs to the public status.
    out$ln_joint <- list(resolution = ln_compact_fields(fit$ln_joint$resolution, c(
      "strict_inference_route", "strict_route", "enabled_inference_route",
      "prospective_route", "strict_inference_failure_reason",
      "strict_failure_reason", "enabled_inference_failure_reason",
      "prospective_failure_reason", "selected_root_certified", "would_nonveto",
      "global_veto_classes", "policy_version", "mode"
    )))
    out$mode <- ln_compact_fields(fit$mode, c(
      "method", "status", "success", "max_score_residual",
      "max_newton_displacement", "alpha", "lambda", "parameterization"
    ))
    out$debug <- ln_compact_fields(fit$debug, "rescue")
    if (!is.null(fit$lrt)) {
      out$lrt <- lapply(fit$lrt, ln_compact_optional,
        keep_details = keep_lrt_details)
    }
    fit <- out
  }
  if (!is.null(recovery)) fit$precomputed_recovery <- recovery
  fit
}
