# Native qualification uses the estimator's existing optimizer/curvature and
# coefficient inference rules. It never evaluates a recovery certificate.
working_native_inference <- function(model, native, ctx) {
  if (!is.null(native$working_inference)) return(native$working_inference)
  if (model == "PGMM") {
    accepted <- pmg_optimizer_converged(native)
    return(list(
      model = pmg_coefficient_wald(native$theta$beta, diag(native$naive_cov), accepted)$success,
      robust = pmg_coefficient_wald(native$theta$beta, diag(native$robust_cov), accepted)$success))
  }
  p <- ncol(ctx$pred)
  robust <- rep(FALSE, p)
  if (isTRUE(native$point_success) && isTRUE(native$covariance_success)) {
    beta <- as.numeric(ctx$M_X %*% native$beta)
    V <- ctx$M_X %*% native$covariance$V_beta %*% t(ctx$M_X)
    robust <- vapply(seq_len(p), function(j) {
      R <- matrix(0, 1L, p); R[1L,j] <- 1
      isTRUE(hl_fo_robust_wald(beta,V,R,0,n_subjects=length(unique(ctx$id)))$success)
    }, logical(1))
  }
  model_variance <- if(isTRUE(native$point_success) && isTRUE(native$model_covariance_success))
    diag(ctx$M_X %*% native$model_covariance$V_beta %*% t(ctx$M_X)) else rep(NA_real_,p)
  beta <- if(isTRUE(native$point_success))as.numeric(ctx$M_X %*% native$beta) else rep(NA_real_,p)
  list(model=pmg_coefficient_wald(beta,model_variance,isTRUE(native$point_success))$success,
    robust=robust)
}

working_native_qualification <- function(model, native, ctx, unresolved) {
  inference <- working_native_inference(model,native,ctx)
  point <- if (model == "HL") isTRUE(native$point_success) && !unresolved else
    pmg_optimizer_converged(native) && !unresolved &&
      all(is.finite(native$theta$beta)) && !isTRUE(native$diagnostics$boundary_any) &&
      isTRUE(native$diagnostics$information_positive_definite) &&
      isTRUE(is.finite(native$diagnostics$information_rcond) &&
        native$diagnostics$information_rcond >= sqrt(.Machine$double.eps))
  p <- ncol(ctx$pred)
  complete_family <- function(x) length(x)==p && all(x %in% TRUE)
  robust <- point && complete_family(inference$robust)
  model_complete <- point && complete_family(inference$model)
  list(point=point, robust_complete=robust, model_complete=model_complete,
    complete=robust && isTRUE(model_complete), model_supported=TRUE)
}

# Complete recovery supplies both coefficient families for every estimable
# coefficient. Never splice covariance or points into a native partial fit.
working_recovery_complete <- function(ans) {
  if (!isTRUE(ans$point_success) || !isTRUE(ans$inference_available)) return(FALSE)
  keep <- ans$estimable %in% TRUE
  if (!length(keep) || !any(keep) || length(ans$beta_report)!=length(keep) ||
      !is.matrix(ans$V_report) || !identical(dim(ans$V_report),rep(length(keep),2L))) return(FALSE)
  beta <- ans$beta_report[keep]; variance <- diag(ans$V_report)[keep]
  robust <- all(is.finite(beta) & is.finite(variance) & variance>0 & is.finite((beta/sqrt(variance))^2))
  if(!robust || !isTRUE(ans$model_inference_available) || !is.matrix(ans$model_V_report) ||
    !identical(dim(ans$model_V_report),rep(length(keep),2L)))return(FALSE)
  variance <- diag(ans$model_V_report)[keep]
  all(is.finite(variance) & variance>0 & is.finite((beta/sqrt(variance))^2))
}

working_native_selected <- function(ans) {
  isTRUE(ans$native_selected) || isTRUE(ans$native_regular)
}
