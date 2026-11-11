# Model covariance calculations consume the selected, already-qualified system.
# These factorizations are not optimizer/root certificates and never refit.
scrope_model_inverse <- function(information) {
  fail <- function(reason) list(success=FALSE, covariance=NULL, reason=reason)
  if (!is.matrix(information) || nrow(information)!=ncol(information) ||
      !nrow(information) || any(!is.finite(information)))
    return(fail("model_information_unavailable"))
  # Information is a mathematical Hessian. Remove roundoff antisymmetry, as
  # in the native JLN inverse-curvature convention; never use augmented HL A.
  H <- 0.5 * (information + t(information))
  V <- tryCatch(chol2inv(chol(H)),error=function(e)NULL)
  if (is.null(V) || any(!is.finite(V)) || any(diag(V)<=0))
    return(fail("model_information_not_positive_definite_or_solve_failed"))
  list(success=TRUE,covariance=0.5*(V+t(V)),reason=NA_character_)
}

hl_fo_model_covariance <- function(state) {
  # B = sum p_{i,gamma gamma}; modes are profiled out and fitted dispersions
  # are plugged in. This is not the inverse of the augmented sensitivity.
  ans <- scrope_model_inverse(-state$augmented$B)
  if (isTRUE(ans$success)) {
    ans$V_gamma <- ans$covariance
    ans$V_beta <- state$K %*% ans$V_gamma %*% t(state$K)
  }
  ans$coordinate_tag <- "canonical_conditional_beta"
  ans$convention <- "profile_h_information_plugin_dispersion"
  ans
}

scrope_model_convention <- function(model) switch(model,
  PGMM="observed_free_marginal_information",
  LN="joint_LN_free_inverse_curvature",
  HL="profile_h_information_plugin_dispersion")
