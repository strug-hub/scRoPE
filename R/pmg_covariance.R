# Equivalent sandwich algebra without requiring invertible donor-score meat.
# Preserve the historical answer when it agrees numerically with the direct
# solve; a singular meat is not a singular sensitivity or a failed contrast.
pmg_sandwich_covariance <- function(out_at_op, nb, compute_full = TRUE) {
  H <- out_at_op$observed_info
  U <- out_at_op$per_subject_gradients
  if (!is.matrix(H) || nrow(H) != ncol(H) || !is.matrix(U) ||
      nrow(U) != nrow(H) || any(!is.finite(c(H, U))))
    stop("Nonfinite or incompatible PGMM sensitivity/donor scores.")
  if (length(nb) != 1L || !is.finite(nb) || nb < 1L ||
      nb != as.integer(nb) || nb > nrow(H))
    stop("nb must index a nonempty coefficient block of the PGMM sensitivity.")
  influence <- solve(H, U)
  V <- tcrossprod(influence)
  if (any(!is.finite(V))) stop("Nonfinite PGMM influence covariance.")
  legacy <- tryCatch(suppressWarnings(compute_sandwich_variance2(
    out_at_op, nb, compute_full)), error = function(e) NULL)
  beta <- V[seq_len(nb), seq_len(nb), drop = FALSE]
  # Agreement is measured in variance/SD-product units, with no absolute
  # floor: small but positive variances require the same relative accuracy.
  agrees_with <- function(old, direct) {
    if (!is.matrix(old) || !identical(dim(old), dim(direct)) ||
        any(!is.finite(old))) return(FALSE)
    agreement_scale <- outer(sqrt(diag(direct)), sqrt(diag(direct)))
    all((diag(old) > 0) == (diag(direct) > 0)) &&
      all(abs(old - direct) <= 1e-10 * agreement_scale)
  }
  agrees <- !is.null(legacy) && agrees_with(legacy$Var_beta_adjusted, beta) &&
    (!compute_full || agrees_with(legacy$Var_full, V))
  if (agrees) return(legacy)
  ans <- list(Var_beta_adjusted = beta)
  if (compute_full) ans$Var_full <- V
  ans
}
