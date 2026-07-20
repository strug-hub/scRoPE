# Pure inference algebra for certified first-order nested HL states.
#
# This file deliberately contains no fitting, root search, ridge repair, or
# objective-ratio calibration.  All beta quantities are in the tagged
# canonical internal analysis coordinates.

.hl_fo_inference_rank <- function(x) {
  x <- as.matrix(x)
  if (!all(is.finite(x))) {
    stop("Numerical rank requires a finite matrix.")
  }
  info <- hl_fo_svd_info(x, compute_pinv = FALSE)
  list(
    rank = info$rank,
    tolerance = info$threshold,
    singular_values = info$singular_values
  )
}

.hl_fo_reduce_contrast <- function(R, rhs = NULL) {
  if (is.null(dim(R)) && is.numeric(R)) {
    R <- matrix(R, nrow = 1L)
  } else {
    R <- as.matrix(R)
  }
  if (!is.numeric(R) || !all(is.finite(R)) || nrow(R) == 0L || ncol(R) == 0L) {
    stop("R must be a nonempty finite numeric matrix.")
  }
  if (is.null(rhs)) {
    rhs <- numeric(nrow(R))
  }
  rhs <- as.numeric(rhs)
  if (length(rhs) != nrow(R) || any(!is.finite(rhs))) {
    stop("rhs must be finite and have one entry per row of R.")
  }

  reduced <- hl_fo_reduce_linear_rows(R, rhs)
  rank_info <- .hl_fo_inference_rank(R)
  if (!isTRUE(reduced$success)) {
    return(list(
      ok = FALSE,
      status = "contrast_rhs_inconsistent",
      rank = reduced$rank,
      R = NULL,
      rhs = NULL,
      selected_rows = reduced$selected_rows,
      consistency_residual = reduced$consistency$residual,
      consistency_tolerance = reduced$consistency$tolerance,
      rank_info = rank_info,
      utility_reduction = reduced
    ))
  }
  if (reduced$rank == 0L) {
    return(list(
      ok = FALSE,
      status = "contrast_rank_zero",
      rank = 0L,
      R = reduced$matrix,
      rhs = reduced$rhs,
      selected_rows = reduced$selected_rows,
      rank_info = rank_info,
      utility_reduction = reduced
    ))
  }
  list(
    ok = TRUE,
    status = "ok",
    rank = reduced$rank,
    R = reduced$matrix,
    rhs = reduced$rhs,
    selected_rows = reduced$selected_rows,
    consistency_residual = reduced$consistency$residual,
    consistency_tolerance = reduced$consistency$tolerance,
    rank_info = rank_info,
    utility_reduction = reduced
  )
}

.hl_fo_canonical_tag <- function(fit) {
  tag <- fit$coordinate_tag
  if (is.null(tag) && !is.null(fit$coordinates)) {
    tag <- fit$coordinates$tag
  }
  tag
}

hl_fo_unrestricted_wald <- function(
    fit,
    R,
    rhs = NULL,
    expected_coordinate_tag = "canonical_conditional_beta_tau_xi") {
  tag <- .hl_fo_canonical_tag(fit)
  if (is.null(tag) || length(tag) != 1L ||
      !as.character(tag) %in% expected_coordinate_tag) {
    return(list(
      method = "wald",
      stat = NA_real_,
      df = NA_integer_,
      p = NA_real_,
      m = if (!is.null(fit$data$m)) fit$data$m else NA_integer_,
      contrast_rank = NA_integer_,
      covariance_rank = NA_integer_,
      success = FALSE,
      status = "wald_coordinate_mismatch",
      coordinate_tag = tag
    ))
  }

  p <- length(fit$beta)
  constraint <- fit$constraint
  unrestricted <- is.list(constraint) && isTRUE(constraint$success) &&
    identical(as.integer(constraint$rank), 0L) &&
    length(constraint$beta_star) == p &&
    all(is.finite(constraint$beta_star)) &&
    max(abs(constraint$beta_star), 0) <= 10 * .Machine$double.eps &&
    is.matrix(constraint$K) && identical(dim(constraint$K), c(p, p)) &&
    hl_fo_sdisc(constraint$K, diag(p)) <= 10 * .Machine$double.eps
  if (!unrestricted) {
    return(list(
      method = "wald",
      stat = NA_real_,
      df = NA_integer_,
      p = NA_real_,
      m = if (!is.null(fit$data$m)) fit$data$m else NA_integer_,
      contrast_rank = NA_integer_,
      covariance_rank = NA_integer_,
      success = FALSE,
      status = "wald_requires_unrestricted_fit",
      coordinate_tag = tag
    ))
  }
  if (!isTRUE(fit$point_success)) {
    return(list(
      method = "wald",
      stat = NA_real_,
      df = NA_integer_,
      p = NA_real_,
      m = if (!is.null(fit$data$m)) fit$data$m else NA_integer_,
      contrast_rank = NA_integer_,
      covariance_rank = NA_integer_,
      success = FALSE,
      status = "wald_unrestricted_fit_not_certified",
      coordinate_tag = tag
    ))
  }
  if (!isTRUE(fit$covariance_success)) {
    return(list(
      method = "wald",
      stat = NA_real_,
      df = NA_integer_,
      p = NA_real_,
      m = if (!is.null(fit$data$m)) fit$data$m else NA_integer_,
      contrast_rank = NA_integer_,
      covariance_rank = NA_integer_,
      success = FALSE,
      status = "wald_unrestricted_covariance_not_certified",
      coordinate_tag = tag
    ))
  }

  beta <- as.numeric(fit$beta)
  V_beta <- fit$covariance$V_beta
  m <- fit$data$m
  if (!length(beta) || !is.matrix(V_beta) ||
      !identical(dim(V_beta), c(length(beta), length(beta))) ||
      any(!is.finite(c(beta, V_beta))) || length(m) != 1L ||
      !is.finite(m) || m < 1) {
    stop("fit must expose finite beta, covariance$V_beta, and positive data$m.")
  }
  core <- hl_fo_robust_wald(
    beta = beta,
    covariance = V_beta,
    R = R,
    rhs = rhs,
    n_subjects = as.integer(m)
  )
  list(
    method = "wald",
    stat = core$statistic,
    df = core$df,
    p = core$p_value,
    m = as.integer(m),
    contrast_rank = core$contrast_rank,
    covariance_rank = core$covariance_rank,
    success = isTRUE(core$success),
    status = core$status,
    calibration = core$calibration,
    coordinate_tag = tag,
    delta = core$error,
    covariance = core$contrast_covariance,
    contrast = list(
      R = core$contrast,
      rhs = core$rhs,
      selected_rows = core$selected_rows
    ),
    beta_covariance_rank = core$beta_covariance_rank
  )
}

.hl_fo_score_failure <- function(status,
                                 m,
                                 contrast_rank,
                                 covariance_rank = NA_integer_,
                                 extras = list()) {
  c(list(
    method = "score_ee",
    stat = NA_real_,
    df = contrast_rank,
    p = NA_real_,
    m = as.integer(m),
    contrast_rank = contrast_rank,
    covariance_rank = covariance_rank,
    success = FALSE,
    status = status
  ), extras)
}

hl_fo_score_postroot <- function(psi_subject,
                                 A_zeta,
                                 R,
                                 beta,
                                 rhs = NULL,
                                 q = 2L) {
  psi_subject <- as.matrix(psi_subject)
  A_zeta <- as.matrix(A_zeta)
  if (is.null(dim(R)) && is.numeric(R)) {
    R <- matrix(R, nrow = 1L)
  } else {
    R <- as.matrix(R)
  }
  beta <- as.numeric(beta)
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) ||
      q != as.integer(q) || !q %in% 0:2) {
    stop("q must be one of 0, 1, or 2.")
  }
  q <- as.integer(q)
  p <- length(beta)
  d <- 2L * p + q

  if (!length(beta) || any(!is.finite(beta))) {
    stop("beta must be a nonempty finite vector.")
  }
  if (!is.numeric(psi_subject) || nrow(psi_subject) != d ||
      ncol(psi_subject) < 1L || any(!is.finite(psi_subject))) {
    stop("psi_subject must be a finite (2p + q) by m matrix with subjects in columns.")
  }
  if (!is.numeric(A_zeta) || !identical(dim(A_zeta), c(d, d)) ||
      any(!is.finite(A_zeta))) {
    stop("A_zeta must be a finite (2p + q) square sensitivity matrix.")
  }
  if (ncol(R) != p) {
    stop("R must have one column per beta coordinate.")
  }
  m <- ncol(psi_subject)
  reduced <- .hl_fo_reduce_contrast(R, rhs)
  if (!isTRUE(reduced$ok)) {
    return(.hl_fo_score_failure(
      status = paste0("score_", reduced$status),
      m = m,
      contrast_rank = reduced$rank,
      extras = list(contrast = reduced)
    ))
  }

  R_use <- reduced$R
  rhs_use <- reduced$rhs
  r <- reduced$rank
  k <- p - r
  RRt <- tcrossprod(R_use)
  C_R <- t(R_use) %*% solve(RRt)
  K <- .hl_fo_null_basis(R_use)
  Q <- cbind(C_R, K)

  T <- matrix(0, nrow = d, ncol = d)
  T[seq_len(p), seq_len(p)] <- Q
  T[p + seq_len(p), p + seq_len(p)] <- Q
  if (q > 0L) {
    rho_idx <- 2L * p + seq_len(q)
    T[rho_idx, rho_idx] <- diag(q)
  }

  constraint_residual <- as.numeric(R_use %*% beta - rhs_use)
  constraint_scale <- 1 + norm(R_use, type = "I") * max(abs(beta)) +
    max(abs(rhs_use))
  constraint_scaled_residual <- max(abs(constraint_residual)) / constraint_scale
  if (constraint_scaled_residual > 100 * .Machine$double.eps) {
    return(.hl_fo_score_failure(
      status = "score_null_constraint_failure",
      m = m,
      contrast_rank = r,
      extras = list(
        contrast = reduced,
        constraint_residual = constraint_residual,
        constraint_scaled_residual = constraint_scaled_residual
      )
    ))
  }

  u_subject <- crossprod(T, psi_subject)
  A <- crossprod(T, A_zeta %*% T)
  tested_idx <- seq_len(r)
  beta_nuisance_idx <- if (k > 0L) (r + 1L):p else integer(0)
  chi_idx <- p + seq_len(p)
  rho_idx <- if (q > 0L) 2L * p + seq_len(q) else integer(0)
  nuisance_idx <- c(beta_nuisance_idx, chi_idx, rho_idx)

  A_tt <- A[tested_idx, tested_idx, drop = FALSE]
  A_tn <- A[tested_idx, nuisance_idx, drop = FALSE]
  A_nt <- A[nuisance_idx, tested_idx, drop = FALSE]
  A_nn <- A[nuisance_idx, nuisance_idx, drop = FALSE]
  nuisance_rank <- .hl_fo_inference_rank(A_nn)$rank
  if (nuisance_rank < length(nuisance_idx)) {
    return(.hl_fo_score_failure(
      status = "score_nuisance_singular",
      m = m,
      contrast_rank = r,
      extras = list(
        contrast = reduced,
        C_R = C_R,
        K = K,
        Q = Q,
        T = T,
        transformed_sensitivity = A,
        nuisance_rank = nuisance_rank
      )
    ))
  }

  H <- tryCatch(
    t(solve(t(A_nn), t(A_tn))),
    error = function(e) NULL
  )
  if (is.null(H) || any(!is.finite(H))) {
    return(.hl_fo_score_failure(
      status = "score_nuisance_singular",
      m = m,
      contrast_rank = r,
      extras = list(nuisance_rank = nuisance_rank, contrast = reduced)
    ))
  }

  u_tested <- u_subject[tested_idx, , drop = FALSE]
  u_nuisance <- u_subject[nuisance_idx, , drop = FALSE]
  q_subject <- u_tested - H %*% u_nuisance
  adjusted_meat <- tcrossprod(q_subject)

  J <- tcrossprod(u_subject)
  J_tt <- J[tested_idx, tested_idx, drop = FALSE]
  J_tn <- J[tested_idx, nuisance_idx, drop = FALSE]
  J_nt <- J[nuisance_idx, tested_idx, drop = FALSE]
  J_nn <- J[nuisance_idx, nuisance_idx, drop = FALSE]
  adjusted_meat_block <- J_tt - H %*% J_nt - J_tn %*% t(H) +
    H %*% J_nn %*% t(H)
  meat_crosscheck <- max(
    abs(adjusted_meat - adjusted_meat_block) /
      (1 + pmax(abs(adjusted_meat), abs(adjusted_meat_block)))
  )
  if (!is.finite(meat_crosscheck) || meat_crosscheck > 1e-10) {
    return(.hl_fo_score_failure(
      status = "score_meat_crosscheck_failure",
      m = m,
      contrast_rank = r,
      extras = list(
        adjusted_meat = adjusted_meat,
        adjusted_meat_block = adjusted_meat_block,
        meat_crosscheck = meat_crosscheck,
        contrast = reduced
      )
    ))
  }

  adjusted_sensitivity <- A_tt - H %*% A_nt
  adjusted_sensitivity_rank <- .hl_fo_inference_rank(adjusted_sensitivity)$rank
  if (adjusted_sensitivity_rank < r) {
    return(.hl_fo_score_failure(
      status = "score_adjusted_sensitivity_singular",
      m = m,
      contrast_rank = r,
      extras = list(
        adjusted_sensitivity = adjusted_sensitivity,
        adjusted_sensitivity_rank = adjusted_sensitivity_rank,
        contrast = reduced
      )
    ))
  }

  adjusted_meat_rank <- .hl_fo_inference_rank(adjusted_meat)$rank
  if (adjusted_meat_rank < r) {
    return(.hl_fo_score_failure(
      status = "score_meat_rank_deficient",
      m = m,
      contrast_rank = r,
      covariance_rank = adjusted_meat_rank,
      extras = list(
        adjusted_meat = adjusted_meat,
        adjusted_meat_rank = adjusted_meat_rank,
        contrast = reduced
      )
    ))
  }

  U_adjusted <- rowSums(q_subject)
  meat_solve <- tryCatch(solve(adjusted_meat, U_adjusted), error = function(e) NULL)
  if (is.null(meat_solve) || any(!is.finite(meat_solve))) {
    return(.hl_fo_score_failure(
      status = "score_meat_rank_deficient",
      m = m,
      contrast_rank = r,
      covariance_rank = adjusted_meat_rank,
      extras = list(adjusted_meat = adjusted_meat, contrast = reduced)
    ))
  }
  stat <- as.numeric(crossprod(U_adjusted, meat_solve))
  if (stat < 0 && stat >= -1e-10 * max(1, abs(stat))) {
    stat <- 0
  }
  if (!is.finite(stat) || stat < 0) {
    return(.hl_fo_score_failure(
      status = "score_statistic_invalid",
      m = m,
      contrast_rank = r,
      covariance_rank = adjusted_meat_rank,
      extras = list(adjusted_meat = adjusted_meat, contrast = reduced)
    ))
  }

  list(
    method = "score_ee",
    stat = stat,
    df = r,
    p = stats::pchisq(stat, df = r, lower.tail = FALSE),
    m = as.integer(m),
    contrast_rank = r,
    covariance_rank = adjusted_meat_rank,
    success = TRUE,
    status = "ok",
    calibration = "asymptotic_chisq_subject_clusters",
    contrast = reduced,
    constraint_residual = constraint_residual,
    constraint_scaled_residual = constraint_scaled_residual,
    C_R = C_R,
    K = K,
    Q = Q,
    T = T,
    tested_indices = tested_idx,
    nuisance_indices = nuisance_idx,
    transformed_subject = u_subject,
    transformed_sensitivity = A,
    H = H,
    q_subject = q_subject,
    adjusted_score = U_adjusted,
    adjusted_meat = adjusted_meat,
    adjusted_meat_block = adjusted_meat_block,
    meat_crosscheck = meat_crosscheck,
    adjusted_sensitivity = adjusted_sensitivity,
    nuisance_rank = nuisance_rank,
    adjusted_sensitivity_rank = adjusted_sensitivity_rank,
    adjusted_meat_rank = adjusted_meat_rank,
    J_rank = .hl_fo_inference_rank(J)$rank
  )
}

hl_fo_objective_ratio_unsupported <- function(
    raw_difference = NA_real_,
    failure_reason = "beta_constraint_objective_ratio_unsupported",
    wald = NULL) {
  if (length(failure_reason) != 1L || is.na(failure_reason) ||
      !nzchar(failure_reason)) {
    stop("failure_reason must be a nonempty audit reason.")
  }
  if (length(raw_difference) != 1L) {
    stop("raw_difference must be a scalar audit diagnostic.")
  }
  if (!is.null(wald)) {
    if (!is.list(wald)) {
      stop("wald must be NULL or a Wald result list.")
    }
    wald$method <- "wald"
  }
  list(
    method = "unsupported_objective_ratio",
    supported = FALSE,
    success = FALSE,
    status = "objective_ratio_unsupported",
    failure_reason = as.character(failure_reason),
    audit = list(
      raw_difference = as.numeric(raw_difference),
      calibration = "unsupported_audit_only"
    ),
    wald = wald
  )
}
