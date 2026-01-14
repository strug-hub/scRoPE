apply_chain_rule_transform <- function(U,
                                       H,
                                       S = NULL,
                                       intercept_idx,
                                       p,
                                       U_per_subject = NULL,
                                       sigma2 = NULL) {
  dim_expected <- p + 2L
  if (length(intercept_idx) != 1L) {
    stop("intercept_idx must identify exactly one intercept column.")
  }
  if (!is.numeric(intercept_idx) || intercept_idx < 1L || intercept_idx > p) {
    stop("intercept_idx must be an integer in [1, p].")
  }
  if (length(U) != dim_expected) {
    stop("Length of U must equal p + 2 (beta, sigma^2, phi).")
  }
  if (!is.matrix(H) || nrow(H) != dim_expected || ncol(H) != dim_expected) {
    stop("H must be a (p + 2) x (p + 2) matrix.")
  }
  if (!is.null(S) && (!is.matrix(S) || nrow(S) != dim_expected || ncol(S) != dim_expected)) {
    stop("S must be NULL or a (p + 2) x (p + 2) matrix.")
  }
  if (!is.null(U_per_subject) && (!is.matrix(U_per_subject) || nrow(U_per_subject) != dim_expected)) {
    stop("U_per_subject must be NULL or a (p + 2) x k matrix.")
  }

  k_sigma <- p + 1L
  alpha <- 0.5
  factor <- if (is.null(sigma2)) alpha else alpha * sigma2

  ## Score
  U_prime <- as.numeric(U)
  U_prime[k_sigma] <- U[k_sigma] + factor * U[intercept_idx]

  ## Hessian
  H_prime <- H
  col_k_new <- H[, k_sigma] + factor * H[, intercept_idx]
  row_k_new <- H[k_sigma, ] + factor * H[intercept_idx, ]
  diag_k_new <- H[k_sigma, k_sigma] +
    factor * H[intercept_idx, k_sigma] +
    factor * H[k_sigma, intercept_idx] +
    factor * factor * H[intercept_idx, intercept_idx]
  if (!is.null(sigma2)) {
    diag_k_new <- diag_k_new + factor * U[intercept_idx]
  }
  H_prime[, k_sigma] <- col_k_new
  H_prime[k_sigma, ] <- row_k_new
  H_prime[k_sigma, k_sigma] <- diag_k_new
  H_prime <- 0.5 * (H_prime + t(H_prime))

  ## Meat (scalar update or via per-subject gradients)
  if (!is.null(U_per_subject)) {
    U_per_subject_prime <- U_per_subject
    U_per_subject_prime[k_sigma, ] <- U_per_subject[k_sigma, ] + factor * U_per_subject[intercept_idx, ]
  } else {
    U_per_subject_prime <- NULL
  }

  if (!is.null(S)) {
    S_prime <- S
    col_s_new <- S[, k_sigma] + factor * S[, intercept_idx]
    row_s_new <- S[k_sigma, ] + factor * S[intercept_idx, ]
    diag_s_new <- S[k_sigma, k_sigma] +
      factor * S[intercept_idx, k_sigma] +
      factor * S[k_sigma, intercept_idx] +
      factor * factor * S[intercept_idx, intercept_idx]
    S_prime[, k_sigma] <- col_s_new
    S_prime[k_sigma, ] <- row_s_new
    S_prime[k_sigma, k_sigma] <- diag_s_new
    S_prime <- 0.5 * (S_prime + t(S_prime))
  } else if (!is.null(U_per_subject_prime)) {
    S_prime <- U_per_subject_prime %*% t(U_per_subject_prime)
  } else {
    S_prime <- NULL
  }

  if (isTRUE(getOption("scRoPE.debug", FALSE))) {
    delta <- U_prime[k_sigma] - (U[k_sigma] + factor * U[intercept_idx])
    if (abs(delta) > 1e-8 * max(1, abs(U))) {
      stop("Chain-rule score check failed: transformed sigma^2 score mismatch.")
    }
  }

  list(
    U = U_prime,
    H = H_prime,
    S = S_prime,
    U_per_subject = U_per_subject_prime
  )
}
