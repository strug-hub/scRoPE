# Numerical, coordinate, and affine-constraint utilities for the corrected
# first-order nested HL pipeline.  These helpers deliberately have no
# dependency on the legacy PTMG-at-HL-point derivative state.

.hl_fo_assert_numeric_matrix <- function(x, name, square = FALSE) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (square && nrow(x) != ncol(x)) {
    stop(name, " must be square.", call. = FALSE)
  }
  x
}

.hl_fo_assert_scalar <- function(x, name, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(name, " must be one finite numeric value.", call. = FALSE)
  }
  if (positive && x <= 0) {
    stop(name, " must be strictly positive.", call. = FALSE)
  }
  as.numeric(x)
}

.hl_fo_assert_beta <- function(beta, name = "beta") {
  if (!is.numeric(beta) || any(!is.finite(beta))) {
    stop(name, " must be a finite numeric vector.", call. = FALSE)
  }
  as.numeric(beta)
}

.hl_fo_intercept_index <- function(intercept_idx, p) {
  if (!is.numeric(intercept_idx) || length(intercept_idx) != 1L ||
      !is.finite(intercept_idx) || intercept_idx != as.integer(intercept_idx) ||
      intercept_idx < 1L || intercept_idx > p) {
    stop("intercept_idx must be one integer in [1, length(beta)].",
         call. = FALSE)
  }
  as.integer(intercept_idx)
}

.hl_fo_norm_inf <- function(x) {
  if (length(x) == 0L) {
    return(0)
  }
  max(abs(x))
}

.hl_fo_norm2 <- function(x) {
  if (length(x) == 0L) {
    return(0)
  }
  scale <- max(abs(x))
  if (scale == 0) {
    return(0)
  }
  scale * sqrt(sum((x / scale)^2))
}

.hl_fo_blockdiag <- function(...) {
  blocks <- list(...)
  if (length(blocks) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }
  if (any(!vapply(blocks, is.matrix, logical(1)))) {
    stop("Every block must be a matrix.", call. = FALSE)
  }
  nr <- sum(vapply(blocks, nrow, integer(1)))
  nc <- sum(vapply(blocks, ncol, integer(1)))
  out <- matrix(0, nrow = nr, ncol = nc)
  row_at <- 1L
  col_at <- 1L
  for (block in blocks) {
    if (nrow(block) > 0L && ncol(block) > 0L) {
      rr <- row_at:(row_at + nrow(block) - 1L)
      cc <- col_at:(col_at + ncol(block) - 1L)
      out[rr, cc] <- block
    }
    row_at <- row_at + nrow(block)
    col_at <- col_at + ncol(block)
  }
  out
}

# Componentwise scaled discrepancy in Specification Equation (45).
hl_fo_sdisc <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b) || length(a) != length(b) ||
      !identical(dim(a), dim(b))) {
    stop("a and b must be equally shaped numeric objects.", call. = FALSE)
  }
  if (length(a) == 0L) {
    return(0)
  }
  finite <- is.finite(a) & is.finite(b)
  component <- rep(Inf, length(a))
  component[finite] <- abs(a[finite] - b[finite]) /
    (1 + pmax(abs(a[finite]), abs(b[finite])))
  max(component)
}

# SVD rank and Moore--Penrose inverse under the one package-wide convention.
hl_fo_svd_info <- function(x, compute_pinv = TRUE, complete_left = FALSE) {
  x <- .hl_fo_assert_numeric_matrix(x, "x")
  nr <- nrow(x)
  nc <- ncol(x)
  min_dim <- min(nr, nc)

  if (min_dim == 0L) {
    pinv <- if (compute_pinv) matrix(0, nrow = nc, ncol = nr) else NULL
    return(list(
      rank = 0L,
      threshold = 0,
      singular_values = numeric(0),
      retained = logical(0),
      reciprocal_condition = 0,
      left_null_basis = if (complete_left) diag(nr) else NULL,
      pinv = pinv
    ))
  }

  # RHS consistency requests the complete left basis.  Ordinary rank checks do
  # not, because a complete n by n U would be inappropriate for a tall design.
  nu <- if (complete_left) nr else min_dim
  decomp <- base::svd(x, nu = nu, nv = min_dim)
  s1 <- if (length(decomp$d)) decomp$d[[1L]] else 0
  threshold <- max(nr, nc) * .Machine$double.eps * s1
  retained <- decomp$d > threshold
  rank <- sum(retained)

  pinv <- NULL
  if (compute_pinv) {
    if (rank == 0L) {
      pinv <- matrix(0, nrow = nc, ncol = nr)
    } else {
      kept <- which(retained)
      vr <- decomp$v[, kept, drop = FALSE]
      ur <- decomp$u[, kept, drop = FALSE]
      pinv <- vr %*% (diag(1 / decomp$d[kept], nrow = rank) %*% t(ur))
    }
  }

  reciprocal_condition <- if (rank == min_dim && rank > 0L) {
    decomp$d[[rank]] / decomp$d[[1L]]
  } else {
    0
  }
  left_null_basis <- if (complete_left && rank < nr) {
    decomp$u[, (rank + 1L):nr, drop = FALSE]
  } else if (complete_left) {
    matrix(0, nrow = nr, ncol = 0L)
  } else {
    NULL
  }

  list(
    rank = as.integer(rank),
    threshold = threshold,
    singular_values = decomp$d,
    retained = retained,
    reciprocal_condition = reciprocal_condition,
    left_null_basis = left_null_basis,
    pinv = pinv
  )
}

hl_fo_num_rank <- function(x) {
  info <- hl_fo_svd_info(x, compute_pinv = FALSE)
  list(
    rank = info$rank,
    tolerance = info$threshold,
    threshold = info$threshold,
    singular_values = info$singular_values
  )
}

hl_fo_pinv <- function(x) {
  hl_fo_svd_info(x, compute_pinv = TRUE)$pinv
}

hl_fo_rhs_consistency <- function(M, rhs) {
  M <- .hl_fo_assert_numeric_matrix(M, "M")
  if (!is.numeric(rhs) || length(rhs) != nrow(M) || any(!is.finite(rhs))) {
    stop("rhs must be a finite vector with length nrow(M).", call. = FALSE)
  }
  rhs <- as.numeric(rhs)
  info <- hl_fo_svd_info(M, compute_pinv = TRUE, complete_left = TRUE)
  # Once the declared numerical rank equals the number of equations, every
  # right-hand side is in the row space algebraically.  Avoid misclassifying a
  # full-row-rank system because the explicitly formed SVD projector accumulated
  # a few ulps.
  residual <- if (info$rank == nrow(M)) {
    numeric(length(rhs))
  } else if (nrow(M) == 0L) {
    numeric(0)
  } else {
    as.vector(
      info$left_null_basis %*% crossprod(info$left_null_basis, rhs)
    )
  }
  tolerance <- max(nrow(M), ncol(M)) * .Machine$double.eps *
    (1 + .hl_fo_norm_inf(rhs))
  residual_norm <- .hl_fo_norm_inf(residual)
  list(
    consistent = isTRUE(residual_norm <= tolerance),
    residual = residual,
    residual_norm = residual_norm,
    tolerance = tolerance,
    rank = info$rank,
    threshold = info$threshold
  )
}

hl_fo_symmetry_defect <- function(x) {
  x <- .hl_fo_assert_numeric_matrix(x, "x", square = TRUE)
  if (length(x) == 0L) {
    return(0)
  }
  sqrt(sum((x - t(x))^2)) / max(1, sqrt(sum(x^2)))
}

.hl_fo_householder_standard <- function(x) {
  norm_x <- .hl_fo_norm2(x)
  if (norm_x == 0) {
    return(list(v = x, beta = 0))
  }
  # Reflectors are invariant to rescaling v.  Work on the unit scale to avoid
  # overflow in v^T v for otherwise finite design or constraint entries.
  v <- x / norm_x
  target <- if (v[[1L]] >= 0) -1 else 1
  v[[1L]] <- v[[1L]] - target
  list(v = v, beta = 2 / sum(v^2))
}

# A stable Householder vector whose transformed leading entry is positive.
.hl_fo_householder_positive <- function(x) {
  norm_x <- .hl_fo_norm2(x)
  if (norm_x == 0) {
    return(list(v = x, beta = 0))
  }
  x_scaled <- x / norm_x
  if (length(x) == 1L || all(x_scaled[-1L] == 0)) {
    if (x[[1L]] > 0) {
      return(list(v = x * 0, beta = 0))
    }
    v <- x * 0
    v[[1L]] <- 1
    return(list(v = v, beta = 2))
  }

  v <- x_scaled
  if (x_scaled[[1L]] > 0) {
    tail_norm <- .hl_fo_norm2(x_scaled[-1L])
    if (tail_norm <= .Machine$double.eps) {
      return(list(v = x * 0, beta = 0))
    }
    v[[1L]] <- -(tail_norm^2) / (x_scaled[[1L]] + 1)
  } else {
    v[[1L]] <- x_scaled[[1L]] - 1
  }
  list(v = v, beta = 2 / sum(v^2))
}

.hl_fo_cpqr_selected_columns <- function(A, rank) {
  A <- .hl_fo_assert_numeric_matrix(A, "A")
  rank <- as.integer(rank)
  if (rank < 0L || rank > min(nrow(A), ncol(A))) {
    stop("rank is incompatible with A.", call. = FALSE)
  }
  if (rank == 0L) {
    return(integer(0))
  }

  R <- A
  pivot <- seq_len(ncol(A))
  for (kk in seq_len(rank)) {
    remaining <- kk:ncol(R)
    norms <- vapply(remaining, function(jj) {
      .hl_fo_norm2(R[kk:nrow(R), jj])
    }, numeric(1))
    max_norm <- max(norms)
    tied <- remaining[norms == max_norm]
    chosen <- tied[[which.min(pivot[tied])]]
    if (chosen != kk) {
      R[, c(kk, chosen)] <- R[, c(chosen, kk), drop = FALSE]
      pivot[c(kk, chosen)] <- pivot[c(chosen, kk)]
    }

    hh <- .hl_fo_householder_standard(R[kk:nrow(R), kk])
    if (hh$beta != 0) {
      block <- R[kk:nrow(R), kk:ncol(R), drop = FALSE]
      projection <- as.vector(crossprod(hh$v, block))
      R[kk:nrow(R), kk:ncol(R)] <-
        block - hh$beta * tcrossprod(hh$v, projection)
    }
  }
  pivot[seq_len(rank)]
}

.hl_fo_null_basis <- function(C) {
  C <- .hl_fo_assert_numeric_matrix(C, "C")
  r <- nrow(C)
  p <- ncol(C)
  if (r > p) {
    stop("C must have no more independent rows than columns.", call. = FALSE)
  }
  if (r == 0L) {
    return(diag(p))
  }

  R <- t(C)
  Q <- diag(p)
  for (kk in seq_len(r)) {
    hh <- .hl_fo_householder_positive(R[kk:p, kk])
    if (hh$beta != 0) {
      r_block <- R[kk:p, kk:r, drop = FALSE]
      r_projection <- as.vector(crossprod(hh$v, r_block))
      R[kk:p, kk:r] <-
        r_block - hh$beta * tcrossprod(hh$v, r_projection)

      q_block <- Q[, kk:p, drop = FALSE]
      q_projection <- as.vector(q_block %*% hh$v)
      Q[, kk:p] <- q_block - hh$beta * tcrossprod(q_projection, hh$v)
    }
  }

  # Roundoff can leave a signed zero or a tiny negative diagonal.  Flipping a
  # leading Q column enforces the declared positive-diagonal convention and
  # leaves the null-space columns unchanged.
  diagonal <- diag(R[seq_len(r), seq_len(r), drop = FALSE])
  for (kk in which(diagonal < 0)) {
    Q[, kk] <- -Q[, kk]
    R[kk, ] <- -R[kk, ]
  }

  if (r == p) {
    return(matrix(0, nrow = p, ncol = 0L))
  }
  K <- Q[, (r + 1L):p, drop = FALSE]
  for (jj in seq_len(ncol(K))) {
    magnitude <- abs(K[, jj])
    pivot <- which(magnitude == max(magnitude))[[1L]]
    if (K[pivot, jj] < 0) {
      K[, jj] <- -K[, jj]
    }
  }
  K
}

# Deterministically remove redundant rows while retaining original-row labels.
hl_fo_reduce_linear_rows <- function(M, rhs) {
  if (is.null(dim(M)) && is.numeric(M)) {
    M <- matrix(M, nrow = 1L)
  }
  M <- .hl_fo_assert_numeric_matrix(M, "M")
  if (!is.numeric(rhs) || length(rhs) != nrow(M) || any(!is.finite(rhs))) {
    stop("rhs must be a finite vector with length nrow(M).", call. = FALSE)
  }
  rhs <- as.numeric(rhs)
  svd_info <- hl_fo_svd_info(M, compute_pinv = FALSE)
  consistency <- hl_fo_rhs_consistency(M, rhs)
  if (!consistency$consistent) {
    return(list(
      success = FALSE,
      status = "constraint_infeasible",
      rank = svd_info$rank,
      rank_threshold = svd_info$threshold,
      selected_rows = integer(0),
      discarded_rows = seq_len(nrow(M)),
      matrix = matrix(0, nrow = 0L, ncol = ncol(M)),
      rhs = numeric(0),
      original_matrix = M,
      original_rhs = rhs,
      consistency = consistency
    ))
  }

  selected <- .hl_fo_cpqr_selected_columns(t(M), svd_info$rank)
  selected <- sort(selected)
  discarded <- setdiff(seq_len(nrow(M)), selected)
  reduced <- M[selected, , drop = FALSE]
  reduced_rhs <- rhs[selected]
  if (hl_fo_num_rank(reduced)$rank != svd_info$rank) {
    stop("Deterministic row reduction failed to retain the SVD rank.",
         call. = FALSE)
  }

  list(
    success = TRUE,
    status = "ok",
    rank = svd_info$rank,
    rank_threshold = svd_info$threshold,
    selected_rows = selected,
    discarded_rows = discarded,
    matrix = reduced,
    rhs = reduced_rhs,
    original_matrix = M,
    original_rhs = rhs,
    consistency = consistency
  )
}

hl_fo_constraint_residual <- function(C, beta, rhs) {
  if (is.null(dim(C)) && is.numeric(C)) {
    C <- matrix(C, nrow = 1L)
  }
  C <- .hl_fo_assert_numeric_matrix(C, "C")
  beta <- .hl_fo_assert_beta(beta)
  if (ncol(C) != length(beta)) {
    stop("ncol(C) must equal length(beta).", call. = FALSE)
  }
  if (!is.numeric(rhs) || length(rhs) != nrow(C) || any(!is.finite(rhs))) {
    stop("rhs must be finite with length nrow(C).", call. = FALSE)
  }
  rhs <- as.numeric(rhs)
  residual <- as.vector(C %*% beta - rhs)
  c_norm_inf <- if (nrow(C) == 0L) 0 else max(rowSums(abs(C)))
  scaled <- .hl_fo_norm_inf(residual) /
    (1 + c_norm_inf * .hl_fo_norm_inf(beta) + .hl_fo_norm_inf(rhs))
  tolerance <- 100 * .Machine$double.eps
  list(
    residual = residual,
    max_abs = .hl_fo_norm_inf(residual),
    scaled = scaled,
    tolerance = tolerance,
    satisfied = isTRUE(scaled <= tolerance)
  )
}

hl_fo_affine_constraint <- function(C = NULL, rhs = numeric(0), p = NULL) {
  if (is.null(C)) {
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p) ||
        p != as.integer(p) || p < 0L) {
      stop("p must be a nonnegative integer when C is NULL.", call. = FALSE)
    }
    p <- as.integer(p)
    C <- matrix(0, nrow = 0L, ncol = p)
  } else {
    if (is.null(dim(C)) && is.numeric(C)) {
      C <- matrix(C, nrow = 1L)
    }
    C <- .hl_fo_assert_numeric_matrix(C, "C")
    if (!is.null(p) && ncol(C) != p) {
      stop("p does not agree with ncol(C).", call. = FALSE)
    }
    p <- ncol(C)
  }

  if (length(rhs) == 0L && nrow(C) > 0L) {
    rhs <- numeric(nrow(C))
  }

  reduced <- hl_fo_reduce_linear_rows(C, rhs)
  if (!reduced$success) {
    reduced$beta_star <- rep(NA_real_, p)
    reduced$K <- matrix(NA_real_, nrow = p, ncol = 0L)
    reduced$null_basis <- reduced$K
    reduced$free_dim <- NA_integer_
    class(reduced) <- c("hl_fo_affine_constraint", "list")
    return(reduced)
  }

  rank <- reduced$rank
  C_reduced <- reduced$matrix
  beta_star <- if (rank == 0L) {
    numeric(p)
  } else {
    C_pinv <- hl_fo_pinv(C_reduced)
    beta_anchor <- as.vector(C_pinv %*% reduced$rhs)
    # Iterative refinement stays in col(C^T), so it preserves the minimum-norm
    # property while making the true constraint residual as accurate as the
    # stable linear solve permits.
    for (iteration in seq_len(3L)) {
      anchor_residual <- as.vector(C_reduced %*% beta_anchor - reduced$rhs)
      if (.hl_fo_norm_inf(anchor_residual) == 0) {
        break
      }
      beta_anchor <- beta_anchor - as.vector(C_pinv %*% anchor_residual)
    }
    beta_anchor
  }
  K <- .hl_fo_null_basis(C_reduced)
  residual_check <- hl_fo_constraint_residual(C, beta_star, rhs)

  out <- c(reduced, list(
    beta_star = beta_star,
    K = K,
    null_basis = K,
    free_dim = ncol(K),
    residual = residual_check$residual,
    residual_scaled = residual_check$scaled,
    residual_tolerance = residual_check$tolerance,
    residual_satisfied = residual_check$satisfied
  ))
  class(out) <- c("hl_fo_affine_constraint", "list")
  out
}

hl_fo_cov_congruence <- function(covariance, jacobian) {
  covariance <- .hl_fo_assert_numeric_matrix(
    covariance, "covariance", square = TRUE
  )
  jacobian <- .hl_fo_assert_numeric_matrix(jacobian, "jacobian")
  if (ncol(jacobian) != nrow(covariance)) {
    stop("ncol(jacobian) must equal nrow(covariance).", call. = FALSE)
  }
  jacobian %*% covariance %*% t(jacobian)
}

hl_fo_canonical_to_natural <- function(beta, tau, xi) {
  beta <- .hl_fo_assert_beta(beta)
  tau <- .hl_fo_assert_scalar(tau, "tau")
  xi <- .hl_fo_assert_scalar(xi, "xi")
  sigma2 <- exp(tau)
  phi <- exp(xi)
  if (!is.finite(sigma2) || !is.finite(phi)) {
    stop("Exponentiated dispersions must be finite.", call. = FALSE)
  }
  jacobian <- .hl_fo_blockdiag(
    diag(length(beta)), matrix(sigma2, 1L, 1L), matrix(phi, 1L, 1L)
  )
  list(
    beta = beta,
    sigma2 = sigma2,
    phi = phi,
    value = c(beta, sigma2, phi),
    jacobian = jacobian,
    from = "canonical",
    to = "natural_size"
  )
}

hl_fo_natural_to_canonical <- function(beta, sigma2, phi) {
  beta <- .hl_fo_assert_beta(beta)
  sigma2 <- .hl_fo_assert_scalar(sigma2, "sigma2", positive = TRUE)
  phi <- .hl_fo_assert_scalar(phi, "phi", positive = TRUE)
  jacobian <- .hl_fo_blockdiag(
    diag(length(beta)), matrix(1 / sigma2, 1L, 1L),
    matrix(1 / phi, 1L, 1L)
  )
  list(
    beta = beta,
    tau = log(sigma2),
    xi = log(phi),
    value = c(beta, log(sigma2), log(phi)),
    jacobian = jacobian,
    from = "natural_size",
    to = "canonical"
  )
}

hl_fo_canonical_to_public <- function(beta, tau, xi, beta_map = NULL) {
  beta <- .hl_fo_assert_beta(beta)
  tau <- .hl_fo_assert_scalar(tau, "tau")
  xi <- .hl_fo_assert_scalar(xi, "xi")
  p <- length(beta)
  if (is.null(beta_map)) {
    beta_map <- diag(p)
  }
  beta_map <- .hl_fo_assert_numeric_matrix(beta_map, "beta_map")
  if (ncol(beta_map) != p) {
    stop("ncol(beta_map) must equal length(beta).", call. = FALSE)
  }
  sigma2 <- exp(tau)
  delta <- exp(-xi)
  if (!is.finite(sigma2) || !is.finite(delta)) {
    stop("Transformed dispersions must be finite.", call. = FALSE)
  }
  jacobian <- .hl_fo_blockdiag(
    beta_map, matrix(sigma2, 1L, 1L), matrix(-delta, 1L, 1L)
  )
  beta_report <- as.vector(beta_map %*% beta)
  list(
    beta = beta_report,
    sigma2 = sigma2,
    delta = delta,
    value = c(beta_report, sigma2, delta),
    jacobian = jacobian,
    from = "canonical",
    to = "public_reciprocal"
  )
}

hl_fo_public_to_canonical <- function(beta, sigma2, delta,
                                      beta_inverse = NULL) {
  beta <- .hl_fo_assert_beta(beta, "public beta")
  sigma2 <- .hl_fo_assert_scalar(sigma2, "sigma2", positive = TRUE)
  delta <- .hl_fo_assert_scalar(delta, "delta", positive = TRUE)
  p_report <- length(beta)
  if (is.null(beta_inverse)) {
    beta_inverse <- diag(p_report)
  }
  beta_inverse <- .hl_fo_assert_numeric_matrix(
    beta_inverse, "beta_inverse"
  )
  if (ncol(beta_inverse) != p_report) {
    stop("ncol(beta_inverse) must equal length(beta).", call. = FALSE)
  }
  canonical_beta <- as.vector(beta_inverse %*% beta)
  jacobian <- .hl_fo_blockdiag(
    beta_inverse, matrix(1 / sigma2, 1L, 1L),
    matrix(-1 / delta, 1L, 1L)
  )
  list(
    beta = canonical_beta,
    tau = log(sigma2),
    xi = -log(delta),
    value = c(canonical_beta, log(sigma2), -log(delta)),
    jacobian = jacobian,
    from = "public_reciprocal",
    to = "canonical"
  )
}

hl_fo_canonical_to_mean_one <- function(beta, tau, xi, intercept_idx) {
  natural <- hl_fo_canonical_to_natural(beta, tau, xi)
  p <- length(natural$beta)
  intercept_idx <- .hl_fo_intercept_index(intercept_idx, p)
  beta_mean_one <- natural$beta
  beta_mean_one[intercept_idx] <-
    beta_mean_one[intercept_idx] + natural$sigma2 / 2
  jacobian <- .hl_fo_blockdiag(
    diag(p), matrix(natural$sigma2, 1L, 1L),
    matrix(natural$phi, 1L, 1L)
  )
  jacobian[intercept_idx, p + 1L] <- natural$sigma2 / 2
  list(
    beta = beta_mean_one,
    sigma2 = natural$sigma2,
    phi = natural$phi,
    value = c(beta_mean_one, natural$sigma2, natural$phi),
    jacobian = jacobian,
    from = "canonical",
    to = "mean_one_natural"
  )
}

hl_fo_mean_one_to_canonical <- function(beta, sigma2, phi, intercept_idx) {
  beta <- .hl_fo_assert_beta(beta, "mean-one beta")
  sigma2 <- .hl_fo_assert_scalar(sigma2, "sigma2", positive = TRUE)
  phi <- .hl_fo_assert_scalar(phi, "phi", positive = TRUE)
  p <- length(beta)
  intercept_idx <- .hl_fo_intercept_index(intercept_idx, p)
  canonical_beta <- beta
  canonical_beta[intercept_idx] <-
    canonical_beta[intercept_idx] - sigma2 / 2
  jacobian <- .hl_fo_blockdiag(
    diag(p), matrix(1 / sigma2, 1L, 1L), matrix(1 / phi, 1L, 1L)
  )
  jacobian[intercept_idx, p + 1L] <- -1 / 2
  list(
    beta = canonical_beta,
    tau = log(sigma2),
    xi = log(phi),
    value = c(canonical_beta, log(sigma2), log(phi)),
    jacobian = jacobian,
    from = "mean_one_natural",
    to = "canonical"
  )
}

hl_fo_design_map <- function(column_means, column_scales, intercept_idx,
                             raw_constant) {
  if (!is.numeric(column_means) || !is.numeric(column_scales) ||
      length(column_means) != length(column_scales) ||
      any(!is.finite(column_means)) || any(!is.finite(column_scales))) {
    stop("column_means and column_scales must be equally sized finite vectors.",
         call. = FALSE)
  }
  p <- length(column_means)
  intercept_idx <- .hl_fo_intercept_index(intercept_idx, p)
  raw_constant <- .hl_fo_assert_scalar(raw_constant, "raw_constant")
  if (raw_constant == 0) {
    stop("raw_constant must be nonzero.", call. = FALSE)
  }
  varying <- setdiff(seq_len(p), intercept_idx)
  if (any(column_scales[varying] <= 0)) {
    stop("Every varying-column scale must be strictly positive.",
         call. = FALSE)
  }

  to_original <- matrix(0, nrow = p, ncol = p)
  to_original[intercept_idx, intercept_idx] <- 1 / raw_constant
  if (length(varying)) {
    to_original[cbind(varying, varying)] <- 1 / column_scales[varying]
    to_original[intercept_idx, varying] <-
      -column_means[varying] /
      (raw_constant * column_scales[varying])
  }

  to_internal <- matrix(0, nrow = p, ncol = p)
  to_internal[intercept_idx, intercept_idx] <- raw_constant
  if (length(varying)) {
    to_internal[cbind(varying, varying)] <- column_scales[varying]
    to_internal[intercept_idx, varying] <- column_means[varying]
  }

  list(
    to_original = to_original,
    to_internal = to_internal,
    column_means = as.numeric(column_means),
    column_scales = as.numeric(column_scales),
    intercept_idx = intercept_idx,
    raw_constant = raw_constant
  )
}

hl_fo_internal_to_original_beta <- function(beta, design_map) {
  beta <- .hl_fo_assert_beta(beta)
  if (!is.list(design_map) || is.null(design_map$to_original)) {
    stop("design_map must be returned by hl_fo_design_map().", call. = FALSE)
  }
  M <- .hl_fo_assert_numeric_matrix(design_map$to_original, "design map")
  if (ncol(M) != length(beta)) {
    stop("design_map and beta dimensions do not agree.", call. = FALSE)
  }
  list(beta = as.vector(M %*% beta), jacobian = M)
}

hl_fo_original_to_internal_beta <- function(beta, design_map) {
  beta <- .hl_fo_assert_beta(beta)
  if (!is.list(design_map) || is.null(design_map$to_internal)) {
    stop("design_map must be returned by hl_fo_design_map().", call. = FALSE)
  }
  M <- .hl_fo_assert_numeric_matrix(design_map$to_internal, "design map")
  if (ncol(M) != length(beta)) {
    stop("design_map and beta dimensions do not agree.", call. = FALSE)
  }
  list(beta = as.vector(M %*% beta), jacobian = M)
}

hl_fo_transform_constraint_to_internal <- function(C_report, rhs, beta_map) {
  if (is.null(dim(C_report)) && is.numeric(C_report)) {
    C_report <- matrix(C_report, nrow = 1L)
  }
  C_report <- .hl_fo_assert_numeric_matrix(C_report, "C_report")
  beta_map <- .hl_fo_assert_numeric_matrix(beta_map, "beta_map")
  if (ncol(C_report) != nrow(beta_map)) {
    stop("C_report and beta_map dimensions do not agree.", call. = FALSE)
  }
  if (!is.numeric(rhs) || length(rhs) != nrow(C_report) ||
      any(!is.finite(rhs))) {
    stop("rhs must be finite with length nrow(C_report).", call. = FALSE)
  }
  list(C = C_report %*% beta_map, rhs = as.numeric(rhs))
}

.hl_fo_active_names <- function(active) {
  canonical <- c("tau", "xi")
  if (is.null(active)) {
    active <- canonical
  }
  if (is.logical(active)) {
    if (length(active) != 2L || anyNA(active)) {
      stop("Logical active must have length two and no missing values.",
           call. = FALSE)
    }
    active <- canonical[active]
  } else if (is.numeric(active)) {
    if (any(!is.finite(active)) || any(active != as.integer(active)) ||
        any(!active %in% 1:2)) {
      stop("Numeric active must contain only indices 1 and/or 2.",
           call. = FALSE)
    }
    active <- canonical[as.integer(active)]
  } else if (is.character(active)) {
    if (any(!active %in% canonical)) {
      stop("Character active must contain only 'tau' and/or 'xi'.",
           call. = FALSE)
    }
  } else {
    stop("active must be NULL, logical, numeric, or character.",
         call. = FALSE)
  }
  if (anyDuplicated(active)) {
    stop("active coordinates must not be duplicated.", call. = FALSE)
  }
  canonical[canonical %in% active]
}

hl_fo_active_coordinates <- function(active = c("tau", "xi")) {
  names <- .hl_fo_active_names(active)
  indices <- match(names, c("tau", "xi"))
  selector <- matrix(0, nrow = length(indices), ncol = 2L)
  if (length(indices)) {
    selector[cbind(seq_along(indices), indices)] <- 1
  }
  list(names = names, indices = indices, q = length(indices), P_A = selector)
}

hl_fo_select_active_rho <- function(tau, xi, active = c("tau", "xi")) {
  tau <- .hl_fo_assert_scalar(tau, "tau")
  xi <- .hl_fo_assert_scalar(xi, "xi")
  spec <- hl_fo_active_coordinates(active)
  out <- c(tau = tau, xi = xi)[spec$indices]
  names(out) <- spec$names
  out
}

hl_fo_active_insertion <- function(p, active = c("tau", "xi")) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) ||
      p != as.integer(p) || p < 0L) {
    stop("p must be one nonnegative integer.", call. = FALSE)
  }
  p <- as.integer(p)
  spec <- hl_fo_active_coordinates(active)
  .hl_fo_blockdiag(diag(p), t(spec$P_A))
}

hl_fo_expand_active_covariance <- function(covariance, p,
                                           active = c("tau", "xi")) {
  insertion <- hl_fo_active_insertion(p, active)
  covariance <- .hl_fo_assert_numeric_matrix(
    covariance, "covariance", square = TRUE
  )
  if (nrow(covariance) != ncol(insertion)) {
    stop("covariance dimension does not match p plus active coordinates.",
         call. = FALSE)
  }
  list(
    covariance = hl_fo_cov_congruence(covariance, insertion),
    insertion = insertion,
    active = hl_fo_active_coordinates(active)
  )
}

hl_fo_make_beta_seed <- function(beta,
                                 coordinate = c("canonical_conditional",
                                                "mean_one"),
                                 sigma2 = NULL, phi = NULL, source = NULL) {
  beta <- .hl_fo_assert_beta(beta)
  coordinate <- match.arg(coordinate)
  if (!is.null(sigma2)) {
    sigma2 <- .hl_fo_assert_scalar(sigma2, "sigma2", positive = TRUE)
  }
  if (!is.null(phi)) {
    phi <- .hl_fo_assert_scalar(phi, "phi", positive = TRUE)
  }
  structure(
    list(
      beta = beta,
      beta_coordinate = coordinate,
      sigma2 = sigma2,
      phi = phi,
      source = source
    ),
    class = c("hl_fo_beta_seed", "list")
  )
}

hl_fo_canonicalize_beta_seed <- function(seed, intercept_idx = NULL,
                                         expected_p = NULL) {
  if (!is.list(seed) || is.null(seed$beta_coordinate) ||
      !is.character(seed$beta_coordinate) ||
      length(seed$beta_coordinate) != 1L ||
      !seed$beta_coordinate %in% c("canonical_conditional", "mean_one")) {
    stop("Untagged beta seed: an explicit supported beta_coordinate is required.",
         call. = FALSE)
  }
  beta <- .hl_fo_assert_beta(seed$beta, "seed$beta")
  if (!is.null(expected_p)) {
    if (!is.numeric(expected_p) || length(expected_p) != 1L ||
        !is.finite(expected_p) || expected_p != as.integer(expected_p) ||
        expected_p < 0L) {
      stop("expected_p must be one nonnegative integer.", call. = FALSE)
    }
    if (length(beta) != expected_p) {
      stop("Seed beta dimension does not match expected_p.", call. = FALSE)
    }
  }

  sigma2 <- seed$sigma2
  phi <- seed$phi
  if (!is.null(sigma2)) {
    sigma2 <- .hl_fo_assert_scalar(sigma2, "seed$sigma2", positive = TRUE)
  }
  if (!is.null(phi)) {
    phi <- .hl_fo_assert_scalar(phi, "seed$phi", positive = TRUE)
  }

  if (identical(seed$beta_coordinate, "mean_one")) {
    if (is.null(sigma2)) {
      stop("A mean-one beta seed requires sigma2.", call. = FALSE)
    }
    intercept_idx <- .hl_fo_intercept_index(intercept_idx, length(beta))
    beta[intercept_idx] <- beta[intercept_idx] - sigma2 / 2
  }

  list(
    beta = beta,
    beta_coordinate = "canonical_conditional",
    tau = if (is.null(sigma2)) NULL else log(sigma2),
    xi = if (is.null(phi)) NULL else log(phi),
    sigma2 = sigma2,
    phi = phi,
    source = seed$source
  )
}

hl_fo_robust_wald <- function(beta, covariance, R, rhs = NULL,
                              n_subjects = NA_integer_) {
  beta <- .hl_fo_assert_beta(beta)
  covariance <- .hl_fo_assert_numeric_matrix(
    covariance, "covariance", square = TRUE
  )
  if (nrow(covariance) != length(beta)) {
    stop("covariance dimension must equal length(beta).", call. = FALSE)
  }
  if (is.null(dim(R)) && is.numeric(R)) {
    R <- matrix(R, nrow = 1L)
  }
  R <- .hl_fo_assert_numeric_matrix(R, "R")
  if (ncol(R) != length(beta)) {
    stop("ncol(R) must equal length(beta).", call. = FALSE)
  }
  if (is.null(rhs)) {
    rhs <- numeric(nrow(R))
  }
  reduced <- hl_fo_reduce_linear_rows(R, rhs)
  beta_cov_rank <- hl_fo_num_rank(covariance)$rank

  base <- list(
    statistic = NA_real_,
    p_value = NA_real_,
    df = reduced$rank,
    n_subjects = n_subjects,
    contrast_rank = reduced$rank,
    covariance_rank = NA_integer_,
    beta_covariance_rank = beta_cov_rank,
    calibration = "asymptotic_chisq",
    contrast = reduced$matrix,
    rhs = reduced$rhs,
    selected_rows = reduced$selected_rows,
    success = FALSE,
    ridge = 0
  )
  if (!reduced$success) {
    return(c(base, list(status = "wald_contrast_inconsistent")))
  }
  if (reduced$rank == 0L) {
    return(c(base, list(status = "wald_contrast_rank_zero")))
  }

  symmetry <- hl_fo_symmetry_defect(covariance)
  if (symmetry > 1e-8) {
    return(c(base, list(
      status = "wald_covariance_invalid",
      covariance_symmetry_defect = symmetry
    )))
  }
  covariance_symmetric <- (covariance + t(covariance)) / 2
  covariance_eigenvalues <- eigen(
    covariance_symmetric, symmetric = TRUE, only.values = TRUE
  )$values
  covariance_psd_tolerance <- 1e-8 * max(
    1, if (length(covariance_eigenvalues)) max(covariance_eigenvalues) else 0
  )
  base$beta_covariance_eigenvalues <- covariance_eigenvalues
  base$beta_covariance_psd_tolerance <- covariance_psd_tolerance
  if (any(!is.finite(covariance_eigenvalues)) ||
      min(covariance_eigenvalues) < -covariance_psd_tolerance) {
    return(c(base, list(status = "wald_covariance_invalid")))
  }
  error <- as.vector(reduced$matrix %*% beta - reduced$rhs)
  contrast_covariance <-
    reduced$matrix %*% covariance_symmetric %*% t(reduced$matrix)
  contrast_covariance <- (contrast_covariance + t(contrast_covariance)) / 2
  covariance_rank <- hl_fo_num_rank(contrast_covariance)$rank
  base$covariance_rank <- covariance_rank
  base$error <- error
  base$contrast_covariance <- contrast_covariance
  base$covariance_symmetry_defect <- symmetry
  if (covariance_rank != reduced$rank) {
    return(c(base, list(status = "wald_covariance_rank_deficient")))
  }

  eigenvalues <- eigen(contrast_covariance, symmetric = TRUE,
                       only.values = TRUE)$values
  if (any(!is.finite(eigenvalues)) || min(eigenvalues) <= 0) {
    return(c(base, list(
      status = "wald_covariance_invalid",
      contrast_covariance_eigenvalues = eigenvalues
    )))
  }
  solved <- tryCatch(
    solve(contrast_covariance, error),
    error = function(e) NULL
  )
  if (is.null(solved)) {
    return(c(base, list(status = "wald_covariance_rank_deficient")))
  }
  statistic <- as.numeric(crossprod(error, solved))
  if (!is.finite(statistic) || statistic < -1e-12) {
    return(c(base, list(status = "wald_statistic_invalid")))
  }
  statistic <- max(0, statistic)
  base$statistic <- statistic
  base$p_value <- stats::pchisq(statistic, df = reduced$rank,
                               lower.tail = FALSE)
  base$success <- TRUE
  c(base, list(
    status = "ok",
    contrast_covariance_eigenvalues = eigenvalues
  ))
}
