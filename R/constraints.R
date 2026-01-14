constraint_reparameterization <- function(L, b) {
  if (is.null(dim(L))) {
    L <- matrix(L, nrow = 1)
  }
  if (!is.matrix(L)) {
    stop("L must be a matrix or a vector.")
  }
  r <- nrow(L)
  p <- ncol(L)
  if (length(b) != r) {
    stop("Length of b must match number of rows of L.")
  }
  qrL <- qr(L)
  if (qrL$rank != r) {
    stop("L must have full row rank.")
  }
  LLt <- L %*% t(L)
  beta_star <- as.vector(t(L) %*% solve(LLt, b))
  qrt <- qr(t(L))
  Q <- qr.Q(qrt, complete = TRUE)
  if (r < p) {
    K <- Q[, (r + 1):p, drop = FALSE]
  } else {
    K <- matrix(0, nrow = p, ncol = 0)
  }
  list(beta_star = beta_star, null_basis = K)
}
