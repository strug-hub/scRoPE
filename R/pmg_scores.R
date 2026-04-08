pmg_per_subject_gradients <- function(para,
                                      X,
                                      offset,
                                      Y,
                                      fid,
                                      cumsumy,
                                      posind,
                                      posindy,
                                      nb,
                                      nind,
                                      k) {
  beta <- para[seq_len(nb)]
  sigma <- para[nb + 1L]

  exps <- exp(sigma)
  alpha <- 1 / (exps - 1)
  lambda <- 1 / (sqrt(exps) * (exps - 1))
  alpha_pr <- -exps / ((exps - 1)^2)
  lambda_pr <- (1 - 3 * exps) / (2 * sqrt(exps) * (exps - 1)^2)

  eta <- as.numeric(offset + X %*% beta)
  mu <- exp(eta)
  subject_index <- rep.int(seq_len(k), diff(fid))

  mu_sum <- as.numeric(rowsum(mu, subject_index, reorder = FALSE))
  xmu_sum <- rowsum(sweep(X, 1L, mu, `*`), subject_index, reorder = FALSE)

  yx_sum <- matrix(0, nrow = k, ncol = nb)
  if (length(posindy) > 0L) {
    pos_idx <- as.integer(posindy) + 1L
    subject_pos <- subject_index[pos_idx]
    yx_pos <- sweep(X[pos_idx, , drop = FALSE], 1L, as.numeric(Y), `*`)
    yx_tmp <- rowsum(yx_pos, subject_pos, reorder = FALSE)
    yx_sum[as.integer(rownames(yx_tmp)), ] <- as.matrix(yx_tmp)
  }

  ystar <- as.numeric(cumsumy) + alpha
  mustar <- mu_sum + lambda

  beta_grad <- sweep(xmu_sum, 1L, ystar / mustar, `*`) - yx_sum

  sigma_score <- alpha_pr * (
    log(lambda) +
      Digamma(as.numeric(cumsumy) + alpha) -
      Digamma(alpha) -
      log(mustar)
  ) + lambda_pr * (alpha / lambda - ystar / mustar)
  sigma_grad <- -sigma_score

  out <- rbind(t(beta_grad), sigma_grad)
  rownames(out) <- NULL
  colnames(out) <- NULL
  out
}

pmg_ll_der_hes <- function(para,
                           X,
                           offset,
                           Y,
                           fid,
                           cumsumy,
                           posind,
                           posindy,
                           nb,
                           nind,
                           k) {
  objective <- pmg_ll(
    para = para,
    X = X,
    offset = offset,
    Y = Y,
    fid = fid,
    cumsumy = cumsumy,
    posind = posind,
    posindy = posindy,
    nb = nb,
    nind = nind,
    k = k
  )
  gradient <- pmg_der(
    para = para,
    X = X,
    offset = offset,
    Y = Y,
    fid = fid,
    cumsumy = cumsumy,
    posind = posind,
    posindy = posindy,
    nb = nb,
    nind = nind,
    k = k
  )
  hessian <- pmg_hes(
    para = para,
    X = X,
    offset = offset,
    Y = Y,
    fid = fid,
    cumsumy = cumsumy,
    posind = posind,
    posindy = posindy,
    nb = nb,
    nind = nind,
    k = k
  )
  per_subject_gradients <- pmg_per_subject_gradients(
    para = para,
    X = X,
    offset = offset,
    Y = Y,
    fid = fid,
    cumsumy = cumsumy,
    posind = posind,
    posindy = posindy,
    nb = nb,
    nind = nind,
    k = k
  )

  list(
    objective = objective,
    value = objective,
    loglik = -objective,
    gradient = gradient,
    score = gradient,
    hessian = hessian,
    observed_info = -hessian,
    per_subject_gradients = per_subject_gradients,
    per_subject_scores = -per_subject_gradients
  )
}
