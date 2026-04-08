ptmg_ll <- function(para, X, offset, Y, n_one, ytwo, fam, fid, cumsumy, posind, posindy, nb, nind, k) {
  beta <- para[1:nb]
  sigma <- para[(nb + 1):(nb + 2)]

  exps <- exp(sigma[1])
  alpha <- 1 / (exps - 1)
  gamma <- sigma[2]

  fn <- ptmg_ll_eigen(X, offset, Y, fam - 1, fid - 1, as.double(cumsumy), posind - 1, posindy, nind, k, beta, sigma)
  fn <- fn + (sum(Lgamma(cumsumy[posind] + alpha)) - length(posind) * Lgamma(alpha))
  fn <- fn + (sum(Lgamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Lgamma(gamma) + sum(n_one) * log(gamma) + n_one[2] * log(gamma + 1))
  return(-fn)
}

ptmg_der <- function(para, X, offset, Y, n_one, ytwo, fam, fid, cumsumy, posind, posindy, nb, nind, k) {
  beta <- para[1:nb]
  sigma <- para[(nb + 1):(nb + 2)]

  exps <- exp(sigma[1])
  alpha <- 1 / (exps - 1)
  gamma <- sigma[2]
  alpha_pr <- -exps / ((exps - 1)^2)

  temp <- ptmg_der_eigen(X, offset, Y, fam - 1, fid - 1, as.double(cumsumy), posind - 1, posindy, nb, nind, k, beta, sigma)
  temp[nb + 1] <- temp[nb + 1] - alpha_pr * (sum(Digamma(cumsumy[posind] + alpha)) - length(posind) * Digamma(alpha))
  temp[nb + 2] <- temp[nb + 2] - (sum(Digamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Digamma(gamma) + sum(n_one) / gamma + n_one[2] / (gamma + 1))
  return(temp)
}

ptmg_ll_der <- function(para, X, offset, Y, n_one, ytwo, fid, cumsumy, posind, posindy, nb, nind, k) {
  beta <- para[1:nb]
  sigma <- para[(nb + 1):(nb + 2)]

  exps <- exp(sigma[1])
  alpha <- 1 / (exps - 1)
  gamma <- sigma[2]
  alpha_pr <- -exps / ((exps - 1)^2)

  temp <- ptmg_ll_der_eigen(X, offset, Y, fid - 1, as.double(cumsumy), posind - 1, posindy, nb, nind, k, beta, sigma)
  temp$fn <- temp$fn - (sum(Lgamma(cumsumy[posind] + alpha)) - length(posind) * Lgamma(alpha))
  temp$fn <- temp$fn - (sum(Lgamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Lgamma(gamma) + sum(n_one) * log(gamma) + n_one[2] * log(gamma + 1))
  gr <- temp$gr
  gr[nb + 1] <- gr[nb + 1] - alpha_pr * (sum(Digamma(cumsumy[posind] + alpha)) - length(posind) * Digamma(alpha))
  gr[nb + 2] <- gr[nb + 2] - (sum(Digamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Digamma(gamma) + sum(n_one) / gamma + n_one[2] / (gamma + 1))
  return(list("objective" = temp$fn, "gradient" = gr))
}

ptmg_ll_der2 <- function(para, X, offset, Y, n_one, ytwo, fid, cumsumy, posind, posindy, nb, nind, k) {
  beta <- para[1:nb]
  sigma <- exp(para[(nb + 1):(nb + 2)])

  exps <- exp(sigma[1])
  alpha <- 1 / (exps - 1)
  gamma <- sigma[2]
  alpha_pr <- -exps / ((exps - 1)^2)

  temp <- ptmg_ll_der_eigen(X, offset, Y, fid - 1, as.double(cumsumy), posind - 1, posindy, nb, nind, k, beta, sigma)
  temp$fn <- temp$fn - (sum(Lgamma(cumsumy[posind] + alpha)) - length(posind) * Lgamma(alpha))
  temp$fn <- temp$fn - (sum(Lgamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Lgamma(gamma) + sum(n_one) * log(gamma) + n_one[2] * log(gamma + 1))
  gr <- temp$gr
  gr[nb + 1] <- gr[nb + 1] - alpha_pr * (sum(Digamma(cumsumy[posind] + alpha)) - length(posind) * Digamma(alpha))
  gr[nb + 2] <- gr[nb + 2] - (sum(Digamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Digamma(gamma) + sum(n_one) / gamma + n_one[2] / (gamma + 1))
  gr[nb + 1] <- gr[nb + 1] * sigma[1]
  gr[nb + 2] <- gr[nb + 2] * sigma[2]
  # return(list("objective"=temp$fn,"gradient"=gr))
  objective <- temp$fn
  attr(objective, "gradient") <- gr
  objective
}


ptmg_ll_der_hes2 <- function(para, X, offset, Y, n_one, ytwo, fid, cumsumy, posind, posindy, nb, nind, k) {
  beta <- para[1:nb]
  sigma <- exp(para[(nb + 1):(nb + 2)])

  exps <- exp(sigma[1])
  exps_m <- (exps - 1)^2
  alpha <- 1 / (exps - 1)
  gamma <- sigma[2]
  alpha_pr <- -exps / ((exps - 1)^2)
  alpha_dpr <- 2 * exps * exps / (exps_m * (exps - 1)) - exps / exps_m

  temp <- ptmg_ll_der_hes_eigen(X, offset, Y, fid - 1, as.double(cumsumy), posind - 1, posindy, nb, nind, k, beta, sigma)
  temp$fn <- temp$fn - (sum(Lgamma(cumsumy[posind] + alpha)) - length(posind) * Lgamma(alpha))
  temp$fn <- temp$fn - (sum(Lgamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Lgamma(gamma) + sum(n_one) * log(gamma) + n_one[2] * log(gamma + 1))
  gr <- temp$gr
  grdig <- sum(Digamma(cumsumy[posind] + alpha)) - length(posind) * Digamma(alpha)
  gr[nb + 1] <- gr[nb + 1] - alpha_pr * grdig
  gr[nb + 2] <- gr[nb + 2] - (sum(Digamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Digamma(gamma) + sum(n_one) / gamma + n_one[2] / (gamma + 1))
  gr[nb + 1] <- gr[nb + 1] * sigma[1]
  gr[nb + 2] <- gr[nb + 2] * sigma[2]
  hes <- temp$hes
  hes[nb + 2, nb + 2] <- hes[nb + 2, nb + 2] - (sum(Trigamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Trigamma(gamma) - sum(n_one) / gamma / gamma - n_one[2] / (gamma + 1) / (gamma + 1))
  hes[nb + 1, nb + 1] <- hes[nb + 1, nb + 1] - alpha_dpr * grdig - alpha_pr * alpha_pr * (sum(Trigamma(cumsumy[posind] + alpha)) - length(posind) * Trigamma(alpha))

  hes[nb + 2, nb + 2] <- hes[nb + 2, nb + 2] * sigma[2] * sigma[2] + gr[nb + 2]
  hes[nb + 1, nb + 1] <- hes[nb + 1, nb + 1] * sigma[1] * sigma[1] + gr[nb + 1]
  hes[nb + 1, nb + 2] <- hes[nb + 2, nb + 1] <- hes[nb + 1, nb + 2] * sigma[2] * sigma[1]
  hes[nb + 1, 1:nb] <- hes[1:nb, nb + 1] <- hes[nb + 1, 1:nb] * sigma[1]
  hes[nb + 2, 1:nb] <- hes[1:nb, nb + 2] <- hes[nb + 2, 1:nb] * sigma[2]
  objective <- temp$fn
  attr(objective, "gradient") <- gr
  attr(objective, "hessian") <- hes
  objective
  # return(list("objective"=-objective,"gradient"=-gr, "hessian"=-hes))
}

ptmg_ll_der_hes3 <- function(para, X, offset, Y, n_one, ytwo, fid, cumsumy, posind, posindy, nb, nind, k) {
  beta <- para[1:nb]
  sigma <- exp(para[(nb + 1):(nb + 2)])

  exps <- exp(sigma[1])
  exps_m <- (exps - 1)^2
  alpha <- 1 / (exps - 1)
  gamma <- sigma[2]
  alpha_pr <- -exps / ((exps - 1)^2)
  alpha_dpr <- 2 * exps * exps / (exps_m * (exps - 1)) - exps / exps_m

  temp <- ptmg_ll_der_hes_eigen(X, offset, Y, fid - 1, as.double(cumsumy), posind - 1, posindy, nb, nind, k, beta, sigma)
  temp$fn <- temp$fn - (sum(Lgamma(cumsumy[posind] + alpha)) - length(posind) * Lgamma(alpha))
  temp$fn <- temp$fn - (sum(Lgamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Lgamma(gamma) + sum(n_one) * log(gamma) + n_one[2] * log(gamma + 1))
  gr <- temp$gr
  grdig <- sum(Digamma(cumsumy[posind] + alpha)) - length(posind) * Digamma(alpha)
  gr[nb + 1] <- gr[nb + 1] - alpha_pr * grdig
  gr[nb + 2] <- gr[nb + 2] - (sum(Digamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Digamma(gamma) + sum(n_one) / gamma + n_one[2] / (gamma + 1))
  gr[nb + 1] <- gr[nb + 1] * sigma[1]
  gr[nb + 2] <- gr[nb + 2] * sigma[2]
  hes <- temp$hes
  hes[nb + 2, nb + 2] <- hes[nb + 2, nb + 2] - (sum(Trigamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Trigamma(gamma) - sum(n_one) / gamma / gamma - n_one[2] / (gamma + 1) / (gamma + 1))
  hes[nb + 1, nb + 1] <- hes[nb + 1, nb + 1] - alpha_dpr * grdig - alpha_pr * alpha_pr * (sum(Trigamma(cumsumy[posind] + alpha)) - length(posind) * Trigamma(alpha))

  hes[nb + 2, nb + 2] <- hes[nb + 2, nb + 2] * sigma[2] * sigma[2] + gr[nb + 2]
  hes[nb + 1, nb + 1] <- hes[nb + 1, nb + 1] * sigma[1] * sigma[1] + gr[nb + 1]
  hes[nb + 1, nb + 2] <- hes[nb + 2, nb + 1] <- hes[nb + 1, nb + 2] * sigma[2] * sigma[1]
  hes[nb + 1, 1:nb] <- hes[1:nb, nb + 1] <- hes[nb + 1, 1:nb] * sigma[1]
  hes[nb + 2, 1:nb] <- hes[1:nb, nb + 2] <- hes[nb + 2, 1:nb] * sigma[2]
  objective <- temp$fn
  # attr(objective, "gradient") = gr
  # attr(objective, "hessian") = hes
  # objective
  return(list("value" = objective, "gradient" = gr, "hessian" = hes))
}

ptmg_ll_der_hes4 <- function(para, X, offset, Y, n_one, ytwo, fid, cumsumy, posind, posindy, nb, nind, k) {
  beta <- para[1:nb]
  sigma <- exp(para[(nb + 1):(nb + 2)])

  exps <- exp(sigma[1])
  exps_m <- (exps - 1)^2
  alpha <- 1 / (exps - 1)
  gamma <- sigma[2]
  alpha_pr <- -exps / ((exps - 1)^2)
  alpha_dpr <- 2 * exps * exps / (exps_m * (exps - 1)) - exps / exps_m

  temp <- ptmg_ll_der_hes_eigen_per_subject(X, offset, Y, fid - 1, as.double(cumsumy), posind - 1, posindy, nb, nind, k, beta, sigma)
  temp$fn <- temp$fn - (sum(Lgamma(cumsumy[posind] + alpha)) - length(posind) * Lgamma(alpha))
  temp$fn <- temp$fn - (sum(Lgamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Lgamma(gamma) + sum(n_one) * log(gamma) + n_one[2] * log(gamma + 1))
  gr <- temp$gr
  grdig <- sum(Digamma(cumsumy[posind] + alpha)) - length(posind) * Digamma(alpha)
  gr[nb + 1] <- gr[nb + 1] - alpha_pr * grdig
  gr[nb + 2] <- gr[nb + 2] - (sum(Digamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Digamma(gamma) + sum(n_one) / gamma + n_one[2] / (gamma + 1))
  gr[nb + 1] <- gr[nb + 1] * sigma[1]
  gr[nb + 2] <- gr[nb + 2] * sigma[2]
  hes <- temp$hes
  hes[nb + 2, nb + 2] <- hes[nb + 2, nb + 2] - (sum(Trigamma(ytwo + gamma)) - (length(posindy) - sum(n_one)) * Trigamma(gamma) - sum(n_one) / gamma / gamma - n_one[2] / (gamma + 1) / (gamma + 1))
  hes[nb + 1, nb + 1] <- hes[nb + 1, nb + 1] - alpha_dpr * grdig - alpha_pr * alpha_pr * (sum(Trigamma(cumsumy[posind] + alpha)) - length(posind) * Trigamma(alpha))

  hes[nb + 2, nb + 2] <- hes[nb + 2, nb + 2] * sigma[2] * sigma[2] + gr[nb + 2]
  hes[nb + 1, nb + 1] <- hes[nb + 1, nb + 1] * sigma[1] * sigma[1] + gr[nb + 1]
  hes[nb + 1, nb + 2] <- hes[nb + 2, nb + 1] <- hes[nb + 1, nb + 2] * sigma[2] * sigma[1]
  hes[nb + 1, 1:nb] <- hes[1:nb, nb + 1] <- hes[nb + 1, 1:nb] * sigma[1]
  hes[nb + 2, 1:nb] <- hes[1:nb, nb + 2] <- hes[nb + 2, 1:nb] * sigma[2]
  objective <- temp$fn

  # Extract per-subject gradients (without Digamma terms)
  per_subject_gradients <- temp$per_subject_gradients

  # Adjust per-subject gradients for sigma[1] (sigma2)
  digamma_alpha <- Digamma(alpha)
  digamma_cumsumy_alpha <- Digamma(temp$cumsumy_alpha)

  grdig_alpha <- digamma_cumsumy_alpha - digamma_alpha # Vector of length k

  per_subject_gradients[nb + 1, ] <- per_subject_gradients[nb + 1, ] - alpha_pr * grdig_alpha
  per_subject_gradients[nb + 1, ] <- per_subject_gradients[nb + 1, ] * sigma[1] # Multiply by sigma[1]

  # Adjust per-subject gradients for gamma (sigma[2])
  digamma_gamma <- Digamma(gamma)
  digamma_Y_plus_gamma <- Digamma(temp$Y_plus_gamma)

  # Build the mapping from observations to subjects
  obs_to_subject <- rep(1:k, times = diff(fid))

  # Map positive observations to subjects
  positive_obs_indices <- posindy # Indices of positive observations
  subject_indices <- obs_to_subject[positive_obs_indices]

  # Compute n_i (number of positive observations per subject)
  n_i <- tabulate(subject_indices, nbins = k) # Length k

  # Compute sum of digamma_Y_plus_gamma per subject
  sum_digamma_Y_plus_gamma_full <- numeric(k)


  for (i in seq_along(digamma_Y_plus_gamma)) {
    subj <- subject_indices[i]
    sum_digamma_Y_plus_gamma_full[subj] <- sum_digamma_Y_plus_gamma_full[subj] + digamma_Y_plus_gamma[i]
  }

  # Compute grdig_gamma_full
  grdig_gamma_full <- sum_digamma_Y_plus_gamma_full - n_i * digamma_gamma

  # Adjust per-subject gradients for gamma
  per_subject_gradients[nb + 2, ] <- per_subject_gradients[nb + 2, ] - grdig_gamma_full
  per_subject_gradients[nb + 2, ] <- per_subject_gradients[nb + 2, ] * gamma # Adjust for log-transformation

  # Return the objective function value, gradients, and Hessian
  return(list("value" = objective, "gradient" = gr, "hessian" = hes, "per_subject_gradients" = per_subject_gradients))
}

compute_sandwich_variance <- function(out_at_op, nb, return_type = c("full", "beta", "beta_se"), invert_function = NULL) {
  # Match the return_type argument
  return_type <- match.arg(return_type)
  
  # Extract the negative Hessian (observed information matrix)
  H <- if (!is.null(out_at_op$observed_info)) out_at_op$observed_info else out_at_op$hessian
  
  # Extract per-subject gradients
  per_subject_gradients <- out_at_op$per_subject_gradients  # (nb + 2) x k matrix
  
  # Compute S: Variance of per-subject gradients
  S <- per_subject_gradients %*% t(per_subject_gradients)  # Sum over subjects
  
  # Invert the observed information matrix H
  if (is.null(invert_function)) {
    # Use Rfast::spdinv by default
    H_inv <- tryCatch({
      Rfast::spdinv(H)
    }, error = function(e) {
      warning("Rfast::spdinv failed, using solve instead")
      solve(H)
    })
  } else {
    # Use the user-provided inversion function
    H_inv <- invert_function(H)
  }
  
  # Compute the sandwich variance-covariance matrix
  Var <- H_inv %*% S %*% H_inv
  
  # Return based on the specified return_type
  if (return_type == "full") {
    return(Var)
  } else if (return_type == "beta") {
    Var_beta <- Var[1:nb, 1:nb]
    return(Var_beta)
  } else if (return_type == "beta_se") {
    Var_beta <- Var[1:nb, 1:nb]
    SE_beta <- sqrt(diag(Var_beta))
    return(SE_beta)
  } else {
    stop("Invalid return_type specified.")
  }
}

compute_sandwich_variance2 <- function(out_at_op, nb, compute_full = FALSE, invert_function = NULL) {
  # Extract the negative Hessian (observed information matrix)
  H <- if (!is.null(out_at_op$observed_info)) out_at_op$observed_info else out_at_op$hessian

  # Extract per-subject gradients
  U <- out_at_op$per_subject_gradients  # (nb + 2) x k matrix

  # Compute the covariance matrix of per-subject gradients
  S <- U %*% t(U)  # (nb + 2) x (nb + 2)

  # Compute the Godambe information matrix G = H %*% S^{-1} %*% H
  if (is.null(invert_function)) {
    # Invert S to get S_inv
    S_inv <- tryCatch({
      Rfast::spdinv(S)
    }, error = function(e) {
      warning("Rfast::spdinv failed on S, using solve instead")
      solve(S)
    })
  } else {
    S_inv <- invert_function(S)
  }

  G <- H %*% S_inv %*% H

  q <- nrow(G) - nb
  if (q < 0) {
    stop("nb exceeds the parameter dimension in compute_sandwich_variance2().")
  }

  if (q == 0L) {
    G_beta <- G[seq_len(nb), seq_len(nb), drop = FALSE]
  } else {
    idx_beta <- seq_len(nb)
    idx_lambda <- nb + seq_len(q)
    G_bb <- G[idx_beta, idx_beta, drop = FALSE]
    G_bl <- G[idx_beta, idx_lambda, drop = FALSE]
    G_lb <- G[idx_lambda, idx_beta, drop = FALSE]
    G_ll <- G[idx_lambda, idx_lambda, drop = FALSE]

    # Compute the partial Godambe information matrix for beta
    # G_beta = G_bb - G_bl G_ll^{-1} G_lb
    if (is.null(invert_function)) {
      G_ll_inv <- tryCatch({
        Rfast::spdinv(G_ll)
      }, error = function(e) {
        warning("Rfast::spdinv failed on G_ll, using solve instead")
        solve(G_ll)
      })
    } else {
      G_ll_inv <- invert_function(G_ll)
    }
    G_beta <- G_bb - G_bl %*% G_ll_inv %*% G_lb
  }

  # Invert G_beta to get Var_beta_adjusted
  if (is.null(invert_function)) {
    Var_beta_adjusted <- tryCatch({
      Rfast::spdinv(G_beta)
    }, error = function(e) {
      warning("Rfast::spdinv failed on G_beta, using solve instead")
      solve(G_beta)
    })
  } else {
    Var_beta_adjusted <- invert_function(G_beta)
  }

  result <- list(Var_beta_adjusted = Var_beta_adjusted)

  if (compute_full) {
    # Compute the sandwich variance-covariance matrix Var_full
    if (is.null(invert_function)) {
      H_inv <- tryCatch({
        Rfast::spdinv(H)
      }, error = function(e) {
        warning("Rfast::spdinv failed on H, using solve instead")
        solve(H)
      })
    } else {
      H_inv <- invert_function(H)
    }
    Var_full <- H_inv %*% S %*% H_inv
    result$Var_full <- Var_full
  }

  return(result)
}

extract_godambe_components <- function(out_at_op, nb) {
  H <- if (!is.null(out_at_op$observed_info)) out_at_op$observed_info else out_at_op$hessian
  U <- out_at_op$per_subject_gradients

  if (is.null(H) || is.null(U)) {
    return(NULL)
  }

  if (!is.null(out_at_op$S_override)) {
    S <- out_at_op$S_override
  } else {
    S <- U %*% t(U)
  }

  idx_beta <- seq_len(nb)
  idx_lambda <- nb + seq_len(nrow(H) - nb)

  list(
    H = H,
    S = S,
    H_beta = H[idx_beta, idx_beta, drop = FALSE],
    H_beta_lambda = H[idx_beta, idx_lambda, drop = FALSE],
    H_lambda = H[idx_lambda, idx_lambda, drop = FALSE],
    per_subject_gradients = U,
    score = out_at_op$gradient
  )
}
