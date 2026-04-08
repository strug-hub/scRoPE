pmg_initial_beta <- function(posv, ctx) {
  start <- rep(0, ctx$nb)
  start[ctx$intcol] <- log(posv$mct) - ctx$moffset
  start
}

pmg_try_invert <- function(mat) {
  if (length(mat) == 0L) {
    return(matrix(0, nrow = 0, ncol = 0))
  }

  tryCatch(
    {
      solve(mat)
    },
    error = function(e) {
      tryCatch(
        {
          Rfast::spdinv(mat)
        },
        error = function(e2) NULL
      )
    }
  )
}

pmg_reduced_mapping <- function(nb, K) {
  free_dim <- ncol(K)
  map <- matrix(0, nrow = nb + 1L, ncol = free_dim + 1L)
  if (free_dim > 0L) {
    map[seq_len(nb), seq_len(free_dim)] <- K
  }
  map[nb + 1L, free_dim + 1L] <- 1
  map
}

fit_gene_pmg_unconstrained <- function(gene_index, posv, ctx) {
  with(ctx, {
    start <- c(pmg_initial_beta(posv, ctx), 1)
    lower <- c(rep(-100, nb), min[1])
    upper <- c(rep(100, nb), max[1])

    fit <- tryCatch(
      {
        ref <- nlminb(
          start = start,
          objective = pmg_ll,
          gradient = pmg_der,
          posindy = posv$posindy,
          X = pred,
          offset = offset,
          Y = posv$Y,
          fid = fid,
          cumsumy = cumsumy[gene_index, ],
          posind = posind[[gene_index]],
          nb = nb,
          nind = nind,
          k = k,
          lower = lower,
          upper = upper
        )
        list(
          par = ref$par,
          optimizer = "nlminb",
          raw = ref
        )
      },
      error = function(e) {
        ref <- bobyqa(
          par = start,
          fn = pmg_ll,
          posindy = posv$posindy,
          X = pred,
          offset = offset,
          Y = posv$Y,
          fid = fid,
          cumsumy = cumsumy[gene_index, ],
          posind = posind[[gene_index]],
          nb = nb,
          nind = nind,
          k = k,
          lower = lower,
          upper = upper
        )
        list(
          par = ref$par,
          optimizer = "bobyqa",
          raw = ref
        )
      }
    )

    eval_at_opt <- pmg_ll_der_hes(
      para = fit$par,
      X = pred,
      offset = offset,
      Y = posv$Y,
      fid = fid,
      cumsumy = cumsumy[gene_index, ],
      posind = posind[[gene_index]],
      posindy = posv$posindy,
      nb = nb,
      nind = nind,
      k = k
    )

    npar <- nb + 1L
    naive_cov_full <- matrix(NA_real_, nrow = npar, ncol = npar)
    active <- is.finite(diag(eval_at_opt$hessian)) &
      (fit$par != lower) &
      (fit$par != upper)
    if (any(active)) {
      naive_active <- pmg_try_invert(-eval_at_opt$hessian[active, active, drop = FALSE])
      if (!is.null(naive_active)) {
        naive_cov_full[active, active] <- naive_active
      }
    }

    robust_out <- tryCatch(
      {
        compute_sandwich_variance2(
          eval_at_opt,
          nb = nb,
          compute_full = TRUE
        )
      },
      error = function(e) NULL
    )

    convergence <- 1L
    conv_raw <- fit$raw$convergence
    if (!is.null(conv_raw) && conv_raw != 0) {
      convergence <- -50L
    }
    if (all(is.na(naive_cov_full))) {
      convergence <- -25L
    }
    if (is.null(robust_out) && convergence > 0L) {
      convergence <- -27L
    }

    naive_cov_beta <- if (all(is.na(naive_cov_full))) {
      matrix(NA_real_, nrow = nb, ncol = nb)
    } else {
      naive_cov_full[seq_len(nb), seq_len(nb), drop = FALSE]
    }
    robust_cov_beta <- if (is.null(robust_out)) {
      matrix(NA_real_, nrow = nb, ncol = nb)
    } else {
      robust_out$Var_beta_adjusted
    }
    robust_cov_full <- if (is.null(robust_out)) {
      matrix(NA_real_, nrow = nb + 1L, ncol = nb + 1L)
    } else {
      robust_out$Var_full
    }

    godambe_input <- list(
      gradient = eval_at_opt$gradient,
      hessian = eval_at_opt$hessian,
      observed_info = eval_at_opt$observed_info,
      per_subject_gradients = eval_at_opt$per_subject_gradients
    )

    list(
      theta = list(
        beta = fit$par[seq_len(nb)],
        subVar = fit$par[nb + 1L]
      ),
      loglik = eval_at_opt$loglik,
      score = eval_at_opt$gradient,
      hessian = eval_at_opt$hessian,
      observed_info = eval_at_opt$observed_info,
      per_subject_gradients = eval_at_opt$per_subject_gradients,
      naive_cov = naive_cov_beta,
      naive_cov_full = naive_cov_full,
      robust_cov = robust_cov_beta,
      robust_cov_full = robust_cov_full,
      godambe = extract_godambe_components(godambe_input, nb),
      convergence = convergence,
      optimizer = fit$optimizer,
      raw_fit = fit$raw
    )
  })
}

fit_gene_pmg_constrained <- function(gene_index, posv, ctx, L, b, unconstrained_fit = NULL) {
  with(ctx, {
    if (is.null(dim(L))) {
      L <- matrix(L, nrow = 1)
    }
    if (!is.matrix(L)) {
      stop("L must be a matrix or a vector.")
    }
    if (ncol(L) != nb) {
      stop("L must have the same number of columns as the design matrix.")
    }
    if (nrow(L) == 0L) {
      return(fit_gene_pmg_unconstrained(gene_index, posv, ctx))
    }
    if (length(b) != nrow(L)) {
      stop("Length of b must match the number of rows of L.")
    }

    constr <- constraint_reparameterization(L, b)
    beta_star <- constr$beta_star
    K <- constr$null_basis
    if (is.null(K)) {
      K <- matrix(0, nrow = nb, ncol = 0)
    }
    free_dim <- ncol(K)

    compose_beta <- function(gamma) {
      if (free_dim == 0L) {
        beta_star
      } else {
        as.vector(beta_star + K %*% gamma)
      }
    }

    build_full_para <- function(theta_reduced) {
      if (free_dim == 0L) {
        gamma <- numeric(0)
        sigma <- theta_reduced
      } else {
        gamma <- theta_reduced[seq_len(free_dim)]
        sigma <- theta_reduced[free_dim + 1L]
      }
      c(compose_beta(gamma), sigma)
    }

    beta_seed <- pmg_initial_beta(posv, ctx)
    sigma_seed <- 1
    if (!is.null(unconstrained_fit) &&
      !is.null(unconstrained_fit$theta) &&
      !is.null(unconstrained_fit$theta$beta) &&
      !is.null(unconstrained_fit$theta$subVar)) {
      beta_seed <- as.numeric(unconstrained_fit$theta$beta)
      sigma_seed <- as.numeric(unconstrained_fit$theta$subVar)
    }
    sigma_seed <- min(max(sigma_seed, min[1]), max[1])
    gamma_init <- if (free_dim == 0L) {
      numeric(0)
    } else {
      as.numeric(crossprod(K, beta_seed - beta_star))
    }
    theta_init <- c(gamma_init, sigma_seed)
    lower_theta <- c(if (free_dim == 0L) numeric(0) else rep(-100, free_dim), min[1])
    upper_theta <- c(if (free_dim == 0L) numeric(0) else rep(100, free_dim), max[1])

    obj_fn <- function(theta_reduced) {
      pmg_ll(
        para = build_full_para(theta_reduced),
        X = pred,
        offset = offset,
        Y = posv$Y,
        fid = fid,
        cumsumy = cumsumy[gene_index, ],
        posind = posind[[gene_index]],
        posindy = posv$posindy,
        nb = nb,
        nind = nind,
        k = k
      )
    }

    grad_fn <- function(theta_reduced) {
      grad_full <- pmg_der(
        para = build_full_para(theta_reduced),
        X = pred,
        offset = offset,
        Y = posv$Y,
        fid = fid,
        cumsumy = cumsumy[gene_index, ],
        posind = posind[[gene_index]],
        posindy = posv$posindy,
        nb = nb,
        nind = nind,
        k = k
      )
      grad_beta <- grad_full[seq_len(nb)]
      grad_sigma <- grad_full[nb + 1L]
      c(
        if (free_dim == 0L) numeric(0) else as.vector(crossprod(K, grad_beta)),
        grad_sigma
      )
    }

    fit <- tryCatch(
      {
        ref <- nlminb(
          start = theta_init,
          objective = obj_fn,
          gradient = grad_fn,
          lower = lower_theta,
          upper = upper_theta
        )
        list(
          par = ref$par,
          optimizer = "nlminb",
          raw = ref
        )
      },
      error = function(e) {
        ref <- bobyqa(
          par = theta_init,
          fn = obj_fn,
          lower = lower_theta,
          upper = upper_theta
        )
        list(
          par = ref$par,
          optimizer = "bobyqa",
          raw = ref
        )
      }
    )

    full_par <- build_full_para(fit$par)
    eval_at_opt <- pmg_ll_der_hes(
      para = full_par,
      X = pred,
      offset = offset,
      Y = posv$Y,
      fid = fid,
      cumsumy = cumsumy[gene_index, ],
      posind = posind[[gene_index]],
      posindy = posv$posindy,
      nb = nb,
      nind = nind,
      k = k
    )

    map <- pmg_reduced_mapping(nb, K)
    H_red <- crossprod(map, eval_at_opt$observed_info %*% map)
    U_red <- crossprod(map, eval_at_opt$per_subject_gradients)
    score_red <- as.numeric(crossprod(map, eval_at_opt$gradient))

    active_red <- is.finite(diag(H_red)) &
      (fit$par != lower_theta) &
      (fit$par != upper_theta)

    npar_red <- length(fit$par)
    naive_cov_red <- matrix(NA_real_, nrow = npar_red, ncol = npar_red)
    robust_cov_red <- matrix(NA_real_, nrow = npar_red, ncol = npar_red)

    if (any(active_red)) {
      H_red_active <- H_red[active_red, active_red, drop = FALSE]
      H_red_inv <- pmg_try_invert(H_red_active)
      if (!is.null(H_red_inv)) {
        naive_cov_red[active_red, active_red] <- H_red_inv

        U_red_active <- U_red[active_red, , drop = FALSE]
        S_red_active <- U_red_active %*% t(U_red_active)
        robust_cov_active <- H_red_inv %*% S_red_active %*% H_red_inv
        robust_cov_red[active_red, active_red] <- robust_cov_active
      }
    }

    gamma_active <- if (free_dim == 0L) logical(0) else active_red[seq_len(free_dim)]

    if (free_dim == 0L) {
      naive_cov_beta <- matrix(0, nrow = nb, ncol = nb)
      robust_cov_beta <- matrix(0, nrow = nb, ncol = nb)
    } else if (all(gamma_active) && all(is.finite(naive_cov_red[seq_len(free_dim), seq_len(free_dim), drop = FALSE]))) {
      naive_cov_beta <- K %*% naive_cov_red[seq_len(free_dim), seq_len(free_dim), drop = FALSE] %*% t(K)
      robust_cov_beta <- if (all(is.finite(robust_cov_red[seq_len(free_dim), seq_len(free_dim), drop = FALSE]))) {
        K %*% robust_cov_red[seq_len(free_dim), seq_len(free_dim), drop = FALSE] %*% t(K)
      } else {
        matrix(NA_real_, nrow = nb, ncol = nb)
      }
    } else {
      naive_cov_beta <- matrix(NA_real_, nrow = nb, ncol = nb)
      robust_cov_beta <- matrix(NA_real_, nrow = nb, ncol = nb)
    }

    naive_cov_full <- if (all(active_red) && all(is.finite(naive_cov_red))) {
      map %*% naive_cov_red %*% t(map)
    } else {
      matrix(NA_real_, nrow = nb + 1L, ncol = nb + 1L)
    }
    robust_cov_full <- if (all(active_red) && all(is.finite(robust_cov_red))) {
      map %*% robust_cov_red %*% t(map)
    } else {
      matrix(NA_real_, nrow = nb + 1L, ncol = nb + 1L)
    }

    convergence <- 1L
    conv_raw <- fit$raw$convergence
    if (!is.null(conv_raw) && conv_raw != 0) {
      convergence <- -50L
    }
    if (all(is.na(naive_cov_beta))) {
      convergence <- -25L
    }
    if (all(is.na(robust_cov_beta)) && convergence > 0L) {
      convergence <- -27L
    }

    godambe_input <- list(
      gradient = eval_at_opt$gradient,
      hessian = eval_at_opt$hessian,
      observed_info = eval_at_opt$observed_info,
      per_subject_gradients = eval_at_opt$per_subject_gradients
    )

    list(
      theta = list(
        beta = full_par[seq_len(nb)],
        subVar = full_par[nb + 1L]
      ),
      loglik = eval_at_opt$loglik,
      score = eval_at_opt$gradient,
      hessian = eval_at_opt$hessian,
      observed_info = eval_at_opt$observed_info,
      per_subject_gradients = eval_at_opt$per_subject_gradients,
      naive_cov = naive_cov_beta,
      naive_cov_full = naive_cov_full,
      robust_cov = robust_cov_beta,
      robust_cov_full = robust_cov_full,
      godambe = extract_godambe_components(godambe_input, nb),
      convergence = convergence,
      optimizer = fit$optimizer,
      raw_fit = fit$raw,
      constraint = list(
        L = L,
        b = b,
        beta_star = beta_star,
        null_basis = K,
        free_dim = free_dim
      ),
      reduced = list(
        theta = fit$par,
        gamma = if (free_dim == 0L) numeric(0) else fit$par[seq_len(free_dim)],
        score = score_red,
        observed_info = H_red,
        per_subject_gradients = U_red
      )
    )
  })
}
