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

pmg_vector_norms <- function(x) {
  x <- as.numeric(x)
  finite <- is.finite(x)
  list(
    max_abs = if (any(finite)) max(abs(x[finite])) else NA_real_,
    l2 = if (all(finite)) sqrt(sum(x^2)) else NA_real_
  )
}

pmg_reduced_newton_diagnostics <- function(H_red, g_red, theta_red) {
  H_red <- as.matrix(H_red)
  g_red <- as.numeric(g_red)
  theta_red <- as.numeric(theta_red)

  empty <- list(
    H_red_positive_definite = FALSE,
    H_red_rcond = NA_real_,
    newton_step_red = rep(NA_real_, length(g_red)),
    newton_step_red_max_abs = NA_real_,
    max_abs_newton_step_scaled = NA_real_,
    gradient_quadratic_lr_bound = NA_real_
  )

  if (!length(g_red) ||
    !all(dim(H_red) == c(length(g_red), length(g_red))) ||
    length(theta_red) != length(g_red) ||
    !all(is.finite(H_red)) ||
    !all(is.finite(g_red)) ||
    !all(is.finite(theta_red))) {
    return(empty)
  }

  H_red <- 0.5 * (H_red + t(H_red))
  empty$H_red_rcond <- tryCatch(
    rcond(H_red),
    error = function(e) NA_real_
  )

  chol_H <- tryCatch(
    chol(H_red),
    error = function(e) NULL
  )
  if (is.null(chol_H)) {
    return(empty)
  }

  step <- tryCatch(
    backsolve(chol_H, forwardsolve(t(chol_H), g_red)),
    error = function(e) NULL
  )
  if (is.null(step) || !all(is.finite(step))) {
    return(empty)
  }

  quad <- as.numeric(crossprod(g_red, step))
  if (!is.finite(quad)) {
    return(empty)
  }

  list(
    H_red_positive_definite = TRUE,
    H_red_rcond = empty$H_red_rcond,
    newton_step_red = as.numeric(step),
    newton_step_red_max_abs = max(abs(step)),
    max_abs_newton_step_scaled = max(abs(step) / pmax(1, abs(theta_red))),
    gradient_quadratic_lr_bound = quad
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

fit_gene_pmg_unconstrained <- function(gene_index,
                                       posv,
                                       ctx,
                                       start = NULL,
                                       control = NULL,
                                       retry = FALSE,
                                       use_average_objective = FALSE) {
  with(ctx, {
    default_start <- c(pmg_initial_beta(posv, ctx), 1)
    lower <- c(rep(-100, nb), min[1])
    upper <- c(rep(100, nb), max[1])
    objective_scale <- if (isTRUE(use_average_objective)) max(k, 1L) else 1

    sanitize_start <- function(x) {
      x <- as.numeric(x)
      if (length(x) != nb + 1L || any(!is.finite(x))) {
        return(NULL)
      }
      pmin(pmax(x, lower), upper)
    }

    start_list <- list(default_start)
    user_start <- sanitize_start(start)
    if (!is.null(user_start)) {
      start_list <- c(list(user_start), start_list)
    }
    if (isTRUE(retry)) {
      sigma_base <- start_list[[1L]][nb + 1L]
      sigma_grid <- unique(pmin(pmax(c(sigma_base, 0.1, 0.5, 1, 2), min[1]), max[1]))
      for (sigma_val in sigma_grid) {
        candidate <- start_list[[1L]]
        candidate[nb + 1L] <- sigma_val
        start_list[[length(start_list) + 1L]] <- candidate
      }
    }

    obj_fn <- function(para) {
      pmg_ll(
        para = para,
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
      ) / objective_scale
    }

    grad_fn <- function(para) {
      pmg_der(
        para = para,
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
      ) / objective_scale
    }

    run_nlminb <- function(start_value, control_value) {
      tryCatch({
        args <- list(
          start = start_value,
          objective = obj_fn,
          gradient = grad_fn,
          lower = lower,
          upper = upper
        )
        if (!is.null(control_value) && length(control_value)) {
          args$control <- control_value
        }
        ref <- do.call(nlminb, args)
        list(par = ref$par, optimizer = "nlminb", raw = ref)
      }, error = function(e) NULL)
    }

    run_bobyqa <- function(start_value) {
      tryCatch({
        ref <- bobyqa(
          par = start_value,
          fn = obj_fn,
          lower = lower,
          upper = upper
        )
        list(par = ref$par, optimizer = "bobyqa", raw = ref)
      }, error = function(e) NULL)
    }

    control_list <- list(if (is.null(control)) list() else control)
    if (isTRUE(retry)) {
      control_list[[length(control_list) + 1L]] <- modifyList(
        list(eval.max = 500L, iter.max = 500L, rel.tol = 1e-10, x.tol = 1e-8),
        if (is.null(control)) list() else control
      )
    }

    raw_converged <- function(fit_try) {
      conv <- fit_try$raw$convergence
      is.null(conv) || isTRUE(conv == 0L)
    }

    fits <- list()
    retry_used_actual <- FALSE
    first_fit <- run_nlminb(start_list[[1L]], control_list[[1L]])
    if (!is.null(first_fit)) {
      fits[[length(fits) + 1L]] <- first_fit
    }
    if (is.null(first_fit) || (isTRUE(retry) && !raw_converged(first_fit))) {
      retry_used_actual <- TRUE
      for (start_value in start_list) {
        for (control_value in control_list) {
          fit_try <- run_nlminb(start_value, control_value)
          if (!is.null(fit_try)) {
            fits[[length(fits) + 1L]] <- fit_try
          }
        }
      }
      has_converged_gradient_fit <- any(vapply(fits, function(fit_try) {
        identical(fit_try$optimizer, "nlminb") && raw_converged(fit_try)
      }, logical(1)))
      if (isTRUE(retry) && length(start_list) && !has_converged_gradient_fit) {
        fit_try <- run_bobyqa(start_list[[1L]])
        if (!is.null(fit_try)) {
          fits[[length(fits) + 1L]] <- fit_try
        }
      }
    }
    if (!length(fits)) {
      retry_used_actual <- TRUE
      fit_try <- run_bobyqa(default_start)
      if (!is.null(fit_try)) {
        fits[[length(fits) + 1L]] <- fit_try
      }
    }
    if (!length(fits)) {
      stop("PGMM unconstrained optimization failed.")
    }

    fit_objective <- vapply(fits, function(fit_try) {
      pmg_ll(
        para = fit_try$par,
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
    }, numeric(1))
    fit_convergence_rank <- vapply(fits, function(fit_try) {
      conv <- fit_try$raw$convergence
      if (is.null(conv) || isTRUE(conv == 0L)) 0L else 1L
    }, integer(1))
    fit <- fits[[order(fit_objective, fit_convergence_rank)[[1L]]]]

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

    gradient_norms <- pmg_vector_norms(eval_at_opt$gradient)
    raw_convergence <- if (!is.null(fit$raw$convergence)) fit$raw$convergence else NA_integer_
    raw_message <- if (!is.null(fit$raw$message)) fit$raw$message else NA_character_

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
      diagnostics = list(
        raw_convergence = raw_convergence,
        raw_message = raw_message,
        objective = eval_at_opt$objective,
        gradient_max_abs = gradient_norms$max_abs,
        gradient_l2 = gradient_norms$l2,
        active = active,
        at_lower = fit$par <= lower + sqrt(.Machine$double.eps),
        at_upper = fit$par >= upper - sqrt(.Machine$double.eps),
        active_all = all(active),
        boundary_any = any(fit$par <= lower + sqrt(.Machine$double.eps) | fit$par >= upper - sqrt(.Machine$double.eps)),
        retry_used = retry_used_actual,
        use_average_objective = isTRUE(use_average_objective),
        lower = lower,
        upper = upper
      ),
      optimizer = fit$optimizer,
      raw_fit = fit$raw
    )
  })
}

fit_gene_pmg_constrained <- function(gene_index,
                                     posv,
                                     ctx,
                                     L,
                                     b,
                                     unconstrained_fit = NULL,
                                     warm_start = NULL,
                                     control = NULL,
                                     retry = FALSE,
                                     use_average_objective = FALSE) {
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
    objective_scale <- if (isTRUE(use_average_objective)) max(k, 1L) else 1

    reduced_start_from_warm <- function(x) {
      if (is.null(x)) {
        return(NULL)
      }
      if (!is.null(x$reduced) && !is.null(x$reduced$theta)) {
        x <- x$reduced$theta
      } else if (!is.null(x$theta) && !is.null(x$theta$beta) && !is.null(x$theta$subVar)) {
        beta_warm <- as.numeric(x$theta$beta)
        sigma_warm <- as.numeric(x$theta$subVar)
        gamma_warm <- if (free_dim == 0L) {
          numeric(0)
        } else {
          as.numeric(crossprod(K, beta_warm - beta_star))
        }
        x <- c(gamma_warm, sigma_warm)
      }
      x <- as.numeric(x)
      if (length(x) != length(theta_init) || any(!is.finite(x))) {
        return(NULL)
      }
      pmin(pmax(x, lower_theta), upper_theta)
    }

    start_list <- list(theta_init)
    warm_start_reduced <- reduced_start_from_warm(warm_start)
    if (!is.null(warm_start_reduced)) {
      start_list <- c(list(warm_start_reduced), start_list)
    }
    if (isTRUE(retry)) {
      base_start <- start_list[[1L]]
      sigma_idx <- length(base_start)
      sigma_base <- base_start[[sigma_idx]]
      sigma_grid <- unique(pmin(pmax(c(sigma_base, sigma_seed, 0.1, 0.5, 1, 2), min[1]), max[1]))
      for (sigma_val in sigma_grid) {
        candidate <- base_start
        candidate[[sigma_idx]] <- sigma_val
        start_list[[length(start_list) + 1L]] <- candidate
      }
      if (length(base_start) > 1L) {
        jitter_scale <- pmax(abs(base_start), 1) * 1e-4
        for (sign_val in c(-1, 1)) {
          candidate <- base_start + sign_val * jitter_scale
          start_list[[length(start_list) + 1L]] <- pmin(pmax(candidate, lower_theta), upper_theta)
        }
      }
    }

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
      ) / objective_scale
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
      ) / objective_scale
    }

    run_nlminb <- function(start_value, control_value) {
      tryCatch({
        args <- list(
          start = start_value,
          objective = obj_fn,
          gradient = grad_fn,
          lower = lower_theta,
          upper = upper_theta
        )
        if (!is.null(control_value) && length(control_value)) {
          args$control <- control_value
        }
        ref <- do.call(nlminb, args)
        list(par = ref$par, optimizer = "nlminb", raw = ref)
      }, error = function(e) NULL)
    }

    run_bobyqa <- function(start_value) {
      tryCatch({
        ref <- bobyqa(
          par = start_value,
          fn = obj_fn,
          lower = lower_theta,
          upper = upper_theta
        )
        list(par = ref$par, optimizer = "bobyqa", raw = ref)
      }, error = function(e) NULL)
    }

    control_list <- list(if (is.null(control)) list() else control)
    if (isTRUE(retry)) {
      control_list[[length(control_list) + 1L]] <- modifyList(
        list(eval.max = 500L, iter.max = 500L, rel.tol = 1e-10, x.tol = 1e-8),
        if (is.null(control)) list() else control
      )
    }

    raw_converged <- function(fit_try) {
      conv <- fit_try$raw$convergence
      is.null(conv) || isTRUE(conv == 0L)
    }

    fits <- list()
    retry_used_actual <- FALSE
    first_fit <- run_nlminb(start_list[[1L]], control_list[[1L]])
    if (!is.null(first_fit)) {
      fits[[length(fits) + 1L]] <- first_fit
    }
    if (is.null(first_fit) || (isTRUE(retry) && !raw_converged(first_fit))) {
      retry_used_actual <- TRUE
      for (start_value in start_list) {
        for (control_value in control_list) {
          fit_try <- run_nlminb(start_value, control_value)
          if (!is.null(fit_try)) {
            fits[[length(fits) + 1L]] <- fit_try
          }
        }
      }
      has_converged_gradient_fit <- any(vapply(fits, function(fit_try) {
        identical(fit_try$optimizer, "nlminb") && raw_converged(fit_try)
      }, logical(1)))
      if (isTRUE(retry) && length(start_list) && !has_converged_gradient_fit) {
        fit_try <- run_bobyqa(start_list[[1L]])
        if (!is.null(fit_try)) {
          fits[[length(fits) + 1L]] <- fit_try
        }
      }
    }
    if (!length(fits)) {
      retry_used_actual <- TRUE
      fit_try <- run_bobyqa(theta_init)
      if (!is.null(fit_try)) {
        fits[[length(fits) + 1L]] <- fit_try
      }
    }
    if (!length(fits)) {
      stop("PGMM constrained optimization failed.")
    }

    fit_objective <- vapply(fits, function(fit_try) {
      pmg_ll(
        para = build_full_para(fit_try$par),
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
    }, numeric(1))
    fit_objective_finite <- fit_objective[is.finite(fit_objective)]
    fit_objective_min <- if (length(fit_objective_finite)) min(fit_objective_finite) else NA_real_
    fit_objective_max <- if (length(fit_objective_finite)) max(fit_objective_finite) else NA_real_
    fit_objective_spread <- if (length(fit_objective_finite)) {
      fit_objective_max - fit_objective_min
    } else {
      NA_real_
    }
    first_fit_objective <- if (length(fit_objective)) fit_objective[[1L]] else NA_real_
    first_fit_objective_minus_best <- if (is.finite(first_fit_objective) && is.finite(fit_objective_min)) {
      first_fit_objective - fit_objective_min
    } else {
      NA_real_
    }
    fit_convergence_rank <- vapply(fits, function(fit_try) {
      conv <- fit_try$raw$convergence
      if (is.null(conv) || isTRUE(conv == 0L)) 0L else 1L
    }, integer(1))
    fit <- fits[[order(fit_objective, fit_convergence_rank)[[1L]]]]

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
    H_red <- 0.5 * (H_red + t(H_red))
    U_red <- crossprod(map, eval_at_opt$per_subject_gradients)
    score_red <- as.numeric(crossprod(map, eval_at_opt$gradient))
    reduced_newton <- pmg_reduced_newton_diagnostics(
      H_red = H_red,
      g_red = score_red,
      theta_red = fit$par
    )

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

    full_gradient_norms <- pmg_vector_norms(eval_at_opt$gradient)
    reduced_gradient_norms <- pmg_vector_norms(score_red)
    raw_convergence <- if (!is.null(fit$raw$convergence)) fit$raw$convergence else NA_integer_
    raw_message <- if (!is.null(fit$raw$message)) fit$raw$message else NA_character_

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
      lrt_gate = list(
        converged = convergence > 0L,
        interior = all(active_red),
        nuisance_interior = if (length(active_red) > 0L) isTRUE(active_red[length(active_red)]) else FALSE
      ),
      diagnostics = list(
        raw_convergence = raw_convergence,
        raw_message = raw_message,
        objective = eval_at_opt$objective,
        full_gradient_max_abs = full_gradient_norms$max_abs,
        full_gradient_l2 = full_gradient_norms$l2,
        reduced_gradient_max_abs = reduced_gradient_norms$max_abs,
        reduced_gradient_l2 = reduced_gradient_norms$l2,
        H_red_positive_definite = reduced_newton$H_red_positive_definite,
        H_red_rcond = reduced_newton$H_red_rcond,
        newton_step_red = reduced_newton$newton_step_red,
        newton_step_red_max_abs = reduced_newton$newton_step_red_max_abs,
        max_abs_newton_step_scaled = reduced_newton$max_abs_newton_step_scaled,
        gradient_quadratic_lr_bound = reduced_newton$gradient_quadratic_lr_bound,
        active_reduced = active_red,
        at_lower_reduced = fit$par <= lower_theta + sqrt(.Machine$double.eps),
        at_upper_reduced = fit$par >= upper_theta - sqrt(.Machine$double.eps),
        active_reduced_all = all(active_red),
        boundary_reduced_any = any(fit$par <= lower_theta + sqrt(.Machine$double.eps) | fit$par >= upper_theta - sqrt(.Machine$double.eps)),
        retry_used = retry_used_actual,
        optimizer_attempts = length(fits),
        optimizer_objective_min = fit_objective_min,
        optimizer_objective_max = fit_objective_max,
        optimizer_objective_spread = fit_objective_spread,
        optimizer_first_objective = first_fit_objective,
        optimizer_first_objective_minus_best = first_fit_objective_minus_best,
        use_average_objective = isTRUE(use_average_objective),
        lower_reduced = lower_theta,
        upper_reduced = upper_theta
      ),
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
        newton_step = reduced_newton$newton_step_red,
        observed_info = H_red,
        per_subject_gradients = U_red
      )
    )
  })
}
