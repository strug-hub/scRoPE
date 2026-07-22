# Estimator-aligned joint-LN fitting ---------------------------------------

ln_joint_thresholds <- function() {
  list(
    donor_additivity = 1e-10,
    objective_gradient = 1e-5,
    sensitivity_score = 1e-4,
    finite_difference_stability = 1e-4,
    root = 1e-6,
    newton = 1e-6,
    parameter_agreement = 1e-6,
    objective_agreement = 1e-8,
    constraint = 1e-10,
    curvature_min_eigen = 1e-8,
    condition_number = 1e10,
    covariance_identity = 1e-8,
    natural_log_coordinate = 1e-7,
    internal_public_coordinate = 1e-7,
    mode_score = 1e-10,
    mode_newton = 1e-10,
    hard_stop = list(
      donor_additivity = 1e-8,
      objective_gradient = 1e-4,
      sensitivity_score = 1e-3,
      finite_difference_stability = 1e-3,
      root = 1e-4,
      newton = 1e-4,
      parameter_agreement = 1e-4,
      objective_agreement = 1e-6,
      constraint = 1e-8,
      condition_number = 1e12,
      covariance_identity = 1e-7,
      natural_log_coordinate = 1e-6,
      internal_public_coordinate = 1e-6,
      mode_score = 1e-8,
      mode_newton = 1e-8
    )
  )
}

ln_joint_objective_provenance <- function() {
  list(
    id = "joint_ln_ptmg_working_objective_v1",
    estimator = "joint_ln",
    criterion = "ptmg_ll",
    derivative_criterion = "ptmg_ll_der_hes4",
    objective_convention = "minimized_q_equals_negative_working_log_objective",
    score_convention = "maximizing_U_equals_negative_gradient_of_q",
    objective_coordinates = "beta, sigma2, gamma_nb_size",
    derivative_coordinates = "beta, log(sigma2), log(gamma_nb_size)",
    donor_unit = "independent_subject",
    subject_additive = TRUE,
    corrected_zero_based_positive_cell_index = TRUE
  )
}

ln_joint_objective_instance_registry <- local({
  registry <- new.env(parent = emptyenv())
  counter <- 0
  max_entries <- 16L

  matching_key <- function(gene_index, posv, ctx) {
    if (length(gene_index) != 1L || !is.numeric(gene_index) ||
        !is.finite(gene_index) || !is.list(posv) || !is.list(ctx)) {
      return(NULL)
    }
    keys <- ls(registry, all.names = TRUE)
    if (!length(keys)) {
      return(NULL)
    }
    matches <- vapply(keys, function(key) {
      entry <- get(key, envir = registry, inherits = FALSE)
      identical(entry$token$gene_index, as.integer(gene_index)) &&
        identical(entry$objective_context$posv, posv) &&
        identical(entry$objective_context$ctx, ctx)
    }, logical(1))
    keys <- keys[matches]
    if (!length(keys)) {
      return(NULL)
    }
    sequences <- vapply(keys, function(key) {
      get(key, envir = registry, inherits = FALSE)$sequence
    }, numeric(1))
    keys[[which.max(sequences)]]
  }

  function(action = c(
    "register",
    "lookup",
    "store_unrestricted",
    "get_unrestricted",
    "resolve"
  ), gene_index = NULL, posv = NULL, ctx = NULL, fit = NULL,
  instance = NULL) {
    action <- match.arg(action)
    if (identical(action, "resolve")) {
      if (!is.list(instance) || is.null(instance$id)) {
        return(NULL)
      }
      keys <- ls(registry, all.names = TRUE)
      for (key in keys) {
        entry <- get(key, envir = registry, inherits = FALSE)
        if (identical(entry$token, instance)) {
          return(entry$objective_context)
        }
      }
      return(NULL)
    }
    if (length(gene_index) != 1L || !is.numeric(gene_index) ||
        !is.finite(gene_index) || !is.list(posv) || !is.list(ctx)) {
      return(NULL)
    }
    if (identical(action, "register")) {
      counter <<- counter + 1
      key <- paste0("objective_instance_", Sys.getpid(), "_", counter)
      token <- list(
        scheme = "process_local_exact_joint_ln_objective_context_v1",
        id = paste(
          Sys.getpid(),
          format(counter, scientific = FALSE),
          sep = "|"
        ),
        gene_index = as.integer(gene_index)
      )
      assign(
        key,
        list(
          token = token,
          sequence = counter,
          unrestricted_fit = NULL,
          objective_context = list(
            gene_index = as.integer(gene_index),
            posv = posv,
            ctx = ctx
          )
        ),
        envir = registry
      )
      existing <- ls(registry, all.names = TRUE)
      if (length(existing) > max_entries) {
        sequences <- vapply(existing, function(existing_key) {
          get(existing_key, envir = registry, inherits = FALSE)$sequence
        }, numeric(1))
        remove_keys <- existing[order(sequences)][seq_len(
          length(existing) - max_entries
        )]
        rm(list = remove_keys, envir = registry)
      }
      return(token)
    }
    key <- matching_key(gene_index, posv, ctx)
    if (is.null(key)) {
      return(NULL)
    }
    entry <- get(key, envir = registry, inherits = FALSE)
    if (identical(action, "lookup")) {
      return(entry$token)
    }
    if (identical(action, "store_unrestricted")) {
      if (!is.list(fit)) {
        return(NULL)
      }
      entry$unrestricted_fit <- fit
      assign(key, entry, envir = registry)
      return(entry$token)
    }
    entry$unrestricted_fit
  }
})

ln_joint_scaled_max <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (!length(a) || length(a) != length(b) || any(!is.finite(c(a, b)))) {
    return(NA_real_)
  }
  max(abs(a - b) / (1 + pmax(abs(a), abs(b))))
}

ln_joint_relative_frobenius <- function(a, b) {
  if (is.null(a) || is.null(b) || !all(dim(a) == dim(b)) ||
      any(!is.finite(c(a, b)))) {
    return(NA_real_)
  }
  as.numeric(norm(a - b, "F") / max(norm(b, "F"), .Machine$double.eps))
}

ln_joint_map <- function(nb, L = NULL, b = NULL) {
  if (is.null(L)) {
    beta_star <- rep(0, nb)
    K <- diag(nb)
    constrained <- FALSE
  } else {
    constr <- constraint_reparameterization(L, b)
    beta_star <- as.numeric(constr$beta_star)
    K <- constr$null_basis
    if (is.null(K)) {
      K <- matrix(0, nrow = nb, ncol = 0L)
    }
    constrained <- TRUE
  }
  free_dim <- ncol(K)
  compose <- function(zeta) {
    if (!free_dim) beta_star else as.numeric(beta_star + K %*% zeta)
  }
  project <- function(beta) {
    if (!free_dim) {
      numeric(0)
    } else {
      as.numeric(qr.solve(K, as.numeric(beta) - beta_star, tol = 1e-10))
    }
  }
  derivative_map <- matrix(0, nrow = nb + 2L, ncol = free_dim + 2L)
  if (free_dim) {
    derivative_map[seq_len(nb), seq_len(free_dim)] <- K
  }
  derivative_map[nb + 1L, free_dim + 1L] <- 1
  derivative_map[nb + 2L, free_dim + 2L] <- 1
  list(
    constrained = constrained,
    beta_star = beta_star,
    K = K,
    free_dim = free_dim,
    compose = compose,
    project = project,
    derivative_map = derivative_map,
    L = L,
    b = b
  )
}

ln_joint_full_natural <- function(theta_reduced, map) {
  free_dim <- map$free_dim
  zeta <- if (free_dim) theta_reduced[seq_len(free_dim)] else numeric(0)
  dispersion <- theta_reduced[free_dim + seq_len(2L)]
  c(map$compose(zeta), dispersion)
}

ln_joint_objective <- function(theta_reduced, gene_index, posv, ctx, map) {
  full <- ln_joint_full_natural(theta_reduced, map)
  ptmg_ll(
    full,
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    n_one = posv$n_onetwo,
    ytwo = posv$ytwo,
    fam = ctx$id,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[gene_index, ],
    posind = ctx$posind[[gene_index]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k
  )
}

ln_joint_gradient <- function(theta_reduced, gene_index, posv, ctx, map) {
  full <- ln_joint_full_natural(theta_reduced, map)
  gradient_full <- ptmg_der(
    full,
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    n_one = posv$n_onetwo,
    ytwo = posv$ytwo,
    fam = ctx$id,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[gene_index, ],
    posind = ctx$posind[[gene_index]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k
  )
  as.numeric(crossprod(map$derivative_map, gradient_full))
}

ln_joint_derivative_state <- function(theta_reduced_log, gene_index, posv, ctx, map) {
  free_dim <- map$free_dim
  zeta <- if (free_dim) theta_reduced_log[seq_len(free_dim)] else numeric(0)
  beta <- map$compose(zeta)
  full_log <- c(beta, theta_reduced_log[free_dim + seq_len(2L)])
  raw <- ptmg_ll_der_hes4(
    full_log,
    ctx$pred,
    ctx$offset,
    posv$Y,
    posv$n_onetwo,
    posv$ytwo,
    ctx$fid,
    ctx$cumsumy[gene_index, ],
    ctx$posind[[gene_index]],
    posv$posindy,
    ctx$nb,
    ctx$nind,
    ctx$k,
    correct_posindy = TRUE
  )
  T <- map$derivative_map
  A_full <- 0.5 * (raw$hessian + t(raw$hessian))
  U_full <- -as.numeric(raw$gradient)
  U_subject_full <- -as.matrix(raw$per_subject_gradients)
  A_reduced <- crossprod(T, A_full %*% T)
  A_reduced <- 0.5 * (A_reduced + t(A_reduced))
  U_reduced <- as.numeric(crossprod(T, U_full))
  U_subject_reduced <- crossprod(T, U_subject_full)
  list(
    full_log = full_log,
    theta_reduced_log = theta_reduced_log,
    objective_provenance = ln_joint_objective_provenance(),
    beta = beta,
    sigma2 = exp(full_log[ctx$nb + 1L]),
    gamma = exp(full_log[ctx$nb + 2L]),
    q_gradient_full = as.numeric(raw$gradient),
    q_gradient_reduced = -U_reduced,
    U_full = U_full,
    U_reduced = U_reduced,
    U_subject_full = U_subject_full,
    U_subject_reduced = U_subject_reduced,
    A_full = A_full,
    A_reduced = A_reduced,
    raw = raw
  )
}

ln_joint_curvature_diagnostics <- function(A) {
  A <- 0.5 * (A + t(A))
  if (!length(A) || any(!is.finite(A))) {
    return(list(
      scaled_min_eigenvalue = NA_real_,
      scaled_condition_number = Inf,
      positive_definite = FALSE
    ))
  }
  d <- pmax(abs(diag(A)), sqrt(.Machine$double.eps))
  scaled <- A / sqrt(outer(d, d))
  eig <- tryCatch(
    eigen(scaled, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eig) || any(!is.finite(eig))) {
    return(list(
      scaled_min_eigenvalue = NA_real_,
      scaled_condition_number = Inf,
      positive_definite = FALSE
    ))
  }
  min_eig <- min(eig)
  cond <- if (min_eig > 0) max(eig) / min_eig else Inf
  list(
    scaled_min_eigenvalue = min_eig,
    scaled_condition_number = cond,
    positive_definite = min_eig > 0,
    scaled_eigenvalues = eig
  )
}

ln_joint_score_matrix_diagnostics <- function(x,
                                              condition_limit =
                                                ln_joint_thresholds()$condition_number,
                                              symmetry_tolerance =
                                                ln_joint_thresholds()$hard_stop$
                                                  covariance_identity) {
  unavailable <- function(nrow = NA_integer_, ncol = NA_integer_) {
    list(
      nrow = as.integer(nrow),
      ncol = as.integer(ncol),
      finite = FALSE,
      square = FALSE,
      symmetry_error = NA_real_,
      symmetry_acceptable = FALSE,
      rank = 0L,
      rank_tolerance = NA_real_,
      condition_number = Inf,
      scaled_condition_number = Inf,
      eigenvalues = NULL,
      scaled_eigenvalues = NULL,
      positive_definite = FALSE,
      full_rank = FALSE,
      condition_acceptable = FALSE,
      usable = FALSE
    )
  }
  if (is.null(x)) {
    return(unavailable())
  }
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  if (!length(x) || nr != nc || any(!is.finite(x))) {
    out <- unavailable(nr, nc)
    out$finite <- length(x) > 0L && all(is.finite(x))
    out$square <- nr == nc
    return(out)
  }
  symmetry_error <- as.numeric(
    norm(x - t(x), "F") / (1 + norm(x, "F"))
  )
  symmetric_x <- 0.5 * (x + t(x))
  singular_values <- tryCatch(
    svd(symmetric_x, nu = 0L, nv = 0L)$d,
    error = function(e) NULL
  )
  eigenvalues <- tryCatch(
    eigen(symmetric_x, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  diagonal_scale <- sqrt(pmax(
    abs(diag(symmetric_x)),
    .Machine$double.eps
  ))
  scaled_x <- symmetric_x / outer(diagonal_scale, diagonal_scale)
  scaled_singular_values <- tryCatch(
    svd(scaled_x, nu = 0L, nv = 0L)$d,
    error = function(e) NULL
  )
  scaled_eigenvalues <- tryCatch(
    eigen(scaled_x, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(singular_values) || is.null(eigenvalues) ||
      is.null(scaled_singular_values) || is.null(scaled_eigenvalues) ||
      any(!is.finite(c(
        singular_values,
        eigenvalues,
        scaled_singular_values,
        scaled_eigenvalues
      )))) {
    out <- unavailable(nr, nc)
    out$finite <- TRUE
    out$square <- TRUE
    out$symmetry_error <- symmetry_error
    return(out)
  }
  max_singular <- max(singular_values)
  rank_tolerance <- max(nr, nc) * .Machine$double.eps * max_singular
  rank <- if (max_singular > 0) {
    sum(singular_values > rank_tolerance)
  } else {
    0L
  }
  scaled_max_singular <- max(scaled_singular_values)
  scaled_rank_tolerance <- max(nr, nc) * .Machine$double.eps *
    scaled_max_singular
  scaled_rank <- if (scaled_max_singular > 0) {
    sum(scaled_singular_values > scaled_rank_tolerance)
  } else {
    0L
  }
  condition_number <- if (rank == nr && min(singular_values) > 0) {
    max_singular / min(singular_values)
  } else {
    Inf
  }
  scaled_condition_number <- if (scaled_rank == nr &&
      min(scaled_singular_values) > 0) {
    scaled_max_singular / min(scaled_singular_values)
  } else {
    Inf
  }
  eigen_tolerance <- nr * .Machine$double.eps *
    max(abs(eigenvalues))
  positive_definite <- min(eigenvalues) > eigen_tolerance
  full_rank <- rank == nr
  symmetry_acceptable <- is.finite(symmetry_error) &&
    symmetry_error <= symmetry_tolerance
  condition_acceptable <- is.finite(scaled_condition_number) &&
    scaled_condition_number <= condition_limit
  list(
    nrow = as.integer(nr),
    ncol = as.integer(nc),
    finite = TRUE,
    square = TRUE,
    symmetry_error = symmetry_error,
    symmetry_acceptable = symmetry_acceptable,
    rank = as.integer(rank),
    scaled_rank = as.integer(scaled_rank),
    rank_tolerance = rank_tolerance,
    scaled_rank_tolerance = scaled_rank_tolerance,
    condition_number = condition_number,
    scaled_condition_number = scaled_condition_number,
    eigenvalues = eigenvalues,
    scaled_eigenvalues = scaled_eigenvalues,
    positive_definite = positive_definite,
    full_rank = full_rank,
    condition_acceptable = condition_acceptable,
    usable = symmetry_acceptable && full_rank && positive_definite &&
      condition_acceptable
  )
}

ln_joint_score_coordinate_map <- function(L, nb, null_basis = NULL,
                                          interest_basis = NULL) {
  thresholds <- ln_joint_thresholds()
  L <- as.matrix(L)
  r <- nrow(L)
  observed_rank <- if (length(L) && all(is.finite(L))) {
    qr(L)$rank
  } else {
    0L
  }
  unavailable <- function(reason) {
    list(
      success = FALSE,
      failure_reason = reason,
      requested_rows = as.integer(r),
      rank = as.integer(observed_rank)
    )
  }
  if (!length(L) || ncol(L) != nb || r < 1L || r > nb ||
      any(!is.finite(L))) {
    return(unavailable("invalid_hypothesis_matrix"))
  }
  if (observed_rank != r) {
    return(unavailable("hypothesis_rank_deficient"))
  }
  qrt <- qr(t(L))
  if (qrt$rank != r) {
    return(unavailable("hypothesis_rank_deficient"))
  }
  complete_basis <- qr.Q(qrt, complete = TRUE)
  if (is.null(interest_basis)) {
    interest_basis <- complete_basis[, seq_len(r), drop = FALSE]
  } else {
    interest_basis <- as.matrix(interest_basis)
  }
  nuisance_beta_dim <- nb - r
  if (is.null(null_basis)) {
    null_basis <- if (nuisance_beta_dim) {
      complete_basis[, r + seq_len(nuisance_beta_dim), drop = FALSE]
    } else {
      matrix(0, nrow = nb, ncol = 0L)
    }
  } else {
    null_basis <- as.matrix(null_basis)
  }
  if (!all(dim(interest_basis) == c(nb, r)) ||
      any(!is.finite(interest_basis))) {
    return(unavailable("invalid_interest_basis"))
  }
  if (!all(dim(null_basis) == c(nb, nuisance_beta_dim)) ||
      any(!is.finite(null_basis))) {
    return(unavailable("invalid_null_basis"))
  }
  tested_map <- L %*% interest_basis
  if (qr(tested_map)$rank != r) {
    return(unavailable("interest_basis_rank_deficient"))
  }
  null_residual_scaled <- if (nuisance_beta_dim) {
    max(abs(L %*% null_basis)) /
      (1 + norm(L, "F") * norm(null_basis, "F"))
  } else {
    0
  }
  if (!is.finite(null_residual_scaled) ||
      null_residual_scaled > thresholds$hard_stop$constraint) {
    return(unavailable("null_basis_not_in_hypothesis_kernel"))
  }
  beta_map <- cbind(interest_basis, null_basis)
  if (qr(beta_map)$rank != nb) {
    return(unavailable("interest_nuisance_map_rank_deficient"))
  }
  beta_map_condition <- tryCatch(
    kappa(beta_map, exact = TRUE),
    error = function(e) Inf
  )
  if (!is.finite(beta_map_condition) ||
      beta_map_condition > thresholds$hard_stop$condition_number) {
    return(unavailable("interest_nuisance_map_ill_conditioned"))
  }
  full_dim <- nb + 2L
  full_map <- matrix(0, nrow = full_dim, ncol = full_dim)
  full_map[seq_len(nb), seq_len(nb)] <- beta_map
  full_map[nb + 1L, nb + 1L] <- 1
  full_map[nb + 2L, nb + 2L] <- 1
  list(
    success = TRUE,
    failure_reason = NA_character_,
    rank = as.integer(r),
    interest_basis = interest_basis,
    null_basis = null_basis,
    beta = beta_map,
    full = full_map,
    tested_hypothesis_map = tested_map,
    null_residual_scaled = null_residual_scaled,
    beta_map_condition_number = beta_map_condition,
    interest_crossprod = crossprod(interest_basis),
    null_crossprod = crossprod(null_basis)
  )
}

ln_joint_robust_score_diagnostic <- function(constrained_fit,
                                             null_basis = NULL,
                                             interest_basis = NULL) {
  thresholds <- ln_joint_thresholds()
  empty_matrix_diagnostics <- ln_joint_score_matrix_diagnostics(NULL)
  result <- list(
    diagnostic_available = FALSE,
    failure_reason = NA_character_,
    public_available = FALSE,
    public_status = "joint_ln_calibration_pending",
    statistic = NA_real_,
    df = NA_integer_,
    theoretical_df = NA_integer_,
    theoretical_reference = paste(
      "chi-square first-order reference only under smooth interior",
      "independent-donor fixed-dimension regularity"
    ),
    score_sign = "maximizing_U_equals_negative_gradient_of_q",
    no_numerical_repair = TRUE,
    hypothesis = NULL,
    coordinate_map = NULL,
    U_coordinate = NULL,
    U_subject_coordinate = NULL,
    U_eff = NULL,
    U_eff_subject = NULL,
    H_psi_chi = NULL,
    H_chi_chi = NULL,
    J_eff = NULL,
    projection = NULL,
    diagnostics = list(
      constraint_residual_scaled = NA_real_,
      nuisance_root_scaled = NA_real_,
      full_donor_additivity_scaled = NA_real_,
      coordinate_donor_additivity_scaled = NA_real_,
      efficient_donor_additivity_scaled = NA_real_,
      efficient_vs_tested_score_scaled = NA_real_,
      J_eff_projection_identity_scaled = NA_real_,
      score_sign_aggregate_scaled = NA_real_,
      score_sign_donor_scaled = NA_real_,
      nuisance_root_scale = paste(
        "max_j abs(sum_i U_ji) /",
        "(1 + sqrt(sum_i U_ji^2))"
      ),
      H_coordinate_symmetry_scaled = NA_real_,
      H_chi_chi = empty_matrix_diagnostics,
      J_eff = empty_matrix_diagnostics
    )
  )
  fail <- function(reason) {
    result$failure_reason <- reason
    result
  }

  constraint <- constrained_fit$constraint
  L <- constraint$L
  b <- constraint$b
  if (is.null(L) || is.null(b)) {
    return(fail("missing_hypothesis"))
  }
  L <- as.matrix(L)
  b <- as.numeric(b)
  if (!nrow(L) || length(b) != nrow(L) ||
      any(!is.finite(c(L, b)))) {
    return(fail("invalid_hypothesis"))
  }
  hypothesis_rank <- qr(L)$rank
  result$hypothesis <- list(
    L = L,
    b = b,
    requested_rows = as.integer(nrow(L)),
    rank = as.integer(hypothesis_rank)
  )
  if (hypothesis_rank != nrow(L)) {
    return(fail("hypothesis_rank_deficient"))
  }
  result$df <- result$theoretical_df <- as.integer(nrow(L))

  if (!identical(constrained_fit$route, "R0")) {
    reason <- if (identical(constrained_fit$failure_reason,
                            "dispersion_boundary")) {
      "constrained_boundary"
    } else {
      "constrained_route_not_R0"
    }
    return(fail(reason))
  }
  if (!isTRUE(constrained_fit$inference_available) ||
      !isTRUE(constrained_fit$lrt_gate$converged)) {
    return(fail("constrained_endpoint_uncertified"))
  }
  if (!isTRUE(constrained_fit$lrt_gate$interior) ||
      !isTRUE(constrained_fit$lrt_gate$nuisance_interior)) {
    return(fail("constrained_boundary"))
  }
  state <- constrained_fit$ln_joint$state
  if (is.null(state) || is.null(state$beta)) {
    return(fail("missing_constrained_score_state"))
  }
  nb <- length(state$beta)
  full_dim <- nb + 2L
  U_full <- as.numeric(state$U_full)
  U_subject_full <- as.matrix(state$U_subject_full)
  H_full <- as.matrix(state$A_full)
  state_dimensions_ok <- length(U_full) == full_dim &&
    nrow(U_subject_full) == full_dim && ncol(U_subject_full) > 0L &&
    all(dim(H_full) == c(full_dim, full_dim))
  if (!state_dimensions_ok || any(!is.finite(c(
    state$beta,
    U_full,
    U_subject_full,
    H_full
  )))) {
    return(fail("invalid_constrained_score_state"))
  }
  state_full_log <- as.numeric(state$full_log)
  state_coordinates_match <- length(state_full_log) == full_dim &&
    identical(
      state_full_log[seq_len(nb)],
      as.numeric(state$beta)
    ) &&
    length(state$sigma2) == 1L && length(state$gamma) == 1L &&
    identical(
      as.numeric(state$sigma2),
      exp(state_full_log[nb + 1L])
    ) &&
    identical(
      as.numeric(state$gamma),
      exp(state_full_log[nb + 2L])
    )
  if (!isTRUE(state_coordinates_match)) {
    return(fail("constrained_state_coordinate_mismatch"))
  }
  returned_theta <- constrained_fit$theta
  returned_state_matches <- is.list(returned_theta) &&
    identical(as.numeric(returned_theta$beta), as.numeric(state$beta)) &&
    identical(as.numeric(returned_theta$subVar), as.numeric(state$sigma2)) &&
    identical(as.numeric(returned_theta$cellVar), as.numeric(state$gamma)) &&
    identical(
      as.numeric(constrained_fit$beta_constrained),
      as.numeric(state$beta)
    )
  if (!isTRUE(returned_state_matches)) {
    return(fail("constrained_return_state_mismatch"))
  }
  candidate_records <- constrained_fit$ln_joint$candidate_records
  certified_state_match <- length(candidate_records) > 0L && any(vapply(
    candidate_records,
    function(candidate) {
      isTRUE(candidate$certified) &&
        !is.null(candidate$theta_full_log) &&
        identical(
          as.numeric(candidate$theta_full_log),
          state_full_log
        )
    },
    logical(1)
  ))
  if (!certified_state_match) {
    return(fail("constrained_endpoint_certification_missing"))
  }

  map <- ln_joint_score_coordinate_map(
    L,
    nb,
    null_basis = if (is.null(null_basis)) constraint$null_basis else null_basis,
    interest_basis = interest_basis
  )
  if (!isTRUE(map$success)) {
    return(fail(map$failure_reason))
  }
  result$coordinate_map <- map

  constraint_residual <- as.numeric(L %*% state$beta - b)
  constraint_denom <- 1 + max(abs(b)) +
    max(rowSums(abs(L))) * max(abs(state$beta))
  constraint_residual_scaled <- max(abs(constraint_residual)) /
    constraint_denom
  result$diagnostics$constraint_residual_scaled <-
    constraint_residual_scaled
  if (!is.finite(constraint_residual_scaled) ||
      constraint_residual_scaled > thresholds$constraint) {
    return(fail("constraint_not_satisfied"))
  }

  score_sign_aggregate_scaled <- ln_joint_scaled_max(
    U_full,
    -state$q_gradient_full
  )
  score_sign_donor_scaled <- ln_joint_scaled_max(
    U_subject_full,
    -state$raw$per_subject_gradients
  )
  full_donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_subject_full),
    U_full
  )
  result$diagnostics$score_sign_aggregate_scaled <-
    score_sign_aggregate_scaled
  result$diagnostics$score_sign_donor_scaled <- score_sign_donor_scaled
  result$diagnostics$full_donor_additivity_scaled <-
    full_donor_additivity_scaled
  if (!all(is.finite(c(
    score_sign_aggregate_scaled,
    score_sign_donor_scaled
  ))) || max(score_sign_aggregate_scaled, score_sign_donor_scaled) >
      thresholds$donor_additivity) {
    return(fail("maximizing_score_sign_mismatch"))
  }
  if (!is.finite(full_donor_additivity_scaled) ||
      full_donor_additivity_scaled > thresholds$donor_additivity) {
    return(fail("full_score_donor_additivity_failed"))
  }

  T <- map$full
  U_coordinate <- as.numeric(crossprod(T, U_full))
  U_subject_coordinate <- crossprod(T, U_subject_full)
  H_coordinate_raw <- crossprod(T, H_full %*% T)
  H_coordinate_symmetry_scaled <- as.numeric(
    norm(H_coordinate_raw - t(H_coordinate_raw), "F") /
      (1 + norm(H_coordinate_raw, "F"))
  )
  result$diagnostics$H_coordinate_symmetry_scaled <-
    H_coordinate_symmetry_scaled
  if (!is.finite(H_coordinate_symmetry_scaled) ||
      H_coordinate_symmetry_scaled >
        thresholds$hard_stop$covariance_identity) {
    return(fail("transformed_sensitivity_asymmetric"))
  }
  H_coordinate <- H_coordinate_raw
  H_coordinate <- 0.5 * (H_coordinate + t(H_coordinate))
  r <- result$df
  interest <- seq_len(r)
  nuisance <- r + seq_len(full_dim - r)
  H_psi_chi <- H_coordinate[interest, nuisance, drop = FALSE]
  H_chi_chi <- H_coordinate[nuisance, nuisance, drop = FALSE]
  coordinate_donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_subject_coordinate),
    U_coordinate
  )
  nuisance_root_scaled <- ln_joint_root_metric(
    U_subject_coordinate[nuisance, , drop = FALSE],
    U_coordinate[nuisance]
  )
  result$U_coordinate <- U_coordinate
  result$U_subject_coordinate <- U_subject_coordinate
  result$H_psi_chi <- H_psi_chi
  result$H_chi_chi <- H_chi_chi
  result$diagnostics$coordinate_donor_additivity_scaled <-
    coordinate_donor_additivity_scaled
  result$diagnostics$nuisance_root_scaled <- nuisance_root_scaled
  if (!is.finite(coordinate_donor_additivity_scaled) ||
      coordinate_donor_additivity_scaled > thresholds$donor_additivity) {
    return(fail("transformed_score_donor_additivity_failed"))
  }
  if (!is.finite(nuisance_root_scaled) ||
      nuisance_root_scaled > thresholds$root) {
    return(fail("constrained_nuisance_root_not_certified"))
  }

  H_chi_diagnostics <- ln_joint_score_matrix_diagnostics(H_chi_chi)
  result$diagnostics$H_chi_chi <- H_chi_diagnostics
  if (!isTRUE(H_chi_diagnostics$symmetry_acceptable)) {
    return(fail("H_chi_chi_asymmetric"))
  }
  if (!isTRUE(H_chi_diagnostics$full_rank)) {
    return(fail("H_chi_chi_rank_deficient"))
  }
  if (!isTRUE(H_chi_diagnostics$positive_definite)) {
    return(fail("H_chi_chi_not_positive_definite"))
  }
  if (!isTRUE(H_chi_diagnostics$condition_acceptable)) {
    return(fail("H_chi_chi_ill_conditioned"))
  }
  nuisance_adjustment <- tryCatch(
    t(solve(t(H_chi_chi), t(H_psi_chi))),
    error = function(e) NULL
  )
  if (is.null(nuisance_adjustment) ||
      any(!is.finite(nuisance_adjustment))) {
    return(fail("H_chi_chi_solve_failed"))
  }
  projection <- cbind(diag(r), -nuisance_adjustment)
  U_eff <- as.numeric(projection %*% U_coordinate)
  U_eff_subject <- projection %*% U_subject_coordinate
  J_eff <- tcrossprod(U_eff_subject)
  J_coordinate <- tcrossprod(U_subject_coordinate)
  J_eff_projected <- projection %*% J_coordinate %*% t(projection)
  efficient_donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_eff_subject),
    U_eff
  )
  efficient_vs_tested_score_scaled <- ln_joint_scaled_max(
    U_eff,
    U_coordinate[interest]
  )
  J_eff_projection_identity_scaled <- ln_joint_relative_frobenius(
    J_eff,
    J_eff_projected
  )
  result$projection <- projection
  result$U_eff <- U_eff
  result$U_eff_subject <- U_eff_subject
  result$J_eff <- J_eff
  result$diagnostics$efficient_donor_additivity_scaled <-
    efficient_donor_additivity_scaled
  result$diagnostics$efficient_vs_tested_score_scaled <-
    efficient_vs_tested_score_scaled
  result$diagnostics$J_eff_projection_identity_scaled <-
    J_eff_projection_identity_scaled
  if (!is.finite(efficient_donor_additivity_scaled) ||
      efficient_donor_additivity_scaled > thresholds$donor_additivity) {
    return(fail("efficient_score_donor_additivity_failed"))
  }
  if (!is.finite(efficient_vs_tested_score_scaled) ||
      efficient_vs_tested_score_scaled > thresholds$root) {
    return(fail("efficient_score_not_equal_tested_score_at_null"))
  }
  if (!is.finite(J_eff_projection_identity_scaled) ||
      J_eff_projection_identity_scaled > thresholds$covariance_identity) {
    return(fail("J_eff_projection_identity_failed"))
  }

  J_eff_diagnostics <- ln_joint_score_matrix_diagnostics(J_eff)
  result$diagnostics$J_eff <- J_eff_diagnostics
  if (!isTRUE(J_eff_diagnostics$symmetry_acceptable)) {
    return(fail("J_eff_asymmetric"))
  }
  if (!isTRUE(J_eff_diagnostics$full_rank)) {
    return(fail("J_eff_rank_deficient"))
  }
  if (!isTRUE(J_eff_diagnostics$positive_definite)) {
    return(fail("J_eff_not_positive_definite"))
  }
  if (!isTRUE(J_eff_diagnostics$condition_acceptable)) {
    return(fail("J_eff_ill_conditioned"))
  }
  J_eff_solution <- tryCatch(
    solve(J_eff, U_eff),
    error = function(e) NULL
  )
  if (is.null(J_eff_solution) || any(!is.finite(J_eff_solution))) {
    return(fail("J_eff_solve_failed"))
  }
  statistic <- as.numeric(crossprod(U_eff, J_eff_solution))
  if (length(statistic) != 1L || !is.finite(statistic) || statistic < 0) {
    return(fail("nonfinite_or_negative_score_statistic"))
  }
  result$statistic <- statistic
  result$diagnostic_available <- TRUE
  result$failure_reason <- NA_character_
  result
}

ln_joint_ratio_endpoint_audit <- function(fit,
                                          endpoint = c(
                                            "unrestricted",
                                            "constrained"
                                          ),
                                          objective_context = NULL) {
  endpoint <- match.arg(endpoint)
  thresholds <- ln_joint_thresholds()
  fit_is_list <- is.list(fit)
  fit_route <- if (fit_is_list) fit$route else NULL
  fit_debug <- if (fit_is_list && is.list(fit$debug)) fit$debug else list()
  fit_gene_index <- tryCatch(
    suppressWarnings(as.integer(fit_debug$gene_index)),
    error = function(e) NA_integer_
  )
  if (length(fit_gene_index) != 1L || is.na(fit_gene_index)) {
    fit_gene_index <- NA_integer_
  }
  result <- list(
    success = FALSE,
    failure_reason = NA_character_,
    endpoint = endpoint,
    route = fit_route,
    inference_available = fit_is_list && isTRUE(fit$inference_available),
    gene_index = fit_gene_index,
    n_beta = NA_integer_,
    n_donors = NA_integer_,
    full_log = NULL,
    state_objective = NA_real_,
    stored_state_objective = NA_real_,
    reevaluated_objective = NA_real_,
    candidate_objective = NA_real_,
    loglik = NA_real_,
    loglik_ln = NA_real_,
    certified_candidate_index = NA_integer_,
    objective_provenance = NULL,
    objective_instance = NULL,
    diagnostics = list(
      candidate_objective_state_scaled = NA_real_,
      loglik_objective_scaled = NA_real_,
      loglik_ln_objective_scaled = NA_real_,
      score_sign_scaled = NA_real_,
      donor_additivity_scaled = NA_real_,
      objective_recheck_scaled = NA_real_,
      full_log_recheck_scaled = NA_real_,
      U_full_recheck_scaled = NA_real_,
      U_subject_full_recheck_scaled = NA_real_,
      A_full_recheck_scaled = NA_real_,
      U_reduced_recheck_scaled = NA_real_,
      U_subject_reduced_recheck_scaled = NA_real_,
      A_reduced_recheck_scaled = NA_real_,
      root_residual = NA_real_,
      newton_displacement = NA_real_,
      curvature = NULL
    )
  )
  fail <- function(reason) {
    result$failure_reason <- paste0(endpoint, "_", reason)
    result
  }

  if (!fit_is_list) {
    return(fail("fit_invalid"))
  }
  if (!identical(fit$route, "R0") || !isTRUE(fit$inference_available)) {
    return(fail("route_not_R0"))
  }
  expected_method <- if (identical(endpoint, "unrestricted")) {
    "joint_ln"
  } else {
    "joint_ln_constrained"
  }
  if (!identical(fit_debug$method, expected_method)) {
    return(fail("estimator_mismatch"))
  }
  if (!is.list(fit$ln_joint)) {
    return(fail("joint_state_invalid"))
  }
  state <- fit$ln_joint$state
  if (!is.list(state) || !is.list(state$raw) || is.null(state$beta) ||
      is.null(state$raw$value)) {
    return(fail("state_missing"))
  }
  nb <- length(state$beta)
  full_dim <- nb + 2L
  U_full <- as.numeric(state$U_full)
  U_subject_full <- as.matrix(state$U_subject_full)
  A_full <- as.matrix(state$A_full)
  result$n_beta <- as.integer(nb)
  result$n_donors <- as.integer(ncol(U_subject_full))
  state_dimensions_ok <- length(U_full) == full_dim &&
    nrow(U_subject_full) == full_dim && ncol(U_subject_full) > 0L &&
    all(dim(A_full) == c(full_dim, full_dim))
  if (!state_dimensions_ok || any(!is.finite(c(
    state$beta,
    U_full,
    U_subject_full,
    A_full,
    state$raw$value
  )))) {
    return(fail("state_invalid"))
  }
  U_reduced <- as.numeric(state$U_reduced)
  U_subject_reduced <- as.matrix(state$U_subject_reduced)
  A_reduced <- as.matrix(state$A_reduced)
  theta_reduced_log <- as.numeric(state$theta_reduced_log)
  reduced_dim <- length(U_reduced)
  reduced_dimensions_ok <- reduced_dim > 0L &&
    nrow(U_subject_reduced) == reduced_dim &&
    ncol(U_subject_reduced) == ncol(U_subject_full) &&
    all(dim(A_reduced) == c(reduced_dim, reduced_dim)) &&
    length(theta_reduced_log) == reduced_dim
  if (!reduced_dimensions_ok || any(!is.finite(c(
    U_reduced,
    U_subject_reduced,
    A_reduced,
    theta_reduced_log
  )))) {
    return(fail("reduced_state_invalid"))
  }
  result$diagnostics$root_residual <- ln_joint_root_metric(
    U_subject_reduced,
    U_reduced
  )
  result$diagnostics$curvature <- ln_joint_curvature_diagnostics(A_reduced)
  newton_step <- tryCatch(
    solve(A_reduced, U_reduced),
    error = function(e) NULL
  )
  result$diagnostics$newton_displacement <- if (!is.null(newton_step) &&
      all(is.finite(newton_step))) {
    max(abs(newton_step) / (1 + abs(theta_reduced_log)))
  } else {
    NA_real_
  }
  if (!is.finite(result$diagnostics$root_residual) ||
      result$diagnostics$root_residual > thresholds$root) {
    return(fail("selected_root_not_certified"))
  }
  if (!is.finite(result$diagnostics$newton_displacement) ||
      result$diagnostics$newton_displacement > thresholds$newton) {
    return(fail("selected_newton_displacement_not_certified"))
  }
  curvature <- result$diagnostics$curvature
  if (!isTRUE(curvature$positive_definite) ||
      !is.finite(curvature$scaled_min_eigenvalue) ||
      curvature$scaled_min_eigenvalue < thresholds$curvature_min_eigen) {
    return(fail("selected_curvature_not_certified"))
  }
  if (!is.finite(curvature$scaled_condition_number) ||
      curvature$scaled_condition_number > thresholds$condition_number) {
    return(fail("selected_curvature_condition_not_certified"))
  }

  expected_provenance <- ln_joint_objective_provenance()
  fit_provenance <- fit$ln_joint$objective_provenance
  state_provenance <- state$objective_provenance
  if (!identical(fit_provenance, expected_provenance) ||
      !identical(state_provenance, expected_provenance)) {
    return(fail("objective_provenance_mismatch"))
  }
  fit_instance <- fit$ln_joint$objective_instance
  state_instance <- state$objective_instance
  if (is.null(fit_instance) || is.null(state_instance)) {
    return(fail("objective_instance_missing"))
  }
  if (!identical(fit_instance, state_instance)) {
    return(fail("objective_instance_mismatch"))
  }
  if (!identical(
    fit$ln_joint$parameterization,
    expected_provenance$derivative_coordinates
  ) || !identical(fit$ln_joint$score_sign, "maximizing")) {
    return(fail("coordinate_or_score_convention_mismatch"))
  }
  result$objective_provenance <- expected_provenance
  result$objective_instance <- fit_instance

  state_full_log <- as.numeric(state$full_log)
  state_coordinates_match <- length(state_full_log) == full_dim &&
    identical(
      state_full_log[seq_len(nb)],
      as.numeric(state$beta)
    ) &&
    length(state$sigma2) == 1L && length(state$gamma) == 1L &&
    identical(as.numeric(state$sigma2), exp(state_full_log[nb + 1L])) &&
    identical(as.numeric(state$gamma), exp(state_full_log[nb + 2L]))
  if (!isTRUE(state_coordinates_match)) {
    return(fail("state_coordinate_mismatch"))
  }
  result$full_log <- state_full_log
  returned_theta <- fit$theta
  returned_state_matches <- is.list(returned_theta) &&
    identical(as.numeric(returned_theta$beta), as.numeric(state$beta)) &&
    identical(as.numeric(returned_theta$subVar), as.numeric(state$sigma2)) &&
    identical(as.numeric(returned_theta$cellVar), as.numeric(state$gamma))
  if (identical(endpoint, "constrained")) {
    returned_state_matches <- returned_state_matches &&
      identical(as.numeric(fit$beta_constrained), as.numeric(state$beta))
  }
  if (!isTRUE(returned_state_matches)) {
    return(fail("return_state_mismatch"))
  }

  if (!is.list(objective_context) ||
      !is.list(objective_context$ctx) ||
      !is.list(objective_context$posv) ||
      length(objective_context$gene_index) != 1L ||
      !is.finite(objective_context$gene_index) ||
      !identical(
        as.integer(objective_context$gene_index),
        result$gene_index
      )) {
    return(fail("objective_context_missing_or_mismatched"))
  }
  ctx <- objective_context$ctx
  posv <- objective_context$posv
  if (length(ctx$nb) != 1L || !is.finite(ctx$nb) ||
      as.integer(ctx$nb) != nb) {
    return(fail("objective_context_dimension_mismatch"))
  }
  endpoint_map <- tryCatch(
    if (identical(endpoint, "unrestricted")) {
      ln_joint_map(nb)
    } else {
      if (!is.list(fit$constraint) || is.null(fit$constraint$L) ||
          is.null(fit$constraint$b)) {
        stop("constraint map missing")
      }
      ln_joint_map(nb, fit$constraint$L, fit$constraint$b)
    },
    error = function(e) NULL
  )
  if (is.null(endpoint_map) ||
      length(theta_reduced_log) != endpoint_map$free_dim + 2L) {
    return(fail("objective_map_invalid"))
  }
  reevaluated_state <- tryCatch(
    ln_joint_derivative_state(
      theta_reduced_log,
      result$gene_index,
      posv,
      ctx,
      endpoint_map
    ),
    error = function(e) NULL
  )
  if (!is.list(reevaluated_state) || !is.list(reevaluated_state$raw) ||
      !is.finite(reevaluated_state$raw$value) ||
      !identical(
        reevaluated_state$objective_provenance,
        expected_provenance
      )) {
    return(fail("objective_reevaluation_failed"))
  }
  result$reevaluated_objective <- as.numeric(reevaluated_state$raw$value)
  result$stored_state_objective <- as.numeric(state$raw$value)
  result$diagnostics$objective_recheck_scaled <- ln_joint_scaled_max(
    result$stored_state_objective,
    result$reevaluated_objective
  )
  result$diagnostics$full_log_recheck_scaled <- ln_joint_scaled_max(
    state_full_log,
    reevaluated_state$full_log
  )
  result$diagnostics$U_full_recheck_scaled <- ln_joint_scaled_max(
    U_full,
    reevaluated_state$U_full
  )
  result$diagnostics$U_subject_full_recheck_scaled <- ln_joint_scaled_max(
    U_subject_full,
    reevaluated_state$U_subject_full
  )
  result$diagnostics$A_full_recheck_scaled <- ln_joint_scaled_max(
    A_full,
    reevaluated_state$A_full
  )
  result$diagnostics$U_reduced_recheck_scaled <- ln_joint_scaled_max(
    U_reduced,
    reevaluated_state$U_reduced
  )
  result$diagnostics$U_subject_reduced_recheck_scaled <-
    ln_joint_scaled_max(
      U_subject_reduced,
      reevaluated_state$U_subject_reduced
    )
  result$diagnostics$A_reduced_recheck_scaled <- ln_joint_scaled_max(
    A_reduced,
    reevaluated_state$A_reduced
  )
  recheck_diagnostics <- unlist(
    result$diagnostics[c(
      "full_log_recheck_scaled",
      "U_full_recheck_scaled",
      "U_subject_full_recheck_scaled",
      "A_full_recheck_scaled",
      "U_reduced_recheck_scaled",
      "U_subject_reduced_recheck_scaled",
      "A_reduced_recheck_scaled"
    )],
    use.names = FALSE
  )
  if (!is.finite(result$diagnostics$objective_recheck_scaled) ||
      result$diagnostics$objective_recheck_scaled >
        thresholds$objective_agreement) {
    return(fail("objective_reevaluation_mismatch"))
  }
  if (length(recheck_diagnostics) != 7L ||
      any(!is.finite(recheck_diagnostics)) ||
      result$diagnostics$full_log_recheck_scaled >
        thresholds$natural_log_coordinate ||
      any(recheck_diagnostics[-1L] > thresholds$covariance_identity)) {
    return(fail("derivative_state_reevaluation_mismatch"))
  }

  candidate_records <- fit$ln_joint$candidate_records
  if (!is.list(candidate_records) || !length(candidate_records) ||
      any(!vapply(candidate_records, is.list, logical(1)))) {
    return(fail("candidate_records_invalid"))
  }
  state_matches <- which(vapply(
    candidate_records,
    function(candidate) {
      isTRUE(candidate$certified) &&
        isTRUE(candidate$interior) &&
        !is.null(candidate$theta_full_log) &&
        identical(as.numeric(candidate$theta_full_log), state_full_log)
    },
    logical(1)
  ))
  if (!length(state_matches)) {
    return(fail("endpoint_certification_missing"))
  }
  provenance_matches <- state_matches[vapply(
    candidate_records[state_matches],
    function(candidate) {
      identical(candidate$objective_provenance, expected_provenance)
    },
    logical(1)
  )]
  if (!length(provenance_matches)) {
    return(fail("candidate_objective_provenance_mismatch"))
  }

  state_objective <- result$reevaluated_objective
  candidate_errors <- vapply(
    candidate_records[provenance_matches],
    function(candidate) {
      ln_joint_scaled_max(candidate$objective, state_objective)
    },
    numeric(1)
  )
  objective_matches <- provenance_matches[
    is.finite(candidate_errors) &
      candidate_errors <= thresholds$objective_agreement
  ]
  if (!length(objective_matches)) {
    result$diagnostics$candidate_objective_state_scaled <-
      if (any(is.finite(candidate_errors))) {
        min(candidate_errors[is.finite(candidate_errors)])
      } else {
        NA_real_
      }
    return(fail("objective_state_mismatch"))
  }
  candidate_index <- objective_matches[[1L]]
  candidate_objective <- as.numeric(
    candidate_records[[candidate_index]]$objective
  )
  result$state_objective <- state_objective
  result$candidate_objective <- candidate_objective
  result$certified_candidate_index <- as.integer(candidate_index)
  result$diagnostics$candidate_objective_state_scaled <-
    ln_joint_scaled_max(candidate_objective, state_objective)

  loglik <- as.numeric(fit$loglik)
  loglik_ln <- as.numeric(fit$loglik_ln)
  result$loglik <- loglik
  result$loglik_ln <- loglik_ln
  result$diagnostics$loglik_objective_scaled <-
    ln_joint_scaled_max(loglik, -state_objective)
  result$diagnostics$loglik_ln_objective_scaled <-
    ln_joint_scaled_max(loglik_ln, -state_objective)
  if (length(loglik) != 1L || length(loglik_ln) != 1L ||
      any(!is.finite(c(loglik, loglik_ln))) ||
      !is.finite(result$diagnostics$loglik_objective_scaled) ||
      !is.finite(result$diagnostics$loglik_ln_objective_scaled) ||
      result$diagnostics$loglik_objective_scaled >
        thresholds$objective_agreement ||
      result$diagnostics$loglik_ln_objective_scaled >
        thresholds$objective_agreement) {
    return(fail("loglik_objective_mismatch"))
  }

  result$diagnostics$score_sign_scaled <- ln_joint_scaled_max(
    U_full,
    -state$q_gradient_full
  )
  result$diagnostics$donor_additivity_scaled <- ln_joint_scaled_max(
    rowSums(U_subject_full),
    U_full
  )
  if (!all(is.finite(c(
    result$diagnostics$score_sign_scaled,
    result$diagnostics$donor_additivity_scaled
  ))) || result$diagnostics$score_sign_scaled >
      thresholds$donor_additivity) {
    return(fail("maximizing_score_sign_mismatch"))
  }
  if (result$diagnostics$donor_additivity_scaled >
      thresholds$donor_additivity) {
    return(fail("donor_additivity_failed"))
  }

  result$success <- TRUE
  result$failure_reason <- NA_character_
  result
}

ln_joint_pace_ratio_diagnostic <- function(unrestricted_fit,
                                           constrained_fit,
                                           null_basis = NULL,
                                           interest_basis = NULL) {
  empty_matrix_diagnostics <- ln_joint_score_matrix_diagnostics(NULL)
  result <- list(
    diagnostic_available = FALSE,
    failure_reason = NA_character_,
    public_available = FALSE,
    public_status = "joint_ln_calibration_pending",
    statistic = NA_real_,
    adjusted_ratio = NA_real_,
    df = NA_integer_,
    theoretical_df = NA_integer_,
    theoretical_reference = paste(
      "chi-square first-order reference only under smooth interior",
      "independent-donor fixed-dimension regularity"
    ),
    method = "pace_invariant_joint_ln_working_objective_ratio_audit",
    no_numerical_repair = TRUE,
    wP = NA_real_,
    wU = NA_real_,
    qP = NA_real_,
    multiplier = NA_real_,
    scale = NA_real_,
    U_eff = NULL,
    U_eff_subject = NULL,
    H_eff = NULL,
    J_eff = NULL,
    coordinate_map = NULL,
    projection = NULL,
    objective_provenance = NULL,
    objective_instance = NULL,
    endpoints = list(unrestricted = NULL, constrained = NULL),
    score_diagnostic = NULL,
    diagnostics = list(
      same_criterion = FALSE,
      same_objective_instance = FALSE,
      same_objective = FALSE,
      same_gene = FALSE,
      same_dimension = FALSE,
      same_donor_count = FALSE,
      wP_loglik_identity_scaled = NA_real_,
      H_eff_projection_identity_scaled = NA_real_,
      wU_score_identity_scaled = NA_real_,
      wP_minus_qP = NA_real_,
      wP_over_qP = NA_real_,
      adjusted_minus_wU = NA_real_,
      pace_decomposition_identity_scaled = NA_real_,
      H_eff = empty_matrix_diagnostics,
      J_eff = empty_matrix_diagnostics
    )
  )
  fail <- function(reason) {
    result$failure_reason <- reason
    result
  }

  if (!is.list(unrestricted_fit)) {
    return(fail("unrestricted_fit_invalid"))
  }
  if (!length(unrestricted_fit)) {
    return(fail("unrestricted_pair_not_registered"))
  }
  if (!is.list(constrained_fit) || !length(constrained_fit)) {
    return(fail("constrained_fit_invalid"))
  }
  fit_instance <- function(fit) {
    if (!is.list(fit) || !is.list(fit$ln_joint)) {
      return(NULL)
    }
    fit$ln_joint$objective_instance
  }
  unrestricted_instance <- fit_instance(unrestricted_fit)
  constrained_instance <- fit_instance(constrained_fit)
  if (!is.null(unrestricted_instance) && !is.null(constrained_instance) &&
      !identical(unrestricted_instance, constrained_instance)) {
    return(fail("endpoint_objective_instance_mismatch"))
  }
  objective_context <- if (!is.null(unrestricted_instance) &&
      identical(unrestricted_instance, constrained_instance)) {
    ln_joint_objective_instance_registry(
      "resolve",
      instance = unrestricted_instance
    )
  } else {
    NULL
  }

  unrestricted_endpoint <- ln_joint_ratio_endpoint_audit(
    unrestricted_fit,
    "unrestricted",
    objective_context = objective_context
  )
  constrained_endpoint <- ln_joint_ratio_endpoint_audit(
    constrained_fit,
    "constrained",
    objective_context = objective_context
  )
  result$endpoints$unrestricted <- unrestricted_endpoint
  result$endpoints$constrained <- constrained_endpoint
  if (!isTRUE(unrestricted_endpoint$success)) {
    return(fail(unrestricted_endpoint$failure_reason))
  }
  if (!isTRUE(constrained_endpoint$success)) {
    return(fail(constrained_endpoint$failure_reason))
  }

  result$diagnostics$same_criterion <- identical(
    unrestricted_endpoint$objective_provenance,
    constrained_endpoint$objective_provenance
  )
  result$diagnostics$same_objective_instance <- identical(
    unrestricted_endpoint$objective_instance,
    constrained_endpoint$objective_instance
  )
  result$diagnostics$same_objective <-
    isTRUE(result$diagnostics$same_criterion) &&
    isTRUE(result$diagnostics$same_objective_instance)
  result$diagnostics$same_gene <- identical(
    unrestricted_endpoint$gene_index,
    constrained_endpoint$gene_index
  ) && !is.na(unrestricted_endpoint$gene_index)
  result$diagnostics$same_dimension <- identical(
    unrestricted_endpoint$n_beta,
    constrained_endpoint$n_beta
  )
  result$diagnostics$same_donor_count <- identical(
    unrestricted_endpoint$n_donors,
    constrained_endpoint$n_donors
  )
  if (!isTRUE(result$diagnostics$same_criterion)) {
    return(fail("endpoint_criterion_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_objective_instance)) {
    return(fail("endpoint_objective_instance_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_gene)) {
    return(fail("endpoint_gene_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_dimension)) {
    return(fail("endpoint_dimension_mismatch"))
  }
  if (!isTRUE(result$diagnostics$same_donor_count)) {
    return(fail("endpoint_donor_count_mismatch"))
  }
  result$objective_provenance <-
    unrestricted_endpoint$objective_provenance
  result$objective_instance <- unrestricted_endpoint$objective_instance

  wP <- 2 * (
    constrained_endpoint$state_objective -
      unrestricted_endpoint$state_objective
  )
  wP_loglik <- 2 * (
    unrestricted_endpoint$loglik_ln - constrained_endpoint$loglik_ln
  )
  result$wP <- wP
  result$diagnostics$wP_loglik_identity_scaled <- ln_joint_scaled_max(
    wP,
    wP_loglik
  )
  if (!is.finite(wP) ||
      !is.finite(result$diagnostics$wP_loglik_identity_scaled) ||
      result$diagnostics$wP_loglik_identity_scaled >
        ln_joint_thresholds()$objective_agreement) {
    return(fail("working_objective_ratio_identity_failed"))
  }
  if (wP < 0) {
    return(fail("negative_working_objective_ratio"))
  }

  score <- ln_joint_robust_score_diagnostic(
    constrained_fit,
    null_basis = null_basis,
    interest_basis = interest_basis
  )
  result$score_diagnostic <- score
  result$df <- result$theoretical_df <- score$df
  if (!isTRUE(score$diagnostic_available)) {
    return(fail(paste0(
      "constrained_score_",
      if (is.null(score$failure_reason) || is.na(score$failure_reason)) {
        "unavailable"
      } else {
        score$failure_reason
      }
    )))
  }
  result$coordinate_map <- score$coordinate_map
  result$projection <- score$projection
  result$U_eff <- score$U_eff
  result$U_eff_subject <- score$U_eff_subject
  result$J_eff <- score$J_eff
  result$diagnostics$J_eff <- score$diagnostics$J_eff

  state <- constrained_fit$ln_joint$state
  T <- score$coordinate_map$full
  H_coordinate <- crossprod(T, state$A_full %*% T)
  H_coordinate <- 0.5 * (H_coordinate + t(H_coordinate))
  r <- score$df
  interest <- seq_len(r)
  nuisance <- r + seq_len(nrow(H_coordinate) - r)
  H_psi_psi <- H_coordinate[interest, interest, drop = FALSE]
  H_chi_psi <- H_coordinate[nuisance, interest, drop = FALSE]
  H_chi_solution <- tryCatch(
    solve(score$H_chi_chi, H_chi_psi),
    error = function(e) NULL
  )
  if (is.null(H_chi_solution) || any(!is.finite(H_chi_solution))) {
    return(fail("H_chi_chi_solve_failed_for_H_eff"))
  }
  H_eff_schur <- H_psi_psi - score$H_psi_chi %*% H_chi_solution
  H_eff_projected <- score$projection %*% H_coordinate %*%
    t(score$projection)
  result$diagnostics$H_eff_projection_identity_scaled <-
    ln_joint_relative_frobenius(H_eff_schur, H_eff_projected)
  if (!is.finite(
    result$diagnostics$H_eff_projection_identity_scaled
  ) || result$diagnostics$H_eff_projection_identity_scaled >
      ln_joint_thresholds()$covariance_identity) {
    return(fail("H_eff_projection_identity_failed"))
  }
  H_eff <- 0.5 * (H_eff_schur + t(H_eff_schur))
  H_eff_diagnostics <- ln_joint_score_matrix_diagnostics(H_eff)
  result$H_eff <- H_eff
  result$diagnostics$H_eff <- H_eff_diagnostics
  if (!isTRUE(H_eff_diagnostics$symmetry_acceptable)) {
    return(fail("H_eff_asymmetric"))
  }
  if (!isTRUE(H_eff_diagnostics$full_rank)) {
    return(fail("H_eff_rank_deficient"))
  }
  if (!isTRUE(H_eff_diagnostics$positive_definite)) {
    return(fail("H_eff_not_positive_definite"))
  }
  if (!isTRUE(H_eff_diagnostics$condition_acceptable)) {
    return(fail("H_eff_ill_conditioned"))
  }

  qP_solution <- tryCatch(
    solve(H_eff, score$U_eff),
    error = function(e) NULL
  )
  wU_solution <- tryCatch(
    solve(score$J_eff, score$U_eff),
    error = function(e) NULL
  )
  if (is.null(qP_solution) || any(!is.finite(qP_solution))) {
    return(fail("H_eff_solve_failed"))
  }
  if (is.null(wU_solution) || any(!is.finite(wU_solution))) {
    return(fail("J_eff_solve_failed"))
  }
  qP <- as.numeric(crossprod(score$U_eff, qP_solution))
  wU_reconstructed <- as.numeric(crossprod(
    score$U_eff,
    wU_solution
  ))
  wU <- as.numeric(score$statistic)
  result$qP <- qP
  result$wU <- wU
  result$diagnostics$wU_score_identity_scaled <- ln_joint_scaled_max(
    wU,
    wU_reconstructed
  )
  if (!is.finite(result$diagnostics$wU_score_identity_scaled) ||
      result$diagnostics$wU_score_identity_scaled >
        ln_joint_thresholds()$covariance_identity) {
    return(fail("wU_score_identity_failed"))
  }
  if (!is.finite(qP) || qP <= .Machine$double.eps) {
    return(fail("qP_nonpositive_or_numerically_zero"))
  }
  if (!is.finite(wU) || wU <= .Machine$double.eps) {
    return(fail("wU_nonpositive_or_numerically_zero"))
  }

  multiplier <- wU / qP
  adjusted_ratio <- multiplier * wP
  result$diagnostics$wP_minus_qP <- wP - qP
  result$diagnostics$wP_over_qP <- wP / qP
  result$diagnostics$adjusted_minus_wU <- adjusted_ratio - wU
  pace_decomposition <- wU + multiplier * (wP - qP)
  result$diagnostics$pace_decomposition_identity_scaled <-
    ln_joint_scaled_max(adjusted_ratio, pace_decomposition)
  if (!is.finite(multiplier) || multiplier <= 0) {
    return(fail("nonpositive_pace_multiplier"))
  }
  if (!is.finite(adjusted_ratio) || adjusted_ratio < 0) {
    return(fail("invalid_adjusted_working_objective_ratio"))
  }
  if (!is.finite(
    result$diagnostics$pace_decomposition_identity_scaled
  ) || result$diagnostics$pace_decomposition_identity_scaled >
      ln_joint_thresholds()$covariance_identity) {
    return(fail("pace_decomposition_identity_failed"))
  }

  result$multiplier <- result$scale <- multiplier
  result$adjusted_ratio <- result$statistic <- adjusted_ratio
  result$diagnostic_available <- TRUE
  result$failure_reason <- NA_character_
  result
}

ln_joint_root_metric <- function(U_subject, U) {
  if (is.null(U_subject) || !is.matrix(U_subject) || any(!is.finite(U_subject)) ||
      any(!is.finite(U))) {
    return(NA_real_)
  }
  denom <- 1 + sqrt(rowSums(U_subject^2))
  max(abs(U) / denom)
}

ln_joint_fd_gradient <- function(fn, theta, step_multiplier = 1) {
  theta <- as.numeric(theta)
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  vapply(seq_along(theta), function(jj) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    (fn(plus) - fn(minus)) / (2 * step[jj])
  }, numeric(1))
}

ln_joint_fd_jacobian <- function(fn, theta, step_multiplier = 1) {
  theta <- as.numeric(theta)
  base <- as.numeric(fn(theta))
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  out <- matrix(NA_real_, nrow = length(base), ncol = length(theta))
  for (jj in seq_along(theta)) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    out[, jj] <- (as.numeric(fn(plus)) - as.numeric(fn(minus))) /
      (2 * step[jj])
  }
  out
}

ln_joint_donor_objective_components <- function(theta_reduced_log, gene_index,
                                                 posv, ctx, map) {
  free_dim <- map$free_dim
  theta_reduced <- c(
    if (free_dim) theta_reduced_log[seq_len(free_dim)] else numeric(0),
    exp(theta_reduced_log[free_dim + seq_len(2L)])
  )
  full <- ln_joint_full_natural(theta_reduced, map)
  beta <- full[seq_len(ctx$nb)]
  sigma2 <- full[ctx$nb + 1L]
  gamma <- full[ctx$nb + 2L]
  exp_sigma2 <- exp(sigma2)
  alpha <- 1 / (exp_sigma2 - 1)
  lambda <- 1 / (sqrt(exp_sigma2) * (exp_sigma2 - 1))
  components <- lapply(seq_len(ctx$k), function(subject) {
    first_zero <- ctx$fid[subject] - 1L
    last_zero <- ctx$fid[subject + 1L] - 2L
    cell_index <- seq.int(first_zero + 1L, last_zero + 1L)
    positive <- which(
      posv$posindy >= first_zero & posv$posindy <= last_zero
    )
    donor_y <- posv$Y[positive]
    donor_posindy <- posv$posindy[positive] - first_zero
    donor_sum <- sum(donor_y)
    donor_n <- length(cell_index)
    eta <- as.numeric(
      ctx$offset[cell_index] +
        ctx$pred[cell_index, , drop = FALSE] %*% beta
    )
    mu <- exp(eta)
    sum_mu <- sum(mu)
    total <- sum_mu + lambda
    donor_weight <- (donor_sum + alpha) / total
    log_cell_total <- log(donor_weight * mu + gamma)
    fn_components <- c(
      positive_linear = if (length(positive))
        sum(eta[donor_posindy + 1L] * donor_y) else 0,
      donor_log_total = -(donor_sum + alpha) * log(total),
      alpha_log_lambda = alpha * log(lambda),
      cell_gamma_normalizer = donor_n * gamma * log(gamma),
      donor_mean_term = donor_weight * sum_mu,
      cell_log_normalizer = -gamma * sum(log_cell_total),
      positive_cell_log = if (length(positive))
        -sum(donor_y * log_cell_total[donor_posindy + 1L]) else 0,
      donor_lgamma = if (donor_sum > 0)
        lgamma(donor_sum + alpha) - lgamma(alpha) else 0,
      cell_lgamma = if (length(positive))
        sum(lgamma(donor_y + gamma) - lgamma(gamma)) else 0
    )
    -fn_components
  })
  out <- do.call(rbind, components)
  rownames(out) <- paste0("donor_", seq_len(ctx$k))
  out
}

ln_joint_donor_objective_values <- function(theta_reduced_log, gene_index,
                                             posv, ctx, map) {
  rowSums(ln_joint_donor_objective_components(
    theta_reduced_log,
    gene_index,
    posv,
    ctx,
    map
  ))
}

ln_joint_fd_gradient_donor_matrix <- function(fn, theta,
                                              step_multiplier = 1) {
  theta <- as.numeric(theta)
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  donor_count <- nrow(as.matrix(fn(theta)))
  out <- vapply(seq_along(theta), function(jj) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    rowSums(
      (as.matrix(fn(plus)) - as.matrix(fn(minus))) / (2 * step[jj])
    )
  }, numeric(donor_count))
  t(out)
}

ln_joint_fd_gradient_donorwise <- function(fn, theta, step_multiplier = 1) {
  rowSums(ln_joint_fd_gradient_donor_matrix(fn, theta, step_multiplier))
}

ln_joint_fd_sensitivity_donorwise <- function(fn, theta,
                                               step_multiplier = 1) {
  theta <- as.numeric(theta)
  base <- as.matrix(fn(theta))
  step <- step_multiplier * .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  out <- matrix(NA_real_, nrow = nrow(base), ncol = length(theta))
  for (jj in seq_along(theta)) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    donor_jacobian <- (as.matrix(fn(plus)) - as.matrix(fn(minus))) /
      (2 * step[jj])
    out[, jj] <- -rowSums(donor_jacobian)
  }
  out
}

ln_joint_fd_audit <- function(theta_reduced_log, gene_index, posv, ctx, map) {
  thresholds <- ln_joint_thresholds()
  free_dim <- map$free_dim
  to_natural_reduced <- function(z) {
    c(
      if (free_dim) z[seq_len(free_dim)] else numeric(0),
      exp(z[free_dim + seq_len(2L)])
    )
  }
  objective_log <- function(z) {
    ln_joint_objective(to_natural_reduced(z), gene_index, posv, ctx, map)
  }
  q_gradient_log <- function(z) {
    ln_joint_derivative_state(z, gene_index, posv, ctx, map)$q_gradient_reduced
  }
  maximizing_score_log <- function(z) {
    ln_joint_derivative_state(z, gene_index, posv, ctx, map)$U_reduced
  }
  donor_objective_log <- function(z) {
    ln_joint_donor_objective_components(z, gene_index, posv, ctx, map)
  }
  donor_score_log <- function(z) {
    ln_joint_derivative_state(z, gene_index, posv, ctx, map)$U_subject_reduced
  }

  state <- ln_joint_derivative_state(
    theta_reduced_log,
    gene_index,
    posv,
    ctx,
    map
  )
  fd_gradient_global_h <- ln_joint_fd_gradient(
    objective_log,
    theta_reduced_log,
    1
  )
  fd_gradient_global_h2 <- ln_joint_fd_gradient(
    objective_log,
    theta_reduced_log,
    0.5
  )
  fd_donor_gradient_h <- ln_joint_fd_gradient_donor_matrix(
    donor_objective_log,
    theta_reduced_log,
    1
  )
  fd_donor_gradient_h2 <- ln_joint_fd_gradient_donor_matrix(
    donor_objective_log,
    theta_reduced_log,
    0.5
  )
  fd_gradient_h <- rowSums(fd_donor_gradient_h)
  fd_gradient_h2 <- rowSums(fd_donor_gradient_h2)
  fd_sensitivity_global_h <- -ln_joint_fd_jacobian(
    maximizing_score_log,
    theta_reduced_log,
    1
  )
  fd_sensitivity_global_h2 <- -ln_joint_fd_jacobian(
    maximizing_score_log,
    theta_reduced_log,
    0.5
  )
  fd_sensitivity_h <- ln_joint_fd_sensitivity_donorwise(
    donor_score_log,
    theta_reduced_log,
    1
  )
  fd_sensitivity_h2 <- ln_joint_fd_sensitivity_donorwise(
    donor_score_log,
    theta_reduced_log,
    0.5
  )

  full_natural <- c(state$beta, state$sigma2, state$gamma)
  natural_gradient <- function(theta) {
    ptmg_der(
      theta,
      X = ctx$pred,
      offset = ctx$offset,
      Y = posv$Y,
      n_one = posv$n_onetwo,
      ytwo = posv$ytwo,
      fam = ctx$id,
      fid = ctx$fid,
      cumsumy = ctx$cumsumy[gene_index, ],
      posind = ctx$posind[[gene_index]],
      posindy = posv$posindy,
      nb = ctx$nb,
      nind = ctx$nind,
      k = ctx$k
    )
  }
  q_gradient_natural <- natural_gradient(full_natural)
  natural_sensitivity_fd <- ln_joint_fd_jacobian(
    natural_gradient,
    full_natural,
    1
  )
  natural_sensitivity_fd_h2 <- ln_joint_fd_jacobian(
    natural_gradient,
    full_natural,
    0.5
  )
  natural_sensitivity <- state$raw$hessian_natural
  q_gradient_natural_analytic <- state$raw$gradient_natural
  J <- diag(c(rep(1, ctx$nb), state$sigma2, state$gamma))
  q_gradient_log_from_natural <- as.numeric(crossprod(
    J,
    q_gradient_natural
  ))
  q_hessian_log_from_natural <- crossprod(J, natural_sensitivity %*% J)
  q_hessian_log_from_natural[ctx$nb + 1L, ctx$nb + 1L] <-
    q_hessian_log_from_natural[ctx$nb + 1L, ctx$nb + 1L] +
    q_gradient_natural[ctx$nb + 1L] * state$sigma2
  q_hessian_log_from_natural[ctx$nb + 2L, ctx$nb + 2L] <-
    q_hessian_log_from_natural[ctx$nb + 2L, ctx$nb + 2L] +
    q_gradient_natural[ctx$nb + 2L] * state$gamma

  metrics <- c(
    objective_value = abs(objective_log(theta_reduced_log) - state$raw$value) /
      (1 + max(abs(objective_log(theta_reduced_log)), abs(state$raw$value))),
    objective_additivity = abs(
      sum(donor_objective_log(theta_reduced_log)) -
        objective_log(theta_reduced_log)
    ) / (1 + abs(objective_log(theta_reduced_log))),
    objective_gradient = ln_joint_scaled_max(
      state$q_gradient_reduced,
      fd_gradient_h
    ),
    donor_objective_gradient = ln_joint_scaled_max(
      -state$U_subject_reduced,
      fd_donor_gradient_h
    ),
    sensitivity_score = ln_joint_scaled_max(
      state$A_reduced,
      fd_sensitivity_h
    ),
    fd_gradient_stability = ln_joint_scaled_max(fd_gradient_h, fd_gradient_h2),
    fd_donor_gradient_stability = ln_joint_scaled_max(
      fd_donor_gradient_h,
      fd_donor_gradient_h2
    ),
    fd_sensitivity_stability = ln_joint_scaled_max(
      fd_sensitivity_h,
      fd_sensitivity_h2
    ),
    donor_additivity = ln_joint_scaled_max(
      rowSums(state$U_subject_reduced),
      state$U_reduced
    ),
    aggregate_formula_vs_donor_sum = ln_joint_scaled_max(
      state$raw$gradient_aggregate_formula,
      state$q_gradient_full
    ),
    natural_log_gradient = ln_joint_scaled_max(
      state$q_gradient_full,
      q_gradient_log_from_natural
    ),
    natural_log_sensitivity = ln_joint_scaled_max(
      state$A_full,
      q_hessian_log_from_natural
    ),
    natural_gradient_implementation = ln_joint_scaled_max(
      q_gradient_natural_analytic,
      q_gradient_natural
    ),
    natural_sensitivity_score = ln_joint_scaled_max(
      natural_sensitivity,
      natural_sensitivity_fd
    ),
    fd_natural_sensitivity_stability = ln_joint_scaled_max(
      natural_sensitivity_fd,
      natural_sensitivity_fd_h2
    )
  )
  target <- c(
    objective_value = thresholds$objective_agreement,
    objective_additivity = thresholds$objective_agreement,
    objective_gradient = thresholds$objective_gradient,
    donor_objective_gradient = thresholds$objective_gradient,
    sensitivity_score = thresholds$sensitivity_score,
    fd_gradient_stability = thresholds$finite_difference_stability,
    fd_donor_gradient_stability = thresholds$finite_difference_stability,
    fd_sensitivity_stability = thresholds$finite_difference_stability,
    donor_additivity = thresholds$donor_additivity,
    aggregate_formula_vs_donor_sum = thresholds$donor_additivity,
    natural_log_gradient = thresholds$natural_log_coordinate,
    natural_log_sensitivity = thresholds$natural_log_coordinate,
    natural_gradient_implementation = thresholds$objective_gradient,
    natural_sensitivity_score = thresholds$sensitivity_score,
    fd_natural_sensitivity_stability =
      thresholds$finite_difference_stability
  )
  hard <- c(
    objective_value = thresholds$hard_stop$objective_agreement,
    objective_additivity = thresholds$hard_stop$objective_agreement,
    objective_gradient = thresholds$hard_stop$objective_gradient,
    donor_objective_gradient = thresholds$hard_stop$objective_gradient,
    sensitivity_score = thresholds$hard_stop$sensitivity_score,
    fd_gradient_stability = thresholds$hard_stop$finite_difference_stability,
    fd_donor_gradient_stability =
      thresholds$hard_stop$finite_difference_stability,
    fd_sensitivity_stability = thresholds$hard_stop$finite_difference_stability,
    donor_additivity = thresholds$hard_stop$donor_additivity,
    aggregate_formula_vs_donor_sum = thresholds$hard_stop$donor_additivity,
    natural_log_gradient = thresholds$hard_stop$natural_log_coordinate,
    natural_log_sensitivity = thresholds$hard_stop$natural_log_coordinate,
    natural_gradient_implementation = thresholds$hard_stop$objective_gradient,
    natural_sensitivity_score = thresholds$hard_stop$sensitivity_score,
    fd_natural_sensitivity_stability =
      thresholds$hard_stop$finite_difference_stability
  )
  list(
    metrics = metrics,
    roundoff_metrics = c(
      global_objective_gradient = ln_joint_scaled_max(
        state$q_gradient_reduced,
        fd_gradient_global_h
      ),
      global_objective_gradient_stability = ln_joint_scaled_max(
        fd_gradient_global_h,
        fd_gradient_global_h2
      ),
      global_sensitivity_score = ln_joint_scaled_max(
        state$A_reduced,
        fd_sensitivity_global_h
      ),
      global_sensitivity_stability = ln_joint_scaled_max(
        fd_sensitivity_global_h,
        fd_sensitivity_global_h2
      )
    ),
    target = target,
    hard_stop = hard,
    target_pass = is.finite(metrics) & metrics <= target,
    hard_stop_pass = is.finite(metrics) & metrics <= hard,
    fd_gradient_h = fd_gradient_h,
    fd_gradient_h2 = fd_gradient_h2,
    fd_donor_gradient_h = fd_donor_gradient_h,
    fd_donor_gradient_h2 = fd_donor_gradient_h2,
    fd_gradient_global_h = fd_gradient_global_h,
    fd_gradient_global_h2 = fd_gradient_global_h2,
    fd_sensitivity_h = fd_sensitivity_h,
    fd_sensitivity_h2 = fd_sensitivity_h2,
    fd_sensitivity_global_h = fd_sensitivity_global_h,
    fd_sensitivity_global_h2 = fd_sensitivity_global_h2,
    analytic_q_gradient = state$q_gradient_reduced,
    analytic_sensitivity = state$A_reduced,
    natural_q_gradient = q_gradient_natural,
    natural_sensitivity = natural_sensitivity,
    natural_sensitivity_fd = natural_sensitivity_fd,
    natural_sensitivity_fd_h2 = natural_sensitivity_fd_h2,
    q_gradient_log_from_natural = q_gradient_log_from_natural,
    q_hessian_log_from_natural = q_hessian_log_from_natural
  )
}

ln_joint_constraint_metric <- function(beta, map) {
  if (!isTRUE(map$constrained)) {
    return(0)
  }
  L <- as.matrix(map$L)
  b <- as.numeric(map$b)
  residual <- as.numeric(L %*% beta - b)
  denom <- 1 + max(abs(b)) + max(rowSums(abs(L))) * max(abs(beta))
  max(abs(residual)) / denom
}

ln_joint_candidate_state <- function(theta_reduced, objective, optimizer, label,
                                     start_provenance, required, core,
                                     gene_index, posv, ctx, map) {
  thresholds <- ln_joint_thresholds()
  free_dim <- map$free_dim
  finite_endpoint <- length(theta_reduced) == free_dim + 2L &&
    all(is.finite(theta_reduced)) && is.finite(objective) &&
    all(theta_reduced[free_dim + seq_len(2L)] > 0)
  if (!finite_endpoint) {
    return(list(
      label = label,
      start_provenance = start_provenance,
      required = required,
      core = core,
      finite_endpoint = FALSE,
      certified = FALSE,
      objective = objective,
      optimizer = optimizer,
      failure = "nonfinite_endpoint"
    ))
  }
  theta_log <- c(
    if (free_dim) theta_reduced[seq_len(free_dim)] else numeric(0),
    log(theta_reduced[free_dim + seq_len(2L)])
  )
  state <- tryCatch(
    ln_joint_derivative_state(theta_log, gene_index, posv, ctx, map),
    error = function(e) NULL
  )
  if (is.null(state)) {
    return(list(
      label = label,
      start_provenance = start_provenance,
      required = required,
      core = core,
      finite_endpoint = TRUE,
      certified = FALSE,
      theta_reduced = theta_reduced,
      objective = objective,
      optimizer = optimizer,
      failure = "derivative_evaluation_failed"
    ))
  }
  donor_sum <- rowSums(state$U_subject_reduced)
  donor_additivity <- ln_joint_scaled_max(donor_sum, state$U_reduced)
  root <- ln_joint_root_metric(state$U_subject_reduced, state$U_reduced)
  curvature <- ln_joint_curvature_diagnostics(state$A_reduced)
  newton_step <- tryCatch(
    solve(state$A_reduced, state$U_reduced),
    error = function(e) rep(NA_real_, length(state$U_reduced))
  )
  newton <- if (all(is.finite(newton_step))) {
    max(abs(newton_step) / (1 + abs(theta_log)))
  } else {
    NA_real_
  }
  dispersion <- theta_reduced[free_dim + seq_len(2L)]
  boundary_tolerance <- sqrt(.Machine$double.eps) *
    (1 + pmax(abs(ctx$min), abs(ctx$max)))
  interior <- all(
    dispersion > ctx$min + boundary_tolerance &
      dispersion < ctx$max - boundary_tolerance
  )
  constraint <- ln_joint_constraint_metric(state$beta, map)
  checks <- c(
    donor_additivity = is.finite(donor_additivity) &&
      donor_additivity <= thresholds$donor_additivity,
    root = is.finite(root) && root <= thresholds$root,
    newton = is.finite(newton) && newton <= thresholds$newton,
    curvature = isTRUE(curvature$positive_definite) &&
      is.finite(curvature$scaled_min_eigenvalue) &&
      curvature$scaled_min_eigenvalue >= thresholds$curvature_min_eigen,
    condition = is.finite(curvature$scaled_condition_number) &&
      curvature$scaled_condition_number <= thresholds$condition_number,
    interior = interior,
    constraint = is.finite(constraint) && constraint <= thresholds$constraint
  )
  failure <- NA_character_
  if (!checks[["donor_additivity"]]) {
    failure <- if (is.finite(donor_additivity) &&
        donor_additivity <= thresholds$hard_stop$donor_additivity) {
      "donor_additivity_target_failed"
    } else {
      "donor_additivity_hard_stop"
    }
  } else if (!checks[["interior"]]) {
    failure <- "dispersion_boundary"
  } else if (!checks[["root"]]) {
    failure <- if (is.finite(root) && root <= thresholds$hard_stop$root) {
      "root_residual_review_zone"
    } else {
      "root_residual_hard_stop"
    }
  } else if (!checks[["newton"]]) {
    failure <- if (is.finite(newton) && newton <= thresholds$hard_stop$newton) {
      "newton_displacement_review_zone"
    } else {
      "newton_displacement_hard_stop"
    }
  } else if (!isTRUE(curvature$positive_definite) ||
      !is.finite(curvature$scaled_min_eigenvalue) ||
      curvature$scaled_min_eigenvalue <= 0) {
    failure <- "nonpositive_joint_curvature"
  } else if (!checks[["curvature"]]) {
    failure <- "curvature_review_zone"
  } else if (!checks[["condition"]]) {
    failure <- if (is.finite(curvature$scaled_condition_number) &&
        curvature$scaled_condition_number <= thresholds$hard_stop$condition_number) {
      "curvature_condition_review_zone"
    } else {
      "curvature_condition_unusable"
    }
  } else if (!checks[["constraint"]]) {
    failure <- if (is.finite(constraint) &&
        constraint <= thresholds$hard_stop$constraint) {
      "constraint_review_zone"
    } else {
      "constraint_hard_stop"
    }
  }
  list(
    label = label,
    start_provenance = start_provenance,
    required = required,
    core = core,
    finite_endpoint = TRUE,
    certified = all(checks),
    theta_reduced = theta_reduced,
    theta_reduced_log = theta_log,
    theta_full_log = state$full_log,
    objective_provenance = state$objective_provenance,
    objective = objective,
    optimizer = optimizer,
    state = state,
    donor_additivity_scaled = donor_additivity,
    root_residual = root,
    newton_displacement = newton,
    curvature = curvature,
    interior = interior,
    constraint_residual = constraint,
    checks = checks,
    failure = failure
  )
}

ln_joint_polish <- function(theta_reduced, gene_index, posv, ctx, map,
                            lower, upper, max_iter = 25L) {
  free_dim <- map$free_dim
  z <- c(
    if (free_dim) theta_reduced[seq_len(free_dim)] else numeric(0),
    log(theta_reduced[free_dim + seq_len(2L)])
  )
  lower_log <- c(if (free_dim) lower[seq_len(free_dim)] else numeric(0), log(ctx$min))
  upper_log <- c(if (free_dim) upper[seq_len(free_dim)] else numeric(0), log(ctx$max))
  objective_at <- function(zz) {
    natural <- c(
      if (free_dim) zz[seq_len(free_dim)] else numeric(0),
      exp(zz[free_dim + seq_len(2L)])
    )
    ln_joint_objective(natural, gene_index, posv, ctx, map)
  }
  current_obj <- objective_at(z)
  for (iter in seq_len(max_iter)) {
    state <- tryCatch(
      ln_joint_derivative_state(z, gene_index, posv, ctx, map),
      error = function(e) NULL
    )
    if (is.null(state)) break
    metric <- ln_joint_root_metric(state$U_subject_reduced, state$U_reduced)
    if (is.finite(metric) && metric <= ln_joint_thresholds()$root / 10) break
    step <- tryCatch(
      solve(state$A_reduced, state$q_gradient_reduced),
      error = function(e) NULL
    )
    if (is.null(step) || any(!is.finite(step))) break
    accepted <- FALSE
    for (half in 0:20) {
      candidate <- z - step / (2^half)
      if (any(candidate <= lower_log) || any(candidate >= upper_log)) next
      candidate_obj <- tryCatch(objective_at(candidate), error = function(e) NA_real_)
      candidate_state <- tryCatch(
        ln_joint_derivative_state(candidate, gene_index, posv, ctx, map),
        error = function(e) NULL
      )
      if (!is.finite(candidate_obj) || is.null(candidate_state)) next
      candidate_metric <- ln_joint_root_metric(
        candidate_state$U_subject_reduced,
        candidate_state$U_reduced
      )
      objective_ok <- candidate_obj <= current_obj +
        1e-12 * (1 + abs(current_obj))
      if (objective_ok && is.finite(candidate_metric) &&
          (!is.finite(metric) || candidate_metric < metric)) {
        z <- candidate
        current_obj <- candidate_obj
        accepted <- TRUE
        break
      }
    }
    if (!accepted) break
  }
  list(
    theta_reduced = c(
      if (free_dim) z[seq_len(free_dim)] else numeric(0),
      exp(z[free_dim + seq_len(2L)])
    ),
    objective = current_obj
  )
}

ln_joint_clamp <- function(x, lower, upper) {
  gap <- 1e-6 * (upper - lower)
  pmin(pmax(x, lower + gap), upper - gap)
}

ln_joint_pql_start <- function(best, gene_index, posv, ctx, map) {
  full <- ln_joint_full_natural(best$theta_reduced, map)
  beta <- full[seq_len(ctx$nb)]
  sigma2 <- full[ctx$nb + 1L]
  gamma <- full[ctx$nb + 2L]
  gni <- (ctx$nind / ctx$k) * gamma
  cv2p <- if (ctx$ncell > 0L) {
    get_cv(ctx$offset, ctx$pred, beta, ctx$cell_ind, ctx$ncell, ctx$nind)
  } else {
    ctx$cv2
  }
  kappa_obs <- gni / (1 + cv2p)
  pair_trigger <- gni < ctx$cutoff_cell
  sigma_trigger <- !pair_trigger &&
    ((kappa_obs < 20) || ((kappa_obs < ctx$kappa) &&
      (sigma2 < (8 / max(kappa_obs, 1e-8)))))
  if (!pair_trigger && !sigma_trigger) {
    return(list(triggered = FALSE, status = "not_triggered"))
  }
  ord <- if ((posv$mct * (ctx$nind / ctx$k)) < 3) 3 else 1
  common <- list(
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    n_one = posv$n_onetwo,
    ytwo = posv$ytwo,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[gene_index, ],
    posind = ctx$posind[[gene_index]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k,
    betas = beta,
    reml = 0,
    eps = ctx$eps,
    ord = ord,
    intcol = ctx$intcol
  )
  dispersion <- tryCatch({
    if (pair_trigger) {
      fit <- do.call(
        bobyqa,
        c(list(x0 = c(sigma2, gamma), fn = pql_ll,
               lower = ctx$min, upper = ctx$max), common)
      )
      as.numeric(fit$par[seq_len(2L)])
    } else {
      fit <- do.call(
        nlminb,
        c(list(start = sigma2, objective = pql_gamma_ll,
               lower = ctx$min[1L], upper = ctx$max[1L], gamma = gamma), common)
      )
      c(as.numeric(fit$par[[1L]]), gamma)
    }
  }, error = function(e) NULL)
  if (is.null(dispersion) || length(dispersion) != 2L ||
      any(!is.finite(dispersion))) {
    return(list(
      triggered = TRUE,
      status = "generation_failed",
      method = if (pair_trigger) "pql_pair" else "pql_sigma"
    ))
  }
  list(
    triggered = TRUE,
    status = "generated",
    method = if (pair_trigger) "pql_pair" else "pql_sigma",
    start = c(map$project(beta), dispersion),
    gni = gni,
    kappa_obs = kappa_obs,
    cv2p = cv2p
  )
}

ln_joint_candidate_summary <- function(candidates) {
  do.call(rbind, lapply(candidates, function(x) {
    data.frame(
      label = x$label,
      start_provenance = x$start_provenance,
      required = isTRUE(x$required),
      core = isTRUE(x$core),
      finite_endpoint = isTRUE(x$finite_endpoint),
      certified = isTRUE(x$certified),
      optimizer_convergence = if (!is.null(x$optimizer$convergence))
        as.integer(x$optimizer$convergence) else NA_integer_,
      optimizer_message = if (!is.null(x$optimizer$message))
        as.character(x$optimizer$message) else NA_character_,
      start = if (!is.null(x$start))
        paste(format(x$start, digits = 17), collapse = ";") else NA_character_,
      optimizer_endpoint = if (!is.null(x$optimizer_endpoint))
        paste(format(x$optimizer_endpoint, digits = 17), collapse = ";") else
        NA_character_,
      optimizer_objective = if (!is.null(x$optimizer_objective))
        as.numeric(x$optimizer_objective) else NA_real_,
      polished = isTRUE(x$polished),
      theta_reduced = if (!is.null(x$theta_reduced))
        paste(format(x$theta_reduced, digits = 17), collapse = ";") else NA_character_,
      theta_full_log = if (!is.null(x$theta_full_log))
        paste(format(x$theta_full_log, digits = 17), collapse = ";") else NA_character_,
      objective = if (length(x$objective)) as.numeric(x$objective) else NA_real_,
      root_residual = if (!is.null(x$root_residual)) x$root_residual else NA_real_,
      newton_displacement = if (!is.null(x$newton_displacement))
        x$newton_displacement else NA_real_,
      donor_additivity_scaled = if (!is.null(x$donor_additivity_scaled))
        x$donor_additivity_scaled else NA_real_,
      scaled_min_eigenvalue = if (!is.null(x$curvature$scaled_min_eigenvalue))
        x$curvature$scaled_min_eigenvalue else NA_real_,
      scaled_condition_number = if (!is.null(x$curvature$scaled_condition_number))
        x$curvature$scaled_condition_number else NA_real_,
      failure = if (!is.null(x$failure)) x$failure else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}

ln_joint_candidate_records <- function(candidates) {
  lapply(candidates, function(x) {
    list(
      label = x$label,
      start_provenance = x$start_provenance,
      required = isTRUE(x$required),
      core = isTRUE(x$core),
      start = x$start,
      optimizer_endpoint = x$optimizer_endpoint,
      optimizer_objective = x$optimizer_objective,
      polished = isTRUE(x$polished),
      theta_reduced = x$theta_reduced,
      theta_full_log = x$theta_full_log,
      objective_provenance = x$objective_provenance,
      finite_endpoint = isTRUE(x$finite_endpoint),
      certified = isTRUE(x$certified),
      objective = x$objective,
      optimizer_convergence = x$optimizer$convergence,
      optimizer_message = x$optimizer$message,
      donor_additivity_scaled = x$donor_additivity_scaled,
      root_residual = x$root_residual,
      newton_displacement = x$newton_displacement,
      curvature = x$curvature,
      interior = x$interior,
      constraint_residual = x$constraint_residual,
      checks = x$checks,
      failure = x$failure
    )
  })
}

ln_joint_resolve_candidates <- function(candidates) {
  thresholds <- ln_joint_thresholds()
  certified <- which(vapply(candidates, function(x) isTRUE(x$certified), logical(1)))
  finite <- which(vapply(candidates, function(x) {
    isTRUE(x$finite_endpoint) && is.finite(x$objective)
  }, logical(1)))
  best_idx <- if (length(certified)) {
    certified[which.min(vapply(candidates[certified], `[[`, numeric(1), "objective"))]
  } else if (length(finite)) {
    finite[which.min(vapply(candidates[finite], `[[`, numeric(1), "objective"))]
  } else {
    NA_integer_
  }
  failure <- NA_character_
  parameter_agreement <- NA_real_
  objective_agreement <- NA_real_
  competing <- character(0)
  dominated <- character(0)
  required_failed <- vapply(candidates, function(x) {
    isTRUE(x$required) && !isTRUE(x$finite_endpoint)
  }, logical(1))
  if (any(required_failed)) {
    failure <- "required_start_failed"
  } else if (!length(certified)) {
    failure <- candidates[[best_idx]]$failure
    if (is.null(failure) || !length(failure) || is.na(failure[[1L]])) {
      failure <- "no_certified_candidate"
    }
  } else {
    best <- candidates[[best_idx]]
    core_certified <- sum(vapply(candidates, function(x) {
      isTRUE(x$core) && isTRUE(x$certified)
    }, logical(1)))
    if (core_certified < 2L) {
      failure <- "insufficient_certified_core_starts"
    }
    parameter_agreement <- max(vapply(candidates[certified], function(x) {
      ln_joint_scaled_max(x$theta_full_log, best$theta_full_log)
    }, numeric(1)))
    objective_agreement <- max(vapply(candidates[certified], function(x) {
      abs(x$objective - best$objective) /
        (1 + max(abs(x$objective), abs(best$objective)))
    }, numeric(1)))
    if (is.na(failure) && (!is.finite(parameter_agreement) ||
        parameter_agreement > thresholds$parameter_agreement)) {
      failure <- if (is.finite(parameter_agreement) &&
          parameter_agreement <= thresholds$hard_stop$parameter_agreement) {
        "start_parameter_review_zone"
      } else {
        "start_parameter_hard_stop"
      }
    }
    if (is.na(failure) && (!is.finite(objective_agreement) ||
        objective_agreement > thresholds$objective_agreement)) {
      failure <- if (is.finite(objective_agreement) &&
          objective_agreement <= thresholds$hard_stop$objective_agreement) {
        "start_objective_review_zone"
      } else {
        "start_objective_hard_stop"
      }
    }
    uncertified_finite <- which(vapply(candidates, function(x) {
      isTRUE(x$finite_endpoint) && !isTRUE(x$certified)
    }, logical(1)))
    if (length(uncertified_finite)) {
      classifications <- lapply(candidates[uncertified_finite], function(x) {
        relative_gap <- (x$objective - best$objective) /
          (1 + max(abs(x$objective), abs(best$objective)))
        if (!is.finite(relative_gap) ||
            relative_gap <= thresholds$hard_stop$objective_agreement) {
          list(competing = x$label, dominated = NA_character_)
        } else {
          list(competing = NA_character_, dominated = x$label)
        }
      })
      competing <- vapply(classifications, `[[`, character(1), "competing")
      competing <- na.omit(competing)
      dominated <- vapply(classifications, `[[`, character(1), "dominated")
      dominated <- na.omit(dominated)
      if (is.na(failure) && length(competing)) {
        failure <- "competing_uncertified_candidate"
      }
    }
  }
  list(
    route = if (is.na(failure)) "R0" else "R2",
    failure_reason = failure,
    best_index = best_idx,
    parameter_agreement = parameter_agreement,
    objective_agreement = objective_agreement,
    competing_uncertified = as.character(competing),
    dominated_uncertified = as.character(dominated),
    candidates = ln_joint_candidate_summary(candidates)
  )
}

ln_joint_fit_engine <- function(gene_index, posv, ctx, map, extra_starts = NULL) {
  free_dim <- map$free_dim
  lower <- c(if (free_dim) rep(-100, free_dim) else numeric(0), ctx$min)
  upper <- c(if (free_dim) rep(100, free_dim) else numeric(0), ctx$max)
  betae <- rep(0, ctx$nb)
  betae[ctx$intcol] <- log(posv$mct) - ctx$moffset
  beta_starts <- list(
    S_direct = betae,
    S_mid = betae,
    S_up = replace(betae, ctx$intcol, betae[ctx$intcol] + log(2)),
    S_down = replace(betae, ctx$intcol, betae[ctx$intcol] - log(2))
  )
  dispersion_starts <- list(
    S_direct = c(1, ctx$cell_init),
    S_mid = sqrt(ctx$min * ctx$max),
    S_up = c(1, ctx$cell_init),
    S_down = c(1, ctx$cell_init)
  )
  core_starts <- Map(function(beta, dispersion) {
    ln_joint_clamp(c(map$project(beta), dispersion), lower, upper)
  }, beta_starts, dispersion_starts)
  keys <- vapply(core_starts, function(x) paste(format(x, digits = 17), collapse = "|"), character(1))
  core_starts <- core_starts[!duplicated(keys)]

  objective <- function(theta) ln_joint_objective(theta, gene_index, posv, ctx, map)
  gradient <- function(theta) ln_joint_gradient(theta, gene_index, posv, ctx, map)
  controls <- list(eval.max = 500L, iter.max = 300L, rel.tol = 1e-12,
                   x.tol = 1e-10, xf.tol = 1e-12)

  run_start <- function(start, label, start_provenance, required, core,
                        optimizer = "nlminb") {
    fit <- tryCatch({
      if (identical(optimizer, "lbfgs")) {
        ans <- lbfgs(
          start,
          function(theta) list(objective = objective(theta), gradient = gradient(theta)),
          lower = lower,
          upper = upper,
          control = list(ftol_abs = ctx$eps)
        )
        list(
          par = ans$par,
          objective = objective(ans$par),
          convergence = ans$convergence,
          message = ans$message
        )
      } else if (identical(optimizer, "trust")) {
        log_start <- c(
          if (free_dim) start[seq_len(free_dim)] else numeric(0),
          log(start[free_dim + seq_len(2L)])
        )
        ans <- trust(
          objfun = function(z) {
            state <- ln_joint_derivative_state(z, gene_index, posv, ctx, map)
            natural <- c(
              if (free_dim) z[seq_len(free_dim)] else numeric(0),
              exp(z[free_dim + seq_len(2L)])
            )
            list(
              value = objective(natural),
              gradient = state$q_gradient_reduced,
              hessian = state$A_reduced
            )
          },
          parinit = log_start,
          rinit = 2,
          rmax = 100,
          iterlim = 100
        )
        par <- c(
          if (free_dim) ans$argument[seq_len(free_dim)] else numeric(0),
          exp(ans$argument[free_dim + seq_len(2L)])
        )
        list(
          par = par,
          objective = objective(par),
          convergence = if (isTRUE(ans$converged)) 0L else 1L,
          message = "trust"
        )
      } else {
        nlminb(
          start = start,
          objective = objective,
          gradient = gradient,
          lower = lower,
          upper = upper,
          control = controls
        )
      }
    }, error = function(e) {
      list(par = rep(NA_real_, length(start)), objective = NA_real_,
           convergence = NA_integer_, message = conditionMessage(e))
    })
    endpoint <- fit$par
    endpoint_objective <- fit$objective
    optimizer_endpoint <- endpoint
    optimizer_objective <- endpoint_objective
    polished <- FALSE
    if (length(endpoint) == length(start) && all(is.finite(endpoint)) &&
        all(endpoint[free_dim + seq_len(2L)] > 0)) {
      polish <- tryCatch(
        ln_joint_polish(endpoint, gene_index, posv, ctx, map, lower, upper),
        error = function(e) NULL
      )
      if (!is.null(polish) && is.finite(polish$objective)) {
        endpoint <- polish$theta_reduced
        endpoint_objective <- polish$objective
        polished <- !isTRUE(all.equal(
          endpoint,
          optimizer_endpoint,
          tolerance = 0,
          check.attributes = FALSE
        ))
      }
    }
    candidate <- ln_joint_candidate_state(
      endpoint,
      endpoint_objective,
      optimizer = fit,
      label = label,
      start_provenance = start_provenance,
      required = required,
      core = core,
      gene_index = gene_index,
      posv = posv,
      ctx = ctx,
      map = map
    )
    candidate$start <- as.numeric(start)
    candidate$optimizer_endpoint <- optimizer_endpoint
    candidate$optimizer_objective <- optimizer_objective
    candidate$polished <- polished
    candidate
  }

  candidates <- Map(function(start, label) {
    run_start(start, label, "direct", required = TRUE, core = TRUE)
  }, core_starts, names(core_starts))

  primary <- run_start(
    core_starts[[1L]],
    "primary_optimizer",
    "direct",
    required = FALSE,
    core = FALSE,
    optimizer = ctx$opt
  )
  candidates <- c(candidates, list(primary))

  finite_for_pql <- which(vapply(candidates, function(x) {
    isTRUE(x$finite_endpoint) && is.finite(x$objective)
  }, logical(1)))
  pql <- list(triggered = FALSE, status = "not_evaluated")
  if (length(finite_for_pql)) {
    prelim <- candidates[[finite_for_pql[
      which.min(vapply(candidates[finite_for_pql], `[[`, numeric(1), "objective"))
    ]]]
    pql <- ln_joint_pql_start(prelim, gene_index, posv, ctx, map)
    if (isTRUE(pql$triggered) && identical(pql$status, "generated")) {
      candidates <- c(candidates, list(run_start(
        ln_joint_clamp(pql$start, lower, upper),
        "pql_completed_start",
        "pql_start",
        required = TRUE,
        core = FALSE
      )))
    } else if (isTRUE(pql$triggered)) {
      candidates <- c(candidates, list(list(
        label = "pql_completed_start",
        start_provenance = "pql_start",
        required = TRUE,
        core = FALSE,
        finite_endpoint = FALSE,
        certified = FALSE,
        objective = NA_real_,
        optimizer = list(convergence = NA_integer_),
        failure = pql$status
      )))
    }
  }

  if (!is.null(extra_starts)) {
    for (ii in seq_along(extra_starts)) {
      start_full <- as.numeric(extra_starts[[ii]])
      if (length(start_full) != ctx$nb + 2L || any(!is.finite(start_full))) next
      start_reduced <- c(map$project(start_full[seq_len(ctx$nb)]),
                         start_full[ctx$nb + seq_len(2L)])
      candidates <- c(candidates, list(run_start(
        ln_joint_clamp(start_reduced, lower, upper),
        paste0("declared_extra_", ii),
        "other_declared_start",
        required = FALSE,
        core = FALSE
      )))
    }
  }

  resolution <- ln_joint_resolve_candidates(candidates)
  best <- if (!is.na(resolution$best_index)) candidates[[resolution$best_index]] else NULL
  list(
    route = resolution$route,
    failure_reason = resolution$failure_reason,
    best = best,
    start_provenance = if (!is.null(best)) best$start_provenance else NA_character_,
    candidates = candidates,
    diagnostics = resolution,
    derivative_audit = NULL,
    pql = pql,
    map = map
  )
}

ln_joint_covariance <- function(state, nb) {
  A <- state$A_full
  U_subject <- state$U_subject_full
  B <- U_subject %*% t(U_subject)
  A_inv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(A_inv) || any(!is.finite(A_inv))) {
    return(list(success = FALSE, failure_reason = "singular_joint_curvature"))
  }
  model_full <- 0.5 * (A_inv + t(A_inv))
  robust_full <- A_inv %*% B %*% t(A_inv)
  robust_full <- 0.5 * (robust_full + t(robust_full))
  influence <- A_inv %*% U_subject
  direct_influence <- influence %*% t(influence)
  direct_error <- ln_joint_relative_frobenius(robust_full, direct_influence)
  profile_available <- FALSE
  profile_error <- NA_real_
  profile_beta <- NULL
  B_inv <- tryCatch(solve(B), error = function(e) NULL)
  if (!is.null(B_inv) && all(is.finite(B_inv))) {
    G <- t(A) %*% B_inv %*% A
    idx_beta <- seq_len(nb)
    idx_rho <- nb + seq_len(nrow(A) - nb)
    G_rr_inv <- tryCatch(solve(G[idx_rho, idx_rho, drop = FALSE]),
                         error = function(e) NULL)
    if (!is.null(G_rr_inv)) {
      G_beta <- G[idx_beta, idx_beta, drop = FALSE] -
        G[idx_beta, idx_rho, drop = FALSE] %*% G_rr_inv %*%
        G[idx_rho, idx_beta, drop = FALSE]
      profile_beta <- tryCatch(solve(G_beta), error = function(e) NULL)
      if (!is.null(profile_beta) && all(is.finite(profile_beta))) {
        profile_available <- TRUE
        profile_error <- ln_joint_relative_frobenius(
          robust_full[idx_beta, idx_beta, drop = FALSE],
          profile_beta
        )
      }
    }
  }
  thresholds <- ln_joint_thresholds()
  beta_model <- model_full[seq_len(nb), seq_len(nb), drop = FALSE]
  beta_robust <- robust_full[seq_len(nb), seq_len(nb), drop = FALSE]
  success <- is.finite(direct_error) &&
    direct_error <= thresholds$covariance_identity &&
    (!profile_available || (is.finite(profile_error) &&
      profile_error <= thresholds$covariance_identity)) &&
    all(is.finite(c(beta_model, beta_robust))) &&
    all(diag(beta_model) > 0) && all(diag(beta_robust) > 0)
  list(
    success = success,
    failure_reason = if (success) NA_character_ else "covariance_identity_failed",
    A = A,
    B = B,
    model_full = model_full,
    robust_full = robust_full,
    model_beta = beta_model,
    robust_beta = beta_robust,
    influence = influence,
    direct_influence_error = direct_error,
    profile_available = profile_available,
    profile_beta = profile_beta,
    profile_error = profile_error
  )
}

ln_joint_fixed_modes <- function(state, gene_index, posv, ctx) {
  unavailable <- function(status, method = "fixed_beta_pml") {
    list(
      values = rep(NA_real_, ctx$k),
      method = method,
      status = status,
      success = FALSE,
      max_score_residual = NA_real_,
      max_newton_displacement = NA_real_,
      donor_status = rep(status, ctx$k),
      alpha = NA_real_,
      lambda = NA_real_
    )
  }
  if (!isTRUE(ctx$output_re)) {
    return(unavailable("not_requested", method = "none"))
  }
  if (is.null(state)) {
    return(unavailable("fit_unavailable"))
  }
  result <- tryCatch(
    opt_pml_fixed_beta_modes(
      ctx$pred,
      ctx$offset,
      posv$Y,
      ctx$fid - 1L,
      ctx$cumsumy[gene_index, ],
      posv$posindy,
      ctx$nind,
      ctx$k,
      state$beta,
      c(state$sigma2, state$gamma),
      ln_joint_thresholds()$mode_score,
      200L
    ),
    error = function(e) NULL
  )
  if (is.null(result)) {
    return(unavailable("solver_error"))
  }
  max_score <- as.numeric(result$max_score_residual)
  max_newton <- as.numeric(result$max_newton_displacement)
  donor_status <- as.character(result$status)
  success <- length(result$logw) == ctx$k &&
    all(is.finite(result$logw)) &&
    all(donor_status == "ok") &&
    is.finite(max_score) && max_score <= ln_joint_thresholds()$mode_score &&
    is.finite(max_newton) && max_newton <= ln_joint_thresholds()$mode_newton
  list(
    values = if (success) as.numeric(result$logw) else rep(NA_real_, ctx$k),
    method = "fixed_beta_pml",
    status = if (success) "ok" else "gate_failed",
    success = success,
    max_score_residual = max_score,
    max_newton_displacement = max_newton,
    donor_status = donor_status,
    alpha = as.numeric(result$alpha),
    lambda = as.numeric(result$lambda),
    parameterization = result$parameterization
  )
}

fit_gene_joint_ln_unconstrained <- function(gene_index, posv, ctx,
                                             extra_starts = NULL) {
  objective_instance <- ln_joint_objective_instance_registry(
    "register",
    gene_index,
    posv,
    ctx
  )
  map <- ln_joint_map(ctx$nb)
  engine <- ln_joint_fit_engine(gene_index, posv, ctx, map, extra_starts)
  best <- engine$best
  state <- if (!is.null(best)) best$state else NULL
  if (!is.null(state)) {
    state$objective_instance <- objective_instance
  }
  route <- engine$route
  failure <- engine$failure_reason
  covariance <- NULL
  if (identical(route, "R0")) {
    covariance <- ln_joint_covariance(state, ctx$nb)
    if (!isTRUE(covariance$success)) {
      route <- "R2"
      failure <- covariance$failure_reason
    }
  }
  inference_available <- identical(route, "R0")
  mode <- ln_joint_fixed_modes(if (inference_available) state else NULL,
                               gene_index, posv, ctx)
  beta <- if (!is.null(state)) state$beta else rep(NA_real_, ctx$nb)
  sigma2 <- if (!is.null(state)) state$sigma2 else NA_real_
  gamma <- if (!is.null(state)) state$gamma else NA_real_
  model_cov <- if (inference_available) covariance$model_beta else
    matrix(NA_real_, ctx$nb, ctx$nb)
  robust_cov <- if (inference_available) covariance$robust_beta else
    matrix(NA_real_, ctx$nb, ctx$nb)
  loglik <- if (!is.null(best) && is.finite(best$objective)) -best$objective else NA_real_
  conv_code <- if (inference_available) 1L else -70L
  restemp <- c(
    beta,
    sigma2,
    if (is.finite(gamma) && gamma > 0) 1 / gamma else NA_real_,
    diag(model_cov),
    diag(robust_cov),
    conv_code,
    1L,
    loglik
  )
  if (ctx$output_re) {
    restemp <- c(restemp, mode$values)
  }
  derivative_object <- if (!is.null(state)) {
    list(
      gradient = state$U_full,
      hessian = state$A_full,
      observed_info = state$A_full,
      per_subject_gradients = state$U_subject_full
    )
  } else NULL
  godambe <- if (!is.null(derivative_object))
    extract_godambe_components(derivative_object, ctx$nb) else NULL
  theta <- list(beta = beta, subVar = sigma2, cellVar = gamma)
  result <- list(
    vec = restemp,
    godambe = godambe,
    loglik = loglik,
    loglik_ln = loglik,
    theta = theta,
    theta_ln = theta,
    theta_ln_rescued = theta,
    theta_unrescued = theta,
    score = if (!is.null(state)) state$U_full else rep(NA_real_, ctx$nb + 2L),
    naive_cov = model_cov,
    robust_cov = robust_cov,
    route = route,
    start_provenance = engine$start_provenance,
    inference_available = inference_available,
    failure_reason = failure,
    mode_method = mode$method,
    mode_status = mode$status,
    mode = mode,
    ln_joint = list(
      parameterization = "beta, log(sigma2), log(gamma_nb_size)",
      score_sign = "maximizing",
      objective_provenance = if (!is.null(state)) {
        state$objective_provenance
      } else {
        ln_joint_objective_provenance()
      },
      objective_instance = objective_instance,
      state = state,
      covariance = covariance,
      candidates = engine$diagnostics$candidates,
      candidate_records = ln_joint_candidate_records(engine$candidates),
      parameter_agreement = engine$diagnostics$parameter_agreement,
      objective_agreement = engine$diagnostics$objective_agreement,
      derivative_audit = engine$derivative_audit,
      pql_start = engine$pql
    ),
    ln_refinement = list(
      refined = TRUE,
      grad_norm_start = NA_real_,
      grad_norm_final = if (!is.null(best)) best$root_residual else NA_real_,
      extra_attempts = if (is.null(extra_starts)) 0L else length(extra_starts),
      extra_success = 0L
    ),
    debug = list(
      method = "joint_ln",
      gene_index = gene_index,
      fit_code = 1L,
      conv_LN = conv_code,
      conv_code = conv_code,
      route = route,
      failure_reason = failure
    )
  )
  ln_joint_objective_instance_registry(
    "store_unrestricted",
    gene_index,
    posv,
    ctx,
    fit = result
  )
  result
}

fit_gene_joint_ln_constrained <- function(gene_index, posv, ctx, L, b) {
  objective_instance <- ln_joint_objective_instance_registry(
    "lookup",
    gene_index,
    posv,
    ctx
  )
  map <- ln_joint_map(ctx$nb, L, b)
  engine <- ln_joint_fit_engine(gene_index, posv, ctx, map)
  best <- engine$best
  state <- if (!is.null(best)) best$state else NULL
  if (!is.null(state)) {
    state$objective_instance <- objective_instance
  }
  route <- engine$route
  failure <- engine$failure_reason
  beta <- if (!is.null(state)) state$beta else rep(NA_real_, ctx$nb)
  sigma2 <- if (!is.null(state)) state$sigma2 else NA_real_
  gamma <- if (!is.null(state)) state$gamma else NA_real_
  loglik <- if (!is.null(best) && is.finite(best$objective)) -best$objective else NA_real_
  inference_available <- identical(route, "R0")
  derivative_object <- if (!is.null(state)) {
    list(
      gradient = state$U_full,
      hessian = state$A_full,
      observed_info = state$A_full,
      per_subject_gradients = state$U_subject_full
    )
  } else NULL
  godambe <- if (!is.null(derivative_object))
    extract_godambe_components(derivative_object, ctx$nb) else NULL
  theta <- list(beta = beta, subVar = sigma2, cellVar = gamma)
  conv_code <- if (inference_available) 1L else -70L
  restemp <- c(
    beta,
    sigma2,
    if (is.finite(gamma) && gamma > 0) 1 / gamma else NA_real_,
    rep(NA_real_, 2L * ctx$nb),
    conv_code,
    1L,
    loglik
  )
  if (ctx$output_re) restemp <- c(restemp, rep(NA_real_, ctx$k))
  reduced <- if (!is.null(state)) list(
    parameterization = "zeta, log(sigma^2), log(gamma_nb_size)",
    map = map$derivative_map,
    theta = best$theta_reduced_log,
    theta_red = best$theta_reduced_log,
    gamma = if (map$free_dim) best$theta_reduced_log[seq_len(map$free_dim)] else numeric(0),
    score = state$U_reduced,
    g_red = state$U_reduced,
    observed_info = state$A_reduced,
    H_red = state$A_reduced,
    per_subject_gradients = state$U_subject_reduced
  ) else NULL
  result <- list(
    vec = restemp,
    godambe = godambe,
    beta_constrained = beta,
    loglik = loglik,
    loglik_ln = loglik,
    loglik_ln_pre = loglik,
    loglik_ln_rescued = loglik,
    theta = theta,
    theta_ln = theta,
    theta_ln_rescued = theta,
    theta_unrescued = theta,
    constraint = list(
      L = L,
      b = b,
      beta_star = map$beta_star,
      null_basis = map$K,
      free_dim = map$free_dim
    ),
    reduced = reduced,
    diagnostics = if (!is.null(best)) list(
      reduced_gradient_max_abs = max(abs(state$U_reduced)),
      H_red_positive_definite = isTRUE(best$curvature$positive_definite),
      H_red_rcond = if (is.finite(best$curvature$scaled_condition_number))
        1 / best$curvature$scaled_condition_number else 0,
      max_abs_newton_step_scaled = best$newton_displacement,
      constraint_residual_scaled = best$constraint_residual
    ) else list(),
    score = if (!is.null(state)) state$U_full else rep(NA_real_, ctx$nb + 2L),
    robust_cov = matrix(NA_real_, ctx$nb, ctx$nb),
    route = route,
    start_provenance = engine$start_provenance,
    inference_available = inference_available,
    failure_reason = failure,
    lrt_gate = list(
      converged = inference_available,
      interior = if (!is.null(best)) isTRUE(best$interior) else FALSE,
      nuisance_interior = if (!is.null(best)) isTRUE(best$interior) else FALSE,
      constraint_satisfied = if (!is.null(best))
        best$constraint_residual <= ln_joint_thresholds()$constraint else FALSE,
      calibration_available = FALSE,
      calibration_status = "joint_ln_calibration_pending"
    ),
    ln_joint = list(
      parameterization = "beta, log(sigma2), log(gamma_nb_size)",
      score_sign = "maximizing",
      objective_provenance = if (!is.null(state)) {
        state$objective_provenance
      } else {
        ln_joint_objective_provenance()
      },
      objective_instance = objective_instance,
      state = state,
      candidates = engine$diagnostics$candidates,
      candidate_records = ln_joint_candidate_records(engine$candidates),
      parameter_agreement = engine$diagnostics$parameter_agreement,
      objective_agreement = engine$diagnostics$objective_agreement,
      derivative_audit = engine$derivative_audit
    ),
    debug = list(
      method = "joint_ln_constrained",
      gene_index = gene_index,
      fit_code = 1L,
      conv_LN = conv_code,
      conv_code = conv_code,
      free_dim = map$free_dim,
      route = route,
      failure_reason = failure
    )
  )
  result$ln_joint$score_diagnostic <-
    ln_joint_robust_score_diagnostic(result)
  unrestricted_fit <- ln_joint_objective_instance_registry(
    "get_unrestricted",
    gene_index,
    posv,
    ctx
  )
  result$ln_joint$pace_ratio_diagnostic <-
    ln_joint_pace_ratio_diagnostic(
      if (is.list(unrestricted_fit)) unrestricted_fit else list(),
      result
    )
  result
}
