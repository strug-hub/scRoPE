# Standalone fixtures and numerical validation helpers for the corrected
# first-order nested HL pipeline.  These helpers deliberately do not call the
# package fitting path: formula, solver, and integration tests build on them in
# later validation stages.

hl_fo_fixture_ids <- function() {
  c("R1", "C1", "C2", "K0", "Q0", "Q1", "Z1", "D1", "W1", "X0", "N0", "B1", "S1")
}

hl_fo_fixture_definitions <- function() {
  beta_r1 <- c(`(Intercept)` = -1.2, z = 0.35, w = -0.25)
  beta_c2 <- c(`(Intercept)` = -1, z = 0.25, w = -0.2)
  beta_d1 <- c(`(Intercept)` = -1.2, z = 0.2, w = -0.15)

  list(
    R1 = list(
      seed = 104729L, m = 240L, n_i = 5L,
      beta = beta_r1, sigma2 = 0.4, phi = 8,
      active_dispersion = c("tau", "xi"),
      expected = "regular_unrestricted_success"
    ),
    C1 = list(
      source = "R1", seed = 104729L,
      constraint = list(
        C = matrix(c(0, 1, 0), nrow = 1L,
          dimnames = list("z_at_truth", names(beta_r1))),
        rhs = c(z_at_truth = 0.35)
      ),
      active_dispersion = c("tau", "xi"),
      expected = "constrained_and_score_success"
    ),
    C2 = list(
      seed = 104759L, m = 300L, n_i = 5L,
      beta = beta_c2, sigma2 = 0.6, phi = 5,
      constraint = list(
        C = rbind(
          z_at_truth = c(0, 1, 0),
          w_at_truth = c(0, 0, 1)
        ),
        rhs = c(z_at_truth = 0.25, w_at_truth = -0.2)
      ),
      active_dispersion = c("tau", "xi"),
      expected = "multirow_and_alternative_basis_success"
    ),
    K0 = list(
      source = "C2", seed = 104761L,
      constraint = list(
        C = structure(diag(3L), dimnames = list(names(beta_c2), names(beta_c2))),
        rhs = beta_c2
      ),
      active_dispersion = c("tau", "xi"),
      expected = "k_zero_point_and_covariance_success"
    ),
    Q0 = list(
      source = "R1", seed = 104729L,
      active_dispersion = character(0),
      fixed_dispersion = c(sigma2 = 0.4, phi = 8),
      expected = "q_zero_conditional_beta_and_sandwich_success"
    ),
    Q1 = list(
      source = "R1", seed = 104729L,
      active_dispersion = "xi",
      fixed_dispersion = c(sigma2 = 0.4),
      expected = "one_active_dispersion_success"
    ),
    Z1 = list(
      seed = 104773L, m = 300L, n_i = 5L,
      beta = c(`(Intercept)` = -4, z = 0.2, w = -0.2),
      sigma2 = 0.8, phi = 2,
      active_dispersion = c("tau", "xi"),
      expected = "zero_heavy_interior_success"
    ),
    D1 = list(
      seeds = 104779:104782, m = 300L, n_i = 5L,
      beta = beta_d1,
      cases = list(
        low_v_high_phi = c(sigma2 = 0.02, phi = 100),
        low_v_low_phi = c(sigma2 = 0.02, phi = 0.3),
        high_v_high_phi = c(sigma2 = 2, phi = 100),
        high_v_low_phi = c(sigma2 = 2, phi = 0.3)
      ),
      active_dispersion = c("tau", "xi"),
      expected = "four_regular_interior_dispersion_successes"
    ),
    W1 = list(
      seed = 104789L, m = 80L, n_i = 5L,
      beta = beta_r1, sigma2 = 0.4, phi = 8,
      w_rule = "weak",
      display_coordinate = "raw_package_input",
      design_coordinate = "package_centered_scaled_analysis",
      active_dispersion = c("tau", "xi"),
      expected = "weak_curvature_after_rank_is_retained"
    ),
    X0 = list(
      seed = 104801L, m = 80L, n_i = 5L,
      beta = beta_r1, sigma2 = 0.4, phi = 8,
      w_rule = "collinear",
      active_dispersion = c("tau", "xi"),
      expected = "inner_rank_deficient"
    ),
    N0 = list(
      seed = 104803L, m = 30L, n_i = 3L,
      beta = c(`(Intercept)` = NA_real_), sigma2 = 0.4, phi = 8,
      design = "intercept_only", counts = "all_zero", exposure = "one",
      active_dispersion = c("tau", "xi"),
      expected = "inner_no_finite_root"
    ),
    B1 = list(
      seed = 104827L, m = 300L, n_i = 5L,
      beta = beta_d1, sigma2 = 1e-8, phi = 1e8,
      active_dispersion = c("tau", "xi"),
      expected = c("certified_interior_success", "outer_boundary")
    ),
    S1 = list(
      seed = 104831L, m = 6L, n_i = 5L,
      beta = beta_r1, sigma2 = 0.4, phi = 8,
      active_dispersion = c("tau", "xi"),
      expected = "point_fit_with_explicit_small_subject_rank_limits"
    )
  )
}

hl_fo_with_rng <- function(seed, code) {
  stopifnot(is.function(code), length(seed) == 1L, is.finite(seed))
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL

  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(as.integer(seed))
  code()
}

hl_fo_num_rank <- function(x) {
  x <- as.matrix(x)
  if (!all(is.finite(x))) {
    stop("Numerical rank requires a finite matrix.")
  }
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr == 0L || nc == 0L || !length(x) || all(x == 0)) {
    return(list(rank = 0L, tolerance = 0, singular_values = numeric(0)))
  }
  singular_values <- svd(x, nu = 0L, nv = 0L)$d
  tolerance <- max(nr, nc) * .Machine$double.eps * singular_values[[1L]]
  list(
    rank = as.integer(sum(singular_values > tolerance)),
    tolerance = tolerance,
    singular_values = singular_values
  )
}

hl_fo_sdisc <- function(a, b) {
  if (length(a) != length(b) || !identical(dim(a), dim(b))) {
    stop("Scaled discrepancy requires equally shaped objects.")
  }
  if (!length(a)) {
    return(0)
  }
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (!all(is.finite(a)) || !all(is.finite(b))) {
    return(Inf)
  }
  max(abs(a - b) / (1 + pmax(abs(a), abs(b))))
}

hl_fo_step_sequence <- function(z_r, eps = .Machine$double.eps) {
  stopifnot(length(z_r) == 1L, is.finite(z_r), length(eps) == 1L, eps > 0)
  h0 <- eps^(1 / 3) * (1 + abs(z_r))
  h0 / (2^(0:5))
}

hl_fo_package_center_scale <- function(X_raw) {
  X_raw <- as.matrix(X_raw)
  constant <- vapply(seq_len(ncol(X_raw)), function(jj) {
    all(X_raw[, jj] == X_raw[[1L, jj]])
  }, logical(1))
  if (sum(constant) != 1L || X_raw[1L, which(constant)] == 0) {
    stop("Fixture preprocessing requires one nonzero constant column.")
  }
  intcol <- which(constant)
  means <- colMeans(X_raw)
  scales <- sqrt(colMeans(sweep(X_raw, 2L, means, FUN = "-")^2))
  X <- X_raw
  X[, intcol] <- 1
  varying <- setdiff(seq_len(ncol(X_raw)), intcol)
  X[, varying] <- sweep(
    sweep(X_raw[, varying, drop = FALSE], 2L, means[varying], FUN = "-"),
    2L,
    scales[varying],
    FUN = "/"
  )
  list(X = X, means = means, scales = scales, intcol = intcol)
}

hl_fo_fd_trace <- function(lower_order,
                           point,
                           coordinate,
                           analytic = NULL,
                           ad_reference = NULL,
                           eps = .Machine$double.eps) {
  lower_order <- match.fun(lower_order)
  point <- as.numeric(point)
  coordinate <- as.integer(coordinate)
  if (length(coordinate) != 1L || is.na(coordinate) ||
      coordinate < 1L || coordinate > length(point)) {
    stop("coordinate must select one element of point.")
  }

  steps <- hl_fo_step_sequence(point[[coordinate]], eps = eps)
  plus <- minus <- central <- vector("list", length(steps))
  for (ii in seq_along(steps)) {
    xp <- xm <- point
    xp[[coordinate]] <- xp[[coordinate]] + steps[[ii]]
    xm[[coordinate]] <- xm[[coordinate]] - steps[[ii]]
    plus[[ii]] <- lower_order(xp)
    minus[[ii]] <- lower_order(xm)
    if (length(plus[[ii]]) != length(minus[[ii]]) ||
        !identical(dim(plus[[ii]]), dim(minus[[ii]]))) {
      stop("Perturbed lower-order evaluations have inconsistent shapes.")
    }
    central[[ii]] <- (plus[[ii]] - minus[[ii]]) / (2 * steps[[ii]])
  }

  richardson <- lapply(seq_len(5L), function(ii) {
    (4 * central[[ii + 1L]] - central[[ii]]) / 3
  })

  if (!is.null(ad_reference)) {
    discrepancies <- vapply(
      richardson,
      hl_fo_sdisc,
      numeric(1),
      b = ad_reference
    )
    candidate_ells <- 1:3
    selected_ell <- candidate_ells[[which.min(discrepancies[candidate_ells + 1L])]]
    selected_error <- discrepancies[[selected_ell + 1L]]
    left_error <- discrepancies[[selected_ell]]
    right_error <- discrepancies[[selected_ell + 2L]]
    reference <- ad_reference
    selector <- "ad"
  } else {
    discrepancies <- vapply(seq_len(4L), function(ii) {
      hl_fo_sdisc(richardson[[ii]], richardson[[ii + 1L]])
    }, numeric(1))
    candidate_ells <- 1:2
    selected_ell <- candidate_ells[[which.min(discrepancies[candidate_ells + 1L])]]
    selected_error <- discrepancies[[selected_ell + 1L]]
    left_error <- discrepancies[[selected_ell]]
    right_error <- discrepancies[[selected_ell + 2L]]
    reference <- richardson[[selected_ell + 2L]]
    selector <- "finite_difference"
  }

  stable_minimum <- selected_error <= left_error &&
    selected_error <= right_error &&
    left_error <= 4 * max(selected_error, 1e-12) &&
    right_error <= 4 * max(selected_error, 1e-12)

  list(
    coordinate = coordinate,
    point = point,
    initial_step = steps[[1L]],
    steps = steps,
    plus = plus,
    minus = minus,
    central = central,
    richardson = richardson,
    selector = selector,
    candidate_ells = candidate_ells,
    selected_ell = selected_ell,
    selection_discrepancies = discrepancies,
    selected_discrepancy = selected_error,
    stable_minimum = stable_minimum,
    reference = reference,
    analytic_discrepancy = if (is.null(analytic)) NA_real_ else hl_fo_sdisc(analytic, reference)
  )
}

hl_fo_generate_nb_case <- function(id, definition, seed = definition$seed) {
  m <- as.integer(definition$m)
  n_i <- as.integer(definition$n_i)
  n <- m * n_i
  beta <- definition$beta
  sigma2 <- definition$sigma2
  phi <- definition$phi
  w_rule <- if (is.null(definition$w_rule)) "normal" else definition$w_rule

  arrays <- hl_fo_with_rng(seed, function() {
    z_subject <- rbinom(m, size = 1L, prob = 0.5)
    z <- rep(z_subject, each = n_i)
    e <- NULL
    if (identical(w_rule, "normal")) {
      w <- rnorm(n)
    } else if (identical(w_rule, "weak")) {
      e <- rnorm(n)
      w <- 1e-9 * z + 1e-14 * e
    } else if (identical(w_rule, "collinear")) {
      w <- 2 * z
    } else {
      stop("Unknown w-generation rule: ", w_rule)
    }

    X_raw <- cbind(`(Intercept)` = 1, z = z, w = w)
    if (identical(w_rule, "weak")) {
      preprocessing <- hl_fo_package_center_scale(X_raw)
      X <- preprocessing$X
      beta_analysis <- beta
      varying <- setdiff(seq_len(ncol(X)), preprocessing$intcol)
      beta_analysis[preprocessing$intcol] <-
        X_raw[1L, preprocessing$intcol] * beta[preprocessing$intcol] +
        sum(preprocessing$means[varying] * beta[varying])
      beta_analysis[varying] <- preprocessing$scales[varying] * beta[varying]
    } else {
      preprocessing <- list(
        X = X_raw,
        means = setNames(rep(0, ncol(X_raw)), colnames(X_raw)),
        scales = setNames(rep(1, ncol(X_raw)), colnames(X_raw)),
        intcol = 1L
      )
      X <- X_raw
      beta_analysis <- beta
    }
    log_pi <- runif(n, min = -0.3, max = 0.3)
    pi <- exp(log_pi)
    alpha <- 1 / expm1(sigma2)
    lambda <- alpha * exp(-sigma2 / 2)
    omega <- rgamma(m, shape = alpha, rate = lambda)
    # Counts follow the displayed/raw fixture definition.  For W1 the
    # transformed beta below gives the same predictor on package analysis X.
    mu <- pi * exp(drop(X_raw %*% beta)) * rep(omega, each = n_i)
    y <- rnbinom(n, size = phi, mu = mu)

    list(
      z_subject = z_subject,
      z = z,
      w = w,
      e = e,
      X_raw = X_raw,
      X = X,
      design_means = preprocessing$means,
      design_scales = preprocessing$scales,
      intcol = preprocessing$intcol,
      beta_raw = beta,
      beta_analysis = beta_analysis,
      log_pi = log_pi,
      pi = pi,
      alpha = alpha,
      lambda = lambda,
      omega = omega,
      mu = mu,
      y = y
    )
  })

  subject_index <- rep(seq_len(m), each = n_i)
  subject_id <- sprintf("subject_%04d", subject_index)
  count <- matrix(
    arrays$y,
    nrow = 1L,
    dimnames = list(id, sprintf("cell_%05d", seq_len(n)))
  )

  list(
    id = id,
    seed = as.integer(seed),
    rng = c(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection"),
    m = m,
    n_i = n_i,
    n = n,
    p = ncol(arrays$X),
    subject_index = subject_index,
    subject_id = subject_id,
    count = count,
    y = arrays$y,
    X_raw = arrays$X_raw,
    X = arrays$X,
    z_subject = arrays$z_subject,
    z = arrays$z,
    w = arrays$w,
    e = arrays$e,
    log_pi = arrays$log_pi,
    pi = arrays$pi,
    offset = arrays$pi,
    omega = arrays$omega,
    mu = arrays$mu,
    truth = list(
      beta = arrays$beta_analysis,
      beta_raw = arrays$beta_raw,
      sigma2 = sigma2,
      phi = phi,
      tau = log(sigma2),
      xi = log(phi),
      alpha = arrays$alpha,
      lambda = arrays$lambda
    ),
    design_means = arrays$design_means,
    design_scales = arrays$design_scales,
    intcol = arrays$intcol,
    display_coordinate = if (is.null(definition$display_coordinate)) "analysis" else definition$display_coordinate,
    design_coordinate = if (is.null(definition$design_coordinate)) "analysis" else definition$design_coordinate,
    design_rank_raw = hl_fo_num_rank(arrays$X_raw),
    design_rank_analysis = hl_fo_num_rank(arrays$X),
    active_dispersion = definition$active_dispersion,
    fixed_dispersion = definition$fixed_dispersion,
    constraint = definition$constraint,
    expected = definition$expected
  )
}

hl_fo_generate_n0 <- function(definition) {
  m <- as.integer(definition$m)
  n_i <- as.integer(definition$n_i)
  n <- m * n_i
  X <- matrix(1, nrow = n, ncol = 1L,
    dimnames = list(NULL, "(Intercept)"))
  y <- rep(0, n)
  subject_index <- rep(seq_len(m), each = n_i)

  list(
    id = "N0",
    seed = definition$seed,
    rng = c(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection"),
    m = m,
    n_i = n_i,
    n = n,
    p = 1L,
    subject_index = subject_index,
    subject_id = sprintf("subject_%04d", subject_index),
    count = matrix(y, nrow = 1L,
      dimnames = list("N0", sprintf("cell_%05d", seq_len(n)))),
    y = y,
    X_raw = X,
    X = X,
    z_subject = NULL,
    z = NULL,
    w = NULL,
    e = NULL,
    log_pi = rep(0, n),
    pi = rep(1, n),
    offset = rep(1, n),
    omega = NULL,
    mu = NULL,
    truth = list(
      beta = definition$beta,
      sigma2 = definition$sigma2,
      phi = definition$phi,
      tau = log(definition$sigma2),
      xi = log(definition$phi),
      alpha = 1 / expm1(definition$sigma2),
      lambda = exp(-definition$sigma2 / 2) / expm1(definition$sigma2)
    ),
    design_means = c(`(Intercept)` = 0),
    design_scales = c(`(Intercept)` = 1),
    intcol = 1L,
    display_coordinate = "analysis",
    design_coordinate = "analysis",
    design_rank_raw = hl_fo_num_rank(X),
    design_rank_analysis = hl_fo_num_rank(X),
    active_dispersion = definition$active_dispersion,
    fixed_dispersion = NULL,
    constraint = NULL,
    expected = definition$expected
  )
}

.hl_fo_fixture_cache <- new.env(parent = emptyenv())

hl_fo_fixture <- function(id, fresh = FALSE) {
  id <- as.character(id)
  if (length(id) != 1L || !id %in% hl_fo_fixture_ids()) {
    stop("Unknown first-order HL fixture id: ", paste(id, collapse = ", "))
  }
  if (!fresh && exists(id, envir = .hl_fo_fixture_cache, inherits = FALSE)) {
    return(get(id, envir = .hl_fo_fixture_cache, inherits = FALSE))
  }

  definitions <- hl_fo_fixture_definitions()
  definition <- definitions[[id]]

  if (identical(id, "D1")) {
    case_names <- names(definition$cases)
    cases <- lapply(seq_along(case_names), function(ii) {
      case_definition <- definition
      case_definition$seed <- definition$seeds[[ii]]
      case_definition$sigma2 <- unname(definition$cases[[ii]][["sigma2"]])
      case_definition$phi <- unname(definition$cases[[ii]][["phi"]])
      case_definition$cases <- NULL
      case_definition$seeds <- NULL
      generated <- hl_fo_generate_nb_case(
        paste0("D1_", case_names[[ii]]),
        case_definition
      )
      generated$case <- case_names[[ii]]
      generated
    })
    names(cases) <- case_names
    fixture <- list(
      id = "D1",
      seeds = as.integer(definition$seeds),
      cases = cases,
      expected = definition$expected
    )
  } else if (identical(id, "N0")) {
    fixture <- hl_fo_generate_n0(definition)
  } else if (!is.null(definition$source)) {
    fixture <- hl_fo_fixture(definition$source, fresh = fresh)
    fixture$data_seed <- fixture$seed
    fixture$seed <- definition$seed
    fixture$id <- id
    fixture$source <- definition$source
    fixture$active_dispersion <- definition$active_dispersion
    fixture$fixed_dispersion <- definition$fixed_dispersion
    fixture$constraint <- definition$constraint
    fixture$expected <- definition$expected
  } else if (identical(id, "Z1")) {
    actual_seed <- definition$seed
    repeat {
      fixture <- hl_fo_generate_nb_case(id, definition, seed = actual_seed)
      if (any(fixture$y > 0)) {
        break
      }
      actual_seed <- actual_seed + 1L
    }
    fixture$requested_seed <- definition$seed
    fixture$actual_seed <- actual_seed
  } else {
    fixture <- hl_fo_generate_nb_case(id, definition)
  }

  if (!fresh) {
    assign(id, fixture, envir = .hl_fo_fixture_cache)
  }
  fixture
}
