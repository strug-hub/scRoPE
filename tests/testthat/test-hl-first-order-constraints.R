.hl_fo_constraint_get <- function(name) {
  if (exists(name, envir = .GlobalEnv, mode = "function", inherits = FALSE)) {
    return(get(name, envir = .GlobalEnv, inherits = FALSE))
  }
  getFromNamespace(name, "scRoPE")
}

.hl_fo_constraint_case <- function(id) {
  fixture <- hl_fo_fixture(id)
  data <- .hl_fo_constraint_get("hl_fo_make_data")(
    y = fixture$y,
    X = fixture$X,
    log_offset = fixture$log_pi,
    subject = fixture$subject_index,
    intcol = fixture$intcol
  )
  seed <- list(
    beta = fixture$truth$beta,
    rho_full = c(fixture$truth$tau, fixture$truth$xi),
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  list(fixture = fixture, data = data, seed = seed)
}

.hl_fo_constraint_profile <- function(case, affine, K = affine$K,
                                      active = 1:2,
                                      rho_active = case$seed$rho_full[active]) {
  .hl_fo_constraint_get("hl_fo_profile_evaluate")(
    data = case$data,
    rho_active = rho_active,
    rho_template = case$seed$rho_full,
    active = active,
    beta_star = affine$beta_star,
    K = K,
    beta_start = case$seed$beta,
    eta_start = NULL,
    cold = TRUE
  )
}

.hl_fo_constraint_mapped_profile_covariance <- function(state) {
  k_free <- ncol(state$K)
  q <- length(state$active)
  sandwich <- .hl_fo_constraint_get("hl_fo_sandwich")(
    state$augmented$A,
    state$augmented$J
  )
  if (!isTRUE(sandwich$success)) {
    stop("Profile sandwich failed: ", sandwich$status, "/", sandwich$reason)
  }
  scientific_index <- c(
    if (k_free) seq_len(k_free) else integer(0),
    if (q) 2L * k_free + seq_len(q) else integer(0)
  )
  reduced <- sandwich$V[
    scientific_index,
    scientific_index,
    drop = FALSE
  ]
  map <- matrix(0, state$data$p + q, k_free + q)
  if (k_free) {
    map[seq_len(state$data$p), seq_len(k_free)] <- state$K
  }
  if (q) {
    map[state$data$p + seq_len(q), k_free + seq_len(q)] <- diag(q)
  }
  map %*% reduced %*% t(map)
}

.hl_fo_constraint_trace_info <- function(label, trace) {
  paste(
    sprintf(
      "%s: selector=%s selected_ell=%d discrepancy=%.17g stable=%s",
      label,
      trace$selector,
      trace$selected_ell,
      trace$analytic_discrepancy,
      trace$stable_minimum
    ),
    paste(capture.output(dput(trace)), collapse = "\n"),
    sep = "\nfull step-halving trace:\n"
  )
}

test_that("arbitrary intcol and raw constant maps preserve predictors, covariance, and contrasts", {
  design_map_fn <- .hl_fo_constraint_get("hl_fo_design_map")
  constraint_map_fn <- .hl_fo_constraint_get(
    "hl_fo_transform_constraint_to_internal"
  )
  covariance_map_fn <- .hl_fo_constraint_get("hl_fo_cov_congruence")
  public_map_fn <- .hl_fo_constraint_get("hl_fo_canonical_to_public")
  mean_one_fn <- .hl_fo_constraint_get("hl_fo_canonical_to_mean_one")
  canonical_fn <- .hl_fo_constraint_get("hl_fo_mean_one_to_canonical")
  wald_fn <- .hl_fo_constraint_get("hl_fo_robust_wald")

  Z <- cbind(
    z = c(-2, -1, 0, 1, 2, 3),
    raw_intercept = rep(2.5, 6),
    w = c(0.25, -0.5, 1, -1.25, 2, -2.5)
  )
  intcol <- 2L
  varying <- setdiff(seq_len(ncol(Z)), intcol)
  means <- colMeans(Z)
  scales <- sqrt(colMeans(sweep(Z, 2L, means, FUN = "-")^2))
  scales[intcol] <- 1
  X <- Z
  X[, intcol] <- 1
  X[, varying] <- sweep(
    sweep(Z[, varying, drop = FALSE], 2L, means[varying], FUN = "-"),
    2L,
    scales[varying],
    FUN = "/"
  )

  design_map <- design_map_fn(
    column_means = means,
    column_scales = scales,
    intercept_idx = intcol,
    raw_constant = Z[[1L, intcol]]
  )
  beta_internal <- c(z = 0.7, raw_intercept = -1.1, w = 0.25)
  beta_original <- as.vector(design_map$to_original %*% beta_internal)
  expect_equal(
    as.vector(Z %*% beta_original),
    as.vector(X %*% beta_internal),
    tolerance = 1e-13
  )
  expect_equal(
    design_map$to_internal %*% design_map$to_original,
    diag(3),
    tolerance = 1e-14
  )
  expect_equal(design_map$intercept_idx, 2L)
  expect_equal(design_map$raw_constant, 2.5)

  V_beta <- matrix(
    c(
      0.40, 0.08, -0.03,
      0.08, 0.30, 0.04,
      -0.03, 0.04, 0.20
    ),
    nrow = 3,
    byrow = TRUE
  )
  V_original <- covariance_map_fn(V_beta, design_map$to_original)
  expect_lte(
    hl_fo_sdisc(
      V_original,
      design_map$to_original %*% V_beta %*% t(design_map$to_original)
    ),
    1e-14
  )
  expect_lte(
    hl_fo_sdisc(
      covariance_map_fn(V_original, design_map$to_internal),
      V_beta
    ),
    1e-13
  )

  tau <- log(0.4)
  xi <- log(8)
  V_canonical <- matrix(0, 5, 5)
  V_canonical[seq_len(3), seq_len(3)] <- V_beta
  V_canonical[4:5, 4:5] <- matrix(c(0.06, 0.01, 0.01, 0.08), 2, 2)
  V_canonical[1:3, 4:5] <- matrix(
    c(0.01, -0.02, 0.005, 0.01, -0.01, 0.015),
    3,
    2,
    byrow = TRUE
  )
  V_canonical[4:5, 1:3] <- t(V_canonical[1:3, 4:5])
  public <- public_map_fn(
    beta_internal,
    tau,
    xi,
    beta_map = design_map$to_original
  )
  V_public <- covariance_map_fn(V_canonical, public$jacobian)
  expect_lte(
    hl_fo_sdisc(
      V_public,
      public$jacobian %*% V_canonical %*% t(public$jacobian)
    ),
    1e-14
  )
  expect_equal(public$jacobian[4, 4], 0.4, tolerance = 1e-14)
  expect_equal(public$jacobian[5, 5], -1 / 8, tolerance = 1e-14)

  mean_one <- mean_one_fn(beta_internal, tau, xi, intercept_idx = intcol)
  canonical <- canonical_fn(
    mean_one$beta,
    mean_one$sigma2,
    mean_one$phi,
    intercept_idx = intcol
  )
  expect_equal(
    unname(mean_one$beta[intcol]),
    unname(beta_internal[intcol] + 0.4 / 2),
    tolerance = 1e-14
  )
  expect_equal(
    unname(canonical$value),
    unname(c(beta_internal, tau, xi)),
    tolerance = 1e-14
  )
  expect_lte(
    hl_fo_sdisc(canonical$jacobian %*% mean_one$jacobian, diag(5)),
    1e-14
  )
  expect_lte(
    hl_fo_sdisc(
      covariance_map_fn(
        covariance_map_fn(V_canonical, mean_one$jacobian),
        canonical$jacobian
      ),
      V_canonical
    ),
    1e-13
  )

  C_report <- rbind(
    first = c(1, 2, -1),
    intercept = c(0, 1, 0)
  )
  rhs <- c(first = 0.35, intercept = -0.6)
  transformed <- constraint_map_fn(
    C_report,
    rhs,
    design_map$to_original
  )
  expect_equal(
    as.vector(transformed$C %*% beta_internal - transformed$rhs),
    as.vector(C_report %*% beta_original - rhs),
    tolerance = 1e-14
  )
  internal_wald <- wald_fn(
    beta_internal,
    V_beta,
    transformed$C,
    transformed$rhs,
    n_subjects = 40L
  )
  report_wald <- wald_fn(
    beta_original,
    V_original,
    C_report,
    rhs,
    n_subjects = 40L
  )
  expect_true(internal_wald$success)
  expect_true(report_wald$success)
  expect_equal(internal_wald$statistic, report_wald$statistic, tolerance = 1e-12)
  expect_equal(internal_wald$p_value, report_wald$p_value, tolerance = 1e-12)
})

test_that("C1 reduced derivatives pass full-reprofile finite differences", {
  case <- .hl_fo_constraint_case("C1")
  affine <- .hl_fo_constraint_get("hl_fo_affine_constraint")(
    C = case$fixture$constraint$C,
    rhs = case$fixture$constraint$rhs,
    p = case$data$p
  )
  expect_true(affine$success)
  expect_equal(ncol(affine$K), 2L)
  rho <- case$seed$rho_full
  evaluate <- function(rho_active) {
    out <- .hl_fo_constraint_profile(
      case,
      affine,
      active = 1:2,
      rho_active = rho_active
    )
    if (!isTRUE(out$success)) {
      stop("Perturbed constrained profile failed: ", out$status, "/", out$reason)
    }
    out
  }
  state <- evaluate(rho)

  gradient_trace <- hl_fo_fd_trace(
    function(rho_active) evaluate(rho_active)$F,
    point = rho,
    coordinate = 1L,
    analytic = state$g[[1L]]
  )
  expect_true(
    gradient_trace$analytic_discrepancy <= 1e-5,
    info = .hl_fo_constraint_trace_info(
      "C1 reduced total gradient tau coordinate",
      gradient_trace
    )
  )
  expect_length(gradient_trace$steps, 6L)
  expect_length(gradient_trace$richardson, 5L)
  expect_identical(gradient_trace$candidate_ells, 1:2)

  hessian_trace <- hl_fo_fd_trace(
    function(rho_active) evaluate(rho_active)$g,
    point = rho,
    coordinate = 1L,
    analytic = state$H_F[, 1L]
  )
  expect_true(
    hessian_trace$analytic_discrepancy <= 1e-5,
    info = .hl_fo_constraint_trace_info(
      "C1 reduced exact Hessian tau column",
      hessian_trace
    )
  )
  expect_length(hessian_trace$steps, 6L)
  expect_length(hessian_trace$richardson, 5L)
  expect_identical(hessian_trace$candidate_ells, 1:2)
})

test_that("C1 records the true constraint residual at every outer evaluation", {
  case <- .hl_fo_constraint_case("C1")
  affine <- .hl_fo_constraint_get("hl_fo_affine_constraint")(
    case$fixture$constraint$C,
    case$fixture$constraint$rhs,
    p = case$data$p
  )
  evaluator <- .hl_fo_constraint_get("hl_fo_make_outer_evaluator")(
    data = case$data,
    rho_template = case$seed$rho_full,
    active = 1:2,
    beta_star = affine$beta_star,
    K = affine$K,
    seed_beta = case$seed$beta,
    constraint_C = case$fixture$constraint$C,
    constraint_rhs = case$fixture$constraint$rhs
  )
  proposals <- list(
    case$seed$rho_full,
    case$seed$rho_full + c(0.02, -0.01),
    case$seed$rho_full + c(-0.015, 0.02)
  )
  states <- lapply(proposals, evaluator$evaluate, cold = TRUE)
  expect_true(all(vapply(states, `[[`, logical(1), "success")))
  trace <- evaluator$trace()
  expect_length(trace, length(states))
  for (ii in seq_along(trace)) {
    actual <- as.numeric(
      case$fixture$constraint$C %*% states[[ii]]$beta -
        case$fixture$constraint$rhs
    )
    expect_length(trace[[ii]]$constraint_residual, nrow(case$fixture$constraint$C))
    expect_equal(trace[[ii]]$constraint_residual, actual, tolerance = 0)
    expect_lte(trace[[ii]]$constraint_residual_scaled, 100 * .Machine$double.eps)
    expect_length(trace[[ii]]$affine_reconstruction_residual, case$data$p)
  }
})

test_that("K equals I reproduces the unrestricted Q0 fit exactly", {
  case <- .hl_fo_constraint_case("Q0")
  fit_fn <- .hl_fo_constraint_get("hl_fo_fit")
  fit_affine_fn <- .hl_fo_constraint_get("hl_fo_fit_affine")
  fixed <- case$fixture$fixed_dispersion
  unrestricted <- fit_fn(
    data = case$data,
    active = character(0),
    fixed_natural = fixed,
    seed = case$seed
  )
  identity <- fit_affine_fn(
    data = case$data,
    beta_star = numeric(case$data$p),
    K = diag(case$data$p),
    active = integer(0),
    fixed_natural = fixed,
    seed = case$seed
  )
  expect_identical(identity$status, unrestricted$status)
  expect_identical(identity$fit_status, unrestricted$fit_status)
  expect_identical(identity$covariance_status, unrestricted$covariance_status)
  expect_identical(identity$coordinate_tag, unrestricted$coordinate_tag)
  for (field in c("beta", "eta", "rho", "rho_full", "g")) {
    expect_lte(hl_fo_sdisc(identity[[field]], unrestricted[[field]]), 1e-8)
  }
  for (field in c("F", "P")) {
    expect_lte(hl_fo_sdisc(identity[[field]], unrestricted[[field]]), 1e-8)
  }
  for (field in c("H_F")) {
    expect_lte(hl_fo_sdisc(identity[[field]], unrestricted[[field]]), 1e-8)
  }
  for (field in c("A", "J", "psi_subject", "augmented_residual")) {
    expect_lte(
      hl_fo_sdisc(identity$augmented[[field]], unrestricted$augmented[[field]]),
      1e-8
    )
  }
  expect_lte(
    hl_fo_sdisc(
      identity$covariance$V_canonical_full,
      unrestricted$covariance$V_canonical_full
    ),
    1e-8
  )
})

test_that("C1 and C2 are invariant to admissible alternative null bases", {
  for (id in c("C1", "C2")) {
    case <- .hl_fo_constraint_case(id)
    affine <- .hl_fo_constraint_get("hl_fo_affine_constraint")(
      case$fixture$constraint$C,
      case$fixture$constraint$rhs,
      p = case$data$p
    )
    k_free <- ncol(affine$K)
    O <- if (k_free == 2L) {
      matrix(c(0, -1, 1, 0), nrow = 2L)
    } else {
      matrix(-1, nrow = 1L)
    }
    K_alternative <- affine$K %*% O
    expect_equal(crossprod(K_alternative), diag(k_free), tolerance = 1e-14)
    expect_lte(
      max(abs(case$fixture$constraint$C %*% K_alternative)),
      1e-14
    )

    reference <- .hl_fo_constraint_profile(case, affine, K = affine$K)
    alternative <- .hl_fo_constraint_profile(case, affine, K = K_alternative)
    expect_true(reference$success, info = id)
    expect_true(alternative$success, info = id)
    for (field in c("beta", "eta", "F", "g", "H_F")) {
      expect_true(
        hl_fo_sdisc(reference[[field]], alternative[[field]]) <= 1e-7,
        info = paste(id, field)
      )
    }
    expect_true(
      hl_fo_sdisc(
        .hl_fo_constraint_mapped_profile_covariance(reference),
        .hl_fo_constraint_mapped_profile_covariance(alternative)
      ) <= 1e-7,
      info = paste(id, "mapped covariance")
    )
    rm(reference, alternative, case, affine)
    invisible(gc())
  }
})

test_that("C1 and C2 complete active constrained fits in reduced coordinates", {
  fit_fn <- .hl_fo_constraint_get("hl_fo_fit")
  for (id in c("C1", "C2")) {
    case <- .hl_fo_constraint_case(id)
    fit <- fit_fn(
      data = case$data,
      C = case$fixture$constraint$C,
      rhs = case$fixture$constraint$rhs,
      active = c("tau", "xi"),
      seed = case$seed
    )
    expect_true(
      fit$point_success,
      info = paste(id, fit$fit_status, fit$fit_failure_reason)
    )
    expect_true(
      fit$covariance_success,
      info = paste(id, fit$covariance_status, fit$covariance_failure_reason)
    )
    expect_identical(fit$fit_status, "ok", info = id)
    expect_identical(fit$covariance_status, "ok", info = id)
    expect_identical(fit$active, 1:2)
    expect_identical(ncol(fit$K), case$data$p - nrow(case$fixture$constraint$C))
    expect_lte(
      .hl_fo_constraint_get("hl_fo_constraint_residual")(
        case$fixture$constraint$C,
        fit$beta,
        case$fixture$constraint$rhs
      )$scaled,
      100 * .Machine$double.eps
    )
    successful_trace <- Filter(
      function(entry) isTRUE(entry$success),
      fit$outer_trace
    )
    expect_gt(length(successful_trace), 0L)
    expect_lte(
      max(vapply(
        successful_trace,
        `[[`,
        numeric(1),
        "constraint_residual_scaled"
      )),
      100 * .Machine$double.eps
    )
    rm(fit, case)
    invisible(gc())
  }
})

test_that("C1 reduced sandwich maps by blockdiag K and annihilates constraints", {
  case <- .hl_fo_constraint_case("C1")
  fit <- .hl_fo_constraint_get("hl_fo_fit")(
    data = case$data,
    C = case$fixture$constraint$C,
    rhs = case$fixture$constraint$rhs,
    active = character(0),
    fixed_natural = c(
      sigma2 = case$fixture$truth$sigma2,
      phi = case$fixture$truth$phi
    ),
    seed = case$seed
  )
  expect_true(fit$point_success)
  expect_true(fit$covariance_success)
  expect_identical(fit$fit_status, "ok")
  expect_identical(fit$covariance_status, "ok")
  expect_equal(dim(fit$augmented$A), rep(2L * ncol(fit$K), 2L))
  expect_identical(
    .hl_fo_constraint_get("hl_fo_matrix_rank")(fit$augmented$A),
    2L * ncol(fit$K)
  )

  T_K <- fit$K
  expected <- T_K %*% fit$covariance$V_gamma_rho %*% t(T_K)
  expect_lte(hl_fo_sdisc(fit$covariance$V_beta, expected), 1e-12)
  expect_lte(
    max(abs(case$fixture$constraint$C %*% fit$covariance$V_beta)),
    1e-10
  )
  expect_lte(
    .hl_fo_constraint_get("hl_fo_matrix_rank")(fit$covariance$V_beta),
    ncol(fit$K)
  )
  expect_lte(fit$certification$constraint_residual_scaled, 100 * .Machine$double.eps)
})

test_that("K0 covers k zero with active rho and k equals q equals zero", {
  case <- .hl_fo_constraint_case("K0")
  fit_fn <- .hl_fo_constraint_get("hl_fo_fit")
  C <- case$fixture$constraint$C
  rhs <- case$fixture$constraint$rhs

  active_fit <- fit_fn(
    data = case$data,
    C = C,
    rhs = rhs,
    active = c("tau", "xi"),
    seed = case$seed
  )
  expect_true(active_fit$point_success)
  expect_true(active_fit$covariance_success)
  expect_identical(ncol(active_fit$K), 0L)
  expect_equal(unname(active_fit$beta), unname(rhs), tolerance = 1e-14)
  expect_equal(active_fit$covariance$V_beta, matrix(0, 3, 3), tolerance = 0)
  expect_equal(dim(active_fit$augmented$A), c(2L, 2L))
  expect_equal(dim(active_fit$augmented$J), c(2L, 2L))
  expect_equal(dim(active_fit$covariance$V_gamma_rho), c(2L, 2L))
  expect_lte(
    max(vapply(
      Filter(function(entry) isTRUE(entry$success), active_fit$outer_trace),
      `[[`,
      numeric(1),
      "constraint_residual_scaled"
    )),
    100 * .Machine$double.eps
  )

  fixed_fit <- fit_fn(
    data = case$data,
    C = C,
    rhs = rhs,
    active = character(0),
    fixed_natural = c(
      sigma2 = case$fixture$truth$sigma2,
      phi = case$fixture$truth$phi
    ),
    seed = case$seed
  )
  expect_true(fixed_fit$point_success)
  expect_identical(ncol(fixed_fit$K), 0L)
  expect_length(fixed_fit$rho, 0L)
  expect_length(fixed_fit$g, 0L)
  expect_equal(dim(fixed_fit$H_F), c(0L, 0L))
  expect_equal(dim(fixed_fit$augmented$A), c(0L, 0L))
  expect_equal(dim(fixed_fit$augmented$J), c(0L, 0L))
  expect_false(fixed_fit$covariance_success)
  expect_identical(
    fixed_fit$covariance_status,
    "not_applicable_all_parameters_fixed"
  )
  expect_null(fixed_fit$covariance)
})

test_that("Q0 and Q1 retain zero- and one-dimensional outer contracts", {
  fit_fn <- .hl_fo_constraint_get("hl_fo_fit")
  q0 <- .hl_fo_constraint_case("Q0")
  q0_fit <- fit_fn(
    data = q0$data,
    active = character(0),
    fixed_natural = q0$fixture$fixed_dispersion,
    seed = q0$seed
  )
  expect_true(q0_fit$point_success)
  expect_true(q0_fit$covariance_success)
  expect_length(q0_fit$rho, 0L)
  expect_length(q0_fit$g, 0L)
  expect_equal(dim(q0_fit$H_F), c(0L, 0L))

  q1 <- .hl_fo_constraint_case("Q1")
  q1_fit <- fit_fn(
    data = q1$data,
    active = q1$fixture$active_dispersion,
    fixed_natural = q1$fixture$fixed_dispersion,
    seed = q1$seed
  )
  expect_true(
    q1_fit$point_success,
    info = paste(q1_fit$fit_status, q1_fit$fit_failure_reason)
  )
  expect_true(q1_fit$covariance_success)
  expect_identical(q1_fit$active, 2L)
  expect_length(q1_fit$rho, 1L)
  expect_length(q1_fit$g, 1L)
  expect_equal(dim(q1_fit$H_F), c(1L, 1L))
  expect_equal(q1_fit$rho_full[[1L]], q1$fixture$truth$tau, tolerance = 1e-14)
})

test_that("X0, W1, Z1, and D1 retain their declared rank and fixed-rho behavior", {
  fit_fn <- .hl_fo_constraint_get("hl_fo_fit")

  x0 <- .hl_fo_constraint_case("X0")
  x0_fit <- fit_fn(x0$data, seed = x0$seed)
  expect_false(x0_fit$point_success)
  expect_identical(x0_fit$fit_status, "inner_rank_deficient")

  w1 <- .hl_fo_constraint_case("W1")
  expect_identical(hl_fo_num_rank(w1$fixture$X)$rank, 3L)
  w1_fit <- fit_fn(
    w1$data,
    active = character(0),
    fixed_natural = c(sigma2 = 0.4, phi = 8),
    seed = w1$seed
  )
  expect_false(w1_fit$point_success)
  expect_identical(w1_fit$fit_status, "weak_curvature")

  z1 <- .hl_fo_constraint_case("Z1")
  z1_fit <- fit_fn(
    z1$data,
    active = character(0),
    fixed_natural = c(
      sigma2 = z1$fixture$truth$sigma2,
      phi = z1$fixture$truth$phi
    ),
    seed = z1$seed
  )
  expect_true(any(z1$fixture$y > 0))
  expect_true(z1_fit$point_success)
  expect_true(z1_fit$covariance_success)

  d1 <- hl_fo_fixture("D1")
  for (case_name in names(d1$cases)) {
    fixture <- d1$cases[[case_name]]
    data <- .hl_fo_constraint_get("hl_fo_make_data")(
      fixture$y,
      fixture$X,
      fixture$log_pi,
      fixture$subject_index,
      fixture$intcol
    )
    seed <- list(
      beta = fixture$truth$beta,
      rho_full = c(fixture$truth$tau, fixture$truth$xi),
      coordinate_tag = "canonical_conditional_beta_tau_xi"
    )
    fit <- fit_fn(
      data,
      active = character(0),
      fixed_natural = c(
        sigma2 = fixture$truth$sigma2,
        phi = fixture$truth$phi
      ),
      seed = seed
    )
    expect_true(fit$point_success, info = case_name)
    expect_true(fit$covariance_success, info = case_name)
    expect_identical(fit$fit_status, "ok", info = case_name)
  }
})

test_that("B1 and S1 surface explicit numerical or small-sample limitations", {
  fit_fn <- .hl_fo_constraint_get("hl_fo_fit")

  b1 <- .hl_fo_constraint_case("B1")
  standard_lower <- log(c(1e-4, 1e-4))
  standard_upper <- log(c(10, 1000))
  expect_true(any(b1$seed$rho_full < standard_lower))
  expect_true(any(b1$seed$rho_full > standard_upper))
  b1_fixed <- fit_fn(
    b1$data,
    active = character(0),
    fixed_natural = c(
      sigma2 = b1$fixture$truth$sigma2,
      phi = b1$fixture$truth$phi
    ),
    seed = b1$seed
  )
  if (isTRUE(b1_fixed$point_success)) {
    expect_true(all(is.finite(c(b1_fixed$beta, b1_fixed$eta, b1_fixed$F))))
  } else {
    expect_true(
      b1_fixed$fit_status %in% c(
        "mode_failure",
        "inner_no_convergence",
        "inner_no_finite_root",
        "weak_curvature"
      )
    )
    expect_false(isTRUE(b1_fixed$covariance_success))
  }

  s1 <- .hl_fo_constraint_case("S1")
  s1_fit <- fit_fn(s1$data, seed = s1$seed)
  expect_identical(s1_fit$data$m, 6L)
  expect_false(is.null(s1_fit$certification$J_rank_bound))
  expect_identical(
    s1_fit$certification$J_rank_bound,
    min(2L * s1$data$p + 2L, s1$data$m - 1L)
  )
  if (isTRUE(s1_fit$point_success)) {
    expect_lte(
      s1_fit$certification$J_rank,
      min(nrow(s1_fit$augmented$J), s1$data$m - 1L)
    )
    expect_false(isTRUE(s1_fit$covariance$sandwich$ridge > 0))
  } else {
    expect_true(
      s1_fit$fit_status %in% c(
        "weak_curvature",
        "sensitivity_singular",
        "covariance_invalid",
        "outer_no_stationary_root"
      )
    )
    expect_false(isTRUE(s1_fit$covariance_success))
  }
})
