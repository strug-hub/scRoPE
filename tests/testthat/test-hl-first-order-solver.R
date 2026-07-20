.hl_fo_solver <- function(name) {
  getFromNamespace(name, "scRoPE")
}

.hl_fo_solver_trace_info <- function(label, trace) {
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

.hl_fo_solver_expect_trace <- function(trace, tolerance, label) {
  expect_true(
    is.finite(trace$analytic_discrepancy),
    info = .hl_fo_solver_trace_info(label, trace)
  )
  expect_true(
    trace$analytic_discrepancy <= tolerance,
    info = .hl_fo_solver_trace_info(label, trace)
  )
  invisible(trace)
}

.hl_fo_solver_case <- function() {
  generated <- hl_fo_with_rng(271828L, function() {
    m <- 16L
    n_i <- 4L
    n <- m * n_i
    subject <- rep(seq_len(m), each = n_i)
    x <- rep(seq(-1, 1, length.out = m), each = n_i) +
      rep(c(-0.3, -0.1, 0.1, 0.3), m)
    X <- cbind(`(Intercept)` = 1, x = x)
    log_offset <- rep(log(c(0.8, 0.95, 1.1, 1.25)), m)
    beta <- c(`(Intercept)` = -0.55, x = 0.3)
    sigma2 <- 0.45
    phi <- 6
    alpha <- 1 / expm1(sigma2)
    lambda <- alpha * exp(-sigma2 / 2)
    omega <- rgamma(m, shape = alpha, rate = lambda)
    mu <- exp(log_offset + drop(X %*% beta)) * omega[subject]
    y <- rnbinom(n, size = phi, mu = mu)
    list(
      y = y,
      X = X,
      log_offset = log_offset,
      subject = subject,
      beta = beta,
      rho_full = log(c(sigma2, phi))
    )
  })
  generated$data <- .hl_fo_solver("hl_fo_make_data")(
    y = generated$y,
    X = generated$X,
    log_offset = generated$log_offset,
    subject = generated$subject,
    intcol = 1L
  )
  generated
}

.hl_fo_solver_profile <- function(case, rho_active = case$rho_full,
                                  beta_start = c(4, -4), eta_start = NULL,
                                  cold = TRUE) {
  out <- .hl_fo_solver("hl_fo_profile_evaluate")(
    data = case$data,
    rho_active = rho_active,
    rho_template = case$rho_full,
    active = 1:2,
    beta_star = numeric(case$data$p),
    K = diag(case$data$p),
    beta_start = beta_start,
    eta_start = eta_start,
    cold = cold
  )
  if (!isTRUE(out$success)) {
    stop("Small deterministic profile failed: ", out$status, "/", out$reason)
  }
  out
}

# Evaluate the augmented estimating equations at an arbitrary independent
# zeta = (gamma, chi, rho). Modes are re-solved, but neither the ordinary-score
# beta root nor chi is re-solved. Thus finite differences of this function are
# immediately below the analytic augmented Jacobian being tested.
.hl_fo_solver_psi_at <- function(case, zeta, beta_star, K, active,
                                 rho_template) {
  k_free <- ncol(K)
  q <- length(active)
  gamma <- if (k_free) zeta[seq_len(k_free)] else numeric(0)
  chi <- if (k_free) zeta[k_free + seq_len(k_free)] else numeric(0)
  rho_active <- if (q) zeta[2L * k_free + seq_len(q)] else numeric(0)
  rho_full <- as.numeric(rho_template)
  if (q) {
    rho_full[active] <- rho_active
  }
  beta <- as.numeric(beta_star + if (k_free) K %*% gamma else 0)

  modes <- .hl_fo_solver("hl_fo_solve_modes")(
    case$data,
    beta,
    rho_full,
    eta_start = NULL
  )
  if (!isTRUE(modes$success)) {
    stop("Perturbed mode solve failed: ", modes$status, "/", modes$reason)
  }
  transform <- .hl_fo_solver("hl_fo_scientific_transform")(
    case$data$p,
    K,
    active
  )
  gamma_idx <- if (k_free) seq_len(k_free) else integer(0)
  rho_idx <- if (q) k_free + seq_len(q) else integer(0)
  psi_subject <- matrix(0, 2L * k_free + q, case$data$m)

  for (ii in seq_len(case$data$m)) {
    rows <- case$data$subject_rows[[ii]]
    subject <- .hl_fo_solver("hl_fo_subject_raw")(
      y = case$data$y[rows],
      X = case$data$X[rows, , drop = FALSE],
      log_offset = case$data$log_offset[rows],
      beta = beta,
      eta = modes$eta[[ii]],
      rho_full = rho_full,
      include_laplace_constant = TRUE
    )
    if (!isTRUE(subject$ok)) {
      stop("Perturbed subject derivative failed: ", subject$reason)
    }
    subject <- .hl_fo_solver("hl_fo_contract_subject")(subject, transform)
    psi <- numeric(2L * k_free + q)
    if (k_free) {
      B_i <- subject$p_hess[gamma_idx, gamma_idx, drop = FALSE]
      psi[seq_len(k_free)] <- subject$p_grad[gamma_idx]
      psi[k_free + seq_len(k_free)] <-
        as.numeric(t(B_i) %*% chi) - subject$ell_grad[gamma_idx]
    }
    if (q) {
      E_i <- if (k_free) {
        subject$p_hess[gamma_idx, rho_idx, drop = FALSE]
      } else {
        matrix(0, 0, q)
      }
      psi[2L * k_free + seq_len(q)] <- subject$ell_grad[rho_idx] -
        if (k_free) as.numeric(crossprod(E_i, chi)) else 0
    }
    psi_subject[, ii] <- psi
  }
  psi_subject
}

test_that("scalar and vectorized safeguarded mode solves agree after bisection", {
  case <- .hl_fo_solver_case()
  mode_one <- .hl_fo_solver("hl_fo_solve_mode")
  mode_all <- .hl_fo_solver("hl_fo_solve_modes")
  poor_start <- -20

  vectorized <- mode_all(
    case$data,
    case$beta,
    case$rho_full,
    eta_start = rep(poor_start, case$data$m)
  )
  expect_true(vectorized$success)
  scalar <- lapply(seq_len(case$data$m), function(ii) {
    rows <- case$data$subject_rows[[ii]]
    mode_one(
      y = case$data$y[rows],
      log_m_base = as.numeric(
        case$data$log_offset[rows] +
          case$data$X[rows, , drop = FALSE] %*% case$beta
      ),
      tau = case$rho_full[[1L]],
      xi = case$rho_full[[2L]],
      start = poor_start
    )
  })
  expect_true(all(vapply(scalar, `[[`, logical(1), "success")))
  scalar_eta <- vapply(scalar, `[[`, numeric(1), "eta")
  expect_equal(vectorized$eta, scalar_eta, tolerance = 1e-12)
  expect_gt(sum(vapply(scalar, `[[`, integer(1), "bisections")), 0L)
  expect_gt(sum(vapply(vectorized$diagnostics, `[[`, integer(1), "bisections")), 0L)

  scalar_residual <- vapply(scalar, `[[`, numeric(1), "residual_scaled")
  scalar_step <- vapply(scalar, `[[`, numeric(1), "step_scaled")
  vector_residual <- vapply(
    vectorized$diagnostics,
    `[[`,
    numeric(1),
    "residual_scaled"
  )
  vector_step <- vapply(
    vectorized$diagnostics,
    `[[`,
    numeric(1),
    "step_scaled"
  )
  expect_lte(max(c(scalar_residual, vector_residual)), 1e-10)
  expect_lte(max(c(scalar_step, vector_step)), 1e-10)
  expect_true(all(vapply(scalar, function(x) {
    x$eta >= x$bracket[[1L]] && x$eta <= x$bracket[[2L]] &&
      x$h_etaeta < 0
  }, logical(1))))
})

test_that("profiled ordinary-score Newton uses Armijo ascent and certifies its root", {
  case <- .hl_fo_solver_case()
  inner <- .hl_fo_solver("hl_fo_solve_inner")(
    data = case$data,
    rho_full = case$rho_full,
    beta_star = numeric(case$data$p),
    K = diag(case$data$p),
    beta_start = c(4, -4),
    eta_start = NULL,
    cold = TRUE
  )
  expect_true(inner$success)
  expect_identical(inner$status, "ok")
  expect_lte(inner$score_metric, 1e-8)
  expect_lte(inner$step_metric, 1e-10)
  expect_lte(max(abs(inner$S)), 1e-8 * (1 + max(abs(inner$S_beta))))

  trace_P <- vapply(inner$trace, `[[`, numeric(1), "P")
  expect_true(all(diff(trace_P) >= -1e-12 * (1 + abs(trace_P[-1L]))))
  armijo_steps <- vapply(inner$trace, function(entry) {
    if (is.null(entry$armijo_step)) NA_real_ else entry$armijo_step
  }, numeric(1))
  expect_true(any(armijo_steps < 1, na.rm = TRUE))
  expect_true(all(armijo_steps[is.finite(armijo_steps)] > 0))
  expect_true(all(armijo_steps[is.finite(armijo_steps)] <= 1))
  expect_gt(min(eigen(-0.5 * (inner$B + t(inner$B)),
    symmetric = TRUE,
    only.values = TRUE
  )$values), 0)
})

test_that("fully reprofiled total gradient and exact Hessian pass step halving", {
  case <- .hl_fo_solver_case()
  state <- .hl_fo_solver_profile(case)
  evaluate <- function(rho) .hl_fo_solver_profile(case, rho_active = rho)

  expect_identical(state$inner$acceptance_tier_requested, "derivative")
  expect_identical(state$inner$accepted_tier, "derivative")
  expect_lte(
    .hl_fo_solver("hl_fo_asymmetry")(state$H_F),
    1e-10
  )
  for (rr in seq_along(state$rho)) {
    gradient_trace <- hl_fo_fd_trace(
      function(rho) evaluate(rho)$F,
      point = state$rho,
      coordinate = rr,
      analytic = state$g[[rr]]
    )
    .hl_fo_solver_expect_trace(
      gradient_trace,
      1e-6,
      sprintf("fully reprofiled total gradient coordinate %d", rr)
    )

    hessian_trace <- hl_fo_fd_trace(
      function(rho) evaluate(rho)$g,
      point = state$rho,
      coordinate = rr,
      analytic = state$H_F[, rr]
    )
    .hl_fo_solver_expect_trace(
      hessian_trace,
      1e-6,
      sprintf("exact reprofiled Hessian column %d", rr)
    )
  }
})

test_that("complete subject and summed augmented Jacobians pass finite differences", {
  case <- .hl_fo_solver_case()
  state <- .hl_fo_solver_profile(case)
  k_free <- ncol(state$K)
  q <- length(state$active)
  zeta <- c(state$gamma, state$augmented$chi, state$rho)
  psi_at <- function(z) {
    .hl_fo_solver_psi_at(
      case = case,
      zeta = z,
      beta_star = state$beta_star,
      K = state$K,
      active = state$active,
      rho_template = state$rho_full
    )
  }
  baseline_psi <- psi_at(zeta)
  expect_lte(
    hl_fo_sdisc(baseline_psi, state$augmented$psi_subject),
    1e-9
  )

  d_aug <- 2L * k_free + q
  for (rr in seq_len(d_aug)) {
    analytic_subject <- vapply(
      state$augmented$jacobian_subject,
      function(x) x[, rr],
      numeric(d_aug)
    )
    subject_trace <- hl_fo_fd_trace(
      psi_at,
      point = zeta,
      coordinate = rr,
      analytic = analytic_subject
    )
    .hl_fo_solver_expect_trace(
      subject_trace,
      1e-6,
      sprintf("all subject augmented Jacobians, coordinate %d", rr)
    )

    summed_analytic <- rowSums(analytic_subject)
    expect_equal(
      summed_analytic,
      -state$augmented$A[, rr],
      tolerance = 1e-12
    )
    summed_trace <- hl_fo_fd_trace(
      function(z) rowSums(psi_at(z)),
      point = zeta,
      coordinate = rr,
      analytic = summed_analytic
    )
    .hl_fo_solver_expect_trace(
      summed_trace,
      1e-6,
      sprintf("summed augmented Jacobian coordinate %d", rr)
    )
  }
  expect_equal(
    state$augmented$J,
    tcrossprod(state$augmented$psi_subject),
    tolerance = 1e-12
  )
  expect_lte(
    .hl_fo_solver("hl_fo_asymmetry")(state$augmented$J),
    1e-12
  )
})

test_that("solve-based sandwich equals direct A inverse multiplication without 1/m", {
  case <- .hl_fo_solver_case()
  state <- .hl_fo_solver_profile(case)
  sandwich <- .hl_fo_solver("hl_fo_sandwich")(
    state$augmented$A,
    state$augmented$J
  )
  expect_true(sandwich$success)
  A_inverse <- solve(state$augmented$A)
  direct <- A_inverse %*% state$augmented$J %*% t(A_inverse)
  expect_lte(hl_fo_sdisc(sandwich$V, direct), 1e-10)
  expect_gt(hl_fo_sdisc(sandwich$V, direct / case$data$m), 1e-3)
  expect_lte(sandwich$asymmetry, 1e-10)
  expect_gte(min(sandwich$eigenvalues), -1e-8 * max(1, max(sandwich$eigenvalues)))
})

test_that("trust-region secular solve remains bracketed near the hard case", {
  trust_step <- .hl_fo_solver("hl_fo_trust_step_max")

  # The positive-eigenvalue multiplier differs from lambda_max by 1e-14.
  # Adding that delta to lambda_max rounds too coarsely for the former
  # multiplier-scale uniroot bracket, even though the global boundary step is
  # well defined.
  g <- c(1e-12, 0)
  H <- diag(c(1, -1))
  radius <- 100
  step <- trust_step(g, H, radius)
  predicted <- sum(g * step) + 0.5 * sum(step * (H %*% step))

  expect_true(all(is.finite(step)))
  expect_equal(sqrt(sum(step^2)), radius, tolerance = 1e-12)
  expect_equal(step, c(radius, 0), tolerance = 1e-10)
  expect_gt(predicted, 0)

  # Locked zero-heavy search state that originally exposed the same failure
  # with a negative leading eigenvalue and a root between zero and the former
  # absolute multiplier offset.
  fixture_g <- c(-1.0249607376104438e-08, 2.3139250492635245e-08)
  fixture_H <- matrix(c(
    -1.0249586865134153e-08, 1.3174293803383386e-07,
     1.3174293803383386e-07, -8.8116681014273741
  ), nrow = 2L)
  fixture_step <- trust_step(fixture_g, fixture_H, radius = 1)
  fixture_predicted <- sum(fixture_g * fixture_step) +
    0.5 * sum(fixture_step * (fixture_H %*% fixture_step))
  expect_true(all(is.finite(fixture_step)))
  expect_equal(sqrt(sum(fixture_step^2)), 1, tolerance = 1e-14)
  expect_gt(fixture_predicted, 0)

  # With no signal in a repeated leading eigenspace, retain the declared
  # lexicographically smallest global solution.
  hard_step <- trust_step(c(0, 0), diag(2), radius = 1)
  expect_equal(hard_step, c(-1, 0), tolerance = 0)
})

test_that("q zero has no dummy outer coordinate and retains a certified beta fit", {
  case <- .hl_fo_solver_case()
  seed <- list(
    beta = c(4, -4),
    rho_full = case$rho_full,
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  fit <- .hl_fo_solver("hl_fo_fit_affine")(
    data = case$data,
    beta_star = numeric(case$data$p),
    K = diag(case$data$p),
    active = integer(0),
    fixed_natural = c(sigma2 = exp(case$rho_full[[1L]]),
      phi = exp(case$rho_full[[2L]])),
    seed = seed
  )
  expect_true(fit$success)
  expect_true(fit$point_success)
  expect_identical(fit$fit_status, "ok")
  expect_length(fit$rho, 0L)
  expect_length(fit$active, 0L)
  expect_length(fit$g, 0L)
  expect_equal(dim(fit$H_F), c(0L, 0L))
  expect_equal(fit$rho_full, case$rho_full, tolerance = 1e-14)
  expect_equal(dim(fit$augmented$A), rep(2L * case$data$p, 2L))
  expect_equal(dim(fit$augmented$J), rep(2L * case$data$p, 2L))
  expect_equal(fit$cold_replay$beta_discrepancy, 0)
  expect_equal(fit$cold_replay$F_discrepancy, 0)
  expect_true(fit$covariance_success)
  expect_equal(dim(fit$covariance$V_beta), rep(case$data$p, 2L))
})

test_that("rank deficiency and nonattainment retain distinct explicit statuses", {
  make_fixture_data <- function(id) {
    fixture <- hl_fo_fixture(id)
    list(
      fixture = fixture,
      data = .hl_fo_solver("hl_fo_make_data")(
        y = fixture$y,
        X = fixture$X,
        log_offset = fixture$log_pi,
        subject = fixture$subject_index,
        intcol = fixture$intcol
      )
    )
  }

  x0 <- make_fixture_data("X0")
  x0_seed <- list(
    beta = x0$fixture$truth$beta,
    rho_full = c(x0$fixture$truth$tau, x0$fixture$truth$xi),
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  x0_fit <- .hl_fo_solver("hl_fo_fit_affine")(x0$data, seed = x0_seed)
  expect_false(x0_fit$success)
  expect_false(x0_fit$point_success)
  expect_identical(x0_fit$status, "inner_rank_deficient")
  expect_identical(x0_fit$fit_status, "inner_rank_deficient")
  expect_identical(x0_fit$reason, "analysis_design_rank_deficient")

  n0 <- make_fixture_data("N0")
  n0_seed <- list(
    beta = 0,
    rho_full = c(n0$fixture$truth$tau, n0$fixture$truth$xi),
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  n0_fit <- .hl_fo_solver("hl_fo_fit_affine")(n0$data, seed = n0_seed)
  expect_false(n0_fit$success)
  expect_false(n0_fit$point_success)
  expect_identical(n0_fit$status, "inner_no_finite_root")
  expect_identical(n0_fit$fit_status, "inner_no_finite_root")
  expect_identical(n0_fit$reason, "all_zero_intercept_recession")
})

test_that("a separated nonzero-count inner problem is explicit nonattainment", {
  data <- .hl_fo_solver("hl_fo_make_data")(
    y = c(1, 2, 0, 0, 0, 0),
    X = cbind(`(Intercept)` = 1, separated = c(0, 0, 1, 1, 1, 1)),
    log_offset = numeric(6),
    subject = paste0("sep", seq_len(6))
  )
  seed <- list(
    beta = c(0, 0),
    rho_full = log(c(0.4, 8)),
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )
  fit <- .hl_fo_solver("hl_fo_fit_affine")(
    data,
    active = integer(0),
    fixed_natural = c(sigma2 = 0.4, phi = 8),
    seed = seed
  )
  expect_false(fit$success)
  expect_false(fit$point_success)
  expect_identical(fit$status, "inner_no_finite_root")
  expect_identical(fit$reason, "separated_inner_recession")
  expect_true(fit$recession$recession)
})

test_that("K equals I is exact and cold replay reproduces all derivative states", {
  case <- .hl_fo_solver_case()
  reference <- .hl_fo_solver_profile(case, cold = TRUE)
  warm <- .hl_fo_solver_profile(
    case,
    beta_start = reference$beta,
    eta_start = reference$eta,
    cold = FALSE
  )
  cold <- .hl_fo_solver_profile(
    case,
    beta_start = c(4, -4),
    eta_start = NULL,
    cold = TRUE
  )
  expect_lte(hl_fo_sdisc(warm$beta, cold$beta), 1e-10)
  expect_lte(hl_fo_sdisc(warm$eta, cold$eta), 1e-10)
  expect_lte(hl_fo_sdisc(warm$F, cold$F), 1e-10)
  expect_lte(hl_fo_sdisc(warm$g, cold$g), 1e-10)
  expect_lte(hl_fo_sdisc(warm$H_F, cold$H_F), 1e-10)
  expect_lte(hl_fo_sdisc(warm$augmented$A, cold$augmented$A), 1e-10)
  expect_lte(hl_fo_sdisc(warm$augmented$J, cold$augmented$J), 1e-10)

  # With K = I and both dispersions active, the scientific transform itself
  # is exactly I. The contracted route must therefore reproduce every raw
  # subject contribution and its assembled augmented state without tolerance
  # loss from a coordinate change.
  transform <- .hl_fo_solver("hl_fo_scientific_transform")(
    case$data$p,
    diag(case$data$p),
    1:2
  )
  expect_identical(transform, diag(case$data$p + 2L))
  raw_subjects <- lapply(seq_len(case$data$m), function(ii) {
    rows <- case$data$subject_rows[[ii]]
    .hl_fo_solver("hl_fo_subject_raw")(
      y = case$data$y[rows],
      X = case$data$X[rows, , drop = FALSE],
      log_offset = case$data$log_offset[rows],
      beta = cold$beta,
      eta = cold$eta[[ii]],
      rho_full = cold$rho_full,
      include_laplace_constant = TRUE
    )
  })
  fields <- c("p", "ell", "p_grad", "p_hess", "p_third", "ell_grad", "ell_hess")
  for (ii in seq_len(case$data$m)) {
    expect_true(raw_subjects[[ii]]$ok)
    for (field in fields) {
      expect_identical(cold$subjects[[ii]][[field]], raw_subjects[[ii]][[field]])
    }
  }
  direct_augmented <- .hl_fo_solver("hl_fo_augmented_state")(
    raw_subjects,
    k_free = case$data$p,
    q = 2L
  )
  expect_true(direct_augmented$success)
  expect_identical(direct_augmented$g, cold$augmented$g)
  expect_identical(direct_augmented$H_F, cold$augmented$H_F)
  expect_identical(direct_augmented$A, cold$augmented$A)
  expect_identical(direct_augmented$J, cold$augmented$J)
})

test_that("affine fitting rejects untagged or inconsistent coordinates before evaluation", {
  case <- .hl_fo_solver_case()
  fit_affine <- .hl_fo_solver("hl_fo_fit_affine")
  canonical_seed <- list(
    beta = case$beta,
    rho_full = case$rho_full,
    coordinate_tag = "canonical_conditional_beta_tau_xi"
  )

  expect_error(
    fit_affine(case$data, K = diag(c(1, 2)), seed = canonical_seed),
    "orthonormal"
  )
  expect_error(
    fit_affine(
      case$data,
      beta_star = c(0, 0),
      K = matrix(c(0, 1), nrow = 2L),
      constraint_C = matrix(c(1, 0), nrow = 1L),
      constraint_rhs = 1,
      seed = canonical_seed
    ),
    "do not satisfy"
  )
  expect_error(
    fit_affine(
      case$data,
      seed = list(beta = case$beta, rho_full = case$rho_full)
    ),
    "Untagged"
  )
  expect_error(
    fit_affine(
      case$data,
      seed = list(
        beta = case$beta[-1L],
        rho_full = case$rho_full,
        coordinate_tag = "canonical_conditional_beta_tau_xi"
      )
    ),
    "finite beta"
  )
  expect_error(
    fit_affine(
      case$data,
      active = 1L,
      fixed_natural = c(sigma2 = 0.4),
      seed = canonical_seed
    ),
    "both active and fixed"
  )
  expect_error(
    fit_affine(
      case$data,
      active = integer(0),
      fixed_natural = c(bad = 0.4),
      seed = canonical_seed
    ),
    "named sigma2/phi"
  )
})

test_that("terminal boundary labels contain only actual active coordinates", {
  faces <- .hl_fo_solver("hl_fo_terminal_boundary_faces")(
    list(
      reason = "outer_radius_collapse",
      invalid_guard_faces = c("2:upper"),
      state = list(
        success = TRUE,
        rho = c(0, 1 - 1e-8),
        g = c(-1, 1)
      )
    ),
    lower = c(-1, -1),
    upper = c(1, 1)
  )
  expect_identical(faces, "2:upper")
  expect_false(any(faces %in% c(":lower", ":upper")))
})

test_that("recorded outward guard crossings survive a slightly interior terminal state", {
  boundary_faces <- .hl_fo_solver("hl_fo_terminal_boundary_faces")
  lower <- c(-13.210340372, -9.210340372)
  upper <- c(2.302585093, 10.907755279)
  run <- list(
    reason = "outer_radius_collapse",
    invalid_guard_faces = c("1:lower", "2:upper"),
    state = list(
      success = TRUE,
      rho = c(-13.210223, 10.860289),
      g = c(-8.1e-5, 7.3e-5)
    )
  )

  # Both relative distances exceed 1e-6.  The rejected crossings plus the
  # terminal outward gradient invoke the alternative guard-expansion rule.
  expect_true(all(c(
    (run$state$rho[1L] - lower[1L]) / (upper[1L] - lower[1L]),
    (upper[2L] - run$state$rho[2L]) / (upper[2L] - lower[2L])
  ) > 1e-6))
  expect_identical(
    boundary_faces(run, lower = lower, upper = upper),
    c("1:lower", "2:upper")
  )

  run$state$g <- -run$state$g
  expect_identical(
    boundary_faces(run, lower = lower, upper = upper),
    character(0)
  )
})

test_that("B1-scale fixed-dispersion inner root is certified without clipping", {
  fixture <- hl_fo_fixture("B1")
  data <- .hl_fo_solver("hl_fo_make_data")(
    fixture$y,
    fixture$X,
    fixture$log_pi,
    fixture$subject_index,
    fixture$intcol
  )
  affine <- .hl_fo_solver("hl_fo_affine_constraint")(NULL, p = data$p)
  inner <- .hl_fo_solver("hl_fo_solve_inner")(
    data = data,
    rho_full = c(fixture$truth$tau, fixture$truth$xi),
    beta_star = affine$beta_star,
    K = affine$K,
    beta_start = fixture$truth$beta,
    cold = TRUE
  )
  expect_true(inner$success)
  expect_identical(inner$status, "ok")
  expect_lte(inner$score_metric, 1e-10)
  expect_lte(inner$step_metric, 1e-10)
  expect_true(all(vapply(inner$mode_diagnostics, function(diag) {
    diag$residual_scaled <= 1e-10 && diag$step_scaled <= 1e-10
  }, logical(1))))
})

test_that("inner statistical acceptance is specification-gated at the cap", {
  decide <- .hl_fo_solver("hl_fo_inner_acceptance_decision")

  before_cap <- decide(
    score_metric = 5e-10, step_metric = 2e-11,
    iteration = 99L, max_iter = 100L,
    acceptance_tier = "statistical"
  )
  expect_false(before_cap$accepted)
  expect_true(before_cap$statistical_grade)

  at_cap <- decide(
    score_metric = 5e-10, step_metric = 2e-11,
    iteration = 100L, max_iter = 100L,
    acceptance_tier = "statistical"
  )
  expect_true(at_cap$accepted)
  expect_identical(at_cap$accepted_tier, "statistical")
  expect_false(at_cap$derivative_grade)

  derivative <- decide(
    score_metric = 5e-10, step_metric = 2e-11,
    iteration = 100L, max_iter = 100L,
    acceptance_tier = "derivative"
  )
  expect_false(derivative$accepted)

  outside_specification <- decide(
    score_metric = 2e-8, step_metric = 2e-11,
    iteration = 100L, max_iter = 100L,
    acceptance_tier = "statistical"
  )
  expect_false(outside_specification$accepted)

  polished <- decide(
    score_metric = 5e-11, step_metric = 2e-11,
    iteration = 3L, max_iter = 100L,
    acceptance_tier = "statistical"
  )
  expect_true(polished$accepted)
  expect_identical(polished$accepted_tier, "derivative")
})
