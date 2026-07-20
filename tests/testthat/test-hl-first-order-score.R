.hl_fo_score_ns <- function(name) {
  getFromNamespace(name, "scRoPE")
}

.hl_fo_score_c1 <- function() {
  fixture <- hl_fo_fixture("C1")
  data <- .hl_fo_score_ns("hl_fo_make_data")(
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

.hl_fo_score_locked_mc_case <- function(replicate = 75L, m = 250L) {
  generated <- hl_fo_with_rng(8675309L, function() {
    stream <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (replicate > 1L) {
      for (ii in 2:replicate) stream <- parallel::nextRNGStream(stream)
    }
    assign(".Random.seed", stream, envir = .GlobalEnv)

    m_max <- 2000L
    n_i <- 5L
    n <- m_max * n_i
    subject <- rep(seq_len(m_max), each = n_i)
    z_subject <- rbinom(m_max, size = 1L, prob = 0.5)
    z <- rep(z_subject, each = n_i)
    w <- rnorm(n)
    X <- cbind(`(Intercept)` = 1, z = z, w = w)
    log_offset <- runif(n, min = -0.3, max = 0.3)
    beta <- c(`(Intercept)` = -1.2, z = 0.35, w = -0.25)
    sigma2 <- 0.4
    phi <- 8
    alpha <- 1 / expm1(sigma2)
    lambda <- alpha * exp(-sigma2 / 2)
    omega <- rgamma(m_max, shape = alpha, rate = lambda)
    mu <- exp(log_offset + drop(X %*% beta)) * omega[subject]
    y <- rnbinom(n, size = phi, mu = mu)
    list(
      y = y, X = X, log_offset = log_offset, subject = subject,
      beta = beta, sigma2 = sigma2, phi = phi, n_i = n_i
    )
  })
  cells <- seq_len(as.integer(m) * generated$n_i)
  generated$data <- .hl_fo_score_ns("hl_fo_make_data")(
    y = generated$y[cells],
    X = generated$X[cells, , drop = FALSE],
    log_offset = generated$log_offset[cells],
    subject = generated$subject[cells],
    intcol = 1L
  )
  generated$R <- matrix(
    c(0, 1, 0), nrow = 1L,
    dimnames = list("z_at_truth", names(generated$beta))
  )
  generated$rhs <- c(z_at_truth = generated$beta[["z"]])
  generated
}

.hl_fo_score_trace_info <- function(label, trace) {
  paste(
    sprintf(
      "%s: selected_ell=%d analytic_discrepancy=%.17g stable=%s",
      label, trace$selected_ell, trace$analytic_discrepancy,
      trace$stable_minimum
    ),
    paste(capture.output(dput(trace)), collapse = "\n"),
    sep = "\nfull step-halving trace:\n"
  )
}

test_that("full-chi root Jacobian matches fully reprofiled finite differences", {
  case <- .hl_fo_score_c1()
  R <- case$fixture$constraint$C
  rhs <- case$fixture$constraint$rhs
  affine <- .hl_fo_score_ns("hl_fo_affine_constraint")(
    R, rhs, p = case$data$p
  )
  rho <- case$seed$rho_full
  state_at <- .hl_fo_score_ns("hl_fo_score_state_at_rho")
  baseline <- state_at(
    data = case$data,
    rho_active = rho,
    rho_template = rho,
    active = 1:2,
    beta_star = affine$beta_star,
    K = affine$K,
    beta_start = case$seed$beta,
    cold = TRUE,
    constraint_R = R,
    constraint_rhs = rhs
  )
  expect_true(baseline$success)
  expect_equal(dim(baseline$G0), c(2L, 2L))

  perturbed <- list()
  lower_order <- function(proposed_rho) {
    evaluated <- state_at(
      data = case$data,
      rho_active = proposed_rho,
      rho_template = rho,
      active = 1:2,
      beta_star = affine$beta_star,
      K = affine$K,
      beta_start = case$seed$beta,
      cold = TRUE,
      constraint_R = R,
      constraint_rhs = rhs
    )
    if (!isTRUE(evaluated$success)) {
      stop("Perturbed full-chi state failed: ", evaluated$status, "/",
           evaluated$reason)
    }
    perturbed[[length(perturbed) + 1L]] <<- evaluated
    evaluated$g0
  }
  fd_eps <- c(1e-6, 1e-11)
  for (jj in 1:2) {
    trace <- hl_fo_fd_trace(
      lower_order = lower_order,
      point = rho,
      coordinate = jj,
      analytic = baseline$G0[, jj],
      eps = fd_eps[[jj]]
    )
    expect_true(
      trace$stable_minimum,
      info = .hl_fo_score_trace_info(paste0("G0 column ", jj), trace)
    )
    expect_true(
      trace$analytic_discrepancy <= 2e-6,
      info = .hl_fo_score_trace_info(paste0("G0 column ", jj), trace)
    )
  }
  # Each finite-difference evaluation has independently re-solved modes and
  # the reduced ordinary-score root while retaining the true beta constraint.
  expect_true(length(perturbed) >= 24L)
  expect_true(all(vapply(perturbed, function(x) {
    x$inner$score_metric <= 1e-8 && x$inner$step_metric <= 1e-10 &&
      x$constraint$scaled <= 100 * .Machine$double.eps
  }, logical(1))))
})

test_that("q zero has a separate certified full-chi Score null object", {
  case <- .hl_fo_score_c1()
  fixture <- case$fixture
  score <- .hl_fo_score_ns("hl_fo_score_test")(
    data = case$data,
    R = fixture$constraint$C,
    rhs = fixture$constraint$rhs,
    unrestricted_fit = case$seed,
    active = character(0),
    fixed_natural = c(
      sigma2 = fixture$truth$sigma2,
      phi = fixture$truth$phi
    )
  )
  expect_true(score$success)
  expect_identical(score$status, "ok")
  expect_identical(score$constrained_fit$object_type,
                   "reduced_constrained_maximizer")
  expect_identical(score$null_fit$object_type, "full_chi_score_null")
  expect_length(score$constrained_fit$rho, 0L)
  expect_length(score$null_fit$rho, 0L)
  expect_length(score$null_fit$g0, 0L)
  expect_equal(score$null_fit$root_scaled_residual, 0)
  expect_length(score$constrained_fit$augmented$chi, case$data$p - 1L)
  expect_length(score$null_fit$augmented$chi, case$data$p)
  expect_equal(score$statistic, score$stat)
  expect_equal(score$p_value, score$p)
  expect_lte(score$root_certification$constraint$scaled,
             100 * .Machine$double.eps)
  expect_lte(max(score$root_certification$statistical_displacement),
             score$root_certification$statistical_displacement_limit)
})

test_that("locked q-zero Score case uses specification-grade cap acceptance", {
  case <- .hl_fo_score_locked_mc_case(replicate = 75L, m = 250L)
  fixed <- c(sigma2 = case$sigma2, phi = case$phi)
  fit <- .hl_fo_score_ns("hl_fo_fit")(
    data = case$data,
    active = character(0),
    fixed_natural = fixed
  )
  expect_true(fit$point_success)
  expect_identical(fit$fit_status, "ok")

  affine <- .hl_fo_score_ns("hl_fo_affine_constraint")(
    case$R, case$rhs, p = case$data$p
  )
  solve_inner <- .hl_fo_score_ns("hl_fo_solve_inner")
  derivative <- solve_inner(
    data = case$data,
    rho_full = log(fixed),
    beta_star = affine$beta_star,
    K = affine$K,
    beta_start = fit$beta,
    cold = TRUE,
    max_iter = 3L,
    acceptance_tier = "derivative"
  )
  statistical <- solve_inner(
    data = case$data,
    rho_full = log(fixed),
    beta_star = affine$beta_star,
    K = affine$K,
    beta_start = fit$beta,
    cold = TRUE,
    max_iter = 3L,
    acceptance_tier = "statistical"
  )
  expect_false(derivative$success)
  expect_identical(derivative$reason, "inner_iteration_exhausted")
  expect_true(statistical$success)
  expect_identical(statistical$accepted_tier, "statistical")
  expect_equal(statistical$iterations, 3L)
  expect_lte(statistical$score_metric, 1e-8)
  expect_lte(statistical$step_metric, 1e-10)

  score <- .hl_fo_score_ns("hl_fo_score_test")(
    data = case$data,
    R = case$R,
    rhs = case$rhs,
    unrestricted_fit = fit,
    active = character(0),
    fixed_natural = fixed
  )
  expect_true(score$success)
  expect_identical(score$status, "ok")
  expect_identical(score$constrained_fit$inner$accepted_tier, "statistical")
  expect_identical(score$null_fit$inner$accepted_tier, "statistical")
  expect_lte(score$null_fit$inner$score_metric, 1e-8)
  expect_lte(score$null_fit$inner$step_metric, 1e-10)
  expect_lte(score$root_certification$constraint$scaled,
             100 * .Machine$double.eps)
  expect_true(all(vapply(score$null_fit$inner$mode_diagnostics, function(diag) {
    diag$residual_scaled <= 1e-10 && diag$step_scaled <= 1e-10
  }, logical(1))))
  expect_true(is.finite(score$stat) && score$stat >= 0)
})

test_that("root clustering retains the earlier root and rejects multiplicity", {
  select <- .hl_fo_score_ns("hl_fo_score_select_candidates")
  first <- list(start = 1L, state = list(rho = c(0, 0), beta = c(1, 2)))
  same <- list(start = 2L,
               state = list(rho = c(5e-8, -5e-8),
                            beta = c(1 + 5e-8, 2 - 5e-8)))
  other <- list(start = 3L, state = list(rho = c(0.1, 0), beta = c(1, 2)))

  none <- select(list())
  expect_false(none$success)
  expect_identical(none$status, "score_null_no_convergence")
  one <- select(list(first, same))
  expect_true(one$success)
  expect_identical(one$selected$start, 1L)
  multiple <- select(list(first, same, other))
  expect_false(multiple$success)
  expect_identical(multiple$status, "score_null_multiple_roots")
})

test_that("C1 uses a separately fitted full-chi null root and robust adjusted meat", {
  case <- .hl_fo_score_c1()
  fixture <- case$fixture
  constrained <- .hl_fo_score_ns("hl_fo_fit")(
    data = case$data,
    C = fixture$constraint$C,
    rhs = fixture$constraint$rhs,
    seed = case$seed
  )
  expect_true(constrained$point_success)
  expect_identical(constrained$fit_status, "ok")

  score <- .hl_fo_score_ns("hl_fo_score_test")(
    data = case$data,
    R = fixture$constraint$C,
    rhs = fixture$constraint$rhs,
    unrestricted_fit = case$seed,
    constrained_fit = constrained
  )
  expect_true(score$success)
  expect_identical(score$status, "ok")
  expect_identical(score$constrained_fit$object_type,
                   "reduced_constrained_maximizer")
  expect_identical(score$null_fit$object_type, "full_chi_score_null")
  expect_identical(score$null_fit$inner$acceptance_tier_requested, "derivative")
  expect_identical(score$null_fit$inner$accepted_tier, "derivative")
  expect_length(score$constrained_fit$augmented$chi, case$data$p - 1L)
  expect_length(score$null_fit$augmented$chi, case$data$p)
  expect_gt(
    .hl_fo_score_ns("hl_fo_cache_distance")(
      score$constrained_fit$rho,
      score$null_fit$rho
    ),
    1e-4
  )
  expect_lte(score$null_fit$root_scaled_residual, 1e-8)
  expect_lte(score$root_certification$constraint$scaled,
             100 * .Machine$double.eps)
  expect_lte(score$root_certification$auxiliary_backward_error, 1e-10)
  expect_lte(max(score$root_certification$statistical_displacement),
             score$root_certification$statistical_displacement_limit)
  null_state <- score$null_fit
  null_K <- score$root_certification$K
  expect_lte(
    max(abs(crossprod(null_K, null_state$inner$S_beta))),
    1e-8
  )
  expect_lte(
    max(abs(
      t(null_state$augmented$B) %*% null_state$augmented$chi -
        null_state$augmented$L_gamma
    )),
    1e-8
  )
  expect_lte(max(abs(null_state$g0 / null_state$root_scale)), 1e-8)
  expect_equal(score$df, 1L)
  expect_equal(score$contrast_rank, 1L)
  expect_true(is.finite(score$stat) && score$stat >= 0)
  expect_true(is.finite(score$p) && score$p >= 0 && score$p <= 1)
  expect_lte(score$meat_crosscheck, 1e-10)

  # An admissible rescaling of the same tested row leaves the post-root test
  # unchanged and does not require another nuisance fit.
  rescaled <- .hl_fo_score_ns("hl_fo_score_postroot")(
    psi_subject = score$null_fit$augmented$psi_subject,
    A_zeta = score$null_fit$augmented$A,
    R = -2 * fixture$constraint$C,
    beta = score$null_fit$beta,
    rhs = -2 * fixture$constraint$rhs,
    q = 2L
  )
  expect_true(rescaled$success)
  expect_equal(rescaled$stat, score$stat, tolerance = 1e-10)
  expect_equal(rescaled$p, score$p, tolerance = 1e-10)

  # A full orthogonal change of beta/chi coordinates changes both the tested
  # complement and null-space basis but cannot change the adjusted EE test.
  p <- case$data$p
  q <- length(score$null_fit$active)
  O <- qr.Q(qr(matrix(c(
    1, 2, -1,
    -2, 1, 1,
    1, 0.5, 2
  ), p, p)))
  T <- matrix(0, 2L * p + q, 2L * p + q)
  T[seq_len(p), seq_len(p)] <- O
  T[p + seq_len(p), p + seq_len(p)] <- O
  T[2L * p + seq_len(q), 2L * p + seq_len(q)] <- diag(q)
  rotated <- .hl_fo_score_ns("hl_fo_score_postroot")(
    psi_subject = crossprod(T, score$null_fit$augmented$psi_subject),
    A_zeta = crossprod(T, score$null_fit$augmented$A %*% T),
    R = fixture$constraint$C %*% O,
    beta = as.numeric(crossprod(O, score$null_fit$beta)),
    rhs = fixture$constraint$rhs,
    q = q
  )
  expect_true(rotated$success)
  expect_equal(rotated$stat, score$stat, tolerance = 1e-9)
  expect_equal(rotated$p, score$p, tolerance = 1e-9)
})
