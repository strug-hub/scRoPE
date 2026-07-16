ln_fd_gradient <- function(fn, theta, rel_step = 2e-5) {
  step <- rel_step * pmax(1, abs(theta))
  vapply(seq_along(theta), function(jj) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    (fn(plus) - fn(minus)) / (2 * step[jj])
  }, numeric(1))
}

ln_fd_hessian <- function(fn, theta, rel_step = 2e-4) {
  step <- rel_step * pmax(1, abs(theta))
  p <- length(theta)
  out <- matrix(0, p, p)
  f0 <- fn(theta)
  for (jj in seq_len(p)) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    out[jj, jj] <- (fn(plus) - 2 * f0 + fn(minus)) / step[jj]^2
    if (jj < p) {
      for (kk in (jj + 1L):p) {
        pp <- pm <- mp <- mm <- theta
        pp[c(jj, kk)] <- pp[c(jj, kk)] + step[c(jj, kk)]
        pm[jj] <- pm[jj] + step[jj]
        pm[kk] <- pm[kk] - step[kk]
        mp[jj] <- mp[jj] - step[jj]
        mp[kk] <- mp[kk] + step[kk]
        mm[c(jj, kk)] <- mm[c(jj, kk)] - step[c(jj, kk)]
        out[jj, kk] <- out[kk, jj] <-
          (fn(pp) - fn(pm) - fn(mp) + fn(mm)) / (4 * step[jj] * step[kk])
      }
    }
  }
  out
}

ln_small_fixture <- function() {
  set.seed(20260715)
  n_subjects <- 12L
  cells_per_subject <- 20L
  n <- n_subjects * cells_per_subject
  id <- rep(seq_len(n_subjects), each = cells_per_subject)
  x <- rnorm(n)
  z_subject <- rbinom(n_subjects, 1, 0.5)
  z <- rep(z_subject, each = cells_per_subject)
  sigma2 <- 0.4
  phi <- 0.3
  omega <- rgamma(
    n_subjects,
    shape = 1 / (exp(sigma2) - 1),
    rate = 1 / (exp(sigma2 / 2) * (exp(sigma2) - 1))
  )
  cell_effect <- rgamma(n, shape = 1 / phi, rate = 1 / phi)
  mu <- exp(-0.5 + 0.15 * x - 0.1 * z) * rep(omega, each = cells_per_subject) * cell_effect
  y <- rpois(n, mu)
  count <- Matrix::Matrix(matrix(y, nrow = 1L), sparse = TRUE)
  rownames(count) <- "gene1"
  pred <- cbind("(Intercept)" = 1, x = x, z = z)
  list(count = count, id = id, pred = pred)
}

test_that("LN reduced map matches finite differences of the raw profile objective", {
  fx <- ln_small_fixture()
  fit <- suppressWarnings(
    scrope(
      fx$count,
      fx$id,
      pred = fx$pred,
      offset = rep(1, ncol(fx$count)),
      ncore = 1,
      verbose = FALSE,
      cpc = 0,
      mincp = 1,
      additional_tests = c("score", "lrt"),
      lrt_details = TRUE,
      keep_diagnostics = TRUE
    )
  )

  raw_gene <- fit$diagnostics$raw[[1L]]
  fit_c <- raw_gene$lrt$x$constrained
  lrt_x <- raw_gene$lrt$x$raw
  expect_equal(fit_c$reduced$parameterization, "gamma, log(sigma^2), log(phi)")
  expect_equal(fit_c$reduced$theta_red, fit_c$reduced$theta)
  expect_equal(fit_c$reduced$g_red, fit_c$reduced$score)
  expect_equal(fit_c$reduced$H_red, fit_c$reduced$observed_info)
  expect_true(is.finite(lrt_x$wU))
  expect_true(is.finite(lrt_x$qP))
  expect_equal(lrt_x$wP_minus_qP, lrt_x$wP - lrt_x$qP, tolerance = 1e-12)
  expect_equal(lrt_x$wP_over_qP, lrt_x$wP / lrt_x$qP, tolerance = 1e-12)

  count <- as(fx$count, "dgCMatrix")
  nind <- ncol(count)
  pred_centered <- center_m(as.matrix(fx$pred))$pred
  offset_state <- cv_offset(rep(1, nind), 1, nind)
  id_numeric <- as.numeric(factor(as.character(fx$id), levels = unique(as.character(fx$id))))
  k <- length(unique(id_numeric))
  fid <- which(c(1, diff(id_numeric)) == 1)
  fid <- as.integer(c(fid, nind + 1L))
  cumsumy <- matrix(call_cumsumy(count, fid - 1L, k, 1L), nrow = 1L)
  posind <- which(cumsumy[1L, ] > 0)
  count_t <- Matrix::t(count)
  posv <- call_posindy(count_t, 0L, nind)
  K <- fit_c$constraint$null_basis
  beta_star <- fit_c$constraint$beta_star
  free_dim <- ncol(K)

  reduced_objective <- function(theta_red) {
    gamma <- if (free_dim) theta_red[seq_len(free_dim)] else numeric(0)
    beta <- if (free_dim) as.vector(beta_star + K %*% gamma) else beta_star
    sigma2 <- exp(theta_red[free_dim + 1L])
    phi <- exp(theta_red[free_dim + 2L])
    ptmg_ll(
      c(beta, sigma2, phi),
      X = pred_centered,
      offset = offset_state$offset,
      Y = posv$Y,
      n_one = posv$n_onetwo,
      ytwo = posv$ytwo,
      fam = id_numeric,
      fid = fid,
      cumsumy = cumsumy[1L, ],
      posind = posind,
      posindy = posv$posindy,
      nb = ncol(pred_centered),
      nind = nind,
      k = k
    )
  }

  theta_red <- fit_c$reduced$theta_red
  expect_equal(
    fit_c$reduced$g_red,
    ln_fd_gradient(reduced_objective, theta_red),
    tolerance = 2e-4
  )
  expect_equal(
    fit_c$reduced$H_red,
    ln_fd_hessian(reduced_objective, theta_red),
    tolerance = 3e-3
  )
})

test_that("LN reduced Newton diagnostics require finite positive-definite information", {
  diagnostic_fn <- getFromNamespace("reduced_newton_diagnostics", "scRoPE")
  good <- diagnostic_fn(diag(c(2, 3)), c(1e-5, 2e-5), c(0, 0))
  expect_true(good$H_red_positive_definite)
  expect_true(is.finite(good$max_abs_newton_step_scaled))
  expect_true(is.finite(good$gradient_quadratic_lr_bound))

  indefinite <- diagnostic_fn(diag(c(1, -1)), c(1e-5, 2e-5), c(0, 0))
  expect_false(indefinite$H_red_positive_definite)
  expect_true(is.na(indefinite$max_abs_newton_step_scaled))

  nonfinite <- diagnostic_fn(diag(2), c(NA_real_, 0), c(0, 0))
  expect_false(nonfinite$H_red_positive_definite)
  expect_true(is.na(nonfinite$gradient_quadratic_lr_bound))
})

ln_mock_constrained_fit <- function(conv_code = -20L,
                                    legacy_converged = FALSE,
                                    positive_definite = TRUE,
                                    newton = 5e-6,
                                    quadratic = 2e-6) {
  list(
    loglik_ln = -10,
    theta_ln = list(beta = c(1, 0), subVar = 0.5, cellVar = 1),
    lrt_gate = list(
      converged = legacy_converged,
      interior = TRUE,
      nuisance_interior = TRUE
    ),
    debug = list(conv_code = conv_code),
    diagnostics = list(
      reduced_gradient_max_abs = 1,
      H_red_positive_definite = positive_definite,
      H_red_rcond = if (positive_definite) 0.5 else 0,
      max_abs_newton_step_scaled = newton,
      gradient_quadratic_lr_bound = quadratic
    ),
    reduced = list(parameterization = "gamma, log(sigma^2), log(phi)")
  )
}

test_that("LN information-scaled gate uses the documented OR rule", {
  assess_fn <- getFromNamespace("ln_assess_lrt_convergence", "scRoPE")

  newton_only <- assess_fn(ln_mock_constrained_fit(newton = 5e-6, quadratic = 2e-6))
  expect_true(newton_only$accepted)
  expect_equal(newton_only$convergence_basis, "information_scaled")
  expect_true(newton_only$newton_step_scaled_ok)
  expect_false(newton_only$gradient_quadratic_lr_bound_ok)

  quadratic_only <- assess_fn(ln_mock_constrained_fit(newton = 2e-5, quadratic = 5e-7))
  expect_true(quadratic_only$accepted)
  expect_false(quadratic_only$newton_step_scaled_ok)
  expect_true(quadratic_only$gradient_quadratic_lr_bound_ok)

  disabled <- assess_fn(
    ln_mock_constrained_fit(newton = 5e-6, quadratic = 5e-7),
    use_information_scaled_convergence = FALSE
  )
  expect_false(disabled$accepted)
  expect_false(disabled$information_scaled_convergence_ok)
})

test_that("LN information-scaled gate cannot override hard legacy failures", {
  assess_fn <- getFromNamespace("ln_assess_lrt_convergence", "scRoPE")

  curvature <- assess_fn(ln_mock_constrained_fit(conv_code = -25L, quadratic = 5e-7))
  expect_false(curvature$information_scaled_eligible)
  expect_false(curvature$accepted)

  boundary <- assess_fn(ln_mock_constrained_fit(conv_code = -60L, quadratic = 5e-7))
  expect_false(boundary$information_scaled_eligible)
  expect_false(boundary$accepted)

  indefinite <- assess_fn(
    ln_mock_constrained_fit(positive_definite = FALSE, newton = 5e-6, quadratic = 5e-7)
  )
  expect_false(indefinite$information_scaled_convergence_ok)
  expect_false(indefinite$accepted)

  nonfinite <- assess_fn(ln_mock_constrained_fit(newton = NA_real_, quadratic = NA_real_))
  expect_false(nonfinite$diagnostics_finite)
  expect_false(nonfinite$information_scaled_convergence_ok)
  expect_false(nonfinite$accepted)

  nonfinite_rcond_fit <- ln_mock_constrained_fit(newton = 5e-6, quadratic = 5e-7)
  nonfinite_rcond_fit$diagnostics$H_red_rcond <- NA_real_
  nonfinite_rcond <- assess_fn(nonfinite_rcond_fit)
  expect_false(nonfinite_rcond$diagnostics_finite)
  expect_false(nonfinite_rcond$accepted)
})

test_that("LN gate audits the raw-profile constraint without mutating fit convergence", {
  gate_fn <- getFromNamespace("ln_lrt_gate_override", "scRoPE")
  fit <- ln_mock_constrained_fit(newton = 5e-6, quadratic = 2e-6)
  L <- matrix(c(0, 1), nrow = 1)

  gate <- gate_fn(fit, L, 0)
  expect_true(gate$converged)
  expect_true(gate$constraint_satisfied)
  expect_equal(gate$convergence_basis, "information_scaled")
  expect_false(fit$lrt_gate$converged)

  failed_constraint <- gate_fn(fit, L, 1)
  expect_false(failed_constraint$constraint_satisfied)
  expect_gt(failed_constraint$constraint_residual_max_abs, 0)
})
