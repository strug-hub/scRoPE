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
