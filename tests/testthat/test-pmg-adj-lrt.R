test_that("PMG adjusted LRT is null at the exact unconstrained estimate", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(1, 0), nrow = 1)
  res <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = list(list(L = L, b = fit_u$theta$beta[1], label = "intercept")),
    unconstrained = fit_u
  )[[1]]

  expect_equal(res$adj$lr$wP, 0, tolerance = 1e-8)
  expect_equal(res$adj$lr$wP_adjusted, 0, tolerance = 1e-6)
  expect_equal(as.numeric(res$adj$lr$score), 0, tolerance = 1e-6)
  expect_true(isTRUE(res$adj$lr$success))
  expect_true(isTRUE(res$adj$lr$adjusted_available))
  expect_false(isTRUE(res$adj$lr$fallback))
  expect_equal(dim(res$adj$lr$Hpsi), c(1, 1))
  expect_equal(dim(res$adj$lr$Gpsi), c(1, 1))
  expect_equal(dim(res$adj$lr$Jpsi), c(1, 1))
  expect_true(as.numeric(res$adj$lr$Hpsi) > 0)
  expect_true(as.numeric(res$adj$lr$Gpsi) > 0)
})

test_that("PMG adjusted LRT is available off-null when the scale gate passes", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")
  direct_lrt_fn <- getFromNamespace("compute_profile_lrt", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(0, 1), nrow = 1)
  b <- 0

  res <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = list(list(L = L, b = b, label = "slope")),
    unconstrained = fit_u
  )[[1]]

  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = b,
    unconstrained_fit = fit_u
  )
  direct_lrt <- direct_lrt_fn(fit_u, fit_c, L, b)

  expect_gt(res$adj$lr$wP, 0)
  expect_gt(res$adj$lr$wP_adjusted, 0)
  expect_true(isTRUE(res$adj$lr$adjusted_available))
  expect_false(isTRUE(res$adj$lr$fallback))
  expect_true(is.finite(res$p_values$lr))
  expect_equal(res$adj$lr$wP, direct_lrt$wP, tolerance = 1e-8)
  expect_equal(res$adj$lr$wP_adjusted, direct_lrt$wP_adjusted, tolerance = 1e-8)
  expect_equal(res$adj$lr$Gpsi, direct_lrt$Gpsi, tolerance = 1e-8)
  expect_equal(res$adj$lr$Hpsi, direct_lrt$Hpsi, tolerance = 1e-8)
  expect_equal(res$adj$lr$Jpsi, direct_lrt$Jpsi, tolerance = 1e-8)
  expect_equal(res$adj$lr$Vpsi_model, direct_lrt$Vpsi_model, tolerance = 1e-8)
  expect_equal(res$adj$lr$Vpsi_godambe, direct_lrt$Vpsi_godambe, tolerance = 1e-8)
})

test_that("PMG adjusted LRT gates out nonpositive curvature and falls back to Wald", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")
  direct_lrt_fn <- getFromNamespace("compute_profile_lrt", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(1, 0), nrow = 1)
  b <- 0

  res <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = list(list(L = L, b = b, label = "intercept")),
    unconstrained = fit_u
  )[[1]]

  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = b,
    unconstrained_fit = fit_u
  )
  direct_lrt <- direct_lrt_fn(fit_u, fit_c, L, b)

  expect_false(isTRUE(res$adj$lr$adjusted_available))
  expect_true(isTRUE(res$adj$lr$fallback))
  expect_match(res$adj$lr$method, "^wald")
  expect_equal(res$adj$lr$adjusted_failure_reason, "nonpositive_curvature")
  expect_false(res$adj$lr$gate$profiled_curvature_positive)
  expect_false(res$adj$lr$gate$scale_positive)
  expect_gt(res$adj$lr$wP_adjusted, 0)
  expect_equal(res$adj$lr$wP_adjusted, direct_lrt$wP_adjusted, tolerance = 1e-8)
})

test_that("PMG adjusted LRT accepts an explicit information-scaled convergence gate", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  assess_fn <- getFromNamespace("pmg_assess_constrained_convergence", "scRoPE")
  lrt_fn <- getFromNamespace("compute_profile_lrt", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )
  L <- matrix(c(0, 1), nrow = 1)
  b <- 0
  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = b,
    unconstrained_fit = fit_u
  )
  fit_c$convergence <- -50L
  fit_c$lrt_gate$converged <- FALSE
  fit_c$diagnostics$reduced_gradient_max_abs <- 1
  fit_c$diagnostics$H_red_positive_definite <- TRUE
  fit_c$diagnostics$max_abs_newton_step_scaled <- 5e-6
  fit_c$diagnostics$gradient_quadratic_lr_bound <- 5e-7

  assessment <- assess_fn(
    diagnostics = fit_c$diagnostics,
    optimizer_converged = FALSE,
    n_subjects = fx$prep$ctx$k,
    gradient_tol = 1e-6,
    average_gradient_tol = 1e-6
  )
  override <- list(
    converged = assessment$accepted,
    convergence_basis = assessment$convergence_basis,
    optimizer_converged = assessment$optimizer_converged,
    information_scaled_convergence_ok = assessment$information_scaled_convergence_ok,
    H_red_positive_definite = assessment$H_red_positive_definite,
    max_abs_newton_step_scaled = assessment$max_abs_newton_step_scaled,
    gradient_quadratic_lr_bound = assessment$gradient_quadratic_lr_bound,
    constraint_satisfied = TRUE
  )

  rejected <- lrt_fn(fit_u, fit_c, L, b)
  accepted <- lrt_fn(
    fit_u,
    fit_c,
    L,
    b,
    lrt_gate_override = override
  )

  expect_false(rejected$success)
  expect_equal(rejected$adjusted_failure_reason, "nonconverged_constrained")
  expect_true(accepted$success)
  expect_equal(accepted$gate$convergence_basis, "information_scaled")
  expect_true(accepted$gate$information_scaled_convergence_ok)
})

test_that("PMG model-specific LRT gate retains hard validity checks", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  lrt_fn <- getFromNamespace("compute_profile_lrt", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )
  L <- matrix(c(0, 1), nrow = 1)
  b <- 0
  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = b,
    unconstrained_fit = fit_u
  )
  valid_override <- list(
    converged = TRUE,
    convergence_basis = "optimizer",
    optimizer_converged = TRUE,
    constraint_satisfied = TRUE
  )

  constraint_failure <- lrt_fn(
    fit_u,
    fit_c,
    L,
    b,
    lrt_gate_override = modifyList(valid_override, list(constraint_satisfied = FALSE))
  )
  expect_false(constraint_failure$success)
  expect_equal(constraint_failure$adjusted_failure_reason, "constraint_not_satisfied")

  boundary_fit <- fit_c
  boundary_fit$lrt_gate$interior <- FALSE
  boundary_failure <- lrt_fn(
    fit_u,
    boundary_fit,
    L,
    b,
    lrt_gate_override = valid_override
  )
  expect_false(boundary_failure$success)
  expect_equal(boundary_failure$adjusted_failure_reason, "boundary_constrained")

  negative_unconstrained <- fit_u
  negative_unconstrained$loglik <- fit_c$loglik - 1
  negative_failure <- lrt_fn(
    negative_unconstrained,
    fit_c,
    L,
    b,
    lrt_gate_override = valid_override
  )
  expect_false(negative_failure$success)
  expect_equal(negative_failure$failure_reason, "negative_raw_lr")

  zero_variability_fit <- fit_c
  zero_variability_fit$per_subject_gradients[,] <- 0
  variability_failure <- lrt_fn(
    fit_u,
    zero_variability_fit,
    L,
    b,
    lrt_gate_override = valid_override
  )
  expect_false(variability_failure$success)
  expect_false(variability_failure$gate$profiled_variability_positive)
})

test_that("PMG retry wiring preserves successful adjusted-test statistics", {
  fx <- pmg_interior_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  adj_contrasts_fn <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )
  contrasts <- list(list(
    L = matrix(c(0, 1), nrow = 1),
    b = 0,
    label = "slope"
  ))
  without_retry <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = contrasts,
    unconstrained = fit_u,
    retry_constrained = FALSE
  )[[1]]
  with_retry <- adj_contrasts_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    contrasts = contrasts,
    unconstrained = fit_u,
    retry_constrained = TRUE
  )[[1]]

  expect_equal(with_retry$adj$score$stat, without_retry$adj$score$stat, tolerance = 1e-8)
  expect_equal(with_retry$adj$lr$wP, without_retry$adj$lr$wP, tolerance = 1e-8)
  expect_equal(with_retry$adj$lr$wP_adjusted, without_retry$adj$lr$wP_adjusted, tolerance = 1e-8)
  expect_true(with_retry$adj$lr$gate$convergence_basis %in% c(
    "optimizer", "raw_gradient", "information_scaled"
  ))
  expect_true(all(c(
    "retry_triggered",
    "retry_used",
    "optimizer_attempts",
    "constraint_residual_max_abs",
    "unconstrained_repair_attempted"
  ) %in% names(with_retry$adj$lr$gate)))
})
