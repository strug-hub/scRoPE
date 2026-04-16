test_that("PGMM profile helpers return the expected pointwise structure", {
  fx <- pmg_interior_test_fixture()
  fit_uncached_fn <- getFromNamespace("pgmm_fit_unconstrained", "scRoPE")
  point_fn <- getFromNamespace("pgmm_profile_point", "scRoPE")

  fit_u <- fit_uncached_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    psi_index = 2
  )

  expect_true(is.list(fit_u))
  expect_true(is.list(fit_u$fit))
  expect_equal(fit_u$psi_hat, fit_u$fit$theta$beta[2], tolerance = 1e-10)
  expect_equal(fit_u$target$psi_index, 2L)
  expect_equal(as.numeric(fit_u$target$c), c(0, 1))

  point <- point_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    psi_value = 0,
    unconstrained = fit_u
  )

  expect_equal(point$psi, 0)
  expect_equal(point$psi_hat, fit_u$psi_hat, tolerance = 1e-10)
  expect_true(is.finite(point$loglik_profile))
  expect_true(point$diagnostics$constraint_satisfied)
  expect_true(is.list(point$score_blocks))
  expect_true(is.list(point$per_subject_scores))
  expect_true(is.list(point$observed_info_blocks))
  expect_equal(length(point$score_blocks$eta), point$target$free_dim + 1L)
  expect_equal(length(point$per_subject_scores$psi), fx$prep$ctx$k)
  expect_equal(
    dim(point$observed_info_blocks$eta_eta),
    c(point$target$free_dim + 1L, point$target$free_dim + 1L)
  )
  expect_true(all(is.finite(point$profile$subject_scores)))
  expect_true(is.finite(point$profile$score))
  expect_true(is.finite(point$profile$variability))
  expect_true(is.finite(point$profile$sensitivity))
})

test_that("PGMM profile grid returns tidy structure and finite relative likelihoods", {
  fx <- pmg_interior_test_fixture()
  fit_uncached_fn <- getFromNamespace("pgmm_fit_unconstrained", "scRoPE")
  grid_fn <- getFromNamespace("pgmm_relative_profile_grid", "scRoPE")

  fit_u <- fit_uncached_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    psi_index = 2
  )

  grid_vals <- fit_u$psi_hat + c(-0.3, -0.1, 0.1, 0.3)
  grid <- grid_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    psi_values = grid_vals,
    unconstrained = fit_u
  )

  expect_s3_class(grid, "data.frame")
  expect_equal(nrow(grid), length(grid_vals))
  expect_equal(grid$psi, grid_vals, tolerance = 1e-12)
  expect_true(!is.null(attr(grid, "target")))
  expect_true(!is.null(attr(grid, "unconstrained")))

  expected_cols <- c(
    "psi",
    "psi_hat",
    "loglik_profile",
    "wP",
    "relLik",
    "UP",
    "Jpsi",
    "Hpsi",
    "Gpsi",
    "wU",
    "qP",
    "wP_star",
    "relLik_star",
    "relLik_std",
    "relLik_star_std",
    "converged",
    "optimizer",
    "robust_available"
  )
  expect_true(all(expected_cols %in% colnames(grid)))

  expect_true(all(is.finite(grid$relLik[grid$converged])))
  expect_true(all(is.finite(grid$relLik_star[grid$robust_available])))
  expect_true(all(is.finite(grid$wP[grid$converged])))
  expect_true(all(is.finite(grid$wP_star[grid$robust_available])))
  expect_equal(max(grid$relLik_std, na.rm = TRUE), 1, tolerance = 1e-10)
  expect_equal(max(grid$relLik_star_std, na.rm = TRUE), 1, tolerance = 1e-10)
})

test_that("PGMM profile robust adjustment agrees with the current scalar robust LRT", {
  fx <- pmg_interior_test_fixture()
  fit_uncached_fn <- getFromNamespace("pgmm_fit_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")
  point_fn <- getFromNamespace("pgmm_relative_profile_point", "scRoPE")
  direct_lrt_fn <- getFromNamespace("compute_profile_lrt", "scRoPE")

  fit_u <- fit_uncached_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    psi_index = 2
  )

  psi0 <- 0
  rel_point <- point_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    psi_value = psi0,
    unconstrained = fit_u
  )

  L <- matrix(c(0, 1), nrow = 1)
  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = psi0,
    unconstrained_fit = fit_u$fit
  )
  direct_lrt <- direct_lrt_fn(
    unconstrained = fit_u$fit,
    constrained = fit_c,
    L = L,
    b = psi0
  )

  expect_true(rel_point$relative_profile$robust_available)
  expect_equal(rel_point$relative_profile$wP, direct_lrt$wP, tolerance = 1e-8)
  expect_equal(rel_point$relative_profile$wP_star, direct_lrt$wP_adjusted, tolerance = 1e-8)
  expect_equal(rel_point$relative_profile$relLik_star, exp(-0.5 * direct_lrt$wP_adjusted), tolerance = 1e-8)
})

test_that("PGMM coefficient-grid wrapper handles a grid containing psi_hat", {
  fx <- pmg_interior_test_fixture()
  fit_uncached_fn <- getFromNamespace("pgmm_fit_unconstrained", "scRoPE")
  coef_grid_fn <- getFromNamespace("pgmm_relative_profile_coef_grid", "scRoPE")

  fit_u <- fit_uncached_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    psi_index = 2
  )

  grid_vals <- fit_u$psi_hat + c(-0.2, 0, 0.2)
  grid <- coef_grid_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    psi_index = 2,
    psi_values = grid_vals
  )

  idx_hat <- which.min(abs(grid$psi - grid$psi_hat))
  expect_equal(grid$psi[idx_hat], grid$psi_hat[idx_hat], tolerance = 1e-12)
  expect_equal(grid$wP[idx_hat], 0, tolerance = 1e-6)
  expect_equal(grid$relLik[idx_hat], 1, tolerance = 1e-8)
  expect_equal(grid$relLik_std[idx_hat], 1, tolerance = 1e-8)
  expect_equal(attr(grid, "target")$psi_index, 2L)
})

test_that("PGMM profile grid supports a custom scalar contrast c", {
  fx <- pmg_interior_test_fixture()
  fit_uncached_fn <- getFromNamespace("pgmm_fit_unconstrained", "scRoPE")
  grid_fn <- getFromNamespace("pgmm_relative_profile_grid", "scRoPE")

  c_vec <- c(1, 1)
  fit_u <- fit_uncached_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    c = c_vec
  )

  expected_hat <- sum(c_vec * fit_u$fit$theta$beta)
  grid_vals <- fit_u$psi_hat + c(-0.25, 0, 0.25)
  grid <- grid_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    psi_values = grid_vals,
    unconstrained = fit_u
  )

  idx_hat <- which.min(abs(grid$psi - grid$psi_hat))
  expect_true(is.na(attr(grid, "target")$psi_index))
  expect_equal(as.numeric(attr(grid, "target")$c), c_vec, tolerance = 1e-12)
  expect_equal(unique(grid$psi_hat), expected_hat, tolerance = 1e-10)
  expect_equal(grid$wP[idx_hat], 0, tolerance = 1e-6)
  expect_equal(grid$relLik[idx_hat], 1, tolerance = 1e-8)
  expect_true(all(is.finite(grid$relLik[grid$converged])))
  expect_true(all(is.finite(grid$relLik_star[grid$robust_available])))
})

test_that("PGMM-truth simulation gives finite profile curves and peaks near the true scaled effect", {
  fx <- pmg_truth_sim_fixture()
  fit_uncached_fn <- getFromNamespace("pgmm_fit_unconstrained", "scRoPE")
  coef_grid_fn <- getFromNamespace("pgmm_relative_profile_coef_grid", "scRoPE")

  fit_u <- fit_uncached_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    psi_index = 2
  )

  true_scaled <- fx$truth$beta_scaled[["x"]]
  grid_vals <- seq(true_scaled - 0.5, true_scaled + 0.5, length.out = 9)
  grid <- coef_grid_fn(
    ctx = fx$prep$ctx,
    gene_index = fx$gene_index,
    posv = fx$posv,
    psi_index = 2,
    psi_values = grid_vals
  )

  ord_peak <- grid$psi[[which.max(grid$relLik_std)]]
  rob_peak <- grid$psi[[which.max(grid$relLik_star_std)]]
  grid_step <- diff(grid_vals)[1]

  expect_true(all(is.finite(grid$relLik[grid$converged])))
  expect_true(all(is.finite(grid$relLik_star[grid$robust_available])))
  expect_equal(max(grid$relLik_std, na.rm = TRUE), 1, tolerance = 1e-10)
  expect_equal(max(grid$relLik_star_std, na.rm = TRUE), 1, tolerance = 1e-10)
  expect_lte(abs(fit_u$psi_hat - true_scaled), 0.15)
  expect_lte(abs(ord_peak - true_scaled), grid_step)
  expect_lte(abs(rob_peak - true_scaled), grid_step)
})
