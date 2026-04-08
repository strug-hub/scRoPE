test_that("constrained PGMM fit recovers the unconstrained optimum at beta_hat", {
  fx <- pmg_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(0, 1), nrow = 1)
  b <- fit_u$theta$beta[2]

  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = b,
    unconstrained_fit = fit_u
  )

  expect_equal(as.numeric(L %*% fit_c$theta$beta), as.numeric(b), tolerance = 1e-10)
  expect_equal(fit_c$theta$beta, fit_u$theta$beta, tolerance = 1e-6)
  expect_equal(fit_c$theta$subVar, fit_u$theta$subVar, tolerance = 1e-6)
  expect_equal(fit_c$loglik, fit_u$loglik, tolerance = 1e-6)
  expect_equal(
    unname(rowSums(fit_c$per_subject_gradients)),
    fit_c$score,
    tolerance = 1e-6
  )
  expect_equal(
    fit_c$reduced$score[seq_len(fit_c$constraint$free_dim)],
    rep(0, fit_c$constraint$free_dim),
    tolerance = 1e-5
  )
})

test_that("constrained PGMM fit enforces the null and returns constrained beta covariance", {
  fx <- pmg_test_fixture()
  fit_unconstrained_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")
  fit_constrained_fn <- getFromNamespace("fit_gene_pmg_constrained", "scRoPE")

  fit_u <- fit_unconstrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  L <- matrix(c(0, 1), nrow = 1)
  fit_c <- fit_constrained_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx,
    L = L,
    b = 0
  )

  expect_equal(as.numeric(L %*% fit_c$theta$beta), 0, tolerance = 1e-10)
  expect_lte(fit_c$loglik, fit_u$loglik + 1e-8)
  expect_equal(dim(fit_c$godambe$H_beta_lambda), c(fx$prep$ctx$nb, 1))
  expect_true(is.finite(fit_c$naive_cov[1, 1]))
  expect_true(is.finite(fit_c$robust_cov[1, 1]))
  expect_equal(fit_c$naive_cov[2, 2], 0, tolerance = 1e-12)
  expect_equal(fit_c$robust_cov[2, 2], 0, tolerance = 1e-12)
  expect_equal(fit_c$naive_cov[1, 2], 0, tolerance = 1e-12)
  expect_equal(fit_c$robust_cov[1, 2], 0, tolerance = 1e-12)
})
