test_that("PMG per-subject gradients aggregate to pmg_der and match pmg_hes", {
  fx <- pmg_test_fixture()
  ctx <- fx$prep$ctx
  i <- fx$gene_index
  posv <- fx$posv

  eval_fn <- getFromNamespace("pmg_ll_der_hes", "scRoPE")
  der_fn <- getFromNamespace("pmg_der", "scRoPE")
  hes_fn <- getFromNamespace("pmg_hes", "scRoPE")

  theta <- c(0.15, -0.25, 0.8)

  eval_at_theta <- eval_fn(
    para = theta,
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[i, ],
    posind = ctx$posind[[i]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k
  )

  grad_global <- der_fn(
    para = theta,
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[i, ],
    posind = ctx$posind[[i]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k
  )
  hes_global <- hes_fn(
    para = theta,
    X = ctx$pred,
    offset = ctx$offset,
    Y = posv$Y,
    fid = ctx$fid,
    cumsumy = ctx$cumsumy[i, ],
    posind = ctx$posind[[i]],
    posindy = posv$posindy,
    nb = ctx$nb,
    nind = ctx$nind,
    k = ctx$k
  )

  expect_equal(
    unname(rowSums(eval_at_theta$per_subject_gradients)),
    grad_global,
    tolerance = 1e-10
  )
  expect_equal(eval_at_theta$hessian, hes_global, tolerance = 1e-10)
})

test_that("fit_gene_pmg_unconstrained agrees with legacy nebula PMM fit", {
  fx <- pmg_test_fixture()
  fit_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")

  fit <- fit_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  legacy <- nebula(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    model = "PMM",
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1
  )

  beta_cols <- grep("^logFC_", names(legacy$summary), value = TRUE)
  se_cols <- grep("^se_", names(legacy$summary), value = TRUE)
  sds <- fx$prep$design_sds
  sds[fx$prep$ctx$intcol] <- 1

  expect_equal(
    as.numeric(fit$theta$beta),
    as.numeric(legacy$summary[1, beta_cols]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(as.numeric(fit$theta$subVar)),
    unname(as.numeric(legacy$overdispersion[1])),
    tolerance = 1e-6
  )
  expect_equal(
    diag(fit$naive_cov) / (sds^2),
    as.numeric(legacy$summary[1, se_cols])^2,
    tolerance = 1e-6
  )
})

test_that("generalized robust covariance works for PGMM with q = 1", {
  fx <- pmg_test_fixture()
  fit_fn <- getFromNamespace("fit_gene_pmg_unconstrained", "scRoPE")

  fit <- fit_fn(
    gene_index = fx$gene_index,
    posv = fx$posv,
    ctx = fx$prep$ctx
  )

  H <- fit$godambe$H
  S <- tcrossprod(fit$per_subject_gradients)
  G <- H %*% solve(S) %*% H

  nb <- fx$prep$ctx$nb
  G_bb <- G[seq_len(nb), seq_len(nb), drop = FALSE]
  G_bl <- G[seq_len(nb), nb + 1L, drop = FALSE]
  G_lb <- G[nb + 1L, seq_len(nb), drop = FALSE]
  G_ll <- matrix(G[nb + 1L, nb + 1L], nrow = 1)
  G_beta <- G_bb - G_bl %*% solve(G_ll) %*% G_lb
  robust_expected <- solve(G_beta)

  expect_equal(dim(fit$godambe$H_beta_lambda), c(nb, 1))
  expect_true(all(is.finite(fit$robust_cov)))
  expect_equal(fit$robust_cov, robust_expected, tolerance = 1e-8)
})
