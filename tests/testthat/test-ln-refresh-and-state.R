test_that("fit_gene_constrained evaluates derivatives at rescued constrained beta", {
  fit_fn <- getFromNamespace("fit_gene_constrained", "scRoPE")
  captured <- new.env(parent = emptyenv())
  captured$para <- NULL

  local_mocked_bindings(
    constraint_reparameterization = function(L, b) {
      list(
        beta_star = c(5, 0),
        null_basis = matrix(c(1, 0), ncol = 1)
      )
    },
    lbfgs = function(start, objective, lower, upper, control) {
      list(par = c(0.2, 2, 3))
    },
    ptmg_ll_der = function(...) {
      list(objective = 0, gradient = c(0, 0, 0, 0))
    },
    ptmg_ll = function(...) {
      1
    },
    opt_pml_constrained = function(...) {
      list(
        beta = c(9, 0),
        var = diag(2),
        loglik = 7,
        loglikp = 7,
        logw = 0
      )
    },
    check_conv = function(...) {
      -25L
    },
    ptmg_ll_der_hes4 = function(para, ...) {
      captured$para <- para
      list(
        gradient = c(1, 2, 3, 4),
        hessian = diag(4),
        per_subject_gradients = matrix(1, nrow = 4, ncol = 1)
      )
    },
    apply_chain_rule_transform = function(U, H, S, intercept_idx, p, U_per_subject, sigma2) {
      list(U = U, H = H, S = S, U_per_subject = U_per_subject)
    },
    extract_godambe_components = function(out_at_op, nb) {
      list(
        H = out_at_op$hessian,
        S = tcrossprod(out_at_op$per_subject_gradients),
        H_beta = diag(nb),
        H_beta_lambda = matrix(0, nrow = nb, ncol = 2),
        H_lambda = diag(2),
        per_subject_gradients = out_at_op$per_subject_gradients,
        score = out_at_op$gradient
      )
    },
    compute_sandwich_variance2 = function(...) {
      list(Var_beta_adjusted = diag(2))
    },
    .package = "scRoPE"
  )

  ctx <- list(
    allow_per_gene_switch = FALSE,
    opt = "lbfgs",
    intcol = 1L,
    nb = 2L,
    min = c(1e-4, 1e-4),
    max = c(10, 1000),
    cell_init = 1,
    eps = 1e-6,
    pred = matrix(0, nrow = 1, ncol = 2),
    offset = 0,
    id = 1L,
    fid = c(1L, 2L),
    cumsumy = matrix(0, nrow = 1, ncol = 1),
    posind = list(integer()),
    nind = 1L,
    k = 1L,
    output_re = FALSE
  )
  posv <- list(
    mct = 1,
    Y = numeric(),
    n_onetwo = c(0, 0),
    ytwo = numeric(),
    posindy = integer()
  )
  L <- matrix(c(0, 1), nrow = 1)

  out <- fit_fn(
    gene_index = 1L,
    posv = posv,
    ctx = ctx,
    L = L,
    b = 0
  )

  expect_equal(captured$para[1:ctx$nb], c(9, 0))
  expect_equal(out$theta$beta, c(9, 0))
  expect_equal(out$theta_ln$beta, c(5.2, 0))
  expect_equal(out$theta_unrescued$beta, out$theta_ln$beta)
})

test_that("fit_gene_unconstrained uses beta_ln when use_betas='final_betas'", {
  fit_fn <- getFromNamespace("fit_gene_unconstrained", "scRoPE")

  run_with_mode <- function(mode) {
    captured <- new.env(parent = emptyenv())
    captured$betas <- NULL

    local_mocked_bindings(
      lbfgs = function(start, objective, ...) {
        list(par = c(2.0, -0.3, 1.0, 0.1))
      },
      ptmg_ll_der = function(...) {
        list(objective = 0, gradient = c(0, 0, 0, 0))
      },
      ptmg_ll = function(...) {
        0
      },
      ptmg_der = function(...) {
        c(0, 0, 0, 0)
      },
      bobyqa = function(par, fn, ...) {
        list(par = c(1.5, 0.2))
      },
      opt_pml = function(X_c, offset_c, Y_c, fid_c, cumsumy_c, posind_c, posindy_c, nb_c, nind_c, k_c, beta_c, sigma_c, reml, eps, ord) {
        captured$betas <- beta_c
        list(
          beta = beta_c,
          var = diag(2),
          loglik = 1,
          loglikp = 1,
          logw = 0,
          iter = 1,
          damp = 0
        )
      },
      check_conv = function(...) {
        -25L
      },
      ptmg_ll_der_hes4 = function(...) {
        list(
          gradient = c(0, 0, 0, 0),
          hessian = diag(4),
          per_subject_gradients = matrix(0, nrow = 4, ncol = 1)
        )
      },
      extract_godambe_components = function(...) {
        NULL
      },
      compute_sandwich_variance2 = function(...) {
        list(Var_beta_adjusted = diag(2))
      },
      .package = "scRoPE"
    )

    ctx <- list(
      pred = matrix(c(1, 0), nrow = 1),
      offset = 0,
      fid = c(1L, 2L),
      cumsumy = matrix(0, nrow = 1, ncol = 1),
      posind = list(integer()),
      nb = 2L,
      nind = 1L,
      k = 1L,
      intcol = 1L,
      moffset = 0,
      min = c(1e-4, 1e-4),
      max = c(10, 1000),
      allow_per_gene_switch = FALSE,
      use_betas = mode,
      cutoff_cell = 20,
      cell_ind = integer(),
      ncell = 0L,
      cv2 = 0,
      kappa = 800,
      eps = 1e-6,
      opt = "lbfgs",
      cell_init = 1,
      output_re = FALSE,
      id = 1L
    )
    posv <- list(
      mct = 2,
      Y = numeric(),
      n_onetwo = c(0, 0),
      ytwo = numeric(),
      posindy = integer()
    )

    out <- fit_fn(
      gene_index = 1L,
      posv = posv,
      ctx = ctx
    )
    list(betas = captured$betas, out = out)
  }

  re_betae <- run_with_mode("betae")
  re_final <- run_with_mode("final_betas")

  expect_false(isTRUE(all.equal(re_betae$betas, re_final$betas)))
  expect_equal(re_betae$out$theta_ln$beta, re_betae$out$theta_unrescued$beta)
  expect_equal(re_final$out$theta_ln$beta, re_final$out$theta_unrescued$beta)
})
