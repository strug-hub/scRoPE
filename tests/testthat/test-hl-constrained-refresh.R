test_that("hl_refresh_constrained_state maps log(phi) and tau indices correctly", {
  refresh_fn <- getFromNamespace("hl_refresh_constrained_state", "scRoPE")

  nb <- 2L
  phi <- 5

  mock_gradient <- c(10, 20, 30, 40) # beta1, beta2, tau, log(phi)
  mock_hessian <- matrix(
    c(
      1, 2, 3, 4,
      2, 5, 6, 7,
      3, 6, 8, 9,
      4, 7, 9, 10
    ),
    nrow = nb + 2L,
    byrow = TRUE
  )
  mock_per_subject <- matrix(
    c(
      1, 5,
      2, 6,
      3, 7,
      4, 8
    ),
    nrow = nb + 2L,
    byrow = TRUE
  )

  local_mocked_bindings(
    ptmg_ll_der_hes4 = function(...) {
      list(
        gradient = mock_gradient,
        hessian = mock_hessian,
        per_subject_gradients = mock_per_subject
      )
    },
    apply_chain_rule_transform = function(...) {
      list(
        U = mock_gradient,
        H = mock_hessian,
        S = NULL,
        U_per_subject = mock_per_subject
      )
    },
    hl_efficient_components = function(H_bb, H_bl, H_ll, score, nb) {
      list(
        H_beta_beta_schur = diag(nb),
        score_beta_eff = rep(0, nb)
      )
    },
    .package = "scRoPE"
  )

  res <- list(
    beta = c(0.1, -0.2),
    dispersion = c(2, phi),
    hl_state = list()
  )
  ctx <- list(
    nb = nb,
    pred = matrix(0, nrow = 1, ncol = 1),
    offset = 0,
    fid = c(1L, 2L),
    cumsumy = matrix(0, nrow = 1, ncol = 1),
    posind = list(integer()),
    nind = 1L,
    k = 2L,
    intcol = 1L
  )
  posv <- list(
    Y = numeric(),
    n_onetwo = c(0, 0),
    ytwo = numeric(),
    posindy = integer()
  )

  out <- refresh_fn(
    res = res,
    ctx = ctx,
    posv = posv,
    gene_index = 1L
  )

  # [beta, phi, tau] after converting log(phi) -> phi.
  expect_equal(out$hl_state$score, c(10, 20, 40 / phi, 30))
  # Guard against the prior swapped mapping [beta, tau/phi, log(phi)].
  expect_false(isTRUE(all.equal(out$hl_state$score, c(10, 20, 30 / phi, 40))))

  expected_per_subject <- rbind(
    mock_per_subject[1:nb, , drop = FALSE],
    mock_per_subject[nb + 2L, , drop = FALSE] / phi,
    mock_per_subject[nb + 1L, , drop = FALSE]
  )
  expect_equal(out$hl_state$score_subject, expected_per_subject)
  expect_equal(out$hl_state$meat, tcrossprod(expected_per_subject))

  h_scaled <- mock_hessian
  idx_tau <- nb + 1L
  idx_log_phi <- nb + 2L
  h_scaled[idx_log_phi, ] <- h_scaled[idx_log_phi, ] / phi
  h_scaled[, idx_log_phi] <- h_scaled[, idx_log_phi] / phi
  reorder_idx <- c(seq_len(nb), idx_log_phi, idx_tau)
  h_nat <- h_scaled[reorder_idx, reorder_idx, drop = FALSE]
  h_nat <- 0.5 * (h_nat + t(h_nat))

  expect_equal(out$hl_state$H_obs$H_bb, h_nat[1:nb, 1:nb, drop = FALSE])
  expect_equal(out$hl_state$H_obs$H_bphi, matrix(h_nat[1:nb, nb + 1L], ncol = 1))
  expect_equal(out$hl_state$H_obs$H_btau, matrix(h_nat[1:nb, nb + 2L], ncol = 1))
  expect_equal(out$hl_state$H_obs$H_phiphi, matrix(h_nat[nb + 1L, nb + 1L], ncol = 1))
  expect_equal(out$hl_state$H_obs$H_phitau, matrix(h_nat[nb + 1L, nb + 2L], ncol = 1))
  expect_equal(out$hl_state$H_obs$H_tautau, matrix(h_nat[nb + 2L, nb + 2L], ncol = 1))
  expect_equal(out$hl_state$H_blocks$H_bl, cbind(out$hl_state$H_obs$H_bphi, out$hl_state$H_obs$H_btau))
})
