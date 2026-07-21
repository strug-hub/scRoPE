test_that("joint LN reports one certified objective state and aligned covariance", {
  result <- ln_joint_fit_fixture()
  fit <- result$fit
  fx <- result$fx
  raw <- fit$diagnostics$raw[[1L]]
  state <- raw$ln_joint$state
  covariance <- raw$ln_joint$covariance
  prep <- ln_joint_prepare_fixture(fx)

  expect_identical(raw$route, "R0")
  expect_true(raw$inference_available)
  expect_equal(raw$theta, raw$theta_ln, tolerance = 0)
  expect_equal(raw$theta, raw$theta_unrescued, tolerance = 0)
  expect_equal(raw$loglik, raw$loglik_ln, tolerance = 0)
  expect_equal(raw$loglik, -state$raw$value, tolerance = 1e-10)
  expect_equal(raw$score, state$U_full, tolerance = 0)
  expect_equal(rowSums(state$U_subject_full), state$U_full, tolerance = 1e-12)
  expect_gte(length(raw$ln_joint$candidate_records), 4L)
  expect_true(all(vapply(
    raw$ln_joint$candidate_records,
    function(candidate) !is.null(candidate$start),
    logical(1)
  )))
  expect_true(all(c(
    "optimizer_endpoint",
    "optimizer_objective",
    "theta_reduced",
    "theta_full_log",
    "root_residual",
    "newton_displacement"
  ) %in% names(raw$ln_joint$candidates)))

  expect_equal(raw$naive_cov, covariance$model_beta, tolerance = 0)
  expect_equal(raw$robust_cov, covariance$robust_beta, tolerance = 0)
  expect_lte(
    covariance$direct_influence_error,
    ln_joint_thresholds()$covariance_identity
  )
  if (isTRUE(covariance$profile_available)) {
    expect_lte(covariance$profile_error, ln_joint_thresholds()$covariance_identity)
  }

  beta_public <- as.numeric(
    fit$provenance$coefficient_map_internal_to_public %*% state$beta
  )
  expect_equal(as.numeric(fit$summary[1L, seq_len(prep$ctx$nb)]), beta_public)
  expect_equal(
    as.numeric(fit$summary[1L, paste0("se_", colnames(fx$pred))]),
    sqrt(diag(covariance$model_beta)) / prep$public_scale
  )
  expect_equal(
    as.numeric(fit$summary[1L, paste0("se_robust_", colnames(fx$pred))]),
    sqrt(diag(covariance$robust_beta)) / prep$public_scale
  )
})

test_that("joint LN derivative and coordinate identities pass locked targets", {
  result <- ln_joint_fit_fixture()
  fit <- result$fit
  prep <- ln_joint_prepare_fixture(result$fx)
  state <- fit$diagnostics$raw[[1L]]$ln_joint$state
  audit <- ln_joint_fd_audit(
    state$theta_reduced_log,
    1L,
    prep$posv,
    prep$ctx,
    ln_joint_map(prep$ctx$nb)
  )

  expect_true(all(audit$target_pass), info = paste(
    names(audit$target_pass)[!audit$target_pass],
    collapse = ", "
  ))
  expect_true(all(audit$hard_stop_pass))
  expect_lte(audit$metrics[["donor_additivity"]], 1e-10)
  expect_lte(audit$metrics[["objective_gradient"]], 1e-5)
  expect_lte(audit$metrics[["sensitivity_score"]], 1e-4)
  expect_lte(audit$metrics[["natural_log_gradient"]], 1e-7)
  expect_lte(audit$metrics[["natural_log_sensitivity"]], 1e-7)

  beta_internal <- state$beta
  beta_public <- as.numeric(
    fit$provenance$coefficient_map_internal_to_public %*% beta_internal
  )
  public_design <- as.matrix(result$fx$pred)
  varying <- setdiff(seq_len(ncol(public_design)), prep$ctx$intcol)
  public_design[, varying] <- scale(
    public_design[, varying, drop = FALSE],
    center = TRUE,
    scale = FALSE
  )
  public_design[, prep$ctx$intcol] <- 1
  expect_equal(
    as.numeric(prep$ctx$pred %*% beta_internal),
    as.numeric(public_design %*% beta_public),
    tolerance = 1e-10
  )
})

test_that("corrected LN donor allocation includes observation zero and donor boundaries", {
  count <- Matrix::Matrix(
    matrix(c(1, 0, 0, 0, 2, 0, 3, 0), nrow = 1L),
    sparse = TRUE
  )
  id <- rep(1:2, each = 4)
  pred <- cbind("(Intercept)" = 1, x = seq(-1, 1, length.out = 8))
  fx <- list(count = count, id = id, pred = pred, offset = rep(1, 8))
  prep <- ln_joint_prepare_fixture(fx)
  theta <- c(0.1, -0.2, log(0.4), log(1.3))

  corrected <- ptmg_ll_der_hes4(
    theta,
    prep$ctx$pred,
    prep$ctx$offset,
    prep$posv$Y,
    prep$posv$n_onetwo,
    prep$posv$ytwo,
    prep$ctx$fid,
    prep$ctx$cumsumy[1L, ],
    prep$ctx$posind[[1L]],
    prep$posv$posindy,
    prep$ctx$nb,
    prep$ctx$nind,
    prep$ctx$k,
    correct_posindy = TRUE
  )

  expect_identical(prep$posv$posindy[c(1L, 2L)], c(0L, 4L))
  expect_equal(
    rowSums(corrected$per_subject_gradients),
    corrected$gradient,
    tolerance = 1e-12
  )
  expect_equal(
    rowSums(corrected$per_subject_gradients),
    corrected$gradient_aggregate_formula,
    tolerance = 1e-10
  )
  expect_equal(
    corrected$gradient_aggregate_formula_raw +
      corrected$aggregate_formula_roundoff_reconciliation_log,
    corrected$gradient_aggregate_formula,
    tolerance = 0
  )
  expect_lte(
    max(abs(c(
      corrected$base_gradient_roundoff_reconciliation_log,
      corrected$donor_roundoff_reconciliation_log_gamma,
      corrected$gamma_special_roundoff_reconciliation_log,
      corrected$aggregate_formula_roundoff_reconciliation_log
    ))),
    ln_joint_thresholds()$hard_stop$donor_additivity
  )

  donor_objective <- function(z) {
    ln_joint_donor_objective_values(
      z,
      1L,
      prep$posv,
      prep$ctx,
      ln_joint_map(prep$ctx$nb)
    )
  }
  step <- .Machine$double.eps^(1 / 3) * (1 + abs(theta))
  donor_fd <- vapply(seq_along(theta), function(jj) {
    plus <- minus <- theta
    plus[jj] <- plus[jj] + step[jj]
    minus[jj] <- minus[jj] - step[jj]
    (donor_objective(plus) - donor_objective(minus)) / (2 * step[jj])
  }, numeric(prep$ctx$k))
  donor_score_error <- ln_joint_scaled_max(
    corrected$per_subject_gradients,
    t(donor_fd)
  )
  expect_lte(donor_score_error, ln_joint_thresholds()$objective_gradient)
  expect_equal(
    as.numeric(corrected$per_subject_gradients[prep$ctx$nb + 2L, ]),
    as.numeric(t(donor_fd)[prep$ctx$nb + 2L, ]),
    tolerance = 1e-5
  )
})
