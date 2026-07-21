test_that("joint LN constraints use the public-to-internal map and reduced root", {
  base <- ln_joint_fit_fixture()
  beta_public <- as.numeric(base$fit$summary[1L, seq_len(ncol(base$fx$pred))])
  L_mixed <- matrix(c(1, 0, -0.5), nrow = 1L)
  b_mixed <- as.numeric(L_mixed %*% beta_public + 0.02)
  L_multi <- rbind(c(1, 0, 0), c(0, 0, 1))
  b_multi <- as.numeric(L_multi %*% beta_public + c(0.01, -0.01))
  L_intercept <- matrix(c(0, 1, 0), nrow = 1L)
  b_intercept <- beta_public[2L] + 0.01
  contrasts <- list(
    list(name = "mixed_nonzero", L = L_mixed, b = b_mixed),
    list(name = "mixed_rescaled", L = 7 * L_mixed, b = 7 * b_mixed),
    list(name = "multirow", L = L_multi, b = b_multi),
    list(name = "intercept_nonzero", L = L_intercept, b = b_intercept)
  )
  result <- ln_joint_fit_fixture(
    additional_tests = c("score", "lrt"),
    lrt_contrasts = contrasts
  )
  fit <- result$fit
  prep <- ln_joint_prepare_fixture(result$fx)
  M <- fit$provenance$coefficient_map_internal_to_public

  expect_true(all(fit$lrt_details$method == "unavailable_joint_ln"))
  expect_true(all(!fit$lrt_details$success))
  expect_true(all(!fit$lrt_details$score_success))
  expect_true(all(is.na(fit$lrt_details$stat)))
  expect_true(all(is.na(fit$lrt_details$p_value)))
  expect_true(all(is.na(fit$lrt_details$score_stat)))
  expect_true(all(is.na(fit$lrt_details$score_p_value)))
  expect_true(all(
    fit$lrt_details$failure_reason == "joint_ln_calibration_pending"
  ))
  expect_true(all(fit$lrt_details$constrained_route == "R0"))
  expect_true(all(fit$lrt_details$gate_converged))
  expect_true(all(
    fit$lrt_details$convergence_basis == "joint_ln_endpoint_pair"
  ))
  expect_false(any(fit$lrt_details$fallback))

  raw_tests <- fit$diagnostics$raw[[1L]]$lrt
  definitions <- list(
    mixed_nonzero = list(L = L_mixed, b = b_mixed),
    mixed_rescaled = list(L = 7 * L_mixed, b = 7 * b_mixed),
    multirow = list(L = L_multi, b = b_multi),
    intercept_nonzero = list(L = L_intercept, b = b_intercept)
  )
  for (name in names(definitions)) {
    constrained <- raw_tests[[name]]$constrained
    definition <- definitions[[name]]
    L_internal <- definition$L %*% M
    expect_equal(constrained$constraint$L, L_internal, tolerance = 0)
    expect_equal(
      as.numeric(definition$L %*% M %*% constrained$theta$beta),
      definition$b,
      tolerance = 1e-9
    )
    expect_equal(
      constrained$reduced$score,
      as.numeric(crossprod(constrained$reduced$map, constrained$score)),
      tolerance = 1e-12
    )
    expect_lte(
      ln_joint_root_metric(
        constrained$reduced$per_subject_gradients,
        constrained$reduced$score
      ),
      ln_joint_thresholds()$root
    )
    audit <- ln_joint_fd_audit(
      constrained$reduced$theta,
      1L,
      prep$posv,
      prep$ctx,
      ln_joint_map(prep$ctx$nb, L_internal, definition$b)
    )
    expect_true(all(audit$target_pass), info = name)
  }

  expect_equal(
    raw_tests$mixed_nonzero$constrained$theta$beta,
    raw_tests$mixed_rescaled$constrained$theta$beta,
    tolerance = 1e-7
  )
  expect_equal(
    raw_tests$mixed_nonzero$constrained$loglik_ln,
    raw_tests$mixed_rescaled$constrained$loglik_ln,
    tolerance = 1e-8
  )
})

test_that("adjusted score cannot bypass an explicit unavailable gate", {
  blocked <- compute_adj_score(
    constrained = list(),
    L = matrix(c(1, 0), nrow = 1L),
    score_gate = list(
      available = FALSE,
      failure_reason = "joint_ln_calibration_pending"
    )
  )
  expect_false(blocked$success)
  expect_true(is.na(blocked$stat))
  expect_true(is.na(blocked$pval))
  expect_identical(blocked$failure_reason, "joint_ln_calibration_pending")
})
