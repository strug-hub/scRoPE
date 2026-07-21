test_that("legacy_post_pml protects declared version-0.2.3 compatibility fields", {
  data(sample_data)
  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
  fit <- scrope(
    sample_data$count[1L, , drop = FALSE],
    sample_data$sid,
    pred = design,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ln_estimator = "legacy_post_pml"
  )

  # Compact values captured before the joint-LN implementation. Coefficients,
  # dispersions, model covariance, convergence, and architecture are protected;
  # robust SE values intentionally changed after the zero-based donor fix.
  expect_equal(
    as.numeric(fit$summary[1L, 1:4]),
    c(-1.902455437, -0.016755225, -0.097867225, 0.047278197),
    tolerance = 2e-7
  )
  expect_equal(
    as.numeric(fit$summary[1L, 5:8]),
    c(0.06335820, 0.03534659, 0.06449424, 0.06879634),
    tolerance = 2e-7
  )
  expect_equal(
    as.numeric(fit$overdispersion[1L, ]),
    c(0.0812525643, 0.8840821),
    tolerance = 2e-7
  )
  expect_identical(as.numeric(fit$convergence), 1)
  expect_identical(as.numeric(fit$algorithm), 1)
  expect_identical(fit$provenance$estimator, "legacy_post_pml_output_reproduction")
  expect_identical(fit$provenance$per_gene$route, "R1")
  expect_false(fit$provenance$per_gene$inference_available)

  old_robust_se <- c(0.06687127, 0.03155076, 0.04847650, 0.06201601)
  corrected_robust_se <- as.numeric(fit$summary[1L, 9:12])
  expect_gt(max(abs(corrected_robust_se - old_robust_se)), 1e-6)
  expect_equal(
    corrected_robust_se,
    c(0.06686006, 0.03156698, 0.04847424, 0.06205160),
    tolerance = 2e-7
  )
})

test_that("the untouched LN default is explicit legacy output reproduction", {
  fx <- ln_joint_test_fixture()
  default <- scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    output_re = TRUE
  )
  explicit <- scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    output_re = TRUE,
    ln_estimator = "legacy_post_pml"
  )
  expect_equal(default$summary, explicit$summary, tolerance = 0)
  expect_equal(default$overdispersion, explicit$overdispersion, tolerance = 0)
  expect_equal(default$random_effect, explicit$random_effect, tolerance = 0)
  expect_identical(default$provenance$per_gene$route, "R1")
  expect_identical(default$provenance$per_gene$mode_method, "legacy_pml_mode")
})

test_that("legacy LN marginal and conditional residual coordinates are preserved", {
  fx <- ln_joint_test_fixture()
  fit <- scrope(
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    output_re = TRUE,
    ln_estimator = "legacy_post_pml"
  )
  marginal <- nbresidual(
    fit,
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    conditional = FALSE
  )$residuals
  conditional <- nbresidual(
    fit,
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    conditional = TRUE
  )$residuals

  design <- as.matrix(fx$pred)
  intercept <- which(apply(design, 2L, sd) == 0)
  varying <- setdiff(seq_len(ncol(design)), intercept)
  design[, varying] <- scale(
    design[, varying, drop = FALSE],
    center = TRUE,
    scale = FALSE
  )
  design[, intercept] <- 1
  beta <- as.numeric(fit$summary[1L, seq_len(ncol(design))])
  sigma2 <- fit$overdispersion$Subject[[1L]]
  cell_variance <- fit$overdispersion$Cell[[1L]]
  eta_fixed <- log(fx$offset) + as.numeric(design %*% beta)
  mu_marginal <- exp(eta_fixed + sigma2 / 2)
  var_marginal <- mu_marginal + mu_marginal^2 *
    (exp(sigma2) * (1 + cell_variance) - 1)
  donor_id <- as.numeric(factor(
    as.character(fx$id),
    levels = unique(as.character(fx$id))
  ))
  w <- as.numeric(fit$random_effect[1L, ])
  mu_conditional <- exp(eta_fixed + w[donor_id])
  var_conditional <- mu_conditional + mu_conditional^2 * cell_variance
  observed <- as.numeric(fx$count[1L, ])

  expect_equal(
    as.numeric(marginal[1L, ]),
    (observed - mu_marginal) / sqrt(var_marginal),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(conditional[1L, ]),
    (observed - mu_conditional) / sqrt(var_conditional),
    tolerance = 1e-10
  )
})
