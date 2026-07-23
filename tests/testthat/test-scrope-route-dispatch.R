test_that("top-level HL dispatch exactly matches the compatibility wrapper", {
  data(sample_data)
  design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
  args <- list(
    count = sample_data$count[1L, , drop = FALSE],
    id = sample_data$sid,
    pred = design,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    additional_tests = c("score", "lrt"),
    lrt_details = TRUE,
    keep_diagnostics = TRUE
  )

  direct <- do.call(scrope_hl, args)
  routed <- do.call(scrope, c(args, list(model = "HL")))

  expect_identical(routed, direct)
  expect_identical(routed$provenance$route, "HL")
  expect_identical(routed$provenance$requested_model, "HL")
  expect_identical(routed$provenance$effective_model, "HL")
  expect_identical(
    routed$provenance$test_capabilities$test,
    c("wald", "score", "lrt")
  )
  expect_identical(
    routed$provenance$test_capabilities$public_available,
    c(TRUE, TRUE, FALSE)
  )

  residual <- nbresidual(
    routed,
    count = sample_data$count[1L, , drop = FALSE],
    id = sample_data$sid,
    pred = design
  )
  expect_identical(dim(residual$residuals), c(1L, ncol(sample_data$count)))
  expect_true(all(is.finite(residual$residuals)))
})

test_that("adding the route selector preserves the default legacy LN call", {
  fx <- ln_joint_test_fixture()
  args <- list(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    keep_diagnostics = FALSE
  )

  set.seed(20260723)
  implicit <- do.call(scrope, args)
  set.seed(20260723)
  explicit <- do.call(scrope, c(args, list(
    model = "LN",
    ln_estimator = "legacy_post_pml"
  )))

  expect_identical(implicit, explicit)
  expect_identical(implicit$provenance$route, "LN")
  expect_identical(
    implicit$provenance$inference_policy,
    "legacy_compatibility_output"
  )
})

test_that("top-level routes expose honest method capabilities", {
  fx <- ln_joint_test_fixture()
  joint <- scrope(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    keep_diagnostics = FALSE,
    model = "LN",
    ln_estimator = "joint_ln"
  )
  expect_identical(
    joint$provenance$test_capabilities$public_available,
    c(TRUE, FALSE, FALSE)
  )
  expect_identical(
    joint$provenance$test_capabilities$status[2:3],
    rep("joint_ln_calibration_pending", 2L)
  )

  pmg <- pmg_subject_level_test_fixture()
  pgmm <- scrope(
    count = pmg$count,
    id = pmg$id,
    pred = pmg$pred,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    model = "PGMM",
    keep_diagnostics = FALSE
  )
  expect_identical(pgmm$provenance$route, "PGMM")
  expect_identical(pgmm$provenance$schema_version, 2L)
  expect_true(all(pgmm$provenance$test_capabilities$public_available))
})

test_that("HL dispatch rejects LN-only routing controls", {
  data(sample_data)
  base <- list(
    count = sample_data$count[1L, , drop = FALSE],
    id = sample_data$sid,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    model = "HL"
  )

  expect_error(
    do.call(scrope, c(base, list(ln_estimator = "joint_ln"))),
    "ln_estimator applies only"
  )
  expect_error(
    do.call(scrope, c(base, list(
      ln_operational_fallback = "legacy_post_pml"
    ))),
    "applies only to model='LN'"
  )
  expect_error(
    do.call(scrope, c(base, list(allow_per_gene_switch = TRUE))),
    "unavailable for model='HL'"
  )
  expect_error(
    do.call(scrope, c(base, list(cutoff_cell = 10))),
    "controls legacy LN refinement"
  )
  expect_error(
    do.call(scrope, c(base, list(kappa = 100))),
    "controls legacy LN refinement"
  )
})

test_that("PGMM dispatch rejects LN-only controls and warns for no-op outputs", {
  pmg <- pmg_subject_level_test_fixture()
  base <- list(
    count = pmg$count,
    id = pmg$id,
    pred = pmg$pred,
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    model = "PGMM",
    keep_diagnostics = FALSE
  )

  expect_error(
    do.call(scrope, c(base, list(allow_per_gene_switch = TRUE))),
    "unavailable for model='PGMM'"
  )
  expect_error(
    do.call(scrope, c(base, list(cutoff_cell = 10))),
    "controls legacy LN refinement"
  )
  expect_error(
    do.call(scrope, c(base, list(kappa = 100))),
    "controls legacy LN refinement"
  )
  expect_error(
    do.call(scrope, c(base, list(use_betas = "final_betas"))),
    "applies only to model='LN'"
  )

  warnings <- character(0)
  fit <- withCallingHandlers(
    do.call(scrope, c(base, list(opt = "trust", output_re = TRUE))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("opt='trust' is ignored", warnings, fixed = TRUE)))
  expect_true(any(grepl("output_re=TRUE is not available", warnings, fixed = TRUE)))
  expect_null(fit$random_effect)
})

test_that("skipped legacy optional tests warn and record the fit-specific reason", {
  fx <- ln_joint_test_fixture()
  expect_warning(
    fit <- scrope(
      count = fx$count,
      id = fx$id,
      pred = fx$pred,
      offset = fx$offset,
      verbose = FALSE,
      cpc = 0,
      mincp = 1,
      ncore = 1L,
      opt = "trust",
      additional_tests = c("score", "lrt"),
      lrt_details = TRUE,
      keep_diagnostics = FALSE
    ),
    "require opt='lbfgs'"
  )

  capabilities <- fit$provenance$test_capabilities
  expect_identical(
    capabilities$status[capabilities$test %in% c("score", "lrt")],
    rep("not_computed_opt_requires_lbfgs", 2L)
  )
  expect_null(fit$lrt_details)
})
