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

test_that("adding the route selector preserves the default joint LN call", {
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
    ln_estimator = "joint_ln"
  )))

  expect_identical(implicit$summary, explicit$summary)
  expect_identical(implicit$overdispersion, explicit$overdispersion)
  expect_identical(implicit$availability, explicit$availability)
  expect_identical(implicit$provenance$route, "LN")
  expect_identical(
    implicit$provenance$inference_policy,
    "joint_ln_same_basin_wald_calibration_candidate"
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
    c(TRUE, TRUE, TRUE)
  )
  expect_identical(joint$provenance$schema_version, 3L)
  expect_identical(
    joint$provenance$inference_policy,
    "joint_ln_same_basin_wald_calibration_candidate"
  )
  expect_identical(
    joint$provenance$candidate_resolution_policy,
    "joint_ln_candidate_resolution_calibration_v1"
  )
  expect_identical(
    joint$provenance$candidate_resolution_mode,
    "calibration_candidate_enabled"
  )
  expect_identical(
    joint$provenance$pace_status,
    "experimental_asymptotic_opt_in"
  )
  expect_identical(
    joint$provenance$test_capabilities$status[[1L]],
    "available_when_inference_status_available"
  )
  expect_identical(
    joint$provenance$test_capabilities$status[2:3],
    rep("experimental_asymptotic_when_requested_and_eligible", 2L)
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
    do.call(scrope, c(base, list(ln_estimator = "legacy_post_pml"))),
    "retired"
  )
  expect_error(
    do.call(scrope, c(base, list(
      ln_operational_fallback = "legacy_post_pml"
    ))),
    "retired"
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

test_that("optional tests cannot reactivate retired legacy dispatch", {
  for (opt in c("lbfgs", "trust")) {
    expect_error(scrope(ln_estimator = "legacy_post_pml", opt = opt,
                        additional_tests = c("score", "lrt")), "retired")
  }
})

test_that("shared default estimator is inert for explicit non-LN routes", {
  local_mocked_bindings(
    hl_fo_public_scrope_hl = function(...) list(route = "HL"),
    scrope_pmg_internal = function(...) list(route = "PGMM"),
    .package = "scRoPE"
  )
  for (model in c("HL", "PGMM")) {
    expect_identical(scrope(NULL, NULL, model = model),
                     scrope(NULL, NULL, model = model, ln_estimator = "joint_ln"))
  }
})
