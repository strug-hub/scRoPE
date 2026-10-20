test_that("legacy estimator and fallback requests are rejected before fitting", {
  for (model in c("LN", "HL", "PGMM")) {
    expect_error(
      scrope(stop("count must not be evaluated"), stop("id must not be evaluated"),
             model = model, ln_estimator = "legacy_post_pml"),
      "legacy_post_pml.*retired.*joint_ln.*exact preserved scRoPE version"
    )
    expect_error(
      scrope(stop("count must not be evaluated"), stop("id must not be evaluated"),
             model = model, ln_operational_fallback = "legacy_post_pml"),
      "legacy_post_pml.*retired.*none.*exact preserved scRoPE version"
    )
  }
  expect_error(scrope(ln_estimator = "legacy"), "retired")
  expect_error(scrope(ln_operational_fallback = "legacy"), "retired")
})

test_that("default and explicit joint LN preserve estimates and conditional modes", {
  fx <- ln_joint_test_fixture()
  args <- list(count = fx$count, id = fx$id, pred = fx$pred, offset = fx$offset,
               verbose = FALSE, cpc = 0, mincp = 1, output_re = TRUE)
  default <- do.call(scrope, args)
  explicit <- do.call(scrope, c(args, list(ln_estimator = "joint_ln")))
  expect_identical(formals(scrope)$ln_estimator, "joint_ln")
  for (field in c("summary", "overdispersion", "random_effect", "convergence",
                  "algorithm", "availability")) {
    expect_identical(default[[field]], explicit[[field]])
  }
  expect_identical(default$provenance$estimator, "ln_joint_root")
  expect_identical(default$provenance$random_effect, "pml_mode_given_ln_joint_root")
  expect_false("operational_fallback" %in% names(default))
})

test_that("retired exported nebula LN is not silently translated", {
  expect_error(nebula(), "default nebula.*retired.*joint_ln")
  for (model in c("NBGMM", "PMM", "NBLMM")) {
    expect_error(
      nebula(stop("count must not be evaluated"), model = model, method = "LN"),
      "retired.*exact preserved scRoPE version"
    )
  }
  for (method in list(c(method = "LN"), factor("LN"))) {
    expect_error(nebula(stop("count must not be evaluated"), method = method),
                 "retired.*exact preserved scRoPE version")
  }
  # The explicit historical HL API still enters its own validation path.
  expect_error(nebula(method = "HL", reml = 2), "reml|REML")
})
