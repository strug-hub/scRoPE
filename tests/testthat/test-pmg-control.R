test_that("public PGMM budgets are explicit and validated before fitting", {
  resolve <- getFromNamespace("pmg_resolve_control", "scRoPE")
  expect_identical(resolve(), list(iter.max = 500L, eval.max = 500L))
  expect_identical(resolve(list()), resolve())
  expect_identical(resolve(list(iter.max = 25)), list(iter.max = 25L, eval.max = 500L))
  expect_identical(resolve(list(eval.max = 42)), list(iter.max = 500L, eval.max = 42L))
  expect_identical(resolve(list(iter.max = 150, eval.max = 200)),
                   list(iter.max = 150L, eval.max = 200L))
  for (value in list(NA_real_, NaN, Inf, -1, 0, 1.5, TRUE, "500", c(1, 2),
                     list(1), NULL, 1 + 0i, .Machine$integer.max + 1)) {
    expect_error(resolve(list(iter.max = value)), "positive finite integer")
  }
  for (control in list(500, c(iter.max = 500), list(500),
                       list(rel.tol = 1e-10), list(iter.max = 1, iter.max = 2),
                       setNames(list(1), NA_character_), data.frame(iter.max = 1))) {
    expect_error(resolve(control), "named list")
  }
})

test_that("PGMM public controls reach its dispatcher and survive disabled diagnostics", {
  fn <- getFromNamespace("scrope", "scRoPE")
  env <- new.env(parent = environment(fn))
  environment(fn) <- env
  seen <- NULL
  env$scrope_pmg_internal <- function(..., pgmm_control, keep_diagnostics) {
    seen <<- pgmm_control
    expect_false(keep_diagnostics)
    list(summary = data.frame(), diagnostics = NULL)
  }
  result <- fn(count = NULL, id = NULL, model = "PGMM", keep_diagnostics = FALSE)
  expect_identical(seen, list(iter.max = 500L, eval.max = 500L))
  expect_identical(result$provenance$optimizer$main_control, seen)
  expect_identical(result$provenance$optimizer$scope, "main_unconstrained_fit")
  expect_false(result$provenance$optimizer$retry_requested)
  result <- fn(count = NULL, id = NULL, model = "PGMM", keep_diagnostics = FALSE,
               pgmm_control = list(iter.max = 150L, eval.max = 200L))
  expect_identical(seen, list(iter.max = 150L, eval.max = 200L))
  for (model in c("LN", "HL")) {
    expect_error(fn(count = NULL, id = NULL, model = model,
                    pgmm_control = list(iter.max = 500L)), "applies only")
  }
  expect_error(fn(count = NULL, id = NULL, model = "PGMM",
                  pgmm_control = list(typo = 500L)), "named list")
})
