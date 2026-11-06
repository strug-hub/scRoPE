research_fixture <- function() readRDS(test_path("fixtures", "jln-research.rds"))
research_call <- function(x, ...) {
  args <- list(count = x$input$y, pred = x$input$X, id = x$input$id,
    offset = x$input$offset, target = "x", cpc = 0, mincp = 1L,
    ncore = 1L, verbose = FALSE)
  do.call(scrope_jln_research, modifyList(args, list(...)))
}

test_that("the public research routes reproduce the frozen prototype", {
  f <- research_fixture()
  for (name in c("zero", "ordinary", "warning", "corner")) {
    x <- f[[name]]
    fit <- suppressWarnings(research_call(x))
    expect_identical(fit$summary$route, x$expected$route)
    expect_true(fit$summary$available)
    expect_equal(fit$summary$coefficient, x$expected$beta, tolerance = 1e-9)
    expect_equal(fit$summary$se_robust, x$expected$se, tolerance = 1e-9)
    expect_equal(fit$summary$p_robust, x$expected$p, tolerance = 1e-9)
    expect_false(fit$summary$selection_adjusted)
    expect_identical(fit$summary$finite_cluster_multiplier, 1)
  }
  expect_warning(research_call(f$warning), "Opt-in JLN research workflow")
  expect_warning(research_call(f$zero), "not selection-adjusted")
})

test_that("native results and disabled research policies remain explicit", {
  f <- research_fixture()
  fit <- research_call(f$ordinary)
  x <- f$ordinary$input
  native <- scrope(count = x$y, id = x$id, pred = x$X, offset = x$offset,
    model = "LN", cpc = 0, mincp = 1L, ncore = 1L, verbose = FALSE,
    recovery = "none", keep_diagnostics = TRUE)
  expect_recovery_refit_table(fit$primary$summary, native$summary)
  expect_recovery_refit_table(fit$primary$overdispersion, native$overdispersion,
    numerical = c("Subject", "Cell"))
  expect_recovery_refit_numeric(fit$primary$diagnostics$raw[[1]]$ln_joint$state$A_full,
                   native$diagnostics$raw[[1]]$ln_joint$state$A_full)
  expect_false(fit$summary$research_added)
  disabled <- suppressWarnings(research_call(f$zero, boundary = FALSE, dominated_warning = FALSE))
  expect_identical(disabled$summary$route, "withheld")
  expect_false(disabled$summary$available)
  expect_true(all(is.na(disabled$summary[c("coefficient", "se_robust", "p_robust", "ci_lower", "ci_upper")])))
  warn_only <- suppressWarnings(research_call(f$warning, boundary = FALSE))
  boundary_only <- suppressWarnings(research_call(f$warning, dominated_warning = FALSE))
  expect_identical(warn_only$summary$route, "ordinary_jln_with_search_warning")
  expect_identical(boundary_only$summary$route, "withheld")
})

test_that("target selection and diagnostic retention do not change the fit", {
  f <- research_fixture()$zero
  full <- suppressWarnings(research_call(f, target = "z"))
  lean <- suppressWarnings(research_call(f, target = "z", keep_diagnostics = FALSE))
  expect_equal(full$summary$coefficient, unname(f$z_expected["beta"]), tolerance = 1e-9)
  expect_equal(full$summary$se_robust, unname(f$z_expected["se"]), tolerance = 1e-9)
  expect_recovery_refit_table(full$summary, lean$summary)
  expect_null(lean$diagnostics)
  expect_null(lean$primary$diagnostics)
  expect_false(is.null(full$diagnostics))
})

test_that("filtering and duplicate display names use stable gene IDs", {
  f <- research_fixture()$zero
  count <- rbind(rep(0, length(f$input$y)), f$input$y, f$input$y)
  rownames(count) <- c("filtered", "same", "same")
  fit <- suppressWarnings(research_call(f, count = count))
  expect_equal(fit$summary$gene_id, c(2L, 3L))
  expect_identical(fit$summary$gene, c("same", "same"))
  expect_equal(fit$summary$coefficient, rep(f$expected$beta, 2), tolerance = 1e-9)
  expect_equal(fit$summary$se_robust, rep(f$expected$se, 2), tolerance = 1e-9)
})

test_that("unsupported inputs and options fail before fitting", {
  f <- research_fixture()$zero
  expect_error(research_call(f, target = "(Intercept)"), "non-intercept")
  expect_error(research_call(f, target = c("x", "z")), "non-intercept")
  expect_error(research_call(f, boundary = NA), "TRUE or FALSE")
  expect_error(research_call(f, offset = rep(0, length(f$input$y))), "positive exposure")
  expect_error(research_call(f, ncore = 0), "positive integers")
  expect_error(research_call(f, additional_tests = "score"), "unused argument")
  id <- f$input$id; id[c(1, 9)] <- id[c(9, 1)]
  expect_error(research_call(f, id = id), "contiguous")
})

test_that("the finite-gradient trial repair retains the matched sandwich", {
  f <- research_fixture()$failed_comparator
  dat <- scRoPE:::ln_research_context(f$input$y, f$input$X, f$input$id, f$input$offset)
  fit <- scRoPE:::ln_research_fixed_fit(dat, f$fixed, f$starts)
  expect_true(fit$success)
  expect_gt(length(fit$trace), 0L)
  expect_equal(fit$coefficient, f$expected$coefficient, tolerance = 1e-9)
  expect_equal(fit$se, f$expected$se, tolerance = 1e-9)
  st <- fit$state; inv <- solve(st$A)
  expect_equal(fit$covariance, inv %*% tcrossprod(st$U) %*% t(inv), tolerance = 1e-10)
  expect_identical(st$A, st$full_state$A_full[fit$spec$free, fit$spec$free, drop = FALSE])
  expect_equal(rowSums(st$U), -st$g, tolerance = 1e-10)
  expect_true(all(is.finite(c(st$A, st$U))))
})
