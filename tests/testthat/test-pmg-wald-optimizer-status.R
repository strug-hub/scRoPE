pgmm_wald_status_fixture <- function(code = 0L, convergence = 1L,
                                     optimizer = "nlminb") {
  list(
    theta = list(beta = c(0.4, 0.6), subVar = 0.3),
    naive_cov = diag(c(0.2, 0.4)),
    robust_cov = diag(c(0.3, 0.5)),
    convergence = as.integer(convergence),
    optimizer = optimizer,
    raw_fit = list(convergence = code),
    loglik = -10,
    godambe = NULL
  )
}

pgmm_mock_wald_summary <- function(fits) {
  fn <- getFromNamespace("scrope_pmg_internal", "scRoPE")
  env <- new.env(parent = environment(fn))
  environment(fn) <- env
  pred <- cbind(Intercept = 1, x = c(-1, -1, 1, 1))
  count <- matrix(1L, nrow = length(fits), ncol = 4L)
  prep <- list(
    ctx = list(nb = 2L, nind = c(2L, 2L), k = 2L, intcol = 1L,
               pred = pred, fid = c(1L, 3L, 5L)),
    count = count,
    gid = seq_along(fits),
    gene_names = paste0("gene", seq_along(fits)),
    predictor_names = colnames(pred),
    design_sds = c(1, 1)
  )
  env$hl_prepare_data <- function(...) prep
  env$call_posindy <- function(...) NULL
  env$fit_gene_pmg_unconstrained <- function(gene_index, ...) fits[[gene_index]]
  env$registerDoFuture <- function(...) invisible(NULL)
  env$plan <- function(...) invisible(NULL)
  env$foreach <- function(gene_idx) gene_idx
  env$`%dorng%` <- function(obj, expr) {
    block <- substitute(expr)
    eval_env <- new.env(parent = parent.frame())
    lapply(obj, function(idx) {
      eval_env$gene_idx <- idx
      eval(block, eval_env)
    })
  }
  fn(count = count, id = rep(1:2, each = 2L), pred = pred,
     verbose = FALSE, cpc = 0, mincp = 1, ncore = 1)
}

test_that("public PGMM Wald inference requires native optimizer success", {
  failed <- pgmm_wald_status_fixture(code = 1L, convergence = -50L)
  unknown <- pgmm_wald_status_fixture(code = 0L, optimizer = "unknown")
  missing <- pgmm_wald_status_fixture(code = NA_integer_)
  bobyqa_budget <- pgmm_wald_status_fixture(code = 5L, optimizer = "bobyqa")
  bobyqa_time <- pgmm_wald_status_fixture(code = 6L, optimizer = "bobyqa")
  res <- pgmm_mock_wald_summary(list(failed, unknown, missing,
                                    bobyqa_budget, bobyqa_time))
  sm <- res$summary

  expect_false(any(sm$optimizer_success))
  expect_true(all(is.na(sm$p_naive_x)))
  expect_true(all(is.na(sm$p_robust_x)))
  expect_true(all(is.na(sm$model_wald_stat_x)))
  expect_true(all(is.na(sm$wald_stat_x)))
  expect_false(any(sm$model_wald_success_x))
  expect_false(any(sm$wald_success_x))
  expect_equal(sm$model_wald_failure_x, rep("optimizer_not_converged", 5L))
  expect_equal(sm$wald_failure_x, rep("optimizer_not_converged", 5L))
  expect_equal(sm$logFC_x, rep(0.6, 5L))
  expect_true(all(is.na(sm$se_x)))
  expect_true(all(is.na(sm$se_robust_x)))
  expect_identical(res$diagnostics$raw[[1]]$theta, failed$theta)
  expect_identical(res$diagnostics$raw[[1]]$naive_cov, failed$naive_cov)
  expect_identical(res$diagnostics$raw[[1]]$robust_cov, failed$robust_cov)
  expect_true(all(res$diagnostics$test_status$failure_reason ==
                    "optimizer_not_converged"))
  expect_false(any(res$summary$test_status == "ok"))
})

test_that("model and sandwich Wald availability stay independent", {
  no_sandwich <- pgmm_wald_status_fixture(convergence = -27L)
  no_sandwich$robust_cov[,] <- NA_real_
  no_model <- pgmm_wald_status_fixture(convergence = -25L)
  no_model$naive_cov[,] <- NA_real_
  res <- pgmm_mock_wald_summary(list(no_sandwich, no_model))
  sm <- res$summary

  expect_true(all(sm$optimizer_success))
  expect_equal(sm$model_wald_success_x, c(TRUE, FALSE))
  expect_equal(sm$wald_success_x, c(FALSE, TRUE))
  expect_equal(sm$p_naive_x[1], pchisq(0.6^2 / 0.4, 1, lower.tail = FALSE))
  expect_equal(sm$p_robust_x[2], pchisq(0.6^2 / 0.5, 1, lower.tail = FALSE))
  expect_true(is.na(sm$p_naive_x[2]))
  expect_true(is.na(sm$p_robust_x[1]))
  expect_identical(sm$model_wald_failure_x[2], "singular_cov")
  expect_identical(sm$wald_failure_x[1], "singular_cov")
  expect_equal(res$convergence, c(-27L, -25L))
})

test_that("successful optimizer codes retain public Wald values", {
  fits <- c(list(pgmm_wald_status_fixture()), lapply(1:4, function(code) {
    pgmm_wald_status_fixture(code = code, optimizer = "bobyqa")
  }))
  sm <- pgmm_mock_wald_summary(fits)$summary
  expect_true(all(sm$optimizer_success))
  expect_true(all(sm$model_wald_success_x))
  expect_true(all(sm$wald_success_x))
  expect_equal(sm$model_wald_stat_x, rep(0.6^2 / 0.4, 5L))
  expect_equal(sm$wald_stat_x, rep(0.6^2 / 0.5, 5L))
  expect_true(all(is.na(sm$model_wald_failure_x)))
  expect_true(all(is.na(sm$wald_failure_x)))
  expect_true(all(sm$test_status == "ok"))
})

test_that("coefficient Wald rejects invalid variance and nonfinite statistics", {
  fn <- getFromNamespace("pmg_coefficient_wald", "scRoPE")
  res <- fn(beta = c(1, 1, Inf, 1), variance = c(0, -1, 1, NA_real_),
            optimizer_converged = rep(TRUE, 4L))
  expect_false(any(res$success))
  expect_true(all(is.na(res$stat)))
  expect_true(all(is.na(res$p_value)))
  expect_equal(res$failure_reason,
               c("singular_cov", "singular_cov", "undefined_stat", "singular_cov"))
})

test_that("contrast Wald preserves diagnostics while gating failed optimizers", {
  fn <- getFromNamespace("pmg_robust_wald", "scRoPE")
  L <- matrix(c(0, 1), nrow = 1L)
  fit <- pgmm_wald_status_fixture(code = 1L, convergence = -50L)
  res <- fn(fit, L, b = 0)
  expect_false(res$success)
  expect_true(is.na(res$stat))
  expect_identical(res$failure_reason, "optimizer_not_converged")
  expect_equal(res$delta, matrix(0.6))
  expect_equal(res$variance, matrix(0.5))

  fit$raw_fit$convergence <- 0L
  fit$convergence <- -25L
  res_ok <- fn(fit, L, b = 0)
  expect_true(res_ok$success)
  expect_equal(res_ok$stat, 0.6^2 / 0.5)
  expect_true(is.na(res_ok$failure_reason))
})
