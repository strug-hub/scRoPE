test_that("PGMM callers preserve raw status and covariance failures independently", {
  run_case <- function(caller, backend, raw_code, covariance) {
    env <- new.env(parent = baseenv())
    for (name in c("pmg_optimizer_converged", "pmg_run_bobyqa", "pmg_initial_beta",
                   "pmg_vector_norms", "pmg_reduced_mapping",
                   "fit_gene_pmg_unconstrained", "fit_gene_pmg_constrained")) {
      fun <- getFromNamespace(name, "scRoPE")
      environment(fun) <- env
      env[[name]] <- fun
    }
    observed <- new.env(parent = emptyenv())
    observed$nlminb_calls <- 0L
    observed$bobyqa_calls <- 0L
    observed$objective_stub_calls <- 0L
    expected_lower <- if (caller == "unconstrained") c(-100, -100, 0.01) else c(-100, 0.01)
    expected_upper <- if (caller == "unconstrained") c(100, 100, 10) else c(100, 10)
    expected_start <- if (caller == "unconstrained") c(0, 0, 1) else c(0, 1)
    make_raw <- function(par) list(
      par = par, convergence = raw_code, message = "mock termination detail",
      iterations = 11L, evaluations = setNames(c(13L, 7L), c("function", "gradient")),
      sentinel = "raw return must survive unchanged"
    )
    check_arguments <- function(start, lower, upper) {
      stopifnot(identical(start, expected_start),
                identical(lower, expected_lower), identical(upper, expected_upper))
    }
    env$nlminb <- function(start, objective, gradient, lower, upper, control = NULL) {
      observed$nlminb_calls <- observed$nlminb_calls + 1L
      check_arguments(start, lower, upper)
      stopifnot(is.function(objective), is.function(gradient), is.null(control))
      if (backend == "bobyqa") stop("Mock nlminb unavailable; exercise fallback")
      observed$expected_raw <- make_raw(start)
      observed$expected_raw
    }
    mock_bobyqa <- function(x0, fn, lower, upper) {
      observed$bobyqa_calls <- observed$bobyqa_calls + 1L
      check_arguments(x0, lower, upper)
      stopifnot(is.function(fn), backend == "bobyqa")
      observed$expected_raw <- make_raw(x0)
      observed$expected_raw
    }
    source_adapter <- env$pmg_run_bobyqa
    env$pmg_run_bobyqa <- function(start, objective, lower, upper) {
      source_adapter(start, objective, lower, upper, optimizer = mock_bobyqa)
    }
    env$pmg_ll <- function(...) {
      observed$objective_stub_calls <- observed$objective_stub_calls + 1L
      0
    }
    env$pmg_der <- function(...) stop("No optimizer may evaluate a gradient")
    env$pmg_ll_der_hes <- function(...) list(
      loglik = 0, objective = 0, gradient = rep(0, 3L),
      hessian = -diag(3L), observed_info = diag(3L),
      per_subject_gradients = if (covariance == "robust_failure") {
        matrix(NA_real_, nrow = 3L, ncol = 3L)
      } else diag(3L)
    )
    env$pmg_try_invert <- function(mat) {
      if (covariance == "naive_failure") NULL else diag(nrow(mat))
    }
    env$compute_sandwich_variance2 <- function(...) {
      if (covariance == "robust_failure") stop("Mock sandwich failure")
      list(Var_beta_adjusted = diag(2L), Var_full = diag(3L))
    }
    env$extract_godambe_components <- function(...) list(mocked = TRUE)
    env$constraint_reparameterization <- function(L, b) {
      stopifnot(identical(L, matrix(c(0, 1), nrow = 1L)), identical(b, 0))
      list(beta_star = c(0, 0), null_basis = matrix(c(1, 0), ncol = 1L))
    }
    env$pmg_reduced_newton_diagnostics <- function(H_red, g_red, theta_red) list(
      H_red_positive_definite = TRUE, H_red_rcond = 1,
      newton_step_red = rep(0, length(theta_red)),
      newton_step_red_max_abs = 0, max_abs_newton_step_scaled = 0,
      gradient_quadratic_lr_bound = 0
    )
    ctx <- list(nb = 2L, intcol = 1L, moffset = 0, min = 0.01, max = 10, k = 3L)
    posv <- list(mct = 1)
    result <- if (caller == "unconstrained") {
      env$fit_gene_pmg_unconstrained(1L, posv, ctx, retry = FALSE)
    } else {
      env$fit_gene_pmg_constrained(
        1L, posv, ctx, L = matrix(c(0, 1), nrow = 1L), b = 0, retry = FALSE
      )
    }
    expected_convergence <- if ((backend == "nlminb" && raw_code == 0L) ||
                               (backend == "bobyqa" && raw_code %in% 1:4)) 1L else -50L
    if (covariance == "naive_failure") expected_convergence <- -25L
    if (covariance == "robust_failure" && expected_convergence > 0L) expected_convergence <- -27L
    stopifnot(identical(result$optimizer, backend),
              identical(result$raw_fit, observed$expected_raw),
              identical(result$diagnostics$raw_convergence, raw_code),
              identical(result$diagnostics$raw_message, "mock termination detail"),
              identical(result$convergence, expected_convergence),
              identical(result$diagnostics$retry_used, backend == "bobyqa"),
              identical(result$diagnostics$use_average_objective, FALSE),
              identical(result$theta$beta, c(0, 0)), identical(result$theta$subVar, 1),
              observed$objective_stub_calls == 1L,
              identical(dim(result$naive_cov), c(2L, 2L)),
              identical(dim(result$robust_cov), c(2L, 2L)))
    if (backend == "nlminb") {
      # A returned nonzero code with retry=FALSE must not activate fallback.
      stopifnot(observed$nlminb_calls == 1L, observed$bobyqa_calls == 0L)
    } else {
      # Preserve the pre-existing no-result path: first call, repeated call,
      # then BOBYQA. No successful optimizer result was retried here.
      stopifnot(observed$nlminb_calls == 2L, observed$bobyqa_calls == 1L)
    }
    if (covariance == "naive_failure") stopifnot(all(is.na(result$naive_cov)))
    else stopifnot(all(is.finite(result$naive_cov)))
    if (covariance == "robust_failure") stopifnot(all(is.na(result$robust_cov)))
    if (covariance == "present") stopifnot(all(is.finite(result$robust_cov)))
    if (caller == "constrained") {
      stopifnot(identical(result$lrt_gate$converged, expected_convergence > 0L),
                result$diagnostics$optimizer_attempts == 1L,
                identical(result$reduced$theta, expected_start))
    }
    data.frame(caller = caller, backend = backend, raw_code = raw_code,
               covariance = covariance, aggregate_code = result$convergence,
               nlminb_mock_calls = observed$nlminb_calls,
               bobyqa_mock_calls = observed$bobyqa_calls,
               real_optimizer_calls = 0L, statistical_fits = 0L, generators = 0L)
  }

  rows <- list()
  for (caller in c("unconstrained", "constrained")) {
    for (backend in c("nlminb", "bobyqa")) {
      codes <- if (backend == "nlminb") c(0L, 1L) else c(1L, 2L, 3L, 4L, 5L, 6L, -1L)
      for (raw_code in codes) {
        for (covariance in c("present", "naive_failure", "robust_failure")) {
          rows[[length(rows) + 1L]] <- run_case(caller, backend, raw_code, covariance)
        }
      }
    }
  }
  results <- do.call(rbind, rows)
  expect_equal(nrow(results), 54L)
  expect_true(all(results$real_optimizer_calls == 0L))
})
