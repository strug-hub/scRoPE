test_that("PGMM honors backend-specific optimizer codes and the BOBYQA x0 interface", {
  env <- new.env(parent = baseenv())
  for (name in c("pmg_optimizer_converged", "pmg_run_bobyqa")) {
    env[[name]] <- getFromNamespace(name, "scRoPE")
  }
  fixture <- function(optimizer, code) {
    list(optimizer = optimizer, raw = list(convergence = code))
  }
  for (code in c(-5:-1, 0:8)) {
    stopifnot(identical(env$pmg_optimizer_converged(fixture("nlminb", code)), code == 0))
    stopifnot(identical(env$pmg_optimizer_converged(fixture("bobyqa", code)), code %in% 1:4))
  }
  for (code in list(NULL, NA_real_, Inf, NaN, numeric(), c(0, 1), "0", 0.5, list(0), 1+0i)) {
    stopifnot(!env$pmg_optimizer_converged(fixture("nlminb", code)))
    stopifnot(!env$pmg_optimizer_converged(fixture("bobyqa", code)))
  }
  stopifnot(!env$pmg_optimizer_converged(NULL))
  stopifnot(!env$pmg_optimizer_converged(list(raw = 0)))
  stopifnot(!env$pmg_optimizer_converged(fixture("unknown", 0)))

  # Completed gene fits use raw_fit; exact field access must fail closed.
  for (code in c(0L, 1L, NA_integer_)) {
    expect_identical(env$pmg_optimizer_converged(list(
      optimizer = "nlminb", raw_fit = list(convergence = code))), isTRUE(code == 0L))
  }
  expect_false(env$pmg_optimizer_converged(list(
    optimizer = "nlminb", raw_fitting = list(convergence = 0L))))

  expect_false(env$pmg_optimizer_converged(list(
    optimizer_backend = "nlminb", raw_fit = list(convergence = 0L))))

  # Requiring x0 catches the old par= argument mismatch. The objective must
  # never be evaluated; this checks interface plumbing, not optimization.
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  mock_bobyqa <- function(x0, fn, lower, upper) {
    observed$calls <- observed$calls + 1L
    observed$start <- x0
    observed$lower <- lower
    observed$upper <- upper
    list(par = x0, convergence = 4L, iter = 7L, message = "NLOPT_XTOL_REACHED")
  }
  result <- env$pmg_run_bobyqa(
    start = c(0, 1),
    objective = function(...) stop("No objective evaluation is allowed"),
    lower = c(-1, 0.01), upper = c(1, 10), optimizer = mock_bobyqa
  )
  stopifnot(observed$calls == 1L)
  stopifnot(identical(observed$start, c(0, 1)))
  stopifnot(identical(observed$lower, c(-1, 0.01)))
  stopifnot(identical(observed$upper, c(1, 10)))
  stopifnot(identical(result$par, observed$start))
  stopifnot(identical(result$optimizer, "bobyqa"))
  stopifnot(identical(result$raw$convergence, 4L))
  stopifnot(env$pmg_optimizer_converged(result))
  for (code in c(5L, 6L)) {
    result$raw$convergence <- code
    stopifnot(!env$pmg_optimizer_converged(result))
  }
  expect_true(TRUE)
})
