hia_get <- function(name) getFromNamespace(name, "scRoPE")
hia_fixture <- function() readRDS(test_path("fixtures", "hl-inner-event106.rds"))

hia_independent_score <- function(input, beta) {
  # Bracket each scalar root, then evaluate the score and Schur curvature
  # directly there. No production mode solver, evaluator or primitive is used.
  data <- input$data
  rho <- input$rho_full
  base <- data$log_offset + as.numeric(data$X %*% beta)
  v <- exp(rho[1]); phi <- exp(rho[2]); alpha <- 1 / expm1(v)
  eta <- vapply(data$subject_rows, function(rows) {
    equation <- function(z) {
      lm <- base[rows] + z
      sum(data$y[rows] * plogis(rho[2] - lm) - phi * plogis(lm - rho[2])) -
        alpha * expm1(z - v / 2)
    }
    uniroot(equation, c(-20, 20), tol = 8 * .Machine$double.eps,
            maxiter = 1000L)$root
  }, numeric(1))
  subjects <- lapply(seq_along(data$subject_rows), function(j) {
    rows <- data$subject_rows[[j]]
    X <- data$X[rows, , drop = FALSE]
    lm <- base[rows] + eta[j]
    r <- plogis(lm - rho[2]); one_minus_r <- plogis(rho[2] - lm)
    score <- data$y[rows] * one_minus_r - phi * r
    curvature <- -(data$y[rows] + phi) * r * one_minus_r
    mixed <- as.numeric(crossprod(X, curvature))
    donor_curvature <- sum(curvature) - alpha * exp(eta[j] - v / 2)
    list(score = as.numeric(crossprod(X, score)),
         curvature = crossprod(X, X * curvature) -
           tcrossprod(mixed) / donor_curvature)
  })
  S_beta <- Reduce(`+`, lapply(subjects, `[[`, "score"))
  B_beta <- Reduce(`+`, lapply(subjects, `[[`, "curvature"))
  S <- as.numeric(crossprod(input$K, S_beta))
  B <- crossprod(input$K, B_beta %*% input$K)
  list(metric = sqrt(sum(S * solve(-B, S))), eta = eta, S = S, B = B)
}

test_that("captured event106 cold roots satisfy independent profile scores", {
  input <- hia_fixture()$input
  for (reference in c(FALSE, TRUE)) local({
    if (reference) local_mocked_bindings(
      hl_fo_inner_quantities = hia_get("hl_fo_inner_quantities_r"), .package = "scRoPE"
    )
    detailed <- do.call(hia_get("hl_fo_solve_inner"), input)
    input$keep_diagnostics <- FALSE
    compact <- do.call(hia_get("hl_fo_solve_inner"), input)
    expect_true(detailed$success)
    expect_identical(detailed$accepted_tier, "derivative")
    expect_lte(detailed$iterations, 10L)
    expect_lte(detailed$inner_acceptance$score_metric, 1e-10)
    expect_lte(detailed$inner_acceptance$step_metric, 1e-10)
    expect_lte(hia_independent_score(input, detailed$beta)$metric, 1e-10)
    expect_identical(compact$beta, detailed$beta)
    expect_identical(compact$eta, detailed$eta)
    expect_identical(compact$mode_accuracy, detailed$mode_accuracy)
    expect_length(compact$trace, 0L)
    expect_true(any(vapply(detailed$trace, function(t) t$mode_polish_rounds > 0,
                           logical(1))))
  })
})

test_that("donor error cannot turn a nonroot into accepted beta convergence", {
  for (reference in c(FALSE, TRUE)) local({
    if (reference) local_mocked_bindings(
      hl_fo_inner_quantities = hia_get("hl_fo_inner_quantities_r"), .package = "scRoPE"
    )
    fixture <- hia_fixture(); input <- fixture$input
    false <- fixture$false_convergence
    raw <- hia_get("hl_fo_inner_quantities")(
      input$data, false$beta, input$rho_full, false$eta
    )
    geometry <- hia_get("hl_fo_inner_geometry")(raw, input$K, input$beta_star)
    expect_lt(geometry$score_metric, 1e-10)
    expect_true(all(vapply(raw$mode_diagnostics, function(d)
      d$residual_scaled <= 1e-10 && d$step_scaled <= 1e-10, logical(1))))
    accuracy <- hia_get("hl_fo_inner_mode_accuracy")(
      input$data, raw, input$rho_full, input$K, geometry
    )
    expect_gt(accuracy$score_error, 1e-10)
    expect_gt(hia_independent_score(input, false$beta)$metric, 2e-10)
    input$beta_start <- false$beta; input$eta_start <- false$eta; input$cold <- FALSE
    result <- do.call(hia_get("hl_fo_solve_inner"), input)
    expect_true(result$success)
    expect_gt(result$iterations, 1L)
    expect_lte(hia_independent_score(input, result$beta)$metric, 1e-10)
    expect_lte(result$inner_acceptance$score_metric, 1e-10)
    expect_lte(result$inner_acceptance$step_metric, 1e-10)
  })
})

test_that("donor polishing has a strict budget and refuses invalid states", {
  fixture <- hia_fixture(); x <- fixture$input
  state <- hia_get("hl_fo_inner_quantities")(
    x$data, fixture$false_convergence$beta, x$rho_full, fixture$false_convergence$eta
  )
  polish <- function(s, max_polish = 3L) hia_get("hl_fo_polish_inner_modes")(
    x$data, s, x$rho_full, x$K, x$beta_star, 1e-12, 1e-12, max_polish
  )
  none <- polish(state, 0L)
  expect_false(none$success)
  expect_identical(none$reason, "inner_mode_accuracy_exhausted")
  expect_identical(none$beta, state$beta)
  expect_identical(none$eta, state$eta)
  local({
    calls <- 0L
    local_mocked_bindings(hl_fo_inner_quantities = function(...) {
      calls <<- calls + 1L; state
    }, .package = "scRoPE")
    exhausted <- polish(state)
    expect_false(exhausted$success)
    expect_identical(exhausted$reason, "inner_mode_accuracy_exhausted")
    expect_identical(calls, 3L)
    expect_equal(exhausted$mode_polish_rounds, 3L)
  })
  local({
    original <- hia_get("hl_fo_inner_mode_accuracy")
    local_mocked_bindings(hl_fo_inner_mode_accuracy = function(...) {
      out <- original(...); out$correction[] <- 0; out
    }, .package = "scRoPE")
    stalled <- polish(state)
    expect_false(stalled$success)
    expect_identical(stalled$reason, "inner_mode_accuracy_stalled")
  })
  for (kind in c("missing", "nonfinite", "positive")) {
    invalid <- state
    invalid$mode_diagnostics[[1L]]$h_etaeta <- switch(kind,
      missing = NULL, nonfinite = NA_real_, positive = 1)
    result <- polish(invalid)
    expect_false(result$success, info = kind)
    expect_identical(result$reason, "inner_mode_accuracy_invalid", info = kind)
  }
  invalid <- state; invalid$B_beta <- diag(ncol(x$K))
  result <- polish(invalid)
  expect_false(result$success)
  expect_identical(result$status, "weak_curvature")
  expect_identical(result$reason, "inner_curvature_not_positive_definite")
  invalid <- state; invalid$S_beta[] <- 1e200
  result <- polish(invalid)
  expect_false(result$success)
  expect_identical(result$reason, "inner_nonfinite_newton_metric")
  invalid <- state; invalid$success <- FALSE; invalid$status <- "mode_failure"
  expect_identical(polish(invalid), invalid)
})

test_that("positive objective roundoff cannot override accurate Armijo evidence", {
  state <- function(P, S, error = 0) list(success = TRUE, P = P, S_beta = S,
    mode_accuracy = list(success = TRUE, component_error = error))
  decide <- hia_get("hl_fo_inner_armijo")
  K <- matrix(1, 1L, 1L); displacement <- 1e-8
  current <- state(1000, 1e-8)
  trial <- state(1000 + 2e-13, -5e-8)
  undecided <- decide(current, trial, K, displacement)
  expect_false(undecided$accepted)
  expect_true(undecided$needs_integral)
  negative <- decide(current, trial, K, displacement, state(1000, -2e-8))
  expect_lt(negative$integrated_gain, 0)
  expect_false(negative$accepted)
  positive <- decide(current, state(1000, 0), K, displacement, state(1000, 5e-9))
  expect_true(positive$accepted)
  expect_identical(positive$method, "profile_score_integral")
  uncertain <- decide(current, state(1000, 0, 2e-8), K, displacement,
                       state(1000, 5e-9, 2e-8))
  expect_false(uncertain$accepted)
  expect_false(decide(current, trial, K, 0, state(1000, 5e-9))$accepted)
  expect_false(decide(current, trial, K, -displacement)$accepted)
  invalid <- state(1000, 5e-9); invalid$success <- FALSE
  expect_false(decide(current, trial, K, displacement, invalid)$accepted)
  direct <- decide(current, state(1001, 0), K, displacement)
  expect_true(direct$accepted)
  expect_false(direct$needs_integral)
  expect_identical(direct$method, "direct_profile_criterion")
})

test_that("inner iteration limits and invalid roots are not relaxed by polishing", {
  input <- hia_fixture()$input
  limited <- do.call(hia_get("hl_fo_solve_inner"), c(input, list(max_iter = 1L)))
  expect_false(limited$success)
  expect_identical(limited$reason, "inner_iteration_exhausted")
  expect_length(limited$trace, 1L)
  input$data$y[] <- 0
  invalid <- do.call(hia_get("hl_fo_solve_inner"), input)
  expect_false(invalid$success)
  expect_identical(invalid$status, "inner_no_finite_root")
})

test_that("a genuine Armijo failure retains the finite backtracking budget", {
  input <- hia_fixture()$input
  evaluate <- hia_get("hl_fo_inner_quantities")
  initial <- do.call(evaluate, list(data = input$data, beta = input$beta_start,
                                    rho_full = input$rho_full))
  calls <- 0L
  local_mocked_bindings(hl_fo_inner_quantities = function(...) {
    calls <<- calls + 1L
    result <- evaluate(...)
    if (calls > 1L) result$P <- initial$P - 1
    result
  }, .package = "scRoPE")
  result <- do.call(hia_get("hl_fo_solve_inner"), input)
  expect_false(result$success)
  expect_identical(result$reason, "inner_armijo_failure")
  expect_identical(calls, 27L) # Initial state plus full step and 25 halvings.
  expect_length(result$trace, 1L)
})

test_that("an earlier line-search refusal retains an established weak-curvature diagnosis", {
  f <- hl_fo_fixture("W1")
  data <- hia_get("hl_fo_make_data")(f$y, f$X, f$log_pi, f$subject_index, f$intcol)
  result <- hia_get("hl_fo_fit")(
    data, active = character(0),
    fixed_natural = c(sigma2 = f$truth$sigma2, phi = f$truth$phi),
    seed = list(beta = f$truth$beta, rho_full = c(f$truth$tau, f$truth$xi),
                coordinate_tag = "canonical_conditional_beta_tau_xi")
  )
  expect_false(result$point_success)
  expect_identical(result$fit_status, "weak_curvature")
  if (identical(result$reason, "inner_armijo_failure")) {
    expect_lt(result$B_rcond, sqrt(.Machine$double.eps))
  }
})

test_that("polished roots retain transformed-design and representative fit accuracy", {
  input <- hia_fixture()$input
  original <- do.call(hia_get("hl_fo_solve_inner"), input)
  factors <- c(-2, 5, -3)
  input$data$X <- sweep(input$data$X, 2L, factors, `*`)
  input$beta_start <- input$beta_start / factors
  transformed <- do.call(hia_get("hl_fo_solve_inner"), input)
  expect_true(transformed$success)
  expect_equal(transformed$beta * factors, original$beta, tolerance = 1e-10)
  expect_lte(hia_independent_score(input, transformed$beta)$metric, 1e-10)
  for (id in c("R1", "B1")) {
    f <- hl_fo_fixture(id)
    d <- hia_get("hl_fo_make_data")(f$y, f$X, f$log_pi, f$subject_index, f$intcol)
    a <- list(data = d, rho_full = c(f$truth$tau, f$truth$xi),
              beta_star = numeric(d$p), K = diag(d$p), beta_start = f$truth$beta,
              cold = TRUE)
    result <- do.call(hia_get("hl_fo_solve_inner"), a)
    expect_true(result$success, info = id)
    expect_lte(result$inner_acceptance$score_metric, 1e-10, label = id)
    expect_lte(result$inner_acceptance$step_metric, 1e-10, label = id)
    expect_lte(hia_independent_score(a, result$beta)$metric, 1e-10, label = id)
  }
})
