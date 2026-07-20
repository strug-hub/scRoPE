.hl_fo_cpp_functions <- function() {
  list(
    make_data = getFromNamespace("hl_fo_make_data", "scRoPE"),
    inner = getFromNamespace("hl_fo_inner_quantities_r", "scRoPE"),
    inner_dispatch = getFromNamespace("hl_fo_inner_quantities", "scRoPE"),
    cpp = getFromNamespace("hl_fo_inner_eval_cpp", "scRoPE")
  )
}

.hl_fo_cpp_sdisc <- function(a, b) {
  expect_identical(dim(a), dim(b))
  expect_identical(length(a), length(b))
  if (!length(a)) {
    return(0)
  }
  max(abs(as.numeric(a) - as.numeric(b)) /
    (1 + pmax(abs(as.numeric(a)), abs(as.numeric(b)))))
}

.hl_fo_cpp_compare_success <- function(reference, compiled, tolerance = 1e-10) {
  expect_true(reference$success)
  expect_true(compiled$success)
  expect_identical(names(compiled), names(reference))
  expect_identical(compiled$status, reference$status)
  expect_identical(compiled$beta, reference$beta)
  expect_lte(.hl_fo_cpp_sdisc(compiled$eta, reference$eta), tolerance)
  expect_lte(.hl_fo_cpp_sdisc(compiled$P, reference$P), tolerance)
  expect_lte(.hl_fo_cpp_sdisc(compiled$S_beta, reference$S_beta), tolerance)
  expect_lte(.hl_fo_cpp_sdisc(compiled$B_beta, reference$B_beta), tolerance)
  expect_identical(dimnames(compiled$B_beta), dimnames(reference$B_beta))

  expect_identical(length(compiled$mode_diagnostics),
    length(reference$mode_diagnostics))
  for (ii in seq_along(reference$mode_diagnostics)) {
    rdiag <- reference$mode_diagnostics[[ii]]
    cdiag <- compiled$mode_diagnostics[[ii]]
    expect_identical(names(cdiag), names(rdiag))
    expect_identical(cdiag$success, rdiag$success)
    expect_identical(cdiag$status, rdiag$status)
    expect_identical(cdiag$reason, rdiag$reason)
    expect_identical(cdiag$iterations, rdiag$iterations)
    expect_identical(cdiag$bracket_expansions, rdiag$bracket_expansions)
    expect_identical(cdiag$bisections, rdiag$bisections)
    expect_identical(cdiag$start_type, rdiag$start_type)
    numeric_fields <- c(
      "eta", "h_eta", "h_etaeta", "residual_scaled", "step_scaled"
    )
    expect_lte(
      .hl_fo_cpp_sdisc(
        unlist(cdiag[numeric_fields], use.names = FALSE),
        unlist(rdiag[numeric_fields], use.names = FALSE)
      ),
      tolerance
    )
    expect_lte(.hl_fo_cpp_sdisc(cdiag$bracket, rdiag$bracket), tolerance)
  }

  expect_identical(length(compiled$subject), length(reference$subject))
  for (ii in seq_along(reference$subject)) {
    rsubject <- reference$subject[[ii]]
    csubject <- compiled$subject[[ii]]
    expect_identical(names(csubject), names(rsubject))
    expect_lte(.hl_fo_cpp_sdisc(csubject$p, rsubject$p), tolerance)
    expect_lte(
      .hl_fo_cpp_sdisc(csubject$h_beta, rsubject$h_beta),
      tolerance
    )
    expect_identical(csubject$p_beta_beta, rsubject$p_beta_beta)
    expect_lte(.hl_fo_cpp_sdisc(csubject$a, rsubject$a), tolerance)
  }
}

.hl_fo_cpp_evaluate <- function(case, eta_start = NULL, use_warm = FALSE) {
  fn <- .hl_fo_cpp_functions()
  data <- fn$make_data(
    case$y,
    case$X,
    case$log_offset,
    case$subject
  )
  reference <- fn$inner(
    data,
    case$beta,
    case$rho,
    eta_start = if (isTRUE(use_warm)) eta_start else NULL
  )
  compiled <- fn$cpp(
    y = case$y,
    X = case$X,
    log_offset = case$log_offset,
    subject_zero = as.integer(data$subject - 1L),
    beta = case$beta,
    tau = case$rho[[1L]],
    xi = case$rho[[2L]],
    eta_start = if (is.null(eta_start)) numeric(0) else eta_start,
    use_warm = use_warm
  )
  list(reference = reference, compiled = compiled, data = data)
}

.hl_fo_cpp_regular_case <- function() {
  list(
    y = c(0, 1, 3, 2, 0, 4, 1, 1, 0, 5, 2, 0),
    X = cbind(
      `(Intercept)` = 1,
      z = c(-1, 0, 1, -0.5, 0.5, 1.5, -1.5, -0.25, 0.75, -0.8, 0.2, 1.2)
    ),
    log_offset = log(c(
      0.8, 1, 1.2, 0.9, 1.1, 1.3, 0.7, 0.95, 1.4, 0.85, 1.05, 1.25
    )),
    subject = rep(paste0("s", 1:4), each = 3),
    beta = c(`(Intercept)` = -0.7, z = 0.25),
    rho = c(tau = log(0.4), xi = log(8))
  )
}

.hl_fo_cpp_zero_heavy_case <- function() {
  n_by_subject <- c(2L, 4L, 3L, 5L, 2L)
  n <- sum(n_by_subject)
  list(
    y = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1),
    X = cbind(
      `(Intercept)` = 1,
      z = seq(-1.4, 1.2, length.out = n),
      w = rep(c(-0.7, 0.3), length.out = n)
    ),
    log_offset = seq(-0.35, 0.3, length.out = n),
    subject = rep(paste0("s", seq_along(n_by_subject)), n_by_subject),
    beta = c(`(Intercept)` = -3.8, z = 0.2, w = -0.1),
    rho = c(tau = log(0.8), xi = log(2))
  )
}

test_that("compiled first-order inner wrapper has the fixed transaction signature", {
  cpp <- .hl_fo_cpp_functions()$cpp
  expect_identical(
    names(formals(cpp)),
    c(
      "y", "X", "log_offset", "subject_zero", "beta", "tau", "xi",
      "eta_start", "use_warm"
    )
  )
})

test_that("installed first-order inner dispatch uses the compiled transaction", {
  case <- .hl_fo_cpp_regular_case()
  fn <- .hl_fo_cpp_functions()
  data <- fn$make_data(case$y, case$X, case$log_offset, case$subject)
  dispatched <- fn$inner_dispatch(data, case$beta, case$rho)
  compiled <- fn$cpp(
    case$y, case$X, case$log_offset, as.integer(data$subject - 1L),
    case$beta, case$rho[1L], case$rho[2L], numeric(0), FALSE
  )
  .hl_fo_cpp_compare_success(dispatched, compiled, tolerance = 0)
})

test_that("compiled inner transaction matches the regular R reference", {
  result <- .hl_fo_cpp_evaluate(.hl_fo_cpp_regular_case())
  .hl_fo_cpp_compare_success(result$reference, result$compiled)
})

test_that("compiled inner transaction matches the zero-heavy R reference", {
  result <- .hl_fo_cpp_evaluate(.hl_fo_cpp_zero_heavy_case())
  .hl_fo_cpp_compare_success(result$reference, result$compiled)
})

test_that("compiled mode path matches difficult dispersion and warm-start cases", {
  case <- .hl_fo_cpp_zero_heavy_case()
  dispersion_cases <- list(
    low_variance_high_size = c(tau = log(0.02), xi = log(100)),
    high_variance_low_size = c(tau = log(2), xi = log(0.3))
  )
  for (rho in dispersion_cases) {
    case$rho <- rho
    result <- .hl_fo_cpp_evaluate(case)
    .hl_fo_cpp_compare_success(result$reference, result$compiled)
  }

  regular <- .hl_fo_cpp_regular_case()
  cold <- .hl_fo_cpp_evaluate(regular)
  eta_start <- cold$reference$eta + c(3, -2, 0.7, -4)
  warm <- .hl_fo_cpp_evaluate(regular, eta_start = eta_start, use_warm = TRUE)
  .hl_fo_cpp_compare_success(warm$reference, warm$compiled)
  expect_true(any(vapply(
    warm$compiled$mode_diagnostics,
    `[[`,
    integer(1),
    "bisections"
  ) > 0L))
  expect_true(all(vapply(
    warm$compiled$mode_diagnostics,
    `[[`,
    character(1),
    "start_type"
  ) == "warm"))
})

test_that("compiled mode failure status matches the R reference exactly", {
  case <- .hl_fo_cpp_regular_case()
  case$rho <- c(tau = 1000, xi = log(8))
  result <- suppressWarnings(.hl_fo_cpp_evaluate(case))

  expect_false(result$reference$success)
  expect_false(result$compiled$success)
  expect_identical(names(result$compiled), names(result$reference))
  expect_identical(result$compiled$status, result$reference$status)
  expect_identical(result$compiled$reason, result$reference$reason)
  expect_identical(
    result$compiled$failed_subject,
    result$reference$failed_subject
  )
  expect_identical(
    result$compiled$diagnostics$unresolved_subjects,
    result$reference$diagnostics$unresolved_subjects
  )
})
