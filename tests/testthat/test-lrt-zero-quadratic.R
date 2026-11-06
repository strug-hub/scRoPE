# Diagnostic fixtures isolate the shared adjustment from optimizer tolerances.
# H and donor meat are full rank, and donor scores sum to the supplied score.
lrt_zero_fixture <- function(direction = c(0, 0), raw_lr = 0) {
  H <- diag(c(2, 3, 4))
  score <- c(c(2, 3) * direction, 0)
  D <- diag(sqrt(c(5, 11, 7) / 2))
  U <- cbind(D, -D, score)
  fit <- list(theta = list(beta = c(0, 0), lambda = 0), loglik = 0,
    score = score, observed_info = H, per_subject_gradients = U,
    robust_cov = diag(c(5/4, 11/9)), convergence = 1L,
    optimizer = "nlminb", raw_fit = list(convergence = 0L), diagnostics = list(),
    lrt_gate = list(converged = TRUE, interior = TRUE, nuisance_interior = TRUE))
  fit$godambe <- getFromNamespace("extract_godambe_components", "scRoPE")(
    list(observed_info = H, per_subject_gradients = U, gradient = score), 2L)
  unrestricted <- fit
  unrestricted$loglik <- raw_lr/2
  list(u = unrestricted, c = fit, L = diag(2), b = c(0, 0))
}

lrt_zero_call <- function(fx, ...) {
  getFromNamespace("compute_profile_lrt", "scRoPE")(
    fx$u, fx$c, fx$L, fx$b, ...)
}

test_that("joint exact zero extends the statistic without inventing a multiplier", {
  fx <- lrt_zero_fixture()
  expect_identical(rowSums(fx$c$per_subject_gradients), fx$c$score)
  out <- lrt_zero_call(fx)
  expect_true(out$success)
  expect_true(out$adjusted_available)
  expect_identical(out$method, "lrt")
  expect_identical(out$wP_adjusted, 0)
  expect_identical(out$qP, 0)
  expect_true(is.na(out$scale))
  expect_false(out$fallback)
  expect_false(out$gate$scale_positive)
  expect_true(out$gate$adjustment_defined)
  expect_identical(out$gate$small_quadratic$mode, "coincident_zero_extension")
  expect_true(is.na(out$failure_reason))
  expect_equal(pchisq(out$wP_adjusted, 2, lower.tail = FALSE), 1)
})

test_that("tiny joint directions retain their Rayleigh ratio including underflow", {
  for (magnitude in c(1e-9, 1e-50, 1e-200)) {
    for (direction in list(c(1, 0), c(0, 1), c(1, -2))) {
      fx <- lrt_zero_fixture(magnitude * direction, raw_lr = 2)
      out <- lrt_zero_call(fx)
      # Closed-form diagonal matrices at this precision: A=(1/2,1/3),
      # B=(5/4,11/9). The scale depends on direction, never its magnitude.
      expected <- sum(direction^2 * c(4/5, 9/11)) /
        sum(direction^2 * c(2, 3))
      expect_true(out$success)
      expect_identical(out$gate$small_quadratic$mode, "normalized_direction")
      expect_equal(out$scale, expected, tolerance = 2e-14)
      expect_equal(out$wP_adjusted, 2 * expected, tolerance = 2e-14)
      expect_false(out$fallback)
      if (magnitude == 1e-200) expect_identical(out$qP, 0)
    }
  }
})

test_that("joint adjustment is invariant to well-conditioned contrast basis changes", {
  for (direction in list(c(0, 0), c(1e-100, -2e-100))) {
    fx <- lrt_zero_fixture(direction, raw_lr = if (all(direction == 0)) 0 else 2)
    expected <- lrt_zero_call(fx)
    for (C in list(matrix(c(1, 0.4, -0.7, 2), 2),
                   diag(c(1e-4, 1e4)), 1e-100 * diag(2))) {
      transformed <- fx
      transformed$L <- C %*% fx$L
      transformed$b <- as.numeric(C %*% fx$b)
      out <- lrt_zero_call(transformed)
      expect_true(out$success)
      expect_equal(out$wP_adjusted, expected$wP_adjusted, tolerance = 1e-12)
      expect_equal(out$scale, expected$scale, tolerance = 1e-12)
    }
  }
})

test_that("undefined zero direction cannot become a finite fallback p-value", {
  cases <- list()
  cases$positive_lr <- lrt_zero_fixture(raw_lr = 1)
  cases$negative_clamp <- lrt_zero_fixture(raw_lr = -1e-8)
  cases$negative_lr <- lrt_zero_fixture(raw_lr = -1)
  cases$missing_lr <- lrt_zero_fixture(raw_lr = NA_real_)
  cases$beta_mismatch <- lrt_zero_fixture()
  cases$beta_mismatch$u$theta$beta[1] <- 1e-12
  cases$nuisance_mismatch <- lrt_zero_fixture()
  cases$nuisance_mismatch$u$theta$lambda <- 1e-12
  cases$missing_nuisance <- lrt_zero_fixture()
  cases$missing_nuisance$u$theta$lambda <- NULL
  cases$both_missing_nuisance <- cases$missing_nuisance
  cases$both_missing_nuisance$c$theta$lambda <- NULL
  cases$renamed_parameter <- lrt_zero_fixture()
  names(cases$renamed_parameter$u$theta$beta) <- c("a", "b")
  cases$objective_state_mismatch <- lrt_zero_fixture()
  cases$objective_state_mismatch$u$loglik_ln <- 0
  cases$objective_state_mismatch$c$loglik_ln <- 0
  cases$objective_state_mismatch$u$theta_ln <- list(beta = c(0, 0), lambda = 0)
  cases$objective_state_mismatch$c$theta_ln <- list(beta = c(1, 0), lambda = 0)
  cases$missing_objective_state <- cases$objective_state_mismatch
  cases$missing_objective_state$c$theta_ln <- NULL
  cases$objective_state_differs_from_final <- cases$objective_state_mismatch
  cases$objective_state_differs_from_final$u$theta_ln <-
    cases$objective_state_differs_from_final$c$theta_ln
  cases$nonfinite_parameter <- lrt_zero_fixture()
  cases$nonfinite_parameter$u$theta$lambda <- Inf
  cases$misses_null <- lrt_zero_fixture()
  cases$misses_null$b[1] <- 1e-12
  for (field in c("converged", "interior", "nuisance_interior")) {
    for (endpoint in c("u", "c")) {
      fx <- lrt_zero_fixture()
      fx[[endpoint]]$lrt_gate[[field]] <- FALSE
      cases[[paste(endpoint, field)]] <- fx
    }
  }
  for (name in names(cases)) {
    out <- lrt_zero_call(cases[[name]])
    expect_false(out$success, info = name)
    expect_false(out$adjusted_available, info = name)
    expect_identical(out$method, "lrt_fail", info = name)
    expect_true(is.na(out$wP_adjusted), info = name)
    expect_false(out$fallback, info = name)
    expect_false(is.na(out$failure_reason), info = name)
  }
  out <- lrt_zero_call(cases$negative_clamp)
  expect_identical(out$wP, 0)
  expect_identical(out$gate$small_quadratic$raw_lr, -1e-8)
  expect_identical(out$gate$small_quadratic$reason, "invalid_unclamped_lr")
  override <- list(constraint_satisfied = FALSE)
  expect_false(lrt_zero_call(lrt_zero_fixture(), lrt_gate_override = override)$success)
})

test_that("coincident objective-specific states and tiny positive LR are explicit", {
  fx <- lrt_zero_fixture()
  fx$u$theta_ln <- fx$c$theta_ln <- fx$c$theta
  fx$u$loglik_ln <- fx$c$loglik_ln <- 0
  expect_true(lrt_zero_call(fx)$success)
  for (size in c(1e-9, 1e-100)) {
    fx <- lrt_zero_fixture(c(size, 0), raw_lr = 2 * size^2)
    out <- lrt_zero_call(fx)
    expect_true(out$success)
    expect_equal(out$wP_adjusted / (2 * size^2), 2/5, tolerance = 1e-14)
    expect_identical(out$gate$small_quadratic$mode, "normalized_direction")
  }
})

test_that("an untested indefinite beta direction cannot qualify a zero extension", {
  H <- diag(c(2, 3, -4, 5)); U <- cbind(diag(4), -diag(4))
  fit <- lrt_zero_fixture()$c
  fit$theta$beta <- rep(0, 3)
  fit$score <- rep(0, 4)
  fit$observed_info <- H
  fit$per_subject_gradients <- U
  fit$robust_cov <- diag(3)
  fit$godambe <- getFromNamespace("extract_godambe_components", "scRoPE")(
    list(observed_info = H, per_subject_gradients = U, gradient = fit$score), 3L)
  fx <- list(u = fit, c = fit, L = cbind(diag(2), 0), b = c(0, 0))
  out <- lrt_zero_call(fx)
  expect_false(out$success)
  expect_true(is.na(out$wP_adjusted))
  expect_identical(out$gate$small_quadratic$reason, "nonpositive_profile_curvature")
})

test_that("new joint handling requires raw full rank and positive curvature", {
  for (direction in list(c(0, 0), c(1e-100, -2e-100))) {
    cases <- list()
    fx <- lrt_zero_fixture(direction)
    cases$singular_meat <- fx
    cases$singular_meat$c$godambe$S[1:2, 1:2] <- matrix(1, 2, 2)
    cases$indefinite_meat <- fx
    cases$indefinite_meat$c$godambe$S[1:2, 1:2] <- matrix(c(1, 2, 2, 1), 2)
    cases$missing_meat <- fx
    cases$missing_meat$c$godambe$S[1, 1] <- NA_real_
    cases$indefinite_nuisance <- fx
    cases$indefinite_nuisance$c$godambe$H_lambda[,] <- -1
    cases$indefinite_profile <- fx
    cases$indefinite_profile$c$godambe$H_beta[1, 1] <- -1
    cases$singular_profile <- fx
    cases$singular_profile$c$godambe$H_beta[,] <- 1
    cases$redundant_contrast <- fx
    cases$redundant_contrast$L[,] <- 1
    for (name in names(cases)) {
      # Historical inversion diagnostics can warn before the new rank gate.
      out <- suppressWarnings(lrt_zero_call(cases[[name]]))
      expect_false(out$success, info = name)
      expect_false(out$adjusted_available, info = name)
      expect_true(is.na(out$wP_adjusted), info = name)
      expect_false(out$fallback, info = name)
    }
  }
})

test_that("scalar zero and positive joint diagnostic outcomes retain their formula", {
  for (direction in list(c(0, 0), c(0.1, -0.2))) {
    fx <- lrt_zero_fixture(direction, raw_lr = 2)
    fx$L <- matrix(c(1, 0), 1)
    fx$b <- 0
    scalar <- lrt_zero_call(fx)
    expect_true(scalar$success)
    expect_equal(scalar$scale, as.numeric(scalar$Hpsi / scalar$Jpsi), tolerance = 1e-12)
    expect_equal(scalar$wP_adjusted, 2 * scalar$scale, tolerance = 1e-12)
    expect_null(scalar$gate$small_quadratic)
    fx$c$observed_info <- NULL # general, one-contrast compatibility path
    legacy_scalar <- lrt_zero_call(fx)
    expect_true(legacy_scalar$success)
    expect_equal(legacy_scalar$scale, scalar$scale, tolerance = 1e-12)
    expect_null(legacy_scalar$gate$small_quadratic)
  }
  fx <- lrt_zero_fixture(c(0.1, -0.2), raw_lr = 2)
  out <- lrt_zero_call(fx)
  A <- diag(c(1/2, 1/3))
  B <- A %*% fx$c$godambe$S[1:2, 1:2] %*% A
  u <- c(0.1, -0.2)
  expected <- drop(t(u) %*% solve(B, u)) / drop(t(u) %*% solve(A, u))
  expect_true(out$success)
  expect_equal(out$scale, expected, tolerance = 1e-12)
  expect_null(out$gate$small_quadratic)
})

test_that("PGMM adapter exposes zero p=1 and withholds undefined joint LRT", {
  f <- getFromNamespace("pmg_adjusted_contrasts", "scRoPE")
  env <- new.env(parent = environment(f))
  environment(f) <- env
  # Isolate the actual adjustment and p-value adapter from fitting. The shared
  # compute_profile_lrt is deliberately not mocked.
  env$pmg_assess_constrained_convergence <- function(...) list(accepted = TRUE)
  env$compute_adj_score <- function(...) list(stat = NA_real_, success = FALSE)
  for (raw_lr in c(0, 1, -1e-8)) {
    fx <- lrt_zero_fixture(raw_lr = raw_lr)
    env$pmg_fit_constrained_pair <- function(...) list(unconstrained = fx$u,
      constrained = fx$c, diagnostics = list(raw_lr_final = raw_lr))
    out <- f(list(k = 7), 1L, NULL,
      list(list(L = fx$L, b = fx$b, label = "joint")),
      unconstrained = fx$u, retry_constrained = FALSE)[[1]]
    expect_identical(out$adj$lr$success, raw_lr == 0)
    if (raw_lr == 0) expect_identical(out$p_values$lr, 1)
    else expect_true(is.na(out$p_values$lr))
    expect_true(out$adj$wald$success) # separately qualified Wald remains usable
  }
})
