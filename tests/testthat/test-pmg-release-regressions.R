pgmm_release_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    fn <- function(nm) getFromNamespace(nm, "scRoPE")
    n <- 16L
    id <- rep(seq_len(n), each = 4L)
    x <- 4 + 2.7 * rep(seq(-1.5, 1.5, length.out = n), each = 4L)
    z <- -1 + 0.8 * rep(c(-1.2, -.4, .4, 1.2), n)
    pred <- cbind(x = x, constant = 2, z = z)
    omega <- rep(qgamma((((seq_len(n) * 7) %% 17) + .5) / 18,
                        shape = 2, rate = 1.7), each = 4L)
    mu <- exp(1 + .23 * (x - 4) - .2 * (z + 1)) * omega
    count <- matrix(qpois((((seq_along(id) * 29) %% 67) + .5) / 68, mu), 1L)
    rownames(count) <- "fixed_toy"
    prep <- fn("hl_prepare_data")(count, id, pred = pred, verbose = FALSE,
                                   cpc = 0, mincp = 1)
    posv <- fn("call_posindy")(prep$count, 0L, prep$ctx$nind)
    fit <- fn("fit_gene_pmg_unconstrained")(1L, posv, prep$ctx)
    cached <<- list(count = count, id = id, pred = pred, prep = prep,
                    posv = posv, fit = fit)
    cached
  }
})

test_that("PGMM public contrasts use reported centered coefficients", {
  fx <- pgmm_release_fixture()
  fn <- function(nm) getFromNamespace(nm, "scRoPE")
  scales <- apply(fx$pred, 2L, function(x) sqrt(mean((x - mean(x))^2)))
  scales[2L] <- 1
  D <- diag(1 / scales)
  beta_reported <- as.numeric(fx$fit$theta$beta / scales)
  defs <- list(
    list(name = "nonzero_slope", L = matrix(c(1, 0, 0), 1L), b = .23),
    list(name = "mixed", L = matrix(c(1, .2, -.3), 1L),
         b = sum(c(1, .2, -.3) * beta_reported) + .1),
    list(name = "joint", L = rbind(c(1, 0, .3), c(0, 1, .2)),
         b = as.numeric(rbind(c(1, 0, .3), c(0, 1, .2)) %*% beta_reported) + c(.1, -.1)),
    list(name = "intercept", L = matrix(c(0, 1, 0), 1L), b = beta_reported[2L] + .1)
  )
  result <- suppressWarnings(fn("scrope_pmg_internal")(
    fx$count, fx$id, pred = fx$pred, cpc = 0, mincp = 1,
    verbose = FALSE, ncore = 1, additional_tests = c("score", "lrt"),
    lrt_details = TRUE, lrt_contrasts = defs
  ))
  raw <- result$diagnostics$raw[[1L]]
  expect_true(fn("pmg_optimizer_converged")(raw))
  # The established public intercept is at centered predictors, with constant
  # column normalized to one. It is not the original-design intercept beta/2.
  expect_equal(as.numeric(result$summary[1L, paste0("logFC_", colnames(fx$pred))]),
               as.numeric(raw$theta$beta / scales), tolerance = 1e-10)
  expect_equal(result$summary$logFC_constant, raw$theta$beta[2L], tolerance = 1e-10)
  expect_equal(as.numeric(result$summary[1L, paste0("se_robust_", colnames(fx$pred))]),
               as.numeric(sqrt(diag(raw$robust_cov)) / scales), tolerance = 1e-10)
  expect_equal(result$diagnostics$coefficient_coordinates$D, D, tolerance = 1e-12)

  for (def in defs) {
    # Independent row-wise coefficient substitution: reported beta_j = beta_j/s_j.
    L_internal <- sweep(def$L, 2L, scales, "/")
    expected <- fn("pmg_adjusted_contrasts")(
      fx$prep$ctx, 1L, fx$posv,
      list(list(L = L_internal, b = def$b, label = def$name)), unconstrained = raw
    )[[1L]]
    actual <- raw$lrt[[def$name]]
    expect_equal(as.numeric(L_internal %*% actual$constrained$theta$beta),
                 as.numeric(def$b), tolerance = 1e-8)
    expect_equal(unname(actual$contrast_reported), unname(def$L))
    expect_equal(unname(actual$contrast_internal), unname(L_internal), tolerance = 1e-12)
    expect_equal(actual$rhs, def$b)
    expect_equal(actual$stat, expected$adj$lr$wP_adjusted, tolerance = 1e-8)
    expect_equal(actual$score_stat, expected$adj$score$stat, tolerance = 1e-8)
    expect_equal(actual$wald_stat, expected$adj$wald$stat, tolerance = 1e-8)
    V_reported <- raw$robust_cov / outer(scales, scales)
    delta <- def$L %*% (raw$theta$beta / scales) - def$b
    V_contrast <- def$L %*% V_reported %*% t(def$L)
    expect_equal(actual$wald_stat, as.numeric(crossprod(delta, solve(V_contrast, delta))),
                 tolerance = 1e-8)
  }

  # Existing zero-null single-slope statistics are invariant to row rescaling.
  for (j in c(1L, 3L)) {
    L <- matrix(0, 1L, 3L)
    L[1L, j] <- 1
    legacy_null <- fn("pmg_adjusted_contrasts")(
      fx$prep$ctx, 1L, fx$posv, list(list(L = L, b = 0)), unconstrained = raw
    )[[1L]]
    actual <- raw$lrt[[colnames(fx$pred)[j]]]
    expect_equal(actual$stat, legacy_null$adj$lr$wP_adjusted, tolerance = 1e-8)
    expect_equal(actual$score_stat, legacy_null$adj$score$stat, tolerance = 1e-8)
    expect_equal(actual$wald_stat, legacy_null$adj$wald$stat, tolerance = 1e-8)
  }
})

test_that("PGMM LRT rejects an unrestricted endpoint even with positive raw LR", {
  fx <- pgmm_release_fixture()
  fn <- function(nm) getFromNamespace(nm, "scRoPE")
  bad <- fn("fit_gene_pmg_unconstrained")(
    1L, fx$posv, fx$prep$ctx, control = list(iter.max = 1L, eval.max = 2L)
  )
  expect_false(fn("pmg_optimizer_converged")(bad))
  result <- fn("pmg_adjusted_contrasts")(
    fx$prep$ctx, 1L, fx$posv,
    list(list(L = matrix(c(0, 0, 1), 1L), b = .4)),
    unconstrained = bad, retry_constrained = FALSE
  )[[1L]]
  expect_gt(result$adj$lr$wP, 0)
  expect_false(result$fit_pair_diagnostics$unconstrained_repair_applied)
  expect_false(result$adj$lr$success)
  expect_false(result$adj$lr$adjusted_available)
  expect_false(result$adj$lr$gate$unconstrained_accepted)
  expect_false(result$adj$lr$fallback)
  expect_identical(result$adj$lr$method, "lrt_fail")
  expect_identical(result$adj$lr$failure_reason, "nonconverged_unconstrained")
  expect_true(is.na(result$adj$lr$wP_adjusted))
  expect_true(is.na(result$p_values$lr))
  expect_false(result$adj$wald$success)
  expect_true(is.na(result$p_values$wald))
  # Score uses the accepted null endpoint, independently of the main fit.
  expect_true(result$adj$score$success)
})

test_that("PGMM Score requires an accepted interior constrained endpoint", {
  fx <- pgmm_release_fixture()
  fn <- function(nm) getFromNamespace(nm, "scRoPE")
  L <- matrix(c(1, 0, 0), 1L)
  bad <- fn("fit_gene_pmg_constrained")(
    1L, fx$posv, fx$prep$ctx, L, 0, fx$fit,
    control = list(iter.max = 1L, eval.max = 2L), retry = FALSE
  )
  assessment <- fn("pmg_assess_constrained_convergence")(
    bad$diagnostics, optimizer_converged = isTRUE(bad$convergence > 0L),
    n_subjects = fx$prep$ctx$k
  )
  expect_false(assessment$accepted)
  expect_true(is.finite(fn("compute_adj_score")(bad, L)$stat))
  evaluate_endpoint <- function(endpoint) {
    wrapped <- fn("pmg_adjusted_contrasts")
    env <- new.env(parent = environment(wrapped))
    environment(wrapped) <- env
    env$pmg_fit_constrained_pair <- function(...) list(
      unconstrained = fx$fit, constrained = endpoint,
      diagnostics = list(raw_lr_initial = 1, raw_lr_final = 1,
                         unconstrained_repair_attempted = FALSE,
                         unconstrained_repair_applied = FALSE)
    )
    wrapped(fx$prep$ctx, 1L, fx$posv, list(list(L = L, b = 0)),
            unconstrained = fx$fit, retry_constrained = FALSE)[[1L]]
  }
  rejected <- evaluate_endpoint(bad)
  expect_false(rejected$adj$score$success)
  expect_true(is.na(rejected$p_values$score))
  expect_identical(rejected$adj$score$failure_reason, "nonconverged_constrained")
  expect_false(rejected$adj$lr$success)

  good <- fn("fit_gene_pmg_constrained")(1L, fx$posv, fx$prep$ctx, L, 0, fx$fit)
  expect_true(evaluate_endpoint(good)$adj$score$success)
  # Reuse genuine finite derivatives; only the endpoint validity flag changes.
  boundary <- good
  boundary$lrt_gate$interior <- FALSE
  boundary$lrt_gate$nuisance_interior <- FALSE
  rejected_boundary <- evaluate_endpoint(boundary)
  expect_false(rejected_boundary$adj$score$success)
  expect_true(is.na(rejected_boundary$p_values$score))
  expect_identical(rejected_boundary$adj$score$failure_reason, "boundary_constrained")
  off_constraint <- good
  off_constraint$theta$beta[1L] <- .1
  rejected_constraint <- evaluate_endpoint(off_constraint)
  expect_false(rejected_constraint$adj$score$success)
  expect_true(is.na(rejected_constraint$p_values$score))
  expect_identical(rejected_constraint$adj$score$failure_reason, "constraint_not_satisfied")
})

test_that("PGMM internal profiles gate a rejected unrestricted reference", {
  fx <- pgmm_release_fixture()
  fn <- function(nm) getFromNamespace(nm, "scRoPE")
  bad <- fn("fit_gene_pmg_unconstrained")(
    1L, fx$posv, fx$prep$ctx, control = list(iter.max = 1L, eval.max = 2L)
  )
  point <- fn("pgmm_relative_profile_point")(
    fx$prep$ctx, 1L, fx$posv, psi_value = .4, psi_index = 3L,
    unconstrained = bad, retry_failed_profiles = FALSE
  )
  rel <- point$relative_profile
  expect_gt(rel$raw_lr, 0)
  expect_false(rel$unconstrained_accepted)
  expect_false(rel$robust_available)
  expect_identical(rel$failure_reason, "nonconverged_unconstrained")
  expect_true(is.na(rel$relLik))
  expect_true(is.na(rel$wP_star))
  expect_true(is.na(rel$relLik_star))
  grid <- fn("pgmm_relative_profile_grid")(
    fx$prep$ctx, 1L, fx$posv, psi_values = c(.4, .6), psi_index = 3L,
    unconstrained = bad, retry_failed_profiles = FALSE
  )
  expect_true(all(is.na(grid$relLik_std)))
  expect_true(all(is.na(grid$relLik_star_std)))
  expect_false(any(grid$robust_available))
})
