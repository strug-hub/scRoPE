# Fixed software fixtures, not a finite-sample calibration experiment.
.ln_optional_public_case <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    fx <- ln_joint_test_fixture()
    fit <- function(tests = character(0), keep = TRUE, contrasts = NULL) {
      scrope(fx$count, fx$id, pred = fx$pred, offset = fx$offset,
             verbose = FALSE, cpc = 0, mincp = 1, ncore = 1L,
             ln_estimator = "joint_ln", ln_operational_fallback = "none",
             additional_tests = tests, lrt_details = "lrt" %in% tests,
             lrt_contrasts = contrasts, keep_diagnostics = keep)
    }
    cached <<- list(
      fx = fx,
      baseline = fit(),
      score = fit("score"),
      lrt = fit("lrt"),
      both = fit(c("score", "lrt")),
      lean = fit(c("score", "lrt"), keep = FALSE),
      scaled = fit(c("score", "lrt"), contrasts = list(
        list(name = "scaled_x", L = matrix(c(-3.5, 0, 0), 1L), b = 0)
      ))
    )
    cached
  }
})

.ln_optional_assemble <- function(case, unrestricted = NULL,
                                  constrained = NULL,
                                  L = matrix(c(1, 0, 0), 1L), b = 0) {
  raw <- case$both$diagnostics$raw[[1L]]
  if (is.null(unrestricted)) unrestricted <- raw
  if (is.null(constrained)) constrained <- raw$lrt$x$constrained
  ln_joint_public_optional_tests(
    unrestricted, constrained, L, b, intercept_indices = 2L,
    coef_internal_to_public =
      case$both$provenance$coefficient_map_internal_to_public
  )
}

.ln_optional_expect_withheld <- function(result, reason) {
  for (test in c("score", "lrt")) {
    expect_false(result[[test]]$success, info = test)
    expect_true(is.na(result[[test]]$stat), info = test)
    expect_true(is.na(result[[test]]$pval), info = test)
    expect_identical(result[[test]]$failure_reason, reason, info = test)
  }
}

test_that("optional JLN tests use robust Score and Pace tails, not the raw ratio", {
  case <- .ln_optional_public_case()
  raw <- case$both$diagnostics$raw[[1L]]
  for (name in c("x", "z")) {
    entry <- raw$lrt[[name]]
    null <- entry$constrained
    score <- null$ln_joint$score_diagnostic
    ratio <- null$ln_joint$pace_ratio_diagnostic
    # Reconstruct the statistics from native components independently of the
    # public assembler and use the equivalent normal tail for chi-square(1).
    score_stat <- as.numeric(crossprod(score$U_eff,
                                      solve(score$J_eff, score$U_eff)))
    raw_ratio <- 2 * (raw$loglik_ln - null$loglik_ln)
    qP <- as.numeric(crossprod(ratio$U_eff,
                              solve(ratio$H_eff, ratio$U_eff)))
    wU <- as.numeric(crossprod(ratio$U_eff,
                              solve(ratio$J_eff, ratio$U_eff)))
    pace_stat <- raw_ratio * wU / qP
    expect_true(entry$score_success, info = name)
    expect_true(entry$success, info = name)
    expect_equal(entry$score_stat, score_stat, tolerance = 1e-12, info = name)
    expect_equal(entry$stat, pace_stat, tolerance = 1e-10, info = name)
    expect_equal(entry$score_pval, 2 * pnorm(-sqrt(score_stat)),
                 tolerance = 1e-14, info = name)
    expect_equal(entry$pval, 2 * pnorm(-sqrt(pace_stat)),
                 tolerance = 1e-12, info = name)
    expect_gt(abs(entry$stat - raw_ratio), 1e-6)
    expect_identical(entry$score_test$method, "joint_ln_robust_score")
    expect_identical(entry$method, "joint_ln_pace_adjusted_ratio")
    expect_identical(case$both$summary[[paste0("score_status_", name)]],
                     "experimental_asymptotic")
    expect_identical(case$both$summary[[paste0("lrt_status_", name)]],
                     "experimental_asymptotic")
    expect_false(score$public_available)
    expect_false(ratio$public_available)
    expect_identical(score$public_status, "diagnostic_only")
    expect_identical(ratio$public_status, "diagnostic_only")
    expect_false(any(c("pval", "p_value") %in% names(score)))
    expect_false(any(c("pval", "p_value") %in% names(ratio)))
  }
})

test_that("existing additional_tests is the only opt-in and preserves Wald", {
  case <- .ln_optional_public_case()
  expect_false(any(grepl("^(score_|p_score_|lrt_|p_lrt_)",
                         names(case$baseline$summary))))
  expect_null(case$baseline$lrt_details)
  expect_null(case$baseline$diagnostics$raw[[1L]]$lrt)
  expect_false(any(grepl("^lrt_|^p_lrt_", names(case$score$summary))))
  expect_false(any(grepl("^score_|^p_score_", names(case$lrt$summary))))
  primary <- setdiff(names(case$baseline$summary),
                     c("test_status", "test_failure_detail"))
  for (mode in c("score", "lrt", "both", "lean", "scaled")) {
    fit <- case[[mode]]
    expect_equal(fit$summary[primary], case$baseline$summary[primary],
                 tolerance = 1e-12, info = mode)
    expect_equal(fit$overdispersion, case$baseline$overdispersion,
                 tolerance = 1e-12, info = mode)
    expect_equal(fit$convergence, case$baseline$convergence,
                 tolerance = 0, info = mode)
    expect_equal(fit$algorithm, case$baseline$algorithm,
                 tolerance = 0, info = mode)
  }
  for (name in c("x", "z")) {
    for (prefix in c("score_adj_", "p_score_adj_", "score_success_")) {
      expect_equal(case$score$summary[[paste0(prefix, name)]],
                   case$both$summary[[paste0(prefix, name)]], tolerance = 0)
    }
    for (prefix in c("lrt_adj_", "p_lrt_adj_", "lrt_success_")) {
      expect_equal(case$lrt$summary[[paste0(prefix, name)]],
                   case$both$summary[[paste0(prefix, name)]], tolerance = 0)
    }
  }
})

test_that("public optional status and outputs survive dropping diagnostics", {
  case <- .ln_optional_public_case()
  expect_null(case$lean$diagnostics)
  expect_equal(case$lean$summary, case$both$summary, tolerance = 0)
  expect_equal(case$lean$lrt_details, case$both$lrt_details, tolerance = 0)
  expect_true(is.data.frame(case$lean$test_status))
  expect_equal(case$lean$test_status, case$both$test_status, tolerance = 0)
  optional <- subset(case$lean$test_status, test %in% c("score", "lrt"))
  expect_equal(nrow(optional), 4L)
  expect_true(all(optional$success))
  expect_true(all(optional$status == "experimental_asymptotic"))
  expect_setequal(unique(optional$method),
                  c("joint_ln_robust_score", "joint_ln_pace_adjusted_ratio"))
})

test_that("only exact scalar zero-null nonintercept selectors are public", {
  case <- .ln_optional_public_case()
  expect_true(.ln_optional_assemble(case)$score$success)
  expect_true(.ln_optional_assemble(case)$lrt$success)
  for (hypothesis in list(
    list(L = matrix(c(1, 0, 0), 1L), b = 0.1),
    list(L = matrix(c(1, 0, 0), 1L), b = .Machine$double.eps),
    list(L = matrix(c(0, 1, 0), 1L), b = 0),
    list(L = matrix(c(1, 0, -0.5), 1L), b = 0),
    list(L = matrix(c(1, 0, .Machine$double.eps), 1L), b = 0),
    list(L = rbind(c(1, 0, 0), c(0, 0, 1)), b = c(0, 0)),
    list(L = rbind(c(1, 0, 0), c(2, 0, 0)), b = c(0, 0)),
    list(L = matrix(c(0, 0, 0), 1L), b = 0),
    list(L = matrix(c(NA_real_, 0, 0), 1L), b = 0)
  )) {
    result <- .ln_optional_assemble(case, L = hypothesis$L, b = hypothesis$b)
    .ln_optional_expect_withheld(result, "joint_ln_optional_scope_unsupported")
  }
  scaled <- case$scaled$diagnostics$raw[[1L]]$lrt$scaled_x
  reference <- case$both$diagnostics$raw[[1L]]$lrt$x
  expect_true(scaled$success)
  expect_true(scaled$score_success)
  expect_equal(scaled$stat, reference$stat, tolerance = 1e-7)
  expect_equal(scaled$pval, reference$pval, tolerance = 1e-7)
  expect_equal(scaled$score_stat, reference$score_stat, tolerance = 1e-7)
  expect_equal(scaled$score_pval, reference$score_pval, tolerance = 1e-7)
})

test_that("finite stale diagnostics cannot bypass primary or null acceptance", {
  case <- .ln_optional_public_case()
  raw <- case$both$diagnostics$raw[[1L]]
  null <- raw$lrt$x$constrained
  expect_true(is.finite(null$ln_joint$score_diagnostic$statistic))
  expect_true(is.finite(null$ln_joint$pace_ratio_diagnostic$statistic))
  for (endpoint in c("unrestricted", "constrained")) {
    for (field in c("route", "inference_available")) {
      failed <- if (endpoint == "unrestricted") raw else null
      failed[[field]] <- if (field == "route") "R2" else FALSE
      failed$failure_reason <- "competing_uncertified_candidate"
      args <- list(case = case)
      args[[endpoint]] <- failed
      .ln_optional_expect_withheld(do.call(.ln_optional_assemble, args),
                                   "competing_uncertified_candidate")
    }
  }
})

test_that("Pace curvature and Score variability rejections remain test-specific", {
  case <- .ln_optional_public_case()
  null <- case$both$diagnostics$raw[[1L]]$lrt$x$constrained
  failed_pace <- null
  failed_pace$ln_joint$pace_ratio_diagnostic$diagnostic_available <- FALSE
  failed_pace$ln_joint$pace_ratio_diagnostic$failure_reason <-
    "H_eff_not_positive_definite"
  failed_pace$ln_joint$pace_ratio_diagnostic$diagnostics$H_eff$positive_definite <-
    FALSE
  # Preserve finite components deliberately: availability, never finiteness,
  # is the contract supplied by the independently tested native curvature gate.
  result <- .ln_optional_assemble(case, constrained = failed_pace)
  expect_true(result$score$success)
  expect_true(is.finite(result$score$pval))
  expect_false(result$lrt$success)
  expect_true(is.na(result$lrt$stat))
  expect_true(is.na(result$lrt$pval))
  expect_identical(result$lrt$failure_reason, "H_eff_not_positive_definite")
  expect_false(result$lrt$gate$profiled_curvature_positive)
  missing_curvature <- failed_pace
  missing_curvature$ln_joint$pace_ratio_diagnostic$H_eff <- NULL
  missing_curvature$ln_joint$pace_ratio_diagnostic$diagnostics$H_eff <- NULL
  missing_result <- .ln_optional_assemble(case, constrained = missing_curvature)
  expect_true(is.na(missing_result$lrt$gate$profiled_curvature_positive))
  expect_true(missing_result$score$success)
  expect_false(missing_result$lrt$success)

  failed_score <- null
  failed_score$ln_joint$score_diagnostic$diagnostic_available <- FALSE
  failed_score$ln_joint$score_diagnostic$failure_reason <-
    "J_eff_not_positive_definite"
  failed_score$ln_joint$pace_ratio_diagnostic$diagnostic_available <- FALSE
  failed_score$ln_joint$pace_ratio_diagnostic$failure_reason <-
    "constrained_score_J_eff_not_positive_definite"
  result <- .ln_optional_assemble(case, constrained = failed_score)
  expect_false(result$score$success)
  expect_true(is.na(result$score$stat))
  expect_true(is.na(result$score$pval))
  expect_identical(result$score$failure_reason, "J_eff_not_positive_definite")
  expect_false(result$lrt$success)
  expect_true(is.na(result$lrt$pval))
  expect_identical(result$lrt$failure_reason,
                   "constrained_score_J_eff_not_positive_definite")
})

test_that("objective, public-coordinate and reference metadata cannot be bypassed", {
  case <- .ln_optional_public_case()
  null <- case$both$diagnostics$raw[[1L]]$lrt$x$constrained
  changed <- null
  changed$ln_joint$objective_instance <- list(unrelated = TRUE)
  .ln_optional_expect_withheld(
    .ln_optional_assemble(case, constrained = changed),
    "endpoint_objective_instance_mismatch"
  )
  changed <- null
  changed$ln_joint$objective_provenance <- list(unrelated = TRUE)
  .ln_optional_expect_withheld(
    .ln_optional_assemble(case, constrained = changed),
    "endpoint_criterion_mismatch"
  )
  changed <- null
  changed$constraint$L[1L, 1L] <- changed$constraint$L[1L, 1L] + 0.1
  .ln_optional_expect_withheld(
    .ln_optional_assemble(case, constrained = changed),
    "joint_ln_public_constraint_mismatch"
  )
  for (field in c("df", "theoretical_df", "no_numerical_repair")) {
    changed <- null
    changed$ln_joint$score_diagnostic[[field]] <-
      if (field == "no_numerical_repair") FALSE else 2L
    result <- .ln_optional_assemble(case, constrained = changed)
    expect_false(result$score$success, info = field)
    expect_true(is.na(result$score$pval), info = field)
    expect_false(result$lrt$success, info = field)
    expect_true(is.na(result$lrt$pval), info = field)
  }
})

test_that("a rejected primary fit has no public optional test tail", {
  fit <- suppressWarnings(ln_joint_fit_fixture(
    additional_tests = c("score", "lrt"), max = c(0.05, 1000)
  )$fit)
  raw <- fit$diagnostics$raw[[1L]]
  expect_identical(raw$route, "R2")
  expect_false(raw$inference_available)
  details <- fit$lrt_details
  expect_gt(nrow(details), 0L)
  expect_true(all(!details$success))
  expect_true(all(!details$score_success))
  expect_true(all(is.na(details$stat)))
  expect_true(all(is.na(details$p_value)))
  expect_true(all(is.na(details$score_stat)))
  expect_true(all(is.na(details$score_p_value)))
  expect_true(all(details$failure_reason == raw$failure_reason))
  expect_true(all(details$score_failure_reason == raw$failure_reason))
  optional <- subset(fit$test_status, test %in% c("score", "lrt"))
  expect_true(all(optional$failure_stage == "primary"))
  expect_true(all(optional$status == "unavailable"))
})
