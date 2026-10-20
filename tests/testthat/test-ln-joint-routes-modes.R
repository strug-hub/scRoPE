test_that("joint LN R2 removes every analytic covariance and Wald field", {
  expect_warning(
    result <- ln_joint_fit_fixture(
      max = c(0.01, 1000),
      output_re = TRUE,
      additional_tests = c("score", "lrt")
    ),
    "Joint-LN inference was unavailable for 1 gene"
  )
  fit <- result$fit
  provenance <- fit$provenance$per_gene

  expect_identical(provenance$route, "R2")
  expect_false(provenance$inference_available)
  expect_identical(
    provenance$failure_reason,
    "competitive_hard_stop_candidate"
  )
  expect_identical(provenance$strict_failure_reason, "dispersion_boundary")
  expect_identical(provenance$resolution_status, "hard_stop")
  expect_identical(provenance$mode_status, "fit_unavailable")

  covariance_columns <- grep("^se_|^se_robust_", names(fit$summary), value = TRUE)
  wald_columns <- grep(
    "^p_naive_|^p_robust_|^wald_stat_",
    names(fit$summary),
    value = TRUE
  )
  expect_true(all(is.na(unlist(fit$summary[covariance_columns]))))
  expect_true(all(is.na(unlist(fit$summary[wald_columns]))))
  expect_true(all(!unlist(fit$summary[grep(
    "^wald_success_",
    names(fit$summary)
  )])))
  expect_true(all(unlist(fit$summary[grep(
    "^wald_failure_",
    names(fit$summary)
  )]) == "competitive_hard_stop_candidate"))
  expect_equal(dim(fit$random_effect), c(1L, length(unique(result$fx$id))))
  expect_true(all(is.na(as.matrix(fit$random_effect))))
  expect_true(all(!fit$lrt_details$gate_converged))
  expect_true(all(
    fit$lrt_details$convergence_basis == "joint_ln_endpoint_pair"
  ))
  expect_true(all(fit$lrt_details$unrestricted_route == "R2"))
  expect_true(all(!fit$lrt_details$success))
  expect_true(all(!fit$lrt_details$score_success))
})

test_that("legacy incomplete resolver records fail closed", {
  candidate <- function(label, objective, certified = TRUE,
                        finite = TRUE, required = TRUE, core = TRUE) {
    list(
      label = label,
      start_provenance = "direct",
      required = required,
      core = core,
      finite_endpoint = finite,
      certified = certified,
      objective = objective,
      optimizer = list(convergence = 0L, message = "test"),
      theta_full_log = c(0, 0, 0),
      failure = if (certified) NA_character_ else "root_residual_hard_stop"
    )
  }
  certified <- list(
    candidate("S_direct", 10),
    candidate("S_mid", 10 + 1e-10)
  )
  dominated <- c(
    certified,
    list(candidate("S_up", 10.1, certified = FALSE))
  )
  dominated_resolution <- ln_joint_resolve_candidates(dominated)
  expect_identical(dominated_resolution$route, "R2")
  expect_identical(
    dominated_resolution$failure_reason,
    "selected_candidate_certificate_inconsistent"
  )
  expect_identical(dominated_resolution$resolution$prospective_route, "R2")
  expect_false(dominated_resolution$resolution$selected_root_certified)
  expect_true(all(vapply(
    dominated_resolution$resolution$candidate_resolutions,
    function(record) {
      identical(record$class, "competitive_hard_stop") && isTRUE(record$veto)
    },
    logical(1)
  )))

  competing <- c(
    certified,
    list(candidate("S_up", 10 + 1e-7, certified = FALSE))
  )
  competing_resolution <- ln_joint_resolve_candidates(competing)
  expect_identical(competing_resolution$route, "R2")
  expect_identical(
    competing_resolution$failure_reason,
    "selected_candidate_certificate_inconsistent"
  )
  expect_identical(competing_resolution$resolution$prospective_route, "R2")
  expect_false(competing_resolution$resolution$selected_root_certified)
  expect_true(all(vapply(
    competing_resolution$resolution$candidate_resolutions,
    function(record) {
      identical(record$class, "competitive_hard_stop") && isTRUE(record$veto)
    },
    logical(1)
  )))

  missing_required <- c(
    certified,
    list(candidate(
      "S_up",
      NA_real_,
      certified = FALSE,
      finite = FALSE
    ))
  )
  missing_resolution <- ln_joint_resolve_candidates(missing_required)
  expect_identical(missing_resolution$route, "R2")
  expect_identical(missing_resolution$failure_reason, "required_start_failed")
})

test_that("a completed joint LN fit releases its objective registry entry", {
  ln_joint_objective_instance_registry("clear")
  on.exit(ln_joint_objective_instance_registry("clear"), add = TRUE)
  expect_identical(
    ln_joint_objective_instance_registry("size"),
    0L
  )

  result <- ln_joint_fit_fixture()

  expect_true(is.list(result$fit))
  expect_identical(
    ln_joint_objective_instance_registry("size"),
    0L
  )
})

test_that("fixed-beta donor modes are inferentially invariant and residual-aligned", {
  without_modes <- ln_joint_fit_fixture(output_re = FALSE)
  with_modes <- ln_joint_fit_fixture(output_re = TRUE)
  fit0 <- without_modes$fit
  fit1 <- with_modes$fit
  fx <- with_modes$fx
  prep <- ln_joint_prepare_fixture(fx)
  state <- fit1$diagnostics$raw[[1L]]$ln_joint$state

  expect_equal(fit0$summary, fit1$summary, tolerance = 1e-10)
  expect_equal(fit0$overdispersion, fit1$overdispersion, tolerance = 1e-10)
  expect_equal(fit0$convergence, fit1$convergence)
  expect_equal(fit0$algorithm, fit1$algorithm)
  expect_identical(
    fit1$provenance$random_effect,
    "pml_mode_given_ln_joint_root"
  )
  mode <- fit1$provenance$mode_diagnostics[[1L]]
  expect_true(mode$success)
  expect_identical(mode$status, "ok")
  expect_lte(mode$max_score_residual, ln_joint_thresholds()$mode_score)
  expect_lte(mode$max_newton_displacement, ln_joint_thresholds()$mode_newton)
  sigma2 <- fit1$diagnostics$raw[[1L]]$theta$subVar
  alpha <- 1 / expm1(sigma2)
  expect_equal(mode$alpha, alpha, tolerance = 1e-13)
  expect_equal(mode$lambda, alpha * exp(-sigma2 / 2), tolerance = 1e-13)
  expect_equal(ncol(fit1$random_effect), length(unique(fx$id)))

  w <- as.numeric(fit1$random_effect[1L, ])
  y <- as.numeric(fx$count[1L, ])
  eta_fixed_internal <- prep$ctx$offset +
    as.numeric(prep$ctx$pred %*% state$beta)
  gamma <- state$gamma
  lambda <- alpha * exp(-sigma2 / 2)
  mode_equation <- function(subject, value) {
    cells <- prep$ctx$fid[subject]:(prep$ctx$fid[subject + 1L] - 1L)
    ratio <- plogis(eta_fixed_internal[cells] + value - log(gamma))
    score <- alpha - lambda * exp(value) + sum(
      y[cells] - (gamma + y[cells]) * ratio
    )
    derivative <- -lambda * exp(value) - sum(
      (gamma + y[cells]) * ratio * (1 - ratio)
    )
    c(score = score, derivative = derivative)
  }
  independent_mode_state <- vapply(seq_len(prep$ctx$k), function(subject) {
    mode_equation(subject, w[subject])
  }, numeric(2L))
  expect_lte(
    max(abs(independent_mode_state["score", ])),
    ln_joint_thresholds()$mode_score
  )
  expect_lte(
    max(abs(
      independent_mode_state["score", ] /
        independent_mode_state["derivative", ]
    )),
    ln_joint_thresholds()$mode_newton
  )
  reference_modes <- vapply(seq_len(prep$ctx$k), function(subject) {
    uniroot(
      function(value) mode_equation(subject, value)[["score"]],
      interval = c(-40, 40),
      extendInt = "yes",
      tol = 1e-12
    )$root
  }, numeric(1L))
  expect_equal(w, reference_modes, tolerance = 1e-9)

  prep$ctx$output_re <- TRUE
  invalid_state <- state
  invalid_state$gamma <- -1
  failed_mode <- ln_joint_fixed_modes(
    invalid_state,
    1L,
    prep$posv,
    prep$ctx
  )
  expect_false(failed_mode$success)
  expect_identical(failed_mode$status, "solver_error")
  expect_length(failed_mode$values, prep$ctx$k)
  expect_true(all(is.na(failed_mode$values)))

  marginal <- nbresidual(
    fit1,
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    conditional = FALSE
  )$residuals
  conditional <- nbresidual(
    fit1,
    fx$count,
    fx$id,
    pred = fx$pred,
    offset = fx$offset,
    conditional = TRUE
  )$residuals

  p <- ncol(fx$pred)
  public_design <- as.matrix(fx$pred)
  intcol <- which(apply(public_design, 2L, sd) == 0)
  varying <- setdiff(seq_len(p), intcol)
  public_design[, varying] <- scale(
    public_design[, varying, drop = FALSE],
    center = TRUE,
    scale = FALSE
  )
  public_design[, intcol] <- 1
  beta <- as.numeric(fit1$summary[1L, seq_len(p)])
  sigma2 <- fit1$overdispersion$Subject[[1L]]
  cell_variance <- fit1$overdispersion$Cell[[1L]]
  donor_id <- as.numeric(factor(
    as.character(fx$id),
    levels = unique(as.character(fx$id))
  ))
  observed <- as.numeric(fx$count[1L, ])
  eta_fixed <- log(fx$offset) + as.numeric(public_design %*% beta)
  mu_marginal <- exp(eta_fixed + sigma2 / 2)
  var_marginal <- mu_marginal + mu_marginal^2 *
    (exp(sigma2) * (1 + cell_variance) - 1)
  mu_conditional <- exp(eta_fixed + w[donor_id])
  var_conditional <- mu_conditional + mu_conditional^2 * cell_variance

  expect_equal(
    as.numeric(marginal[1L, ]),
    (observed - mu_marginal) / sqrt(var_marginal),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(conditional[1L, ]),
    (observed - mu_conditional) / sqrt(var_conditional),
    tolerance = 1e-10
  )
})

test_that("joint LN API rejects unsupported switching", {
  fx <- ln_joint_test_fixture()
  expect_error(
    scrope(
      fx$count,
      fx$id,
      pred = fx$pred,
      offset = fx$offset,
      verbose = FALSE,
      cpc = 0,
      mincp = 1,
      allow_per_gene_switch = TRUE,
      ln_estimator = "joint_ln"
    ),
    "does not support allow_per_gene_switch"
  )
})

test_that("joint LN single-core and multicore results agree", {
  skip_if(parallelly::availableCores() < 2L)
  one <- ln_joint_fit_fixture(ncore = 1L)$fit
  two <- ln_joint_fit_fixture(ncore = 2L)$fit
  expect_equal(one$summary, two$summary, tolerance = 1e-8)
  expect_equal(one$overdispersion, two$overdispersion, tolerance = 1e-8)
  expect_equal(one$convergence, two$convergence)
  expect_equal(one$provenance$per_gene, two$provenance$per_gene, tolerance = 1e-8)
})
