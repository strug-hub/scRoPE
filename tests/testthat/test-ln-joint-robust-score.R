.ln_score_focused_case <- local({
  cached <- NULL

  function() {
    if (!is.null(cached)) {
      return(cached)
    }

    base <- ln_joint_fit_fixture()
    p <- ncol(base$fx$pred)
    beta_public <- as.numeric(base$fit$summary[1L, seq_len(p)])

    L_scalar <- matrix(c(1, 0, -0.5), nrow = 1L)
    b_scalar <- as.numeric(L_scalar %*% beta_public + 0.02)
    L_rank2 <- rbind(c(1, 0, 0), c(0, 0, 1))
    b_rank2 <- as.numeric(L_rank2 %*% beta_public + c(0.015, -0.01))
    row_permutation <- rbind(c(0, 1), c(1, 0))
    row_mixing <- matrix(c(1, 0.35, -0.2, 1.15), nrow = 2L,
                         byrow = TRUE)

    definitions <- list(
      scalar = list(L = L_scalar, b = b_scalar),
      scalar_scaled = list(L = -3.5 * L_scalar, b = -3.5 * b_scalar),
      rank2 = list(L = L_rank2, b = b_rank2),
      rank2_permuted = list(
        L = row_permutation %*% L_rank2,
        b = as.numeric(row_permutation %*% b_rank2)
      ),
      rank2_mixed = list(
        L = row_mixing %*% L_rank2,
        b = as.numeric(row_mixing %*% b_rank2)
      )
    )
    contrasts <- Map(function(definition, name) {
      list(name = name, L = definition$L, b = definition$b)
    }, definitions, names(definitions))

    tested <- ln_joint_fit_fixture(
      additional_tests = "score",
      lrt_contrasts = contrasts
    )
    raw_tests <- tested$fit$diagnostics$raw[[1L]]$lrt

    cached <<- list(
      base = base,
      tested = tested,
      prep = ln_joint_prepare_fixture(tested$fx),
      definitions = definitions,
      raw_tests = raw_tests,
      constrained = lapply(raw_tests, `[[`, "constrained")
    )
    cached
  }
})

.ln_score_required_fields <- c(
  "diagnostic_available", "failure_reason", "public_available",
  "public_status", "statistic", "df", "U_eff", "U_eff_subject",
  "H_psi_chi", "H_chi_chi", "J_eff", "projection", "coordinate_map",
  "diagnostics"
)

.ln_score_required_diagnostics <- c(
  "constraint_residual_scaled", "nuisance_root_scaled",
  "full_donor_additivity_scaled", "efficient_donor_additivity_scaled",
  "efficient_vs_tested_score_scaled", "J_eff_projection_identity_scaled",
  "score_sign_aggregate_scaled",
  "H_chi_chi", "J_eff"
)

test_that("joint-LN robust Score diagnostics retain certified null state contracts", {
  case <- .ln_score_focused_case()
  thresholds <- ln_joint_thresholds()

  for (name in c("scalar", "rank2")) {
    constrained <- case$constrained[[name]]
    score <- constrained$ln_joint$score_diagnostic
    r <- nrow(constrained$constraint$L)
    p <- length(constrained$theta$beta)
    d <- p + 2L
    k <- ncol(constrained$ln_joint$state$U_subject_full)

    expect_true(all(.ln_score_required_fields %in% names(score)), info = name)
    expect_true(
      all(.ln_score_required_diagnostics %in% names(score$diagnostics)),
      info = name
    )
    expect_true(score$diagnostic_available, info = name)
    expect_true(is.na(score$failure_reason), info = name)
    expect_false(score$public_available, info = name)
    expect_identical(score$public_status, "joint_ln_calibration_pending")
    expect_identical(score$df, as.integer(r))
    expect_true(is.finite(score$statistic) && score$statistic >= 0, info = name)
    expect_length(score$U_eff, r)
    expect_identical(dim(score$U_eff_subject), c(r, k))
    expect_identical(dim(score$H_psi_chi), c(r, d - r))
    expect_identical(dim(score$H_chi_chi), c(d - r, d - r))
    expect_identical(dim(score$J_eff), c(r, r))
    expect_identical(dim(score$projection), c(r, d))
    expect_identical(dim(score$coordinate_map$full), c(d, d))
    expect_identical(dim(score$coordinate_map$interest_basis), c(p, r))
    expect_identical(dim(score$coordinate_map$null_basis), c(p, p - r))

    constraint_residual <- as.numeric(
      constrained$constraint$L %*% constrained$theta$beta -
        constrained$constraint$b
    )
    expect_lte(max(abs(constraint_residual)), 1e-9)
    expect_lte(
      score$diagnostics$constraint_residual_scaled,
      thresholds$constraint
    )
    expect_lte(
      score$diagnostics$nuisance_root_scaled,
      thresholds$root
    )
    expect_lte(
      score$diagnostics$full_donor_additivity_scaled,
      thresholds$donor_additivity
    )
    expect_lte(
      score$diagnostics$efficient_donor_additivity_scaled,
      thresholds$donor_additivity
    )
    expect_lte(
      score$diagnostics$efficient_vs_tested_score_scaled,
      thresholds$root
    )
    expect_lte(
      score$diagnostics$J_eff_projection_identity_scaled,
      thresholds$covariance_identity
    )
    expect_lte(
      score$diagnostics$score_sign_aggregate_scaled,
      thresholds$donor_additivity
    )
    expect_identical(score$diagnostics$H_chi_chi$rank, as.integer(d - r))
    expect_identical(score$diagnostics$J_eff$rank, as.integer(r))
    expect_true(score$diagnostics$H_chi_chi$positive_definite, info = name)
    expect_true(score$diagnostics$H_chi_chi$condition_acceptable, info = name)
    expect_true(score$diagnostics$H_chi_chi$usable, info = name)
    expect_true(score$diagnostics$J_eff$positive_definite, info = name)
    expect_true(score$diagnostics$J_eff$condition_acceptable, info = name)
    expect_true(score$diagnostics$J_eff$usable, info = name)

    reduced_root <- ln_joint_root_metric(
      constrained$reduced$per_subject_gradients,
      constrained$reduced$score
    )
    expect_lte(reduced_root, thresholds$root)
  }
})

test_that("full constrained LN score and selected Hessian columns match finite differences", {
  case <- .ln_score_focused_case()
  prep <- case$prep
  nb <- prep$ctx$nb
  full_map <- ln_joint_map(nb)

  objective_log <- function(z) {
    ln_joint_objective(
      c(z[seq_len(nb)], exp(z[nb + seq_len(2L)])),
      1L,
      prep$posv,
      prep$ctx,
      full_map
    )
  }
  maximizing_score <- function(z) {
    ln_joint_derivative_state(
      z,
      1L,
      prep$posv,
      prep$ctx,
      full_map
    )$U_full
  }

  for (name in c("scalar", "rank2")) {
    state <- case$constrained[[name]]$ln_joint$state
    theta <- state$full_log
    fd_q_gradient <- ln_joint_fd_gradient(objective_log, theta)
    expect_lte(
      ln_joint_scaled_max(state$U_full, -fd_q_gradient),
      ln_joint_thresholds()$objective_gradient
    )

    selected <- unique(c(1L, nb, nb + 1L, nb + 2L))
    fd_score_jacobian <- ln_joint_fd_jacobian(maximizing_score, theta)
    expect_lte(
      ln_joint_scaled_max(
        state$A_full[, selected, drop = FALSE],
        -fd_score_jacobian[, selected, drop = FALSE]
      ),
      ln_joint_thresholds()$sensitivity_score
    )
    expect_lte(
      ln_joint_scaled_max(
        rowSums(state$U_subject_full),
        state$U_full
      ),
      ln_joint_thresholds()$donor_additivity
    )
    expect_lte(
      ln_joint_scaled_max(state$U_full, -state$q_gradient_full),
      ln_joint_thresholds()$donor_additivity
    )
  }
})

test_that("efficient-score projection and variability identities reconstruct independently", {
  case <- .ln_score_focused_case()

  for (name in c("scalar", "rank2")) {
    constrained <- case$constrained[[name]]
    score <- constrained$ln_joint$score_diagnostic
    state <- constrained$ln_joint$state
    r <- score$df
    p <- length(state$beta)
    d <- length(state$U_full)
    beta_map <- cbind(
      score$coordinate_map$interest_basis,
      score$coordinate_map$null_basis
    )
    T <- diag(d)
    T[seq_len(p), seq_len(p)] <- beta_map
    expect_equal(score$coordinate_map$full, T, tolerance = 0, info = name)
    nuisance <- (r + 1L):d

    U <- as.numeric(crossprod(T, state$U_full))
    U_subject <- crossprod(T, state$U_subject_full)
    H <- crossprod(T, state$A_full %*% T)
    H_psi_chi <- H[seq_len(r), nuisance, drop = FALSE]
    H_chi_chi <- H[nuisance, nuisance, drop = FALSE]
    nuisance_projection <- t(solve(t(H_chi_chi), t(H_psi_chi)))
    projection <- cbind(diag(r), -nuisance_projection)
    U_eff <- as.numeric(projection %*% U)
    U_eff_subject <- projection %*% U_subject
    J_eff <- tcrossprod(U_eff_subject)
    J_full <- tcrossprod(U_subject)
    J_eff_projected <- projection %*% J_full %*% t(projection)
    statistic <- as.numeric(crossprod(U_eff, solve(J_eff, U_eff)))

    expect_equal(score$H_psi_chi, H_psi_chi, tolerance = 1e-12,
                 info = name)
    expect_equal(score$H_chi_chi, H_chi_chi, tolerance = 1e-12,
                 info = name)
    expect_equal(score$projection, projection, tolerance = 1e-12,
                 info = name)
    expect_equal(score$U_eff, U_eff, tolerance = 1e-12, info = name)
    expect_equal(score$U_eff_subject, U_eff_subject, tolerance = 1e-12,
                 info = name)
    expect_equal(rowSums(U_eff_subject), U_eff, tolerance = 1e-12,
                 info = name)
    expect_equal(score$J_eff, J_eff, tolerance = 1e-12, info = name)
    expect_equal(J_eff, J_eff_projected, tolerance = 1e-12, info = name)
    expect_equal(score$statistic, statistic, tolerance = 1e-12,
                 info = name)
  }
})

test_that("audit-only Score path does not calculate ratio or public inference", {
  case <- .ln_score_focused_case()
  forbidden <- c(
    "wP", "qP", "pace", "multiplier", "compute_profile_lrt",
    "compute_adj_score", "pchisq"
  )
  helper_text <- paste(
    deparse(body(ln_joint_robust_score_diagnostic)),
    collapse = "\n"
  )
  helper_text <- tolower(helper_text)
  for (token in forbidden) {
    expect_false(
      grepl(tolower(token), helper_text, fixed = TRUE),
      info = token
    )
  }
  for (name in names(case$raw_tests)) {
    entry <- case$raw_tests[[name]]
    expect_null(entry$raw, info = name)
    expect_true(is.na(entry$stat), info = name)
    expect_true(is.na(entry$pval), info = name)
    expect_true(all(is.na(c(entry$wP, entry$wU, entry$qP))), info = name)
    expect_false(entry$score_success, info = name)
    expect_true(is.na(entry$score_stat), info = name)
    expect_true(is.na(entry$score_pval), info = name)
    expect_identical(
      entry$score_failure,
      "joint_ln_calibration_pending"
    )
    expect_identical(entry$score_test$method, "unavailable_joint_ln")
    expect_true(is.na(entry$score_test$stat), info = name)
    expect_true(is.na(entry$score_test$pval), info = name)
  }
})

test_that("robust Score statistic is invariant to equivalent hypothesis representations", {
  case <- .ln_score_focused_case()
  score <- function(name) {
    case$constrained[[name]]$ln_joint$score_diagnostic$statistic
  }

  expect_equal(score("scalar"), score("scalar_scaled"), tolerance = 1e-8)
  expect_equal(score("rank2"), score("rank2_permuted"), tolerance = 1e-8)
  expect_equal(score("rank2"), score("rank2_mixed"), tolerance = 1e-8)
  expect_equal(
    case$constrained$scalar$theta$beta,
    case$constrained$scalar_scaled$theta$beta,
    tolerance = 1e-7
  )
  expect_equal(
    case$constrained$rank2$theta$beta,
    case$constrained$rank2_permuted$theta$beta,
    tolerance = 1e-7
  )
  expect_equal(
    case$constrained$rank2$theta$beta,
    case$constrained$rank2_mixed$theta$beta,
    tolerance = 1e-7
  )
})

test_that("alternative null and interest bases leave the robust Score unchanged", {
  case <- .ln_score_focused_case()
  scalar <- case$constrained$scalar
  scalar_reference <- scalar$ln_joint$score_diagnostic
  K <- scalar_reference$coordinate_map$null_basis
  K_transform <- matrix(c(1, 0.3, -0.2, 1.15), nrow = 2L,
                        byrow = TRUE)
  K_alternative <- K %*% K_transform
  expect_lte(max(abs(scalar$constraint$L %*% K_alternative)), 1e-12)
  expect_identical(qr(K_alternative)$rank, ncol(K_alternative))
  alternative_K <- ln_joint_robust_score_diagnostic(
    scalar,
    null_basis = K_alternative
  )

  rank2 <- case$constrained$rank2
  rank2_reference <- rank2$ln_joint$score_diagnostic
  C <- rank2_reference$coordinate_map$interest_basis
  tested_transform <- matrix(c(1.25, 0.2, -0.15, 0.85), nrow = 2L,
                             byrow = TRUE)
  C_alternative <- C %*% tested_transform
  expect_identical(qr(C_alternative)$rank, ncol(C_alternative))
  alternative_C <- ln_joint_robust_score_diagnostic(
    rank2,
    interest_basis = C_alternative
  )

  alternatives <- list(
    list(candidate = alternative_K, reference = scalar_reference),
    list(candidate = alternative_C, reference = rank2_reference)
  )
  for (entry in alternatives) {
    candidate <- entry$candidate
    reference <- entry$reference
    expect_true(candidate$diagnostic_available)
    expect_equal(candidate$statistic, reference$statistic, tolerance = 1e-9)
    expect_identical(candidate$df, reference$df)
    expect_false(candidate$public_available)
    expect_identical(candidate$public_status, "joint_ln_calibration_pending")
  }
})

test_that("rank-deficient tested variability is explicitly unavailable without repair", {
  case <- .ln_score_focused_case()
  constrained <- case$constrained$rank2
  mocked <- constrained
  state <- mocked$ln_joint$state
  d <- length(state$U_full)
  k <- ncol(state$U_subject_full)

  state$U_full <- rep(0, d)
  state$U_subject_full <- matrix(0, nrow = d, ncol = k)
  state$q_gradient_full <- rep(0, d)
  state$U_reduced[] <- 0
  state$U_subject_reduced[,] <- 0
  state$q_gradient_reduced[] <- 0
  state$raw$gradient[] <- 0
  state$raw$per_subject_gradients[,] <- 0
  state$raw$gradient_aggregate_formula[] <- 0
  state$raw$gradient_aggregate_formula_raw[] <- 0
  mocked$ln_joint$state <- state
  mocked$score <- state$U_full
  mocked$reduced$score[] <- 0
  mocked$reduced$g_red[] <- 0
  mocked$reduced$per_subject_gradients[,] <- 0

  score <- ln_joint_robust_score_diagnostic(mocked)
  expect_false(score$diagnostic_available)
  expect_true(is.na(score$statistic))
  expect_false(score$public_available)
  expect_identical(score$public_status, "joint_ln_calibration_pending")
  expect_match(score$failure_reason, "J_eff|variability|rank|singular",
               ignore.case = TRUE)
  expect_true(is.list(score$diagnostics$J_eff))
  expect_lt(score$diagnostics$J_eff$rank, score$df)
  expect_false(any(c("ridge", "pseudoinverse", "eigenvalue_clipping") %in%
                     names(score)))
})

test_that("Score diagnostics fail closed for invalid sensitivity and provenance", {
  case <- .ln_score_focused_case()
  reference <- case$constrained$scalar
  d <- nrow(reference$ln_joint$state$A_full)

  rank_deficient_H <- reference
  rank_deficient_H$ln_joint$state$A_full <- matrix(0, d, d)
  rank_failure <- ln_joint_robust_score_diagnostic(rank_deficient_H)
  expect_false(rank_failure$diagnostic_available)
  expect_identical(rank_failure$failure_reason,
                   "H_chi_chi_rank_deficient")

  nonpositive_H <- reference
  nonpositive_H$ln_joint$state$A_full <- -diag(d)
  curvature_failure <- ln_joint_robust_score_diagnostic(nonpositive_H)
  expect_false(curvature_failure$diagnostic_available)
  expect_identical(curvature_failure$failure_reason,
                   "H_chi_chi_not_positive_definite")

  asymmetric_H <- reference
  asymmetric_H$ln_joint$state$A_full[1L, 2L] <-
    asymmetric_H$ln_joint$state$A_full[1L, 2L] + 1
  symmetry_failure <- ln_joint_robust_score_diagnostic(asymmetric_H)
  expect_false(symmetry_failure$diagnostic_available)
  expect_identical(symmetry_failure$failure_reason,
                   "transformed_sensitivity_asymmetric")

  missing_certificate <- reference
  missing_certificate$ln_joint$candidate_records <- lapply(
    missing_certificate$ln_joint$candidate_records,
    function(candidate) {
      candidate$certified <- FALSE
      candidate
    }
  )
  certificate_failure <- ln_joint_robust_score_diagnostic(
    missing_certificate
  )
  expect_false(certificate_failure$diagnostic_available)
  expect_identical(certificate_failure$failure_reason,
                   "constrained_endpoint_certification_missing")

  inconsistent_state <- reference
  inconsistent_state$ln_joint$state$beta[1L] <-
    inconsistent_state$ln_joint$state$beta[1L] + 1e-8
  state_failure <- ln_joint_robust_score_diagnostic(inconsistent_state)
  expect_false(state_failure$diagnostic_available)
  expect_identical(state_failure$failure_reason,
                   "constrained_state_coordinate_mismatch")

  roundtrip_grid <- seq(-20, 20, length.out = 10001L)
  roundtrip_index <- which(vapply(
    roundtrip_grid,
    function(x) !identical(log(exp(x)), x),
    logical(1)
  ))[1L]
  expect_true(is.finite(roundtrip_index))
  roundtrip_log <- roundtrip_grid[roundtrip_index]
  expect_false(identical(log(exp(roundtrip_log)), roundtrip_log))
  roundtrip_state <- reference
  dispersion_index <- length(roundtrip_state$ln_joint$state$full_log)
  roundtrip_state$ln_joint$state$full_log[dispersion_index] <-
    roundtrip_log
  roundtrip_state$ln_joint$state$gamma <- exp(roundtrip_log)
  roundtrip_state$theta$cellVar <- exp(roundtrip_log)
  matching_record <- which(vapply(
    roundtrip_state$ln_joint$candidate_records,
    function(candidate) identical(
      as.numeric(candidate$theta_full_log),
      as.numeric(reference$ln_joint$state$full_log)
    ),
    logical(1)
  ))
  expect_length(matching_record, 1L)
  roundtrip_state$ln_joint$candidate_records[[matching_record]]$
    theta_full_log[dispersion_index] <- roundtrip_log
  roundtrip_result <- ln_joint_robust_score_diagnostic(roundtrip_state)
  expect_false(identical(
    roundtrip_result$failure_reason,
    "constrained_state_coordinate_mismatch"
  ))

  inconsistent_return <- reference
  inconsistent_return$theta$beta[1L] <-
    inconsistent_return$theta$beta[1L] + 1e-8
  return_failure <- ln_joint_robust_score_diagnostic(inconsistent_return)
  expect_false(return_failure$diagnostic_available)
  expect_identical(return_failure$failure_reason,
                   "constrained_return_state_mismatch")

  mismatched_certificate <- reference
  matching_record <- which(vapply(
    mismatched_certificate$ln_joint$candidate_records,
    function(candidate) identical(
      as.numeric(candidate$theta_full_log),
      as.numeric(mismatched_certificate$ln_joint$state$full_log)
    ),
    logical(1)
  ))
  expect_length(matching_record, 1L)
  mismatched_certificate$ln_joint$candidate_records[[matching_record]]$
    theta_full_log[1L] <-
      mismatched_certificate$ln_joint$candidate_records[[matching_record]]$
        theta_full_log[1L] + 1e-12
  exact_certificate_failure <- ln_joint_robust_score_diagnostic(
    mismatched_certificate
  )
  expect_false(exact_certificate_failure$diagnostic_available)
  expect_identical(
    exact_certificate_failure$failure_reason,
    "constrained_endpoint_certification_missing"
  )

  redundant_hypothesis <- reference
  redundant_hypothesis$constraint$L <- rbind(
    reference$constraint$L,
    reference$constraint$L
  )
  redundant_hypothesis$constraint$b <- rep(reference$constraint$b, 2L)
  rank_failure <- ln_joint_robust_score_diagnostic(redundant_hypothesis)
  expect_false(rank_failure$diagnostic_available)
  expect_identical(rank_failure$failure_reason,
                   "hypothesis_rank_deficient")
  expect_identical(rank_failure$hypothesis$rank, 1L)
  expect_true(is.na(rank_failure$df))

  invalid_K <- reference$constraint$null_basis
  invalid_K[, 1L] <- invalid_K[, 1L] +
    reference$ln_joint$score_diagnostic$coordinate_map$interest_basis[, 1L]
  kernel_failure <- ln_joint_robust_score_diagnostic(
    reference,
    null_basis = invalid_K
  )
  expect_false(kernel_failure$diagnostic_available)
  expect_identical(kernel_failure$failure_reason,
                   "null_basis_not_in_hypothesis_kernel")

  nearly_singular <- matrix(
    c(1, 1 - 1e-11, 1 - 1e-11, 1),
    nrow = 2L
  )
  conditioning <- ln_joint_score_matrix_diagnostics(nearly_singular)
  expect_true(conditioning$full_rank)
  expect_true(conditioning$positive_definite)
  expect_false(conditioning$condition_acceptable)
  expect_false(conditioning$usable)
})

test_that("boundary, non-R0, and uncertified constrained fits stay unavailable", {
  expect_warning(
    boundary <- ln_joint_fit_fixture(
      max = c(0.01, 1000),
      additional_tests = "score",
      lrt_contrasts = list(list(
        name = "boundary",
        L = matrix(c(1, 0, 0), nrow = 1L),
        b = 0
      ))
    )$fit,
    "Joint-LN inference was unavailable for 1 gene"
  )
  constrained_boundary <- boundary$diagnostics$raw[[1L]]$lrt$boundary$constrained
  boundary_score <- constrained_boundary$ln_joint$score_diagnostic

  expect_identical(constrained_boundary$route, "R2")
  expect_false(boundary_score$diagnostic_available)
  expect_true(is.na(boundary_score$statistic))
  expect_false(boundary_score$public_available)
  expect_identical(boundary_score$public_status,
                   "joint_ln_calibration_pending")
  boundary_score_success <- grep(
    "^score_success_",
    names(boundary$summary),
    value = TRUE
  )
  boundary_score_stat <- grep(
    "^score_adj_",
    names(boundary$summary),
    value = TRUE
  )
  boundary_test_count <- length(
    boundary$diagnostics$raw[[1L]]$lrt
  )
  expect_gt(boundary_test_count, 0L)
  expect_length(boundary_score_success, boundary_test_count)
  expect_length(boundary_score_stat, boundary_test_count)
  expect_true(all(!unlist(boundary$summary[boundary_score_success])))
  expect_true(all(is.na(unlist(boundary$summary[boundary_score_stat]))))

  case <- .ln_score_focused_case()
  reference <- case$constrained$scalar
  non_r0 <- reference
  non_r0$route <- "R2"
  non_r0$inference_available <- FALSE
  non_r0$failure_reason <- "root_residual_hard_stop"
  non_r0$lrt_gate$converged <- FALSE
  non_r0_score <- ln_joint_robust_score_diagnostic(non_r0)
  expect_false(non_r0_score$diagnostic_available)
  expect_true(is.na(non_r0_score$statistic))

  uncertified <- reference
  uncertified$inference_available <- FALSE
  uncertified$lrt_gate$converged <- FALSE
  uncertified$ln_joint$candidate_records <- lapply(
    uncertified$ln_joint$candidate_records,
    function(candidate) {
      candidate$certified <- FALSE
      candidate
    }
  )
  uncertified_score <- ln_joint_robust_score_diagnostic(uncertified)
  expect_false(uncertified_score$diagnostic_available)
  expect_true(is.na(uncertified_score$statistic))
  expect_false(uncertified_score$public_available)
  expect_identical(uncertified_score$public_status,
                   "joint_ln_calibration_pending")
})

test_that("audit-only Score diagnostics do not change joint Wald or public Score fields", {
  case <- .ln_score_focused_case()
  baseline <- case$base$fit
  tested <- case$tested$fit
  raw_baseline <- baseline$diagnostics$raw[[1L]]
  raw_tested <- tested$diagnostics$raw[[1L]]

  primary_summary_columns <- setdiff(
    names(baseline$summary),
    c("test_status", "test_failure_detail")
  )
  expect_equal(
    tested$summary[, primary_summary_columns, drop = FALSE],
    baseline$summary[, primary_summary_columns, drop = FALSE],
    tolerance = 1e-10
  )
  expect_equal(tested$overdispersion, baseline$overdispersion,
               tolerance = 1e-10)
  expect_equal(tested$convergence, baseline$convergence, tolerance = 0)
  expect_equal(tested$algorithm, baseline$algorithm, tolerance = 0)
  expect_equal(raw_tested$theta, raw_baseline$theta, tolerance = 0)
  expect_equal(raw_tested$naive_cov, raw_baseline$naive_cov, tolerance = 0)
  expect_equal(raw_tested$robust_cov, raw_baseline$robust_cov, tolerance = 0)
  expect_equal(raw_tested$score, raw_baseline$score, tolerance = 0)
  score_success_columns <- grep(
    "^score_success_",
    names(tested$summary),
    value = TRUE
  )
  score_stat_columns <- grep(
    "^score_adj_",
    names(tested$summary),
    value = TRUE
  )
  score_p_columns <- grep(
    "^p_score_adj_",
    names(tested$summary),
    value = TRUE
  )
  score_failure_columns <- grep(
    "^score_failure_",
    names(tested$summary),
    value = TRUE
  )
  expected_score_columns <- length(raw_tested$lrt)
  expect_gt(expected_score_columns, 0L)
  expect_length(score_success_columns, expected_score_columns)
  expect_length(score_stat_columns, expected_score_columns)
  expect_length(score_p_columns, expected_score_columns)
  expect_length(score_failure_columns, expected_score_columns)
  expect_true(all(!unlist(tested$summary[score_success_columns])))
  expect_true(all(is.na(unlist(tested$summary[score_stat_columns]))))
  expect_true(all(is.na(unlist(tested$summary[score_p_columns]))))
  expect_true(all(
    unlist(tested$summary[score_failure_columns]) ==
      "joint_ln_calibration_pending"
  ))
})

test_that("requesting Score does not restore the retired legacy estimator", {
  expect_error(scrope(ln_estimator = "legacy_post_pml",
                      additional_tests = "score"), "retired")
})