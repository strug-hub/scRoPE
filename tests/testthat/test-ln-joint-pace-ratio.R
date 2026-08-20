.ln_ratio_fit <- function(fx, additional_tests = character(0),
                          lrt_contrasts = NULL) {
  args <- list(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    offset = fx$offset,
    min = c(1e-4, 1e-4),
    max = c(10, 1000),
    verbose = FALSE,
    cpc = 0,
    mincp = 1,
    ncore = 1L,
    output_re = FALSE,
    additional_tests = additional_tests,
    lrt_contrasts = lrt_contrasts,
    lrt_details = "lrt" %in% additional_tests,
    keep_diagnostics = TRUE,
    ln_estimator = "joint_ln"
  )
  fallback_supported <- "ln_operational_fallback" %in%
    names(formals(scrope))
  if (fallback_supported) {
    args$ln_operational_fallback <- "none"
  }
  list(
    fit = suppressWarnings(do.call(scrope, args)),
    fallback_setting = if (fallback_supported) {
      "none"
    } else {
      "absent_at_ce623d1"
    }
  )
}

.ln_ratio_focused_case <- local({
  cached <- NULL

  function() {
    if (!is.null(cached)) {
      return(cached)
    }

    original_kind <- RNGkind()
    had_original_seed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    original_seed <- if (had_original_seed) .Random.seed else NULL
    on.exit({
      do.call(RNGkind, as.list(original_kind))
      if (had_original_seed) {
        assign(".Random.seed", original_seed, envir = .GlobalEnv)
      } else if (exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    RNGkind(
      kind = "Mersenne-Twister",
      normal.kind = "Inversion",
      sample.kind = "Rejection"
    )
    set.seed(20260723)

    fx <- ln_joint_test_fixture()
    base <- .ln_ratio_fit(fx)
    p <- ncol(fx$pred)
    beta_public <- as.numeric(base$fit$summary[1L, seq_len(p)])

    L_scalar <- matrix(c(1, 0, -0.5), nrow = 1L)
    b_scalar <- as.numeric(L_scalar %*% beta_public + 0.02)
    L_nonintercept <- matrix(c(0, 1, 0), nrow = 1L)
    b_nonintercept <- as.numeric(
      L_nonintercept %*% beta_public + 0.012
    )
    L_rank2 <- rbind(c(1, 0, 0), c(0, 0, 1))
    b_rank2 <- as.numeric(L_rank2 %*% beta_public + c(0.015, -0.01))
    row_permutation <- rbind(c(0, 1), c(1, 0))
    row_mixing <- matrix(
      c(1, 0.35, -0.2, 1.15),
      nrow = 2L,
      byrow = TRUE
    )

    definitions <- list(
      scalar = list(L = L_scalar, b = b_scalar),
      scalar_scaled = list(L = -3.5 * L_scalar, b = -3.5 * b_scalar),
      nonintercept = list(L = L_nonintercept, b = b_nonintercept),
      nonintercept_scaled = list(
        L = -2.25 * L_nonintercept,
        b = -2.25 * b_nonintercept
      ),
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
    tested <- .ln_ratio_fit(
      fx,
      additional_tests = c("score", "lrt"),
      lrt_contrasts = contrasts
    )
    raw <- tested$fit$diagnostics$raw[[1L]]
    constrained <- lapply(raw$lrt, `[[`, "constrained")
    diagnostics <- lapply(constrained, function(constrained_fit) {
      constrained_fit$ln_joint$pace_ratio_diagnostic
    })

    cached <<- list(
      fx = fx,
      prep = ln_joint_prepare_fixture(fx),
      base = base,
      tested = tested,
      raw = raw,
      definitions = definitions,
      constrained = constrained,
      diagnostics = diagnostics
    )
    cached
  }
})

.ln_ratio_objective_context <- function(case, gene_index = 1L) {
  list(
    gene_index = as.integer(gene_index),
    posv = case$prep$posv,
    ctx = case$prep$ctx
  )
}

.ln_ratio_replay <- function(case, constrained_fit,
                             unrestricted_fit = case$raw,
                             null_basis = NULL,
                             interest_basis = NULL) {
  ln_joint_pace_ratio_diagnostic(
    unrestricted_fit,
    constrained_fit,
    null_basis = null_basis,
    interest_basis = interest_basis,
    objective_context = .ln_ratio_objective_context(case)
  )
}

.ln_ratio_required_fields <- c(
  "diagnostic_available", "failure_reason", "public_available",
  "public_status", "statistic", "adjusted_ratio", "df", "wP", "wU",
  "qP", "multiplier", "scale", "U_eff", "U_eff_subject", "H_eff",
  "J_eff", "coordinate_map", "projection", "objective_provenance",
  "objective_instance", "endpoints", "score_diagnostic", "diagnostics"
)

test_that("Pace ratio audit retains certified same-objective endpoint contracts", {
  case <- .ln_ratio_focused_case()
  thresholds <- ln_joint_thresholds()

  for (name in c("scalar", "nonintercept", "rank2")) {
    ratio <- case$diagnostics[[name]]
    constrained <- case$constrained[[name]]
    r <- nrow(constrained$constraint$L)
    d <- length(constrained$ln_joint$state$U_full)
    k <- ncol(constrained$ln_joint$state$U_subject_full)

    expect_true(all(.ln_ratio_required_fields %in% names(ratio)),
                info = name)
    expect_true(ratio$diagnostic_available, info = name)
    expect_equal(
      ratio,
      .ln_ratio_replay(case, constrained),
      tolerance = 0,
      info = name
    )
    expect_true(is.na(ratio$failure_reason), info = name)
    expect_false(ratio$public_available, info = name)
    expect_identical(ratio$public_status, "joint_ln_calibration_pending")
    expect_true(ratio$no_numerical_repair, info = name)
    expect_identical(ratio$df, as.integer(r))
    expect_identical(ratio$theoretical_df, as.integer(r))
    expect_false("p_value" %in% names(ratio))
    expect_false("pval" %in% names(ratio))

    expect_true(all(vapply(
      ratio$endpoints,
      function(endpoint) isTRUE(endpoint$success),
      logical(1)
    )), info = name)
    for (endpoint in ratio$endpoints) {
      expect_equal(
        endpoint$state_objective,
        endpoint$reevaluated_objective,
        tolerance = 0,
        info = name
      )
      expect_true(
        endpoint$diagnostics$objective_recheck_scaled <=
          thresholds$objective_agreement,
        info = name
      )
      recheck_fields <- grep(
        "_recheck_scaled$",
        names(endpoint$diagnostics),
        value = TRUE
      )
      expect_true(all(is.finite(unlist(
        endpoint$diagnostics[recheck_fields],
        use.names = FALSE
      ))), info = name)
    }
    expect_identical(
      ratio$endpoints$unrestricted$objective_provenance,
      ratio$endpoints$constrained$objective_provenance
    )
    expect_identical(
      ratio$objective_provenance,
      ln_joint_objective_provenance()
    )
    expect_identical(
      ratio$endpoints$unrestricted$objective_instance,
      ratio$endpoints$constrained$objective_instance
    )
    expect_identical(
      ratio$objective_instance,
      ratio$endpoints$unrestricted$objective_instance
    )
    expect_true(ratio$diagnostics$same_criterion, info = name)
    expect_true(ratio$diagnostics$same_objective_instance, info = name)
    expect_true(ratio$diagnostics$same_objective, info = name)
    expect_true(ratio$diagnostics$same_gene, info = name)
    expect_true(ratio$diagnostics$same_dimension, info = name)
    expect_true(ratio$diagnostics$same_donor_count, info = name)

    expect_true(all(is.finite(c(
      ratio$wP,
      ratio$wU,
      ratio$qP,
      ratio$multiplier,
      ratio$statistic
    ))), info = name)
    expect_true(ratio$wP >= 0, info = name)
    expect_true(ratio$wU > 0, info = name)
    expect_true(ratio$qP > 0, info = name)
    expect_true(ratio$multiplier > 0, info = name)
    expect_equal(ratio$scale, ratio$multiplier, tolerance = 0,
                 info = name)
    expect_equal(
      ratio$statistic,
      ratio$wP * ratio$wU / ratio$qP,
      tolerance = 1e-12,
      info = name
    )
    expect_equal(ratio$adjusted_ratio, ratio$statistic, tolerance = 0,
                 info = name)
    expect_equal(
      ratio$wU,
      constrained$ln_joint$score_diagnostic$statistic,
      tolerance = 1e-12,
      info = name
    )

    expect_identical(dim(ratio$H_eff), c(r, r))
    expect_identical(dim(ratio$J_eff), c(r, r))
    expect_length(ratio$U_eff, r)
    expect_identical(dim(ratio$U_eff_subject), c(r, k))
    expect_identical(dim(ratio$projection), c(r, d))
    expect_true(ratio$diagnostics$H_eff$usable, info = name)
    expect_true(ratio$diagnostics$J_eff$usable, info = name)
    expect_true(
      ratio$diagnostics$wP_loglik_identity_scaled <=
        thresholds$objective_agreement,
      info = name
    )
    expect_true(
      ratio$diagnostics$H_eff_projection_identity_scaled <=
        thresholds$covariance_identity,
      info = name
    )
    expect_true(
      ratio$diagnostics$wU_score_identity_scaled <=
        thresholds$covariance_identity,
      info = name
    )
    expect_true(
      ratio$diagnostics$pace_decomposition_identity_scaled <=
        thresholds$covariance_identity,
      info = name
    )
  }
})

test_that("Pace ingredients reconstruct from exact objectives and Schur blocks", {
  case <- .ln_ratio_focused_case()
  prep <- case$prep
  nb <- prep$ctx$nb
  thresholds <- ln_joint_thresholds()
  unrestricted_state <- case$raw$ln_joint$state
  unrestricted_map <- ln_joint_map(nb)
  unrestricted_natural <- c(
    unrestricted_state$beta,
    unrestricted_state$sigma2,
    unrestricted_state$gamma
  )
  unrestricted_objective <- ln_joint_objective(
    unrestricted_natural,
    1L,
    prep$posv,
    prep$ctx,
    unrestricted_map
  )
  expect_equal(
    unrestricted_objective,
    unrestricted_state$raw$value,
    tolerance = 1e-12
  )
  expect_equal(unrestricted_objective, -case$raw$loglik_ln,
               tolerance = 1e-12)
  expect_equal(
    sum(ln_joint_donor_objective_values(
      unrestricted_state$theta_reduced_log,
      1L,
      prep$posv,
      prep$ctx,
      unrestricted_map
    )),
    unrestricted_objective,
    tolerance = 1e-10
  )
  unrestricted_fd <- ln_joint_fd_audit(
    unrestricted_state$theta_reduced_log,
    1L,
    prep$posv,
    prep$ctx,
    unrestricted_map
  )
  expect_true(all(unrestricted_fd$target_pass))

  for (name in c("scalar", "nonintercept", "rank2")) {
    constrained <- case$constrained[[name]]
    ratio <- case$diagnostics[[name]]
    state <- constrained$ln_joint$state
    definition <- case$definitions[[name]]
    L_internal <- constrained$constraint$L
    constrained_map <- ln_joint_map(nb, L_internal, definition$b)
    free_dim <- constrained_map$free_dim
    constrained_natural <- c(
      if (free_dim) {
        state$theta_reduced_log[seq_len(free_dim)]
      } else {
        numeric(0)
      },
      exp(state$theta_reduced_log[free_dim + seq_len(2L)])
    )
    constrained_objective <- ln_joint_objective(
      constrained_natural,
      1L,
      prep$posv,
      prep$ctx,
      constrained_map
    )
    expect_equal(constrained_objective, state$raw$value, tolerance = 1e-12,
                 info = name)
    expect_equal(constrained_objective, -constrained$loglik_ln,
                 tolerance = 1e-12, info = name)
    expect_equal(
      sum(ln_joint_donor_objective_values(
        state$theta_reduced_log,
        1L,
        prep$posv,
        prep$ctx,
        constrained_map
      )),
      constrained_objective,
      tolerance = 1e-10,
      info = name
    )
    constrained_fd <- ln_joint_fd_audit(
      state$theta_reduced_log,
      1L,
      prep$posv,
      prep$ctx,
      constrained_map
    )
    expect_true(all(constrained_fd$target_pass), info = name)

    T <- ratio$coordinate_map$full
    H_coordinate <- crossprod(T, state$A_full %*% T)
    H_coordinate <- 0.5 * (H_coordinate + t(H_coordinate))
    r <- ratio$df
    nuisance <- r + seq_len(nrow(H_coordinate) - r)
    H_eff <- H_coordinate[seq_len(r), seq_len(r), drop = FALSE] -
      H_coordinate[seq_len(r), nuisance, drop = FALSE] %*%
        solve(H_coordinate[nuisance, nuisance, drop = FALSE]) %*%
        H_coordinate[nuisance, seq_len(r), drop = FALSE]
    U_coordinate <- as.numeric(crossprod(T, state$U_full))
    U_subject_coordinate <- crossprod(T, state$U_subject_full)
    projection <- cbind(
      diag(r),
      -H_coordinate[seq_len(r), nuisance, drop = FALSE] %*%
        solve(H_coordinate[nuisance, nuisance, drop = FALSE])
    )
    U_eff <- as.numeric(projection %*% U_coordinate)
    U_eff_subject <- projection %*% U_subject_coordinate
    J_eff <- tcrossprod(U_eff_subject)
    qP <- as.numeric(crossprod(U_eff, solve(H_eff, U_eff)))
    wU <- as.numeric(crossprod(U_eff, solve(J_eff, U_eff)))
    wP <- 2 * (constrained_objective - unrestricted_objective)

    expect_equal(ratio$H_eff, H_eff, tolerance = 1e-12, info = name)
    expect_equal(ratio$J_eff, J_eff, tolerance = 1e-12, info = name)
    expect_equal(ratio$U_eff, U_eff, tolerance = 1e-12, info = name)
    expect_equal(ratio$qP, qP, tolerance = 1e-12, info = name)
    expect_equal(ratio$wU, wU, tolerance = 1e-12, info = name)
    expect_equal(ratio$wP, wP, tolerance = 1e-10, info = name)

    L <- constrained$constraint$L
    lambda <- solve(L %*% t(L), L %*% state$U_full[seq_len(nb)])
    kkt_score <- c(as.numeric(t(L) %*% lambda), 0, 0)
    expect_true(
      ln_joint_scaled_max(state$U_full, kkt_score) <= thresholds$root,
      info = name
    )
  }
})

test_that("Pace ratio is invariant to equivalent constraints and bases", {
  case <- .ln_ratio_focused_case()
  ratio <- function(name) case$diagnostics[[name]]
  fields <- c("wP", "wU", "qP", "multiplier", "statistic")

  for (field in fields) {
    expect_equal(
      ratio("scalar")[[field]],
      ratio("scalar_scaled")[[field]],
      tolerance = 1e-8,
      info = paste("scalar", field)
    )
    expect_equal(
      ratio("rank2")[[field]],
      ratio("rank2_permuted")[[field]],
      tolerance = 1e-8,
      info = paste("rank2 permutation", field)
    )
    expect_equal(
      ratio("rank2")[[field]],
      ratio("rank2_mixed")[[field]],
      tolerance = 1e-8,
      info = paste("rank2 mixing", field)
    )
    expect_equal(
      ratio("nonintercept")[[field]],
      ratio("nonintercept_scaled")[[field]],
      tolerance = 1e-8,
      info = paste("nonintercept", field)
    )
  }

  scalar <- case$constrained$scalar
  scalar_reference <- ratio("scalar")
  K <- scalar_reference$coordinate_map$null_basis
  K_transform <- matrix(c(1, 0.3, -0.2, 1.15), nrow = 2L,
                        byrow = TRUE)
  alternative_K <- .ln_ratio_replay(
    case,
    scalar,
    null_basis = K %*% K_transform
  )

  rank2 <- case$constrained$rank2
  rank2_reference <- ratio("rank2")
  C <- rank2_reference$coordinate_map$interest_basis
  C_transform <- matrix(c(1.25, 0.2, -0.15, 0.85), nrow = 2L,
                        byrow = TRUE)
  alternative_C <- .ln_ratio_replay(
    case,
    rank2,
    interest_basis = C %*% C_transform
  )

  for (alternative in list(alternative_K, alternative_C)) {
    expect_true(alternative$diagnostic_available)
    expect_false(alternative$public_available)
    expect_identical(
      alternative$public_status,
      "joint_ln_calibration_pending"
    )
  }
  for (field in fields) {
    expect_equal(alternative_K[[field]], scalar_reference[[field]],
                 tolerance = 1e-9, info = paste("alternative K", field))
    expect_equal(alternative_C[[field]], rank2_reference[[field]],
                 tolerance = 1e-9, info = paste("alternative C", field))
  }
})

test_that("Pace ratio attaches across deterministic multi-gene execution", {
  fx <- ln_joint_test_fixture()
  fx$count <- rbind(fx$count, fx$count)
  base <- .ln_ratio_fit(fx)$fit
  p <- ncol(fx$pred)
  beta_public <- as.numeric(base$summary[1L, seq_len(p)])
  L <- matrix(c(0, 1, 0), nrow = 1L)
  b <- as.numeric(L %*% beta_public + 0.012)
  tested <- .ln_ratio_fit(
    fx,
    additional_tests = c("score", "lrt"),
    lrt_contrasts = list(list(name = "serial_pair", L = L, b = b))
  )$fit
  raw <- tested$diagnostics$raw

  expect_length(raw, 2L)
  instances <- vector("list", length(raw))
  for (ii in seq_along(raw)) {
    constrained <- raw[[ii]]$lrt$serial_pair$constrained
    ratio <- constrained$ln_joint$pace_ratio_diagnostic
    expect_identical(raw[[ii]]$route, "R0")
    expect_identical(constrained$route, "R0")
    expect_true(ratio$diagnostic_available)
    expect_identical(
      ratio$objective_instance,
      raw[[ii]]$ln_joint$objective_instance
    )
    expect_identical(
      ratio$objective_instance,
      constrained$ln_joint$objective_instance
    )
    instances[[ii]] <- ratio$objective_instance
  }
  expect_false(identical(instances[[1L]], instances[[2L]]))
})

test_that("Pace ratio audit fails closed for endpoint and matrix corruption", {
  case <- .ln_ratio_focused_case()
  unrestricted <- case$raw
  constrained <- case$constrained$scalar

  no_context_result <- ln_joint_pace_ratio_diagnostic(
    unrestricted,
    constrained
  )
  expect_false(no_context_result$diagnostic_available)
  expect_identical(
    no_context_result$failure_reason,
    "unrestricted_objective_context_missing_or_mismatched"
  )

  non_r0 <- unrestricted
  non_r0$route <- "R2"
  non_r0$inference_available <- FALSE
  non_r0_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = non_r0
  )
  expect_false(non_r0_result$diagnostic_available)
  expect_identical(non_r0_result$failure_reason,
                   "unrestricted_route_not_R0")

  missing_certificate <- unrestricted
  missing_certificate$ln_joint$candidate_records <- lapply(
    missing_certificate$ln_joint$candidate_records,
    function(candidate) {
      candidate$certified <- FALSE
      candidate
    }
  )
  certificate_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = missing_certificate
  )
  expect_false(certificate_result$diagnostic_available)
  expect_identical(
    certificate_result$failure_reason,
    "unrestricted_endpoint_certification_missing"
  )

  malformed_records <- unrestricted
  malformed_records$ln_joint$candidate_records[[1L]] <- 1
  malformed_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = malformed_records
  )
  expect_false(malformed_result$diagnostic_available)
  expect_identical(
    malformed_result$failure_reason,
    "unrestricted_candidate_records_invalid"
  )

  atomic_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = 1
  )
  expect_false(atomic_result$diagnostic_available)
  expect_identical(
    atomic_result$failure_reason,
    "unrestricted_fit_invalid"
  )
  missing_pair_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = list()
  )
  expect_false(missing_pair_result$diagnostic_available)
  expect_identical(
    missing_pair_result$failure_reason,
    "unrestricted_pair_not_registered"
  )
  malformed_nested <- unrestricted
  malformed_nested$ln_joint <- 1
  nested_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = malformed_nested
  )
  expect_false(nested_result$diagnostic_available)
  expect_identical(
    nested_result$failure_reason,
    "unrestricted_joint_state_invalid"
  )

  corrupted_root <- unrestricted
  corrupted_root$ln_joint$state$U_reduced[1L] <-
    corrupted_root$ln_joint$state$U_reduced[1L] + 1
  root_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = corrupted_root
  )
  expect_false(root_result$diagnostic_available)
  expect_identical(
    root_result$failure_reason,
    "unrestricted_selected_root_not_certified"
  )

  corrupted_curvature <- unrestricted
  corrupted_curvature$ln_joint$state$A_reduced[1L, 1L] <-
    -abs(corrupted_curvature$ln_joint$state$A_reduced[1L, 1L])
  curvature_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = corrupted_curvature
  )
  expect_false(curvature_result$diagnostic_available)
  expect_identical(
    curvature_result$failure_reason,
    "unrestricted_selected_curvature_not_certified"
  )

  objective_mismatch <- unrestricted
  matching_record <- which(vapply(
    objective_mismatch$ln_joint$candidate_records,
    function(candidate) {
      isTRUE(candidate$certified) && identical(
        as.numeric(candidate$theta_full_log),
        as.numeric(objective_mismatch$ln_joint$state$full_log)
      )
    },
    logical(1)
  ))
  expect_length(matching_record, 1L)
  matching_objective <-
    objective_mismatch$ln_joint$candidate_records[[matching_record]]$objective
  objective_mismatch$ln_joint$candidate_records[[matching_record]]$
    objective <- matching_objective + 1
  objective_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = objective_mismatch
  )
  expect_false(objective_result$diagnostic_available)
  expect_identical(objective_result$failure_reason,
                   "unrestricted_objective_state_mismatch")

  provenance_mismatch <- unrestricted
  provenance_mismatch$ln_joint$objective_provenance$id <- "wrong_objective"
  provenance_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = provenance_mismatch
  )
  expect_false(provenance_result$diagnostic_available)
  expect_identical(
    provenance_result$failure_reason,
    "unrestricted_objective_provenance_mismatch"
  )

  wrong_gene <- unrestricted
  wrong_gene$debug$gene_index <- wrong_gene$debug$gene_index + 1L
  gene_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = wrong_gene
  )
  expect_false(gene_result$diagnostic_available)
  expect_identical(
    gene_result$failure_reason,
    "unrestricted_objective_context_missing_or_mismatched"
  )

  constrained_non_r0 <- constrained
  constrained_non_r0$route <- "R2"
  constrained_non_r0$inference_available <- FALSE
  constrained_route_result <- .ln_ratio_replay(
    case,
    constrained_non_r0,
    unrestricted_fit = unrestricted
  )
  expect_false(constrained_route_result$diagnostic_available)
  expect_identical(
    constrained_route_result$failure_reason,
    "constrained_route_not_R0"
  )

  constrained_loglik <- constrained
  constrained_loglik$loglik_ln <- constrained_loglik$loglik_ln + 1
  constrained_loglik_result <- .ln_ratio_replay(
    case,
    constrained_loglik,
    unrestricted_fit = unrestricted
  )
  expect_false(constrained_loglik_result$diagnostic_available)
  expect_identical(
    constrained_loglik_result$failure_reason,
    "constrained_loglik_objective_mismatch"
  )

  other_fx <- case$fx
  other_fx$offset[1L] <- other_fx$offset[1L] * 1.01
  other_base <- .ln_ratio_fit(other_fx)
  other_beta <- as.numeric(other_base$fit$summary[1L, seq_len(ncol(
    other_fx$pred
  ))])
  other_L <- matrix(c(1, 0, -0.5), nrow = 1L)
  other_b <- as.numeric(other_L %*% other_beta + 0.02)
  other_tested <- .ln_ratio_fit(
    other_fx,
    additional_tests = c("score", "lrt"),
    lrt_contrasts = list(list(
      name = "other",
      L = other_L,
      b = other_b
    ))
  )$fit
  other_raw <- other_tested$diagnostics$raw[[1L]]
  other_constrained <- other_raw$lrt$other$constrained
  expect_identical(other_raw$route, "R0")
  expect_identical(other_constrained$route, "R0")
  expect_false(identical(
    unrestricted$ln_joint$objective_instance,
    other_constrained$ln_joint$objective_instance
  ))
  instance_result <- .ln_ratio_replay(
    case,
    other_constrained,
    unrestricted_fit = unrestricted
  )
  expect_false(instance_result$diagnostic_available)
  expect_identical(
    instance_result$failure_reason,
    "endpoint_objective_instance_mismatch"
  )

  nonpositive_H <- constrained
  reference_ratio <- case$diagnostics$scalar
  T <- reference_ratio$coordinate_map$full
  T_inv <- solve(T)
  H_coordinate <- crossprod(
    T,
    nonpositive_H$ln_joint$state$A_full %*% T
  )
  H_coordinate[1L, 1L] <- -abs(H_coordinate[1L, 1L])
  nonpositive_H$ln_joint$state$A_full <-
    t(T_inv) %*% H_coordinate %*% T_inv
  H_result <- .ln_ratio_replay(
    case,
    nonpositive_H,
    unrestricted_fit = unrestricted
  )
  expect_false(H_result$diagnostic_available)
  expect_identical(
    H_result$failure_reason,
    "constrained_derivative_state_reevaluation_mismatch"
  )

  reported_objective_corruption <- unrestricted
  constrained_objective <- constrained$ln_joint$state$raw$value
  objective_gap <- constrained_objective -
    unrestricted$ln_joint$state$raw$value
  expect_gt(objective_gap, 0)
  corrupted_objective <- unrestricted$ln_joint$state$raw$value +
    objective_gap / 4
  reported_objective_corruption$ln_joint$state$raw$value <-
    corrupted_objective
  reported_objective_corruption$loglik <-
    reported_objective_corruption$loglik_ln <- -corrupted_objective
  matching_record <- which(vapply(
    reported_objective_corruption$ln_joint$candidate_records,
    function(candidate) {
      isTRUE(candidate$certified) && identical(
        as.numeric(candidate$theta_full_log),
        as.numeric(reported_objective_corruption$ln_joint$state$full_log)
      )
    },
    logical(1)
  ))
  expect_length(matching_record, 1L)
  corrupted_records <-
    reported_objective_corruption$ln_joint$candidate_records
  corrupted_records[[matching_record]]$objective <- corrupted_objective
  reported_objective_corruption$ln_joint$candidate_records <-
    corrupted_records
  reported_objective_result <- .ln_ratio_replay(
    case,
    constrained,
    unrestricted_fit = reported_objective_corruption
  )
  expect_false(reported_objective_result$diagnostic_available)
  expect_identical(
    reported_objective_result$failure_reason,
    "unrestricted_objective_reevaluation_mismatch"
  )

  for (result in list(
    no_context_result,
    non_r0_result,
    certificate_result,
    malformed_result,
    atomic_result,
    missing_pair_result,
    nested_result,
    root_result,
    curvature_result,
    objective_result,
    provenance_result,
    gene_result,
    constrained_route_result,
    constrained_loglik_result,
    instance_result,
    H_result,
    reported_objective_result
  )) {
    expect_false(result$public_available)
    expect_identical(result$public_status, "joint_ln_calibration_pending")
    expect_true(is.na(result$statistic))
  }
})

test_that("audit-only Pace ratio leaves public Score LRT and Wald fields unchanged", {
  case <- .ln_ratio_focused_case()
  baseline <- case$base$fit
  tested <- case$tested$fit
  raw_baseline <- baseline$diagnostics$raw[[1L]]
  raw_tested <- tested$diagnostics$raw[[1L]]

  expect_true(case$tested$fallback_setting %in%
                c("none", "absent_at_ce623d1"))
  expect_identical(raw_tested$route, "R0")
  expect_identical(tested$provenance$per_gene$route, "R0")
  primary_summary_columns <- setdiff(
    names(baseline$summary),
    c("test_status", "test_failure_detail")
  )
  expect_equal(
    tested$summary[, primary_summary_columns, drop = FALSE],
    baseline$summary[, primary_summary_columns, drop = FALSE],
    tolerance = 1e-10
  )
  expect_equal(raw_tested$theta, raw_baseline$theta, tolerance = 0)
  expect_equal(raw_tested$naive_cov, raw_baseline$naive_cov, tolerance = 0)
  expect_equal(raw_tested$robust_cov, raw_baseline$robust_cov, tolerance = 0)
  expect_equal(raw_tested$score, raw_baseline$score, tolerance = 0)

  for (name in names(case$definitions)) {
    entry <- raw_tested$lrt[[name]]
    expect_true(is.list(entry$raw), info = name)
    expect_true(is.finite(entry$raw$wP), info = name)
    expect_true(all(is.na(c(entry$raw$wU, entry$raw$qP,
                            entry$raw$wP_adjusted))), info = name)
    expect_identical(entry$raw$method, "unavailable_joint_ln")
    expect_false(entry$raw$success, info = name)
    expect_false(entry$raw$adjusted_available, info = name)
    expect_true(is.na(entry$stat), info = name)
    expect_true(is.na(entry$pval), info = name)
    expect_true(is.finite(entry$wP), info = name)
    expect_true(all(is.na(c(entry$wU, entry$qP))), info = name)
    expect_false(entry$success, info = name)
    expect_false(entry$score_success, info = name)
    expect_identical(entry$failure_reason, "joint_ln_calibration_pending")
    expect_identical(entry$score_failure,
                     "joint_ln_calibration_pending")
    expect_identical(entry$method, "unavailable_joint_ln")
  }

  helper_text <- tolower(paste(
    deparse(body(ln_joint_pace_ratio_diagnostic)),
    collapse = "\n"
  ))
  for (token in c(
    "pchisq", "compute_profile_lrt", "compute_adj_score",
    "pseudoinverse", "eigenvalue_clipping", "wald_fallback"
  )) {
    expect_false(grepl(token, helper_text, fixed = TRUE), info = token)
  }
})

test_that("Pace ratio attachment survives registry eviction in a gene batch", {
  fx <- ln_joint_test_fixture()
  fx$count <- fx$count[rep(1L, 18L), , drop = FALSE]
  reference <- .ln_ratio_focused_case()$base$fit
  p <- ncol(fx$pred)
  beta_public <- as.numeric(reference$summary[1L, seq_len(p)])
  L <- matrix(c(0, 1, 0), nrow = 1L)
  b <- as.numeric(L %*% beta_public + 0.012)
  tested <- .ln_ratio_fit(
    fx,
    additional_tests = c("score", "lrt"),
    lrt_contrasts = list(list(name = "eviction", L = L, b = b))
  )$fit
  raw <- tested$diagnostics$raw

  expect_length(raw, 18L)
  for (ii in seq_along(raw)) {
    constrained <- raw[[ii]]$lrt$eviction$constrained
    ratio <- constrained$ln_joint$pace_ratio_diagnostic
    expect_identical(raw[[ii]]$route, "R0", info = ii)
    expect_identical(constrained$route, "R0", info = ii)
    expect_true(ratio$diagnostic_available, info = ii)
    expect_identical(
      ratio$objective_instance,
      raw[[ii]]$ln_joint$objective_instance,
      info = ii
    )
  }
})
