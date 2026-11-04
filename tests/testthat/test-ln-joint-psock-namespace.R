ln_joint_psock_serial_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- ln_joint_fit_fixture(
        additional_tests = c("score", "lrt"),
        lrt_contrasts = list(list(name = "psock_x", L = matrix(c(1, 0, 0), 1L), b = 0)),
        ncore = 1L
      )$fit
    }
    cached
  }
})

test_that("joint-LN dispatchers expose their fitter dependencies", {
  namespace <- asNamespace("scRoPE")
  expected <- c(
    fit_gene_unconstrained = "fit_gene_joint_ln_unconstrained",
    fit_gene_constrained = "fit_gene_joint_ln_constrained"
  )

  for (name in names(expected)) {
    dispatcher <- get(name, envir = namespace, inherits = FALSE)
    symbols <- all.names(
      body(dispatcher),
      functions = TRUE,
      unique = TRUE
    )
    expect_true(expected[[name]] %in% symbols)
    expect_false("getFromNamespace" %in% symbols)
  }
})

test_that("PSOCK diagnostic comparison tolerates only numerical text differences", {
  baseline <- list(label = "S_direct", required = TRUE, route = "R0",
                   theta_reduced = "0.1;2.6", missing = NA_character_,
                   candidate_key = "joint_ln_candidate_v1|TRUE|TRUE|TRUE|TRUE|263.6|0.1;2.6")
  close <- baseline
  close$theta_reduced <- "0.100000000000001;2.60000000000001"
  close$candidate_key <- "joint_ln_candidate_v1|TRUE|TRUE|TRUE|TRUE|263.6000000000001|0.100000000000001;2.60000000000001"
  expect_length(ln_joint_psock_differences(close, baseline), 0L)
  for (field in c("theta_reduced", "candidate_key", "label", "required", "route")) {
    changed <- baseline
    changed[[field]] <- switch(field,
      theta_reduced = "0.11;2.6",
      candidate_key = "joint_ln_candidate_v1|FALSE|TRUE|TRUE|TRUE|263.6|0.1;2.6",
      label = "S_mid", required = FALSE, route = "R2")
    expect_gt(length(ln_joint_psock_differences(changed, baseline)), 0L)
  }
  changed <- baseline
  changed$candidate_key <- "joint_ln_candidate_v1|TRUE|TRUE|TRUE|TRUE|263.7|0.1;2.6"
  expect_gt(length(ln_joint_psock_differences(changed, baseline)), 0L)
  changed <- baseline
  changed$theta_reduced <- "2.6;0.1"
  expect_gt(length(ln_joint_psock_differences(changed, baseline)), 0L)
  changed$theta_reduced <- c(first = "0.1;2.6")
  expect_gt(length(ln_joint_psock_differences(changed, baseline)), 0L)
  changed$theta_reduced <- c(0.1, 2.6)
  expect_gt(length(ln_joint_psock_differences(changed, baseline)), 0L)
  expect_error(ln_joint_psock_parse_numbers("0.1;typo"), "Malformed")
  expect_error(ln_joint_psock_parse_numbers("0.1;"), "Malformed")
  expect_error(ln_joint_psock_parse_key("joint_ln_candidate_v1|TRUE|0.1"), "Malformed")
  expect_error(ln_joint_psock_parse_key(
    "joint_ln_candidate_v1|TRUE|MAYBE|TRUE|TRUE|263.6|0.1"), "Malformed")
  expect_gt(length(ln_joint_psock_differences(
    list(theta_reduced = NA_character_), list(theta_reduced = "NA"))), 0L)
  expect_gt(length(ln_joint_psock_differences(
    list(theta_reduced = "<missing>"), list(theta_reduced = "<malformed>"))), 0L)
  expect_gt(length(ln_joint_psock_differences(
    list(theta_reduced = "NA"), list(theta_reduced = "NaN"))), 0L)
})

test_that("PSOCK absolute floor applies only to small diagnostic roundoff", {
  fields <- c("score", "g_red", "q_gradient_reduced", "U_reduced",
    "reduced_gradient_max_abs", "max_abs_newton_step_scaled",
    "root_residual", "curvature_compatibility", "nuisance_root_scaled")
  baseline <- c(0, 3.7e-7)
  roundoff <- baseline + 2e-13
  for (field in fields) {
    expect_length(ln_joint_psock_differences(roundoff, baseline, field = field), 0L)
    expect_gt(length(ln_joint_psock_differences(baseline + 2e-11, baseline,
      field = field)), 0L)
  }
  # Coefficients and p-values retain their previous relative-only comparison.
  for (field in c("beta", "p_value")) {
    expect_gt(length(ln_joint_psock_differences(roundoff, baseline,
      field = field)), 0L)
  }
  expect_length(ln_joint_psock_differences(1e6 + 1e-3, 1e6, field = "score"), 0L)
  expect_gt(length(ln_joint_psock_differences(1e6 + 1, 1e6, field = "score")), 0L)
  for (pair in list(list(Inf, -Inf), list(NA_real_, NaN),
                   list(c(0, NA_real_), c(NA_real_, 0)))) {
    expect_gt(length(ln_joint_psock_differences(pair[[1L]], pair[[2L]],
      field = "score")), 0L)
  }
  expect_gt(length(ln_joint_psock_differences(0, 0L, field = "score")), 0L)
  expect_gt(length(ln_joint_psock_differences(c(a = 0), c(b = 0),
    field = "score")), 0L)
  expect_gt(length(ln_joint_psock_differences(list(score = 0, route = "R0"),
    list(score = 0, route = "R2"))), 0L)
})

test_that("joint-LN diagnostic keys are exact identities within one fit", {
  raw <- ln_joint_psock_serial_fixture()$diagnostics$raw
  errors <- ln_joint_psock_key_errors(raw)
  expect_identical(length(errors), 0L, info = paste(errors, collapse = "\n"))
  before <- serialize(raw, NULL)
  expect_length(ln_joint_psock_differences(raw, raw), 0L)
  expect_identical(serialize(raw, NULL), before)
  changed <- raw
  record <- changed[[1L]]$ln_joint$candidate_records[[1L]]
  record$objective <- record$objective + 1e-12
  changed[[1L]]$ln_joint$candidate_records[[1L]] <- record
  expect_gt(length(ln_joint_psock_key_errors(changed)), 0L)
  changed <- raw
  changed[[1L]]$ln_joint$resolution$candidate_resolutions[[1L]]$candidate_key <- "wrong key"
  expect_gt(length(ln_joint_psock_key_errors(changed)), 0L)
})

test_that("Pace witnesses point to the first eligible candidate in their own fit", {
  raw <- ln_joint_psock_serial_fixture()$diagnostics$raw
  before <- serialize(raw, NULL)
  audit <- ln_joint_psock_witness_audit(raw)
  expect_length(audit$errors, 0L)
  expect_length(audit$witnesses, 6L)
  expect_true(all(vapply(audit$witnesses, function(x) x$checked, logical(1))))
  expect_identical(serialize(raw, NULL), before)

  owner <- raw[[1L]]
  endpoint <- owner$lrt$x$constrained$ln_joint$pace_ratio_diagnostic$endpoints$unrestricted
  candidate <- owner$ln_joint$candidate_records[[endpoint$certified_candidate_index]]
  owner$ln_joint$candidate_records <- list(candidate, candidate)
  endpoint$certified_candidate_index <- 1L
  endpoint$candidate_objective <- as.numeric(candidate$objective)
  valid <- ln_joint_psock_endpoint_witness(endpoint, owner)
  expect_length(valid$errors, 0L)
  expect_identical(valid$eligible, 1:2)

  # Two numerically equivalent fits can have different first exact witnesses.
  # Within-fit matching stays exact, including differences as small as 1e-14.
  other_owner <- owner
  other_owner$ln_joint$candidate_records[[1L]]$theta_full_log[[1L]] <-
    candidate$theta_full_log[[1L]] + 1e-14
  other_endpoint <- endpoint
  other_endpoint$certified_candidate_index <- 2L
  other <- ln_joint_psock_endpoint_witness(other_endpoint, other_owner)
  expect_length(other$errors, 0L)
  expect_identical(other$eligible, 2L)
  minimal_raw <- function(u, e) {
    u$lrt <- list(x = list(constrained = list(ln_joint = list(
      pace_ratio_diagnostic = list(endpoints = list(unrestricted = e))))))
    # Add a valid constrained owner with the same synthetic state. This tests
    # the traversal's owner linkage without invoking another optimization.
    c <- u
    c$lrt <- NULL
    c$debug$method <- "joint_ln_constrained"
    constrained_endpoint <- e
    constrained_endpoint$endpoint <- "constrained"
    c$ln_joint$pace_ratio_diagnostic <- list(endpoints = list(
      unrestricted = e, constrained = constrained_endpoint))
    u$lrt$x$constrained <- c
    list(u)
  }
  a <- ln_joint_psock_witness_audit(minimal_raw(owner, endpoint))
  b <- ln_joint_psock_witness_audit(minimal_raw(other_owner, other_endpoint))
  expect_length(a$errors, 0L)
  expect_length(b$errors, 0L)
  expect_length(ln_joint_psock_differences(a$normalized, b$normalized), 0L)
  expect_gt(length(ln_joint_psock_differences(
    minimal_raw(owner, endpoint), minimal_raw(other_owner, other_endpoint))), 0L)
  expect_gt(length(ln_joint_psock_endpoint_witness(
    endpoint, owner, "constrained")$errors), 0L)
  for (field in c("method", "route", "inference_available")) {
    changed <- owner
    if (field == "method") changed$debug$method <- "wrong_method"
    if (field == "route") changed$route <- "R2"
    if (field == "inference_available") changed$inference_available <- FALSE
    expect_gt(length(ln_joint_psock_endpoint_witness(endpoint, changed)$errors), 0L)
  }

  for (index in list(0L, 3L, NA_integer_, c(1L, 2L), 1, 2L)) {
    changed <- endpoint
    changed$certified_candidate_index <- index
    expect_gt(length(ln_joint_psock_endpoint_witness(changed, owner)$errors), 0L)
  }
  # Every corruption below invalidates the *reported* first witness even if
  # another candidate remains eligible. Neither a tolerance nor normalization
  # is allowed to conceal it.
  for (field in c("certified", "interior", "theta_full_log",
                  "objective_provenance", "objective")) {
    changed <- owner
    changed$ln_joint$candidate_records[[1L]][[field]] <- switch(field,
      certified = FALSE, interior = FALSE,
      theta_full_log = candidate$theta_full_log + 1e-14,
      objective_provenance = list(wrong = TRUE),
      objective = candidate$objective + 1)
    expect_gt(length(ln_joint_psock_endpoint_witness(endpoint, changed)$errors), 0L)
    invalid <- ln_joint_psock_witness_audit(minimal_raw(changed, endpoint))
    expect_gt(length(invalid$errors), 0L)
    expect_identical(invalid$normalized[[1L]]$lrt$x$constrained$ln_joint$
      pace_ratio_diagnostic$endpoints$unrestricted$certified_candidate_index, 1L)
  }
  for (field in c("full_log", "candidate_objective", "stored_state_objective",
                  "state_objective", "objective_provenance", "objective_instance",
                  "gene_index")) {
    changed <- endpoint
    changed[[field]] <- switch(field,
      full_log = endpoint$full_log + 1e-14,
      candidate_objective = endpoint$candidate_objective + 1e-12,
      stored_state_objective = endpoint$stored_state_objective + 1e-12,
      state_objective = endpoint$state_objective + 1,
      objective_provenance = list(wrong = TRUE),
      objective_instance = list(wrong = TRUE), gene_index = 2L)
    expect_gt(length(ln_joint_psock_endpoint_witness(changed, owner)$errors), 0L)
  }
})

test_that("PSOCK receipts preserve original synthetic pairs without private metadata", {
  directory <- withr::local_tempdir()
  withr::local_envvar(SCROPE_TEST_DIAGNOSTICS_DIR = directory)
  fit <- ln_joint_psock_serial_fixture()
  before <- serialize(fit, NULL)
  ln_joint_save_psock_pair("optional-tests", fit)
  expect_null(readRDS(file.path(directory, "psock-pair-optional-tests.rds"))$parallel)
  # The receipt is already durable if any later expectation raises an error.
  expect_error({
    ln_joint_save_psock_pair("optional-tests", fit, fit)
    stop("subsequent assertion failure")
  }, "subsequent assertion failure")
  receipt <- readRDS(file.path(directory, "psock-pair-optional-tests.rds"))
  expect_identical(receipt$serial$raw, fit$diagnostics$raw)
  expect_identical(receipt$parallel$raw, fit$diagnostics$raw)
  expect_identical(serialize(fit, NULL), before)
  expect_true(ln_joint_psock_receipt_safe(receipt))
  expect_false(ln_joint_psock_receipt_safe(list(path = "/home/example/private")))
  expect_false(ln_joint_psock_receipt_safe(list(path = "C:\\Users\\example\\private")))
  expect_false(ln_joint_psock_receipt_safe(list(environment = new.env())))
  expect_false(ln_joint_psock_receipt_safe(list(callback = function() NULL)))
})

test_that("owner-aware PSOCK contract preserves scientific and metadata failures", {
  baseline <- ln_joint_psock_serial_fixture()
  before <- serialize(baseline, NULL)
  expect_length(ln_joint_psock_optional_fit_errors(baseline, baseline), 0L)
  expect_identical(serialize(baseline, NULL), before)
  close <- baseline
  close$summary$score_adj_z <- close$summary$score_adj_z + 1e-10
  expect_length(ln_joint_psock_optional_fit_errors(baseline, close), 0L)
  for (column in c("score_adj_z", "lrt_adj_z", "p_score_adj_z", "p_lrt_adj_z",
                   "logFC_x", "se_robust_x")) {
    changed <- baseline
    changed$summary[[column]] <- changed$summary[[column]] + 1e-6
    expect_gt(length(ln_joint_psock_optional_fit_errors(baseline, changed)), 0L)
  }
  for (field in c("schema_version", "candidate_resolution_policy")) {
    changed <- baseline
    changed$provenance[[field]] <- if (field == "schema_version") 4L else "wrong_policy"
    expect_gt(length(ln_joint_psock_optional_fit_errors(baseline, changed)), 0L)
  }
  for (field in c("route", "inference_available", "failure_reason")) {
    changed <- baseline
    changed$diagnostics$raw[[1L]]$lrt$z$constrained[[field]] <- switch(field,
      route = "R2", inference_available = FALSE, failure_reason = "wrong_reason")
    expect_gt(length(ln_joint_psock_optional_fit_errors(baseline, changed)), 0L)
  }
  for (field in c("H_red_rcond", "reduced_gradient_max_abs", "max_abs_newton_step_scaled")) {
    changed <- baseline
    changed$lrt_details[[field]][[2L]] <- changed$lrt_details[[field]][[2L]] + 1e-6
    expect_gt(length(ln_joint_psock_optional_fit_errors(baseline, changed)), 0L)
  }
  owner <- baseline$diagnostics$raw[[1L]]$lrt$z$constrained
  for (field in c("root_residual", "newton_displacement", "certified", "polished")) {
    changed <- owner
    changed$ln_joint$candidate_records[[1L]][[field]] <- switch(field,
      root_residual = 1e-3, newton_displacement = 1e-3,
      certified = FALSE, polished = !changed$ln_joint$candidate_records[[1L]]$polished)
    expect_gt(length(ln_joint_psock_owner_audit(changed)), 0L)
  }
  changed <- owner
  changed$ln_joint$state$full_log[[1L]] <- changed$ln_joint$state$full_log[[1L]] + 1e-14
  expect_gt(length(ln_joint_psock_owner_audit(changed)), 0L)
  changed <- owner
  changed$ln_joint$state$A_reduced <- -diag(nrow(changed$ln_joint$state$A_reduced))
  expect_gt(length(ln_joint_psock_owner_audit(changed)), 0L)
  changed <- owner
  changed$ln_joint$state$U_reduced[] <- 1
  expect_gt(length(ln_joint_psock_owner_audit(changed)), 0L)
  changed <- owner
  changed$ln_joint$candidate_records[[1L]]$objective_instance$id <- "wrong_owner"
  expect_gt(length(ln_joint_psock_owner_audit(changed)), 0L)
  changed <- baseline
  changed$diagnostics$raw[[1L]]$lrt$z$constrained$ln_joint$
    pace_ratio_diagnostic$endpoints$constrained$certified_candidate_index <- 0L
  expect_gt(length(ln_joint_psock_optional_fit_errors(baseline, changed)), 0L)
  changed <- owner
  changed$ln_joint$candidate_records[[1L]]$resolution$candidate_key <- "wrong key"
  expect_gt(length(ln_joint_psock_owner_audit(changed)), 0L)
  for (field in c("godambe", "raw_gradient", "reduced_information", "score_statistic")) {
    changed <- owner
    if (field == "godambe") changed$godambe$S[[1L]] <- changed$godambe$S[[1L]] + 1e-6
    if (field == "raw_gradient") changed$ln_joint$state$raw$gradient[[1L]] <-
      changed$ln_joint$state$raw$gradient[[1L]] + 1e-6
    if (field == "reduced_information") changed$reduced$H_red[[1L]] <-
      changed$reduced$H_red[[1L]] + 1e-6
    if (field == "score_statistic") changed$ln_joint$score_diagnostic$statistic <-
      changed$ln_joint$score_diagnostic$statistic + 1e-6
    expect_gt(length(ln_joint_psock_owner_audit(changed)), 0L)
  }
  changed <- baseline
  changed$diagnostics$raw[[1L]]$lrt$z$constrained$ln_joint$pace_ratio_diagnostic$wP <-
    changed$diagnostics$raw[[1L]]$lrt$z$constrained$ln_joint$pace_ratio_diagnostic$wP + 1e-6
  expect_gt(length(ln_joint_psock_optional_fit_errors(baseline, changed)), 0L)

  # Different, internally truthful polishing flags are allowed only after both
  # runs pass their own audits. Corresponding candidate endpoints still agree.
  owner_a <- owner_b <- owner
  for (which in c("a", "b")) {
    u <- if (which == "a") owner_a else owner_b
    r <- u$ln_joint$candidate_records[[1L]]
    r$optimizer_endpoint <- r$theta_reduced
    if (which == "b") r$optimizer_endpoint[[1L]] <- r$optimizer_endpoint[[1L]] + 1e-15
    r$polished <- which == "b"
    u$ln_joint$candidate_records[[1L]] <- r
    reconstructed <- lapply(u$ln_joint$candidate_records, function(record) {
      record$optimizer <- list(convergence = record$optimizer_convergence,
                               message = record$optimizer_message)
      record
    })
    u$ln_joint$candidates <- ln_joint_candidate_summary(reconstructed,
      u$ln_joint$resolution$candidate_resolutions)
    if (which == "a") owner_a <- u else owner_b <- u
  }
  expect_length(ln_joint_psock_owner_audit(owner_a), 0L)
  expect_length(ln_joint_psock_owner_audit(owner_b), 0L)
  a <- b <- baseline
  a$diagnostics$raw[[1L]]$lrt$z$constrained <- owner_a
  b$diagnostics$raw[[1L]]$lrt$z$constrained <- owner_b
  expect_length(ln_joint_psock_optional_fit_errors(a, b), 0L)
  changed <- b
  changed$diagnostics$raw[[1L]]$lrt$z$constrained$ln_joint$candidates$polished[[1L]] <- FALSE
  expect_gt(length(ln_joint_psock_optional_fit_errors(a, changed)), 0L)
})

test_that("optional mixed budget preserves boundaries, types and missingness", {
  expect_length(ln_joint_psock_optional_differences(0, 5e-9), 0L)
  expect_gt(length(ln_joint_psock_optional_differences(0, 2e-8)), 0L)
  expect_length(ln_joint_psock_optional_differences(1e6, 1e6 + 5e-3), 0L)
  expect_gt(length(ln_joint_psock_optional_differences(1e6, 1e6 + 2e-2)), 0L)
  for (pair in list(list(NA_real_, NaN), list(Inf, -Inf), list(0, NA_real_),
                   list(0, 0L), list(c(a = 0), c(b = 0)),
                   list(TRUE, FALSE), list("R0", "R2"))) {
    expect_gt(length(ln_joint_psock_optional_differences(pair[[1L]], pair[[2L]])), 0L)
  }
})

test_that("installed joint-LN PSOCK fits use the running package", {
  package_path <- normalizePath(
    getNamespaceInfo(asNamespace("scRoPE"), "path"),
    mustWork = FALSE
  )
  skip_if_not(
    file.exists(file.path(package_path, "Meta", "package.rds")),
    "requires an installed package"
  )
  skip_if(parallelly::availableCores() < 2L)

  contrast <- list(list(
    name = "psock_x",
    L = matrix(c(1, 0, 0), nrow = 1L),
    b = 0
  ))
  one <- ln_joint_psock_serial_fixture()
  ln_joint_save_psock_pair("optional-tests", one)

  original_libraries <- .libPaths()
  on.exit(.libPaths(original_libraries), add = TRUE)
  package_library <- normalizePath(dirname(package_path), mustWork = FALSE)
  other_libraries <- original_libraries[
    normalizePath(original_libraries, mustWork = FALSE) != package_library
  ]
  alternate_package <- tryCatch(
    find.package("scRoPE", lib.loc = other_libraries),
    error = function(e) character(0)
  )
  if (length(alternate_package)) {
    # Reproduce a lib.loc install while an older default scRoPE is visible.
    .libPaths(other_libraries)
  }

  two <- ln_joint_fit_fixture(
    additional_tests = c("score", "lrt"),
    lrt_contrasts = contrast,
    ncore = 2L
  )$fit
  # Persist the original, unnormalized fits before *any* comparison can fail.
  ln_joint_save_psock_pair("optional-tests", one, two)

  errors <- ln_joint_psock_optional_fit_errors(one, two)
  expect_identical(length(errors), 0L, info = paste(errors, collapse = "\n"))
  expect_identical(one$provenance$schema_version, 3L)
  expect_identical(
    one$provenance$candidate_resolution_policy,
    "joint_ln_candidate_resolution_calibration_v1"
  )
  expect_identical(
    one$provenance$candidate_resolution_mode,
    "calibration_candidate_enabled"
  )
  expect_true(all(
    one$provenance$per_gene$resolution_mode ==
      "calibration_candidate_enabled"
  ))
  expect_true(all(
    two$provenance$per_gene$resolution_mode ==
      "calibration_candidate_enabled"
  ))
})
