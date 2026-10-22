ln_joint_normalise_psock_objective_ids <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  object_names <- names(x)
  for (index in seq_along(x)) {
    value <- x[[index]]
    field <- if (is.null(object_names)) "" else object_names[[index]]
    if (identical(field, "objective_instance") && is.list(value) &&
        identical(
          value$scheme,
          "process_local_exact_joint_ln_objective_context_v1"
        )) {
      value$id <- paste0("<run-local-id-gene-", value$gene_index, ">")
    }
    x[index] <- list(ln_joint_normalise_psock_objective_ids(value))
  }
  x
}

# Diagnostic text serializes floating-point values at 17 digits. Numerical
# tolerance must apply to those values, not to the characters that encode them.
# Runtime keys remain exact identities within each fit and are never modified.
ln_joint_psock_parse_numbers <- function(text) {
  if (is.na(text)) return(list(kind = "missing_character"))
  if (text %in% c("<missing>", "<malformed>")) {
    return(list(kind = "sentinel", value = text))
  }
  if (!nzchar(text)) return(list(kind = "numeric", value = numeric(0)))
  tokens <- trimws(strsplit(paste0(text, ";<end>"), ";", fixed = TRUE)[[1L]])
  tokens <- tokens[-length(tokens)]
  number <- "^[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?$"
  valid <- grepl(number, tokens) | tokens %in% c("NA", "NaN", "Inf", "-Inf")
  if (!all(valid)) stop("Malformed numerical diagnostic text", call. = FALSE)
  value <- rep(NA_real_, length(tokens))
  value[tokens != "NA"] <- as.numeric(tokens[tokens != "NA"])
  list(kind = "numeric", value = value)
}

ln_joint_psock_parse_key <- function(text) {
  if (is.na(text)) return(list(kind = "missing_character"))
  pieces <- strsplit(text, "|", fixed = TRUE)[[1L]]
  if (identical(pieces, c("joint_ln_candidate_v1", "malformed_nonlist_candidate"))) {
    return(list(kind = "sentinel", value = text))
  }
  if (length(pieces) != 7L || pieces[[1L]] != "joint_ln_candidate_v1" ||
      !all(pieces[2:5] %in% c("TRUE", "FALSE"))) {
    stop("Malformed candidate key", call. = FALSE)
  }
  list(kind = pieces[[1L]], flags = pieces[2:5],
       objective = ln_joint_psock_parse_numbers(pieces[[6L]]),
       theta = ln_joint_psock_parse_numbers(pieces[[7L]]))
}

ln_joint_psock_differences <- function(actual, expected, tolerance = 1e-8,
                                      field = "", path = "raw") {
  fail <- function(reason) paste0(path, ": ", reason)
  if (!identical(typeof(actual), typeof(expected)) ||
      !identical(attributes(actual), attributes(expected))) {
    return(fail("type or attributes differ"))
  }
  if (length(actual) != length(expected)) return(fail("length differs"))
  if (is.list(actual)) {
    return(unlist(lapply(seq_along(actual), function(i) {
      name <- if (is.null(names(actual))) "" else names(actual)[[i]]
      ln_joint_psock_differences(actual[[i]], expected[[i]], tolerance,
                                name, paste0(path, "[[", i, "]]$", name))
    }), use.names = FALSE))
  }
  if (is.character(actual) && field %in%
      c("start", "optimizer_endpoint", "theta_reduced", "theta_full_log",
        "candidate_key", "anchor_candidate_key")) {
    parser <- if (grepl("candidate_key$", field)) {
      ln_joint_psock_parse_key
    } else ln_joint_psock_parse_numbers
    return(unlist(lapply(seq_along(actual), function(i) {
      ln_joint_psock_differences(parser(actual[[i]]), parser(expected[[i]]),
                                tolerance, path = paste0(path, "[", i, "]"))
    }), use.names = FALSE))
  }
  if (is.double(actual)) {
    # Preserve missing/nonfinite locations and kinds; finite values use the
    # same tolerance as the pre-existing numerical PSOCK comparison.
    if (!identical(is.na(actual), is.na(expected)) ||
        !identical(is.nan(actual), is.nan(expected)) ||
        !identical(is.finite(actual), is.finite(expected))) {
      return(fail("missing or nonfinite states differ"))
    }
    comparison <- all.equal(actual, expected, tolerance = tolerance)
    if (!isTRUE(comparison)) return(fail(paste(comparison, collapse = "; ")))
  } else if (!identical(actual, expected)) {
    return(fail("exact values differ"))
  }
  character(0)
}

ln_joint_psock_key_errors <- function(x) {
  if (!is.list(x)) return(character(0))
  errors <- character(0)
  if (all(c("candidate_records", "candidates", "resolution") %in% names(x))) {
    records <- x$candidate_records
    resolved <- x$resolution$candidate_resolutions
    if (length(records) != length(resolved) || nrow(x$candidates) != length(records)) {
      return("candidate record counts differ")
    }
    for (i in seq_along(records)) {
      key <- ln_joint_resolution_candidate_key(records[[i]])
      if (!identical(records[[i]]$resolution$candidate_key, key) ||
          !identical(resolved[[i]]$candidate_key, key)) {
        errors <- c(errors, paste("candidate", i, "has an inconsistent exact key"))
      }
      for (field in c("start", "optimizer_endpoint", "theta_reduced", "theta_full_log")) {
        if (!identical(x$candidates[[field]][[i]],
                       ln_joint_diagnostic_vector_string(records[[i]][[field]]))) {
          errors <- c(errors, paste("candidate", i, field, "text differs from its record"))
        }
      }
    }
  }
  c(errors, unlist(lapply(x, ln_joint_psock_key_errors), use.names = FALSE))
}

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

  expect_equal(one$summary, two$summary, tolerance = 1e-8)
  expect_equal(one$overdispersion, two$overdispersion, tolerance = 1e-8)
  expect_identical(one$convergence, two$convergence)
  expect_equal(
    one$provenance$per_gene,
    two$provenance$per_gene,
    tolerance = 1e-8
  )
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
  expect_equal(one$lrt_details, two$lrt_details, tolerance = 1e-8)
  for (raw in list(one$diagnostics$raw, two$diagnostics$raw)) {
    errors <- ln_joint_psock_key_errors(raw)
    expect_identical(length(errors), 0L, info = paste(errors, collapse = "\n"))
  }
  one_audit <- ln_joint_psock_witness_audit(one$diagnostics$raw)
  two_audit <- ln_joint_psock_witness_audit(two$diagnostics$raw)
  for (audit in list(one_audit, two_audit)) {
    expect_identical(length(audit$errors), 0L, info = paste(audit$errors, collapse = "\n"))
  }
  differences <- ln_joint_psock_differences(
    ln_joint_normalise_psock_objective_ids(one_audit$normalized),
    ln_joint_normalise_psock_objective_ids(two_audit$normalized),
    tolerance = 1e-8
  )
  expect_identical(length(differences), 0L, info = paste(differences, collapse = "\n"))
})
