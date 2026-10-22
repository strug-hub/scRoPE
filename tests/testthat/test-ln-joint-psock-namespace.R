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
  differences <- ln_joint_psock_differences(
    ln_joint_normalise_psock_objective_ids(one$diagnostics$raw),
    ln_joint_normalise_psock_objective_ids(two$diagnostics$raw),
    tolerance = 1e-8
  )
  expect_identical(length(differences), 0L, info = paste(differences, collapse = "\n"))
})
