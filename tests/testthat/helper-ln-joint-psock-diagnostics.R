# These checks inspect saved fits; they never refit, select, or certify a
# candidate. A Pace witness index is local to its owner's candidate list.
ln_joint_psock_endpoint_witness <- function(endpoint, owner,
                                           expected_endpoint = endpoint$endpoint) {
  errors <- character(0)
  fail <- function(message) errors <<- c(errors, message)
  index <- endpoint$certified_candidate_index
  if (identical(index, NA_integer_) && !isTRUE(endpoint$success)) {
    return(list(errors = errors, checked = FALSE, index = index,
                eligible = integer(0), candidates = NULL))
  }
  state <- owner$ln_joint$state
  records <- owner$ln_joint$candidate_records
  provenance <- ln_joint_objective_provenance()
  threshold <- ln_joint_thresholds()$objective_agreement
  if (!identical(endpoint$endpoint, expected_endpoint)) {
    fail("endpoint kind differs from its owner role")
  }
  expected_method <- if (identical(expected_endpoint, "unrestricted")) {
    "joint_ln"
  } else "joint_ln_constrained"
  if (!identical(owner$debug$method, expected_method) ||
      !identical(owner$route, "R0") || !isTRUE(owner$inference_available) ||
      !identical(endpoint$route, owner$route) ||
      !identical(endpoint$inference_available, owner$inference_available)) {
    fail("endpoint owner method, route or inference availability differs")
  }
  if (!is.list(state) || !is.list(records) || !length(records) ||
      any(!vapply(records, is.list, logical(1)))) {
    return(list(errors = "missing owner state or candidate records",
                checked = FALSE, index = index, eligible = integer(0),
                candidates = NULL))
  }
  if (!identical(endpoint$full_log, as.numeric(state$full_log))) {
    fail("endpoint coordinates do not match their owner")
  }
  if (!identical(endpoint$objective_provenance, provenance) ||
      !identical(state$objective_provenance, provenance) ||
      !identical(owner$ln_joint$objective_provenance, provenance)) {
    fail("endpoint/owner objective provenance differs")
  }
  if (!identical(endpoint$objective_instance, owner$ln_joint$objective_instance) ||
      !identical(endpoint$objective_instance, state$objective_instance) ||
      !identical(endpoint$gene_index, as.integer(owner$debug$gene_index))) {
    fail("endpoint objective instance or gene differs from its owner")
  }
  scaled_error <- function(a, b) {
    if (!is.numeric(a) || !is.numeric(b) || !length(a) ||
        length(a) != length(b) || any(!is.finite(c(a, b)))) return(Inf)
    max(abs(a - b) / (1 + pmax(abs(a), abs(b))))
  }
  if (!identical(endpoint$stored_state_objective, as.numeric(state$raw$value)) ||
      !identical(endpoint$state_objective, endpoint$reevaluated_objective) ||
      scaled_error(endpoint$state_objective, state$raw$value) > threshold) {
    fail("endpoint objective does not match its owner's reevaluated state")
  }
  candidates <- do.call(rbind, lapply(seq_along(records), function(i) {
    candidate <- records[[i]]
    data.frame(
      index = as.integer(i), label = candidate$label,
      certified = isTRUE(candidate$certified), interior = isTRUE(candidate$interior),
      exact_state = !is.null(candidate$theta_full_log) &&
        identical(as.numeric(candidate$theta_full_log), as.numeric(state$full_log)),
      provenance = identical(candidate$objective_provenance, provenance),
      state_scaled_error = scaled_error(candidate$theta_full_log, state$full_log),
      objective_scaled_error = scaled_error(candidate$objective,
                                             endpoint$reevaluated_objective),
      stringsAsFactors = FALSE
    )
  }))
  eligible <- with(candidates, index[certified & interior & exact_state &
                                      provenance & objective_scaled_error <= threshold])
  valid_index <- is.integer(index) && length(index) == 1L && !is.na(index) &&
    index >= 1L && index <= length(records)
  if (!valid_index) {
    fail("witness index is not an in-range scalar integer")
  } else {
    if (!length(eligible) || !identical(index, eligible[[1L]])) {
      fail("witness is not the first eligible candidate in its own fit")
    }
    if (!identical(endpoint$candidate_objective,
                   as.numeric(records[[index]]$objective))) {
      fail("witness objective does not match the indexed candidate")
    }
  }
  list(errors = errors, checked = TRUE, index = index, eligible = eligible,
       objective_threshold = threshold, candidates = candidates)
}

ln_joint_psock_witness_audit <- function(raw) {
  errors <- character(0)
  witnesses <- list()
  normalized <- raw
  for (g in seq_along(raw)) {
    gene <- raw[[g]]
    for (contrast in names(gene$lrt)) {
      constrained <- gene$lrt[[contrast]]$constrained
      endpoints <- constrained$ln_joint$pace_ratio_diagnostic$endpoints
      if (is.null(endpoints)) next
      for (kind in c("unrestricted", "constrained")) {
        label <- paste(g, contrast, kind, sep = "/")
        owner <- if (kind == "unrestricted") gene else constrained
        audit <- ln_joint_psock_endpoint_witness(endpoints[[kind]], owner, kind)
        witnesses[[label]] <- audit
        if (length(audit$errors)) errors <- c(errors, paste0(label, ": ", audit$errors))
        if (isTRUE(audit$checked) && !length(audit$errors)) {
          # Replace only an independently validated run-local witness position
          # in a comparison copy. Original diagnostics and every other integer,
          # candidate, endpoint, status and numerical tolerance stay unchanged.
          normalized[[g]]$lrt[[contrast]]$constrained$ln_joint$
            pace_ratio_diagnostic$endpoints[[kind]]$certified_candidate_index <- 0L
        }
      }
    }
  }
  list(errors = errors, witnesses = witnesses, normalized = normalized)
}

# Receipts contain only this small synthetic fixture, public fit fields and
# plain diagnostic values. Exclude sessionInfo/calls/environments and all paths.
ln_joint_psock_receipt_safe <- function(x) {
  if (is.null(x)) return(TRUE)
  if (!(is.atomic(x) || is.list(x)) || is.object(x) && !is.data.frame(x)) return(FALSE)
  if (is.character(x) && any(grepl(
      "(/Users/|/home/|/private/|/tmp/|/var/|[A-Za-z]:[\\\\/]Users[\\\\/]|file://)",
      x, perl = TRUE), na.rm = TRUE)) return(FALSE)
  if (is.list(x) && !all(vapply(x, ln_joint_psock_receipt_safe, logical(1)))) return(FALSE)
  attrs <- attributes(x)
  is.null(attrs) || all(vapply(attrs, ln_joint_psock_receipt_safe, logical(1)))
}

ln_joint_save_psock_pair <- function(name, serial, parallel = NULL) {
  stopifnot(name %in% c("optional-tests", "default-wald"))
  snapshot <- function(fit) {
    if (is.null(fit)) return(NULL)
    selected <- fit[c("summary", "overdispersion", "convergence", "provenance",
                      "lrt_details")]
    selected$raw <- fit$diagnostics$raw
    audit <- tryCatch(ln_joint_psock_witness_audit(selected$raw), error = function(e) {
      list(errors = "witness audit raised an error; inspect the saved raw pair",
           witnesses = list())
    })
    selected$witness_audit <- audit[c("errors", "witnesses")]
    selected
  }
  record <- list(
    purpose = "synthetic_parallel_software_regression_not_calibration",
    schema_version = 1L, test = name,
    fixture = list(seed = 20260720L, genes = 1L, donors = 16L,
                   cells_per_donor = 24L),
    requested_workers = c(serial = 1L, parallel = 2L),
    package_version = as.character(utils::packageVersion("scRoPE")),
    R_version = R.version$version.string, platform = R.version$platform,
    dll_md5 = unname(tools::md5sum(getLoadedDLLs()[["scRoPE"]][["path"]])),
    serial = snapshot(serial), parallel = snapshot(parallel)
  )
  stopifnot(object.size(record) < 8 * 1024^2, ln_joint_psock_receipt_safe(record))
  directory <- Sys.getenv("SCROPE_TEST_DIAGNOSTICS_DIR", unset = "_diagnostics")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(directory, paste0("psock-pair-", name))
  saveRDS(record, paste0(stem, ".rds"))
  lines <- c(paste("Synthetic parallel pair:", name),
             paste("Both fits captured:", !is.null(parallel)))
  for (which in c("serial", "parallel")) {
    side <- record[[which]]
    if (is.null(side)) next
    lines <- c(lines, paste(which, "witness errors:", length(side$witness_audit$errors)),
               side$witness_audit$errors)
    for (label in names(side$witness_audit$witnesses)) {
      witness <- side$witness_audit$witnesses[[label]]
      lines <- c(lines, paste(which, label, "index", witness$index,
                              "eligible", paste(witness$eligible, collapse = ",")))
    }
  }
  writeLines(lines, paste0(stem, ".txt"))
  invisible(record)
}
