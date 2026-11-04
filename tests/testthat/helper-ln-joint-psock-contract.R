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
    # Preserve missing/nonfinite locations and kinds and the existing relative
    # comparison. Near-zero score/probe diagnostics also allow 1e-12 absolute
    # roundoff, well below the 1e-6 root and 1e-8 curvature certificates.
    if (!identical(is.na(actual), is.na(expected)) ||
        !identical(is.nan(actual), is.nan(expected)) ||
        !identical(is.finite(actual), is.finite(expected))) {
      return(fail("missing or nonfinite states differ"))
    }
    comparison <- all.equal(actual, expected, tolerance = tolerance)
    if (!isTRUE(comparison)) {
      absolute_fields <- c(
        "score", "g_red", "q_gradient_reduced", "U_reduced",
        "reduced_gradient_max_abs", "max_abs_newton_step_scaled",
        "root_residual", "curvature_compatibility", "nuisance_root_scaled"
      )
      finite <- is.finite(actual)
      within_floor <- field %in% absolute_fields &&
        identical(actual[!finite], expected[!finite]) &&
        all(abs(actual[finite] - expected[finite]) <=
              1e-12 + tolerance * pmax(abs(actual[finite]), abs(expected[finite])))
      if (!within_floor) return(fail(paste(comparison, collapse = "; ")))
    }
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


# A software regression budget, not an implication of root certification or a
# statistical calibration guarantee. Optional statistics/probabilities and
# selected log coordinates allow 1e-8 * (1 + max(abs(a), abs(b))) elementwise.
# Corresponding optimizer candidates and primary results retain their previous
# comparator. Selected gradients/conditioning are checked against their owner,
# not compared as if independent runs necessarily selected the same endpoint.
ln_joint_psock_optional_differences <- function(a, b, path = "optional", budget = 1e-8) {
  fail <- function(reason) paste0(path, ": ", reason)
  if (!identical(typeof(a), typeof(b)) ||
      !identical(attributes(a), attributes(b)) || length(a) != length(b)) {
    return(fail("type, attributes or length differ"))
  }
  if (is.list(a)) return(unlist(lapply(seq_along(a), function(i) {
    ln_joint_psock_optional_differences(a[[i]], b[[i]],
      paste0(path, "$", if (is.null(names(a))) i else names(a)[[i]]), budget)
  }), use.names = FALSE))
  if (!is.double(a)) return(if (identical(a, b)) character(0) else fail("metadata differ"))
  finite <- is.finite(a)
  if (!identical(finite, is.finite(b)) || !identical(is.na(a), is.na(b)) ||
      !identical(is.nan(a), is.nan(b)) || !identical(a[!finite], b[!finite])) {
    return(fail("missing or nonfinite states differ"))
  }
  if (any(abs(a[finite] - b[finite]) >
          budget * (1 + pmax(abs(a[finite]), abs(b[finite]))))) {
    return(fail("mixed numerical regression budget exceeded"))
  }
  character(0)
}

# Comparison-only structural view: all shapes, types, finite/missing kinds,
# integer policy values, strings and logical decisions remain exact. Numerical
# payloads get their substantive checks in the owner audit, candidate inventory,
# primary comparison and optional public-result comparison below.
ln_joint_psock_structure <- function(x, field = "") {
  if (is.character(x) && field %in% c("start", "optimizer_endpoint",
      "theta_reduced", "theta_full_log", "candidate_key", "anchor_candidate_key")) {
    parser <- if (grepl("candidate_key$", field)) ln_joint_psock_parse_key else
      ln_joint_psock_parse_numbers
    return(lapply(x, function(value) ln_joint_psock_structure(parser(value))))
  }
  if (is.list(x)) {
    for (i in seq_along(x)) x[i] <- list(ln_joint_psock_structure(x[[i]],
      if (is.null(names(x))) "" else names(x)[[i]]))
  } else if (is.double(x)) x[is.finite(x)] <- 0
  x
}

ln_joint_psock_owner_audit <- function(owner) {
  # This gate's synthetic fixture requires R0 owners; it is not a replacement
  # for the production resolver or a new certificate for arbitrary fits.
  errors <- character(0)
  check <- function(ok, message) {
    if (!isTRUE(ok)) errors <<- c(errors, message)
  }
  j <- owner$ln_joint
  records <- j$candidate_records
  state <- j$state
  t <- ln_joint_thresholds()
  check(identical(owner$route, "R0") && identical(owner$inference_available, TRUE),
        "fixture owner must report certified R0 inference")
  if (!is.list(state) || !length(records)) return(c(errors, "owner state/records missing"))
  errors <- c(errors, ln_joint_psock_key_errors(owner))
  check(identical(j$resolution$thresholds, ln_joint_resolution_thresholds()),
        "resolution policy thresholds differ")
  check(identical(j$resolution$selected_root_certified, TRUE),
        "owner resolution does not certify its selected root")
  check(identical(state$objective_instance, j$objective_instance) &&
          identical(state$objective_provenance, j$objective_provenance),
        "selected state objective ownership differs")
  for (i in seq_along(records)) {
    r <- records[[i]]
    dispersion <- tail(r$theta_reduced, 2L)
    margin <- sqrt(.Machine$double.eps) *
      (1 + pmax(abs(r$bounds$lower), abs(r$bounds$upper)))
    interior <- all(dispersion > r$bounds$lower + margin &
                      dispersion < r$bounds$upper - margin)
    finite_below <- function(v, limit) is.double(v) && length(v) == 1L &&
      is.finite(v) && v >= 0 && v <= limit
    checks <- c(
      donor_additivity = finite_below(r$donor_additivity_scaled, t$donor_additivity),
      root = finite_below(r$root_residual, t$root),
      newton = finite_below(r$newton_displacement, t$newton),
      curvature = isTRUE(r$curvature$positive_definite) &&
        is.finite(r$curvature$scaled_min_eigenvalue) &&
        r$curvature$scaled_min_eigenvalue >= t$curvature_min_eigen,
      condition = finite_below(r$curvature$scaled_condition_number, t$condition_number),
      interior = interior,
      constraint = finite_below(r$constraint_residual, t$constraint))
    check(identical(r$interior, interior) && identical(r$checks, checks) &&
            identical(r$certified, all(checks)) &&
            isTRUE(r$finite_endpoint), paste("candidate", i, "certificate inconsistent"))
    check(identical(r$objective_instance, j$objective_instance) &&
            identical(r$objective_provenance, j$objective_provenance),
          paste("candidate", i, "objective ownership differs"))
    polished <- !isTRUE(all.equal(r$theta_reduced, r$optimizer_endpoint,
                                  tolerance = 0, check.attributes = FALSE))
    check(identical(r$polished, polished), paste("candidate", i, "polishing flag inconsistent"))
    check(identical(r$resolution, j$resolution$candidate_resolutions[[i]]),
          paste("candidate", i, "resolution copies differ"))
  }
  reconstructed <- lapply(records, function(r) {
    r$optimizer <- list(convergence = r$optimizer_convergence, message = r$optimizer_message)
    r
  })
  check(identical(j$candidates, ln_joint_candidate_summary(
    reconstructed, j$resolution$candidate_resolutions)), "candidate table differs from its records")
  certified <- which(vapply(records, function(r) isTRUE(r$certified), logical(1)))
  if (!length(certified)) return(c(errors, "selected owner has no certified candidate"))
  selected <- records[[certified[[which.min(vapply(records[certified],
    `[[`, numeric(1), "objective"))]]]]
  check(identical(state$full_log, selected$theta_full_log) &&
          identical(state$theta_reduced_log, selected$theta_reduced_log),
        "selected state is not its own minimum-objective certified candidate")
  root <- ln_joint_root_metric(state$U_subject_reduced, state$U_reduced)
  newton <- max(abs(solve(state$A_reduced, state$U_reduced)) /
                  (1 + abs(state$theta_reduced_log)))
  curvature <- ln_joint_curvature_diagnostics(state$A_reduced)
  check(is.finite(root) && root <= t$root && is.finite(newton) && newton <= t$newton &&
          isTRUE(curvature$positive_definite) &&
          curvature$scaled_min_eigenvalue >= t$curvature_min_eigen &&
          curvature$scaled_condition_number <= t$condition_number,
        "selected state fails recomputed root/Newton/curvature certificate")
  # These are same-run identities, separate from the cross-run regression budget.
  check(isTRUE(all.equal(root, selected$root_residual, tolerance = 1e-12)) &&
          isTRUE(all.equal(newton, selected$newton_displacement, tolerance = 1e-12)) &&
          isTRUE(all.equal(curvature, selected$curvature, tolerance = 1e-12)),
        "selected certificate metrics differ from its state")
  check(identical(owner$score, state$U_full) &&
          identical(state$q_gradient_full, -state$U_full) &&
          identical(state$q_gradient_reduced, -state$U_reduced),
        "selected score/sign copies differ")
  check(ln_joint_scaled_max(rowSums(state$U_subject_reduced), state$U_reduced) <=
          t$donor_additivity, "selected donor scores do not add up")
  own_close <- function(a, b) !length(ln_joint_psock_optional_differences(a, b, budget = 1e-12))
  check(own_close(state$full_log, c(state$beta, log(state$sigma2), log(state$gamma))),
        "selected natural/log coordinates disagree")
  check(identical(state$U_full, -as.numeric(state$raw$gradient)) &&
          identical(state$U_subject_full, -as.matrix(state$raw$per_subject_gradients)) &&
          identical(state$A_full, 0.5 * (state$raw$hessian + t(state$raw$hessian))),
        "selected raw derivatives differ from their state")
  expected_godambe <- extract_godambe_components(list(gradient = state$U_full,
    hessian = state$A_full, per_subject_gradients = state$U_subject_full), length(state$beta))
  check(own_close(owner$godambe, expected_godambe), "Godambe copies differ from their state")
  check(isTRUE(all.equal(-owner$loglik, selected$objective, tolerance = 1e-12)) &&
          isTRUE(all.equal(state$raw$value, selected$objective, tolerance = 1e-12)),
        "selected objective copies differ")
  if (!is.null(owner$constraint)) {
    map <- ln_joint_map(length(state$beta), owner$constraint$L, owner$constraint$b)
    T <- map$derivative_map
    projected <- crossprod(T, state$A_full %*% T)
    check(own_close(state$A_reduced, 0.5 * (projected + t(projected))) &&
            own_close(state$U_reduced, as.numeric(crossprod(T, state$U_full))) &&
            own_close(state$U_subject_reduced, crossprod(T, state$U_subject_full)),
          "selected full/reduced derivatives disagree")
    check(ln_joint_constraint_metric(state$beta, map) <= t$constraint,
          "selected state violates its constraint")
    check(identical(owner$reduced$score, state$U_reduced) &&
            identical(owner$reduced$g_red, state$U_reduced) &&
            identical(owner$reduced$H_red, state$A_reduced) &&
            identical(owner$reduced$per_subject_gradients, state$U_subject_reduced),
          "reduced diagnostic copies differ from their owner")
    check(identical(owner$beta_constrained, state$beta) &&
            identical(owner$theta, list(beta = state$beta, subVar = state$sigma2,
                                        cellVar = state$gamma)) &&
            identical(owner$theta, owner$theta_ln) &&
            identical(owner$theta, owner$theta_ln_rescued) &&
            identical(owner$theta, owner$theta_unrescued) &&
            identical(owner$reduced$theta, state$theta_reduced_log) &&
            identical(owner$reduced$theta_red, state$theta_reduced_log),
          "selected parameter copies differ from their owner")
    expected <- list(reduced_gradient_max_abs = max(abs(state$U_reduced)),
      H_red_positive_definite = curvature$positive_definite,
      H_red_rcond = 1 / curvature$scaled_condition_number,
      max_abs_newton_step_scaled = newton,
      constraint_residual_scaled = ln_joint_constraint_metric(state$beta, map))
    check(isTRUE(all.equal(owner$diagnostics, expected, tolerance = 1e-12)),
          "selected diagnostics differ from their owner")
    check(own_close(j$score_diagnostic, ln_joint_robust_score_diagnostic(owner)),
          "score diagnostic does not reconstruct from its owner")
  }
  errors
}

ln_joint_psock_without_polishing <- function(owner) {
  owner$ln_joint$candidates$polished[] <- FALSE
  for (i in seq_along(owner$ln_joint$candidate_records)) {
    owner$ln_joint$candidate_records[[i]]$polished <- FALSE
  }
  owner
}

ln_joint_psock_optional_fit_errors <- function(one, two) {
  raw <- function(fit) if (is.null(fit$raw)) fit$diagnostics$raw else fit$raw
  errors <- character(0)
  add <- function(x) errors <<- c(errors, x)
  optional <- grep("^(lrt_adj_|p_lrt_adj_|score_adj_|p_score_adj_)",
                   names(one$summary), value = TRUE)
  primary <- setdiff(names(one$summary), optional)
  add(ln_joint_psock_differences(one$summary[primary], two$summary[primary], path = "primary summary"))
  add(ln_joint_psock_optional_differences(one$summary[optional], two$summary[optional], "optional summary"))
  add(ln_joint_psock_differences(one$overdispersion, two$overdispersion, path = "dispersion"))
  add(ln_joint_psock_differences(one$convergence, two$convergence, path = "convergence"))
  add(ln_joint_psock_differences(one$provenance, two$provenance, path = "provenance"))
  sides <- list(raw(one), raw(two))
  audits <- lapply(sides, ln_joint_psock_witness_audit)
  for (audit in audits) add(audit$errors)
  sides <- lapply(audits, function(a) ln_joint_normalise_psock_objective_ids(a$normalized))
  # Diagnostics exposed in lrt_details must agree with that row's own endpoint.
  diagnostic <- c("reduced_gradient_max_abs", "max_abs_newton_step_scaled", "H_red_rcond")
  for (side in seq_along(sides)) {
    fit <- list(one, two)[[side]]
    for (row in seq_len(nrow(fit$lrt_details))) {
      detail <- fit$lrt_details[row, ]
      g <- match(detail$gene_id, fit$summary$gene_id)
      owner <- sides[[side]][[g]]$lrt[[detail$contrast]]$constrained
      for (field in diagnostic) if (!identical(detail[[field]], owner$diagnostics[[field]])) {
        add(paste("lrt_details", row, field, "differs from its owner"))
      }
    }
  }
  add(ln_joint_psock_differences(ln_joint_psock_structure(one$lrt_details),
    ln_joint_psock_structure(two$lrt_details), path = "lrt_details structure"))
  fields <- setdiff(names(one$lrt_details), diagnostic)
  add(ln_joint_psock_optional_differences(one$lrt_details[fields], two$lrt_details[fields], "lrt_details"))
  for (g in seq_along(sides[[1L]])) {
    a <- sides[[1L]][[g]]; b <- sides[[2L]][[g]]
    # Audit originals before creating any comparison-only normalization.
    for (owner in list(raw(one)[[g]], raw(two)[[g]])) add(ln_joint_psock_owner_audit(owner))
    a <- ln_joint_psock_without_polishing(a); b <- ln_joint_psock_without_polishing(b)
    contrasts <- names(a$lrt)
    primary_a <- a; primary_b <- b; primary_a$lrt <- primary_b$lrt <- NULL
    add(ln_joint_psock_differences(primary_a, primary_b, path = "primary raw"))
    add(ln_joint_psock_differences(names(a$lrt), names(b$lrt), path = "contrast names"))
    for (contrast in contrasts) {
      u <- a$lrt[[contrast]]$constrained; v <- b$lrt[[contrast]]$constrained
      for (owner in list(raw(one)[[g]]$lrt[[contrast]]$constrained,
                         raw(two)[[g]]$lrt[[contrast]]$constrained)) {
        add(ln_joint_psock_owner_audit(owner))
      }
      for (original in list(raw(one)[[g]], raw(two)[[g]])) {
        owner <- original$lrt[[contrast]]$constrained
        p <- owner$ln_joint$pace_ratio_diagnostic
        s <- owner$ln_joint$score_diagnostic
        expected <- list(wP = 2 * (original$loglik - owner$loglik),
          wU = s$statistic,
          qP = as.numeric(crossprod(s$U_eff, solve(p$H_eff, s$U_eff))))
        expected$scale <- expected$wU / expected$qP
        expected$statistic <- expected$scale * expected$wP
        add(ln_joint_psock_optional_differences(p[names(expected)], expected,
          paste(contrast, "own Pace identities"), budget = 1e-12))
        if (!identical(p$score_diagnostic, s)) add(paste(contrast, "Pace/score copies differ"))
      }
      u <- ln_joint_psock_without_polishing(u); v <- ln_joint_psock_without_polishing(v)
      wrapper_a <- a$lrt[[contrast]]; wrapper_b <- b$lrt[[contrast]]
      wrapper_a$constrained <- u; wrapper_b$constrained <- v
      add(ln_joint_psock_differences(ln_joint_psock_structure(wrapper_a),
        ln_joint_psock_structure(wrapper_b), path = paste(contrast, "test structure")))
      result_fields <- c("stat", "pval", "scale", "wP", "wU", "qP",
                         "score_stat", "score_pval")
      add(ln_joint_psock_optional_differences(wrapper_a[result_fields],
        wrapper_b[result_fields], paste(contrast, "test results")))
      inventory <- function(owner) lapply(owner$ln_joint$candidate_records, function(r) {
        # Resolution distances refer to the selected anchor; candidate endpoints,
        # certificate metrics/checks and optimizer history remain tightly compared.
        r$resolution <- NULL
        r
      })
      add(ln_joint_psock_differences(inventory(u), inventory(v), path = paste(contrast, "candidates")))
      add(ln_joint_psock_differences(ln_joint_psock_structure(u),
        ln_joint_psock_structure(v), path = paste(contrast, "structure")))
      add(ln_joint_psock_optional_differences(u$ln_joint$state$full_log,
        v$ln_joint$state$full_log, paste(contrast, "selected log coordinates")))
      add(ln_joint_psock_differences(u$loglik, v$loglik, path = paste(contrast, "objective")))
      inference <- function(owner) {
        j <- owner$ln_joint
        list(score = j$score_diagnostic[c("statistic", "df", "public_available")],
          pace = j$pace_ratio_diagnostic[c("statistic", "scale", "wP", "wU", "qP")])
      }
      add(ln_joint_psock_optional_differences(inference(u), inference(v), paste(contrast, "inference")))
    }
  }
  errors
}
