# Reevaluate retained finite observations using the current scalar criterion.
# Native trajectory history is stored in `recovery`, not `bounded_recovery`.
ln_recovery_scalar_evidence <- function(theta, dat) {
  comp <- ln_joint_objective_components(theta, 1L, dat$posv, dat$ctx, dat$map)
  q <- sum(rowSums(comp)); scale <- sum(abs(comp))
  if (any(!is.finite(comp)) || !is.finite(q) || !is.finite(scale))
    stop("Nonfinite retained objective evaluation.")
  eps <- length(comp) * .Machine$double.eps
  list(q = q, scale = scale,
       error = eps / (1 - eps) * scale + .Machine$double.eps * (1 + abs(q)))
}

ln_recovery_known_points <- function(primary, dat) {
  records <- primary$joint$candidate_records
  points <- list(); p <- dat$ctx$nb
  if (!length(records)) return(list(success = FALSE, reason = "no_retained_candidates"))
  instance <- records[[1L]]$objective_instance
  provenance <- ln_joint_objective_provenance()
  add <- function(theta, reported, label, allowance = 0) {
    reported <- as.numeric(reported)
    if (is.null(theta) && !any(is.finite(reported))) return(invisible(NULL))
    valid <- is.numeric(theta) && length(theta) == p + 2L &&
      all(is.finite(theta)) && all(tail(theta, 2L) > 0)
    if (!valid) {
      if (any(is.finite(reported))) stop("Finite objective has no usable retained coordinates.")
      return(invisible(NULL))
    }
    e <- ln_recovery_scalar_evidence(theta, dat)
    discrepancy <- if (any(is.finite(reported)))
      max(abs(reported[is.finite(reported)] - e$q)) else 0
    if (!is.numeric(allowance) || length(allowance) != 1L ||
        !is.finite(allowance) || allowance < 0 ||
        discrepancy > e$error + allowance + 1e-8 * (1 + abs(e$q)))
      stop("Retained scalar objective is inconsistent with its coordinates.")
    d <- tail(theta, 2L)
    tol <- sqrt(.Machine$double.eps) * (1 + pmax(abs(dat$ctx$min), abs(dat$ctx$max)))
    points[[length(points) + 1L]] <<- data.frame(
      label = label, q = e$q, scale = e$scale,
      error = e$error + discrepancy + allowance,
      boundary = any(d <= dat$ctx$min + tol | d >= dat$ctx$max - tol),
      recorded_discrepancy = discrepancy, objective_allowance = allowance)
  }
  same_context <- function(r) {
    if (!identical(r$objective_instance, instance) ||
        !identical(r$objective_provenance, provenance))
      stop("Retained objective provenance does not match.")
  }
  valid <- tryCatch({
    for (r in records) {
      same_context(r)
      add(r$theta_reduced, r$objective, paste0(r$label, ":endpoint"))
      add(r$optimizer_endpoint, r$optimizer_objective, paste0(r$label, ":optimizer"))
      h <- r$recovery
      if (!is.null(h)) {
        allowance <- if (is.null(h$objective_allowance)) 0 else h$objective_allowance
        for (name in c("original", "certificate")) {
          z <- h[[name]]
          if (!is.null(z)) {
            same_context(z)
            add(z$theta_reduced, z$objective, paste0(r$label, ":recovery_", name), allowance)
            add(z$optimizer_endpoint, z$optimizer_objective,
                paste0(r$label, ":recovery_", name, "_optimizer"), allowance)
          }
        }
        add(h$optimizer$par, h$optimizer$value, paste0(r$label, ":recovery_optimizer"), allowance)
        add(h$polish$theta_reduced, h$polish$objective, paste0(r$label, ":recovery_polish"), allowance)
      }
      z <- r$resolution$probe$proposed
      if (!is.null(z)) {
        theta <- z$theta_reduced_log
        if (is.numeric(theta) && length(theta) == p + 2L)
          theta[p + seq_len(2L)] <- exp(theta[p + seq_len(2L)])
        add(theta, z$objective, paste0(r$label, ":basin_probe"))
      }
    }
    TRUE
  }, error = function(e) e)
  if (inherits(valid, "error")) return(list(success = FALSE,
    reason = "retained_objective_audit_failed", error = conditionMessage(valid)))
  if (!length(points)) return(list(success = FALSE, reason = "no_retained_finite_points"))
  list(success = TRUE, observations = do.call(rbind, points))
}

ln_recovery_competitor_audit <- function(primary, dat, q, reference_error = 0) {
  audit <- ln_recovery_known_points(primary, dat)
  if (!isTRUE(audit$success)) return(audit)
  x <- audit$observations
  x$gap <- x$q - q
  x$margin <- x$error + reference_error + 1e-8 * (1 + abs(q))
  x$passes <- x$gap >= -x$margin
  audit$observations <- x
  audit$known_objective_min <- min(x$q)
  audit$success <- all(x$passes)
  audit$reason <- if (audit$success) "no_known_better_full_point" else "known_better_full_point"
  audit
}
