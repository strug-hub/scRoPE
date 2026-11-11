#' Opt-in JLN boundary and search-warning research workflow
#'
#' Runs ordinary JLN first, then optionally reports a corroborated interior root
#' with a search warning or a fixed-positive-dispersion working fit. Existing
#' This compatibility wrapper preserves the earlier opt-in rule and explicitly
#' uses strict native scrope() fits; the main interface now supports recovery.
#' @export
scrope_jln_research <- function(count, id, pred, target, offset = NULL,
                                boundary = TRUE, dominated_warning = TRUE,
                                cpc = 0.005, mincp = 5L, ncore = 1L,
                                keep_diagnostics = TRUE, verbose = TRUE) {
  logical_scalar <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!all(vapply(list(boundary, dominated_warning, keep_diagnostics, verbose),
                  logical_scalar, logical(1)))) {
    stop("Policy and retention flags must each be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.matrix(pred) || !is.numeric(pred) || any(!is.finite(pred)) ||
      ncol(pred) < 2L || any(pred[, 1L] != 1) || is.null(colnames(pred)) ||
      anyNA(colnames(pred)) || any(!nzchar(colnames(pred))) ||
      anyDuplicated(colnames(pred))) {
    stop("pred must be a finite numeric matrix with named columns and an intercept of ones first.", call. = FALSE)
  }
  if (!is.character(target) || length(target) != 1L || is.na(target) ||
      !target %in% colnames(pred)[-1L]) {
    stop("target must name one prespecified non-intercept column of pred.", call. = FALSE)
  }
  target_index <- match(target, colnames(pred))
  if (is.vector(count)) count <- matrix(count, nrow = 1L)
  count <- methods::as(count, "dgCMatrix")
  if (ncol(count) != nrow(pred) || ncol(count) < 2L || nrow(count) < 1L ||
      any(!is.finite(count@x)) || any(count@x < 0)) {
    stop("count must contain finite nonnegative counts aligned with pred.", call. = FALSE)
  }
  if (length(id) != ncol(count) || anyNA(id) ||
      is.unsorted(match(id, unique(id)))) {
    stop("id must be complete, aligned with count, and contiguous by donor.", call. = FALSE)
  }
  if (is.null(offset)) offset <- rep(1, ncol(count))
  if (!is.numeric(offset) || length(offset) != ncol(count) ||
      any(!is.finite(offset)) || any(offset <= 0)) {
    stop("offset must be a finite positive exposure for each cell.", call. = FALSE)
  }
  if (!is.numeric(cpc) || length(cpc) != 1L || !is.finite(cpc) || cpc < 0 ||
      !is.numeric(mincp) || length(mincp) != 1L || !is.finite(mincp) ||
      mincp < 1 || mincp != floor(mincp) ||
      !is.numeric(ncore) || length(ncore) != 1L || !is.finite(ncore) ||
      ncore < 1 || ncore != floor(ncore)) {
    stop("cpc must be nonnegative; mincp and ncore must be positive integers.", call. = FALSE)
  }
  primary_warnings <- character()
  primary <- withCallingHandlers(
    scrope(count = count, id = id, pred = pred, offset = offset,
           model = "LN", ln_estimator = "joint_ln", ln_recovery = "none",
           ln_operational_fallback = "none", allow_per_gene_switch = FALSE,
           opt = "lbfgs", cpc = cpc, mincp = mincp, ncore = ncore,
           additional_tests = character(), output_re = FALSE,
           keep_diagnostics = TRUE, verbose = verbose),
    warning = function(w) {
      primary_warnings <<- c(primary_warnings, conditionMessage(w))
      # This warning describes the native component, before research routing.
      # Other warnings (e.g. rank deficiency) are still emitted immediately.
      if (inherits(w,"scrope_recovery_warning") ||
          startsWith(conditionMessage(w), "Joint-LN inference was unavailable"))
        invokeRestart("muffleWarning")
    }
  )
  rows <- diagnostics <- vector("list", nrow(primary$summary))
  for (j in seq_len(nrow(primary$summary))) {
    sm <- primary$summary[j, , drop = FALSE]
    gid <- sm$gene_id
    pg <- primary$provenance$per_gene
    native_available <- isTRUE(pg$inference_available[match(gid, pg$gene_id)]) &&
      isTRUE(sm[[paste0("wald_success_", target)]]) &&
      is.finite(sm[[paste0("p_robust_", target)]])
    native <- list(success = native_available,
      coefficient = sm[[paste0("logFC_", target)]],
      se = sm[[paste0("se_robust_", target)]],
      p = sm[[paste0("p_robust_", target)]],
      joint = primary$diagnostics$raw[[j]]$ln_joint)
    ans <- ln_research_route(native, as.numeric(count[gid, ]), pred, id,
                             offset, target_index, boundary, dominated_warning)
    f <- ans$result
    available <- isTRUE(f$success)
    fixed <- if (is.null(ans$boundary_spec$fixed)) list() else ans$boundary_spec$fixed
    label <- if (ans$route == "ordinary_jln") "released_JLN_policy" else
      if (ans$route == "ordinary_jln_with_search_warning") "local_intended_branch_with_unresolved_search" else
        if (startsWith(ans$route, "fixed_")) "fixed_mask_epsilon_working_target" else "unavailable"
    caution <- if (ans$route == "ordinary_jln_with_search_warning")
      "A boundary-start trajectory remains unresolved; unvisited competing branches are not excluded." else
      if (startsWith(ans$route, "fixed_"))
        "Fixed-dispersion working target; data-selected route inference is not selectively calibrated." else
        if (!available) "No certified research inference; see reason and primary diagnostics." else ""
    beta <- if (available) f$coefficient else NA_real_
    se <- if (available) f$se else NA_real_
    model_se <- if(ans$route=="ordinary_jln")sm[[paste0("se_",target)]] else
      if(available&&isTRUE(f$model_inference_available))sqrt(f$model_Vbeta[target_index,target_index]) else NA_real_
    model_available <- available&&is.finite(model_se)&&model_se>0&&is.finite((beta/model_se)^2)
    rows[[j]] <- data.frame(gene_id = gid, gene = sm$gene, target = target,
      route = ans$route, available = available, native_available = native_available,
      research_added = available && !native_available,
      model_available=model_available,se_model=if(available&&is.finite(model_se)&&model_se>0)model_se else NA_real_,
      p_model=if(model_available)pchisq((beta/model_se)^2,1,lower.tail=FALSE) else NA_real_,
      model_covariance_convention=scrope_model_convention("LN"),
      coefficient = beta, se_robust = se, p_robust = if (available) f$p else NA_real_,
      ci_lower = beta - stats::qnorm(.975) * se,
      ci_upper = beta + stats::qnorm(.975) * se,
      fixed_donor_variance = if (is.null(fixed$v)) NA_real_ else fixed$v,
      fixed_cell_overdispersion = if (is.null(fixed$phi)) NA_real_ else fixed$phi,
      inferential_scope = label, selection_adjusted = FALSE,
      finite_cluster_multiplier = 1,
      reason = if (available) label else f$reason, warning = caution,
      stringsAsFactors = FALSE)
    if (keep_diagnostics) diagnostics[[j]] <- ans
  }
  summary <- do.call(rbind, rows)
  if (any(summary$research_added) || any(!summary$available)) {
    warning(sprintf(paste0("Opt-in JLN research workflow: %d/%d reported (%d research additions); ",
      "%d withheld. Research-route p-values are not selection-adjusted and empirical undercoverage ",
      "has been observed. Inspect summary$route, warning and inferential_scope."),
      sum(summary$available), nrow(summary), sum(summary$research_added), sum(!summary$available)),
      call. = FALSE)
  }
  if (!keep_diagnostics) primary$diagnostics <- NULL
  list(summary = summary, primary = primary,
       diagnostics = if (keep_diagnostics) diagnostics else NULL,
       primary_warnings = primary_warnings,
       policy = list(version = "jln_research_20261103", boundary = boundary,
         dominated_warning = dominated_warning, donor_epsilon = 1e-4,
         cell_epsilon = 1e-3, target = target, tests = "coefficient_zero_null_Wald_only",
         interpretation = paste0("Fixed-route sandwiches concern a regular fixed-mask working target. ",
           "Warning-route sandwiches concern the intended regular interior branch. ",
           "Automatic selection, scientific-target equality and FDR control are not established. ",
           "The original scrope result is preserved in primary.")))
}

ln_research_route <- function(primary, y, X, id, offset, target,
                               boundary, dominated_warning) {
  result <- primary
  if (!isTRUE(primary$success)) result$reason <- "native_inference_withheld"
  route <- "withheld"
  spec <- normal <- warning_audit <- NULL
  if (isTRUE(primary$success)) route <- "ordinary_jln" else {
    dat <- tryCatch(ln_research_context(y, X, id, offset), error = identity)
    if (inherits(dat, "error")) {
      result <- list(success = FALSE, reason = "input_or_design_failure", error = conditionMessage(dat))
    } else {
      if (dominated_warning) {
        warning_audit <- tryCatch(ln_research_warning_fit(primary, dat, target),
          error = function(e) list(success = FALSE, reason = "warning_evaluation_error", error = conditionMessage(e)))
        if (warning_audit$success) {
          result <- warning_audit
          route <- "ordinary_jln_with_search_warning"
        } else result$reason <- warning_audit$reason
      }
      if (route == "withheld" && boundary) {
        spec <- ln_research_boundary_spec(primary, FALSE)
        if (startsWith(spec$route, "fixed_")) {
          result <- tryCatch(ln_research_fixed_fit(dat, spec$fixed,
            ln_research_starts(primary, dat), target),
            error = function(e) list(success = FALSE, reason = "restricted_evaluation_error", error = conditionMessage(e)))
          if (result$success) {
            normal <- ln_research_boundary_check(result, dat)
            if (normal$success) route <- spec$route else
              result <- list(success = FALSE, reason = normal$reason)
          }
        } else result <- list(success = FALSE, reason = spec$reason)
      }
    }
  }
  list(route = route, result = result, boundary_spec = spec,
       boundary_check = normal, warning_audit = warning_audit)
}
