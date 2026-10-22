# Export-only helper added after the locked fit. No fitting or generation.
# Preserves prespecified gene order, point diagnostics and reported rejection
# status; incomplete/unsupported inference cannot produce SEs, tests or intervals.
extract_toy_results <- function(toy, fit, fit_error = NULL) {
  results <- toy$gene_metadata
  if (!is.data.frame(results) || !"gene" %in% names(results) ||
      anyNA(results$gene) || anyDuplicated(results$gene)) {
    stop("toy$gene_metadata must contain unique, nonmissing gene names.")
  }
  if (!is.numeric(results$truth)) stop("Generating truth must be numeric.")
  results$availability <- if (is.null(fit)) "fit_error" else "not_returned"
  mapping <- c(estimate = "logFC_x", model_se = "se_x", robust_se = "se_robust_x",
               model_p_zero = "p_naive_x", robust_p_zero = "p_robust_x")
  for (nm in names(mapping)) results[[nm]] <- NA_real_
  reason <- rep(if (is.null(fit)) "fit_error_not_reported" else
                  "missing_summary_and_availability", nrow(results))
  if (is.null(fit) && !is.null(fit_error) && length(fit_error)) {
    reason[] <- paste0("fit_error: ", paste(as.character(fit_error), collapse = "; "))
  }
  if (!is.null(fit)) {
    table_or_empty <- function(x, name) {
      if (is.null(x)) return(data.frame(gene = character()))
      if (!is.data.frame(x) || !"gene" %in% names(x)) {
        stop(name, " must be a data frame with gene names.")
      }
      if (anyNA(x$gene) || anyDuplicated(x$gene)) {
        stop(name, " has duplicate or missing gene names.")
      }
      x
    }
    s <- table_or_empty(fit$summary, "fit$summary")
    a <- table_or_empty(fit$availability, "fit$availability")
    j <- match(results$gene, s$gene)
    k <- match(results$gene, a$gene)
    has_summary <- !is.na(j)
    has_availability <- !is.na(k)
    # A present point estimate can be retained even when inference is unsupported.
    for (nm in names(mapping)) {
      column <- mapping[[nm]]
      if (column %in% names(s)) {
        if (!is.numeric(s[[column]])) stop("Summary column must be numeric: ", column)
        results[[nm]] <- s[[column]][j]
      }
    }
    reason[has_summary & !has_availability] <- "missing_availability"
    reason[!has_summary & has_availability] <- "missing_summary"
    both <- has_summary & has_availability
    if ("availability" %in% names(a)) {
      statuses <- as.character(a$availability[k])
    } else {
      statuses <- rep(NA_character_, nrow(results))
    }
    status_present <- !is.na(statuses) & nzchar(statuses)
    results$availability[both & status_present] <- statuses[both & status_present]
    reason[both & !status_present] <- "missing_availability_status"
    # Prefer explicit reason codes, with readable status text as a fallback.
    for (i in which(both & status_present)) {
      candidates <- unlist(lapply(c("reason_group", "technical_code", "status_message"),
        function(nm) if (nm %in% names(a)) as.character(a[[nm]][k[i]]) else NA_character_),
        use.names = FALSE)
      candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
      reason[i] <- if (length(candidates)) candidates[[1L]] else "reason_not_reported"
    }
    supported <- both & status_present & statuses == "inference"
    fields_valid <- is.finite(results$estimate) &
      is.finite(results$model_se) & results$model_se > 0 &
      is.finite(results$robust_se) & results$robust_se > 0 &
      is.finite(results$model_p_zero) & results$model_p_zero >= 0 & results$model_p_zero <= 1 &
      is.finite(results$robust_p_zero) & results$robust_p_zero >= 0 & results$robust_p_zero <= 1
    invalid <- supported & !fields_valid
    results$availability[invalid] <- ifelse(is.finite(results$estimate[invalid]),
                                          "estimate_only", "unavailable")
    reason[invalid] <- "invalid_or_missing_inference_fields"
    supported <- supported & fields_valid
    for (nm in setdiff(names(mapping), "estimate")) results[[nm]][!supported] <- NA_real_
  }
  for (method in c("model", "robust")) {
    lower <- results$estimate - qnorm(.975) * results[[paste0(method, "_se")]]
    upper <- results$estimate + qnorm(.975) * results[[paste0(method, "_se")]]
    results[[paste0(method, "_lower")]] <- lower
    results[[paste0(method, "_upper")]] <- upper
    results[[paste0(method, "_includes_truth")]] <- lower <= results$truth & upper >= results$truth
  }
  # Append the reason, preserving every original public column's position.
  results$reason <- reason
  results
}
