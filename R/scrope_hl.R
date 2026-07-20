#' sc-RoPE HL pipeline
#'
#' Fits the sc-RoPE model using the H-likelihood (APHL) objective for every
#' gene. The output mirrors the structure of `scrope()` while providing the
#' APHL value and subject-level caches required for the adjusted HL tests.
#'
#' @inheritParams scrope
#' @param additional_tests Character vector selecting additional adjusted HL
#'   statistics to compute. Supported values are `"score"` and `"lrt"`. Wald
#'   statistics are always included.
#' @param lrt Logical convenience flag; equivalent to including `"lrt"` in
#'   `additional_tests`.
#' @param lrt_contrasts Optional list of user-specified contrasts. Each entry
#'   may be a matrix `L`, a numeric vector, or a list with elements `L`, `b`,
#'   and an optional `name`.
#' @param lrt_include_scale Logical; include the HL invariant scale factor for
#'   each contrast in the summary output when `"lrt"` is requested.
#' @param lrt_details Logical; if `TRUE`, attach a `hl_adj_details` data frame
#'   with per-contrast diagnostics (statistics, p-values, ridge indicators).
#' @param adj_ridge_factor Baseline ridge multiplier supplied to the adjusted
#'   test solver (see `hl_adjusted_tests()`).
#' @param adj_cond_limit Condition-number threshold for injecting a ridge into
#'   the adjusted test inversions.
#' @param keep_diagnostics Logical; retain per-gene diagnostic objects
#'   (including the HL caches).
#' @return A list matching the structure of `scrope()` with the following HL
#'   additions: the summary gains columns `se_robust_*`, `wald_stat_*`,
#'   `p_robust_*`, and their success/failure flags
#'   (`wald_success_*`, `wald_failure_*`). Optional score/LRT requests append
#'   `hl_score_adj_*`, `p_hl_score_*`, `hl_lrt_adj_*`, `p_hl_lrt_*`, and
#'   `hl_lrt_scale_*`. The tidy `hl_adj_details` frame is populated when
#'   `lrt_details = TRUE`, and `diagnostics$hl_adj` stores the per-contrast
#'   objects when `keep_diagnostics = TRUE`.
#' @export
scrope_hl <- function(
    count,
    id,
    pred = NULL,
    offset = NULL,
    min = c(1e-4, 1e-4),
    max = c(10, 1000),
    opt = "lbfgs",
    verbose = TRUE,
    cpc = 0.005,
    mincp = 5,
    ncore = 1,
    fmaxsize = Inf,
    use_betas = "betae",
    output_re = FALSE,
    additional_tests = character(0),
    lrt = FALSE,
    lrt_contrasts = NULL,
    lrt_include_scale = FALSE,
    lrt_details = FALSE,
    adj_ridge_factor = 1e-8,
    adj_cond_limit = 1e8,
    keep_diagnostics = TRUE) {
  hl_fo_public_scrope_hl(
    count = count,
    id = id,
    pred = pred,
    offset = offset,
    min = min,
    max = max,
    opt = opt,
    verbose = verbose,
    cpc = cpc,
    mincp = mincp,
    ncore = ncore,
    fmaxsize = fmaxsize,
    use_betas = use_betas,
    output_re = output_re,
    additional_tests = additional_tests,
    lrt = lrt,
    lrt_contrasts = lrt_contrasts,
    lrt_include_scale = lrt_include_scale,
    lrt_details = lrt_details,
    adj_ridge_factor = adj_ridge_factor,
    adj_cond_limit = adj_cond_limit,
    keep_diagnostics = keep_diagnostics
  )
}
