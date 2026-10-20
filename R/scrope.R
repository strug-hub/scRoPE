scrope_psock_library_paths <- function() {
  paths <- .libPaths()
  package_path <- tryCatch(
    getNamespaceInfo(asNamespace("scRoPE"), "path"),
    error = function(e) NA_character_
  )
  if (length(package_path) == 1L && !is.na(package_path) &&
      file.exists(file.path(package_path, "Meta", "package.rds"))) {
    paths <- c(dirname(package_path), paths)
  }
  unique(normalizePath(paths, mustWork = FALSE))
}

#' sc-RoPE main entry point
#'
#' Fits sc-RoPE working-likelihood models for multi-subject single-cell data
#' and, where supported by the selected estimator, computes optional adjusted
#' tests for fixed-effect contrasts.
#'
#' @param count Gene-by-cell count matrix (dense matrix or `dgCMatrix`).
#' @param id Vector of subject identifiers aligned with the columns of `count`.
#' @param pred Cell-level design matrix (defaults to intercept only).
#' @param offset Optional positive vector of exposure offsets.
#' @param min Lower bounds for the overdispersion parameters `(sigma^2, phi)`.
#' @param max Upper bounds for the overdispersion parameters `(sigma^2, phi)`.
#' @param opt Optimiser selector (`"lbfgs"` or `"trust"`). It controls the LN
#'   route. The corrected first-order HL route uses its dedicated solver, and
#'   PGMM uses its exact-route optimisers; those routes warn when `"trust"` is
#'   supplied.
#' @param verbose Logical; print progress messages.
#' @param cpc Minimum counts-per-cell threshold for gene filtering.
#' @param mincp Minimum number of non-zero cells required to keep a gene.
#' @param ncore Number of cores for parallel processing.
#' @param fmaxsize Maximum future.globals size (in bytes) when using parallel execution.
#' @param cutoff_cell Threshold (cells per subject × phi) below which legacy LN
#'   refinement is triggered. Non-default values are unavailable outside the
#'   legacy LN route.
#' @param kappa Kappa threshold for adaptive legacy LN refinement. This argument
#'   is unavailable outside the legacy LN route when set to a non-default value.
#' @param allow_per_gene_switch Logical; allow legacy per-gene LN/HL switching.
#'   Defaults to `FALSE` so one requested route is used for every gene. This is
#'   unavailable for `model = "HL"`, `model = "PGMM"`, and joint-LN.
#' @param use_betas Which fixed-effect anchor to use in LN variance-refinement
#'   starts: `"betae"` matches legacy `nebula` intercept-only anchoring, while
#'   `"final_betas"` uses the LN optimizer coefficients. The non-default option
#'   is unavailable for PGMM and ignored with a warning by first-order HL.
#' @param output_re Logical; return subject-level random effects. PGMM does not
#'   currently support this output and warns before returning `NULL`.
#' @param additional_tests Optional character vector selecting extra
#'   contrast-level diagnostics to compute, in addition to the default
#'   adjusted Wald columns (`se_robust_*`, `p_robust_*`). Supported values
#'   are `"score"` and `"lrt"`. Defaults to character(0) (Wald only).
#' @param lrt Logical; legacy flag for adjusted profile LRTs. Equivalent to
#'   including `"lrt"` in `additional_tests`.
#' @param lrt_contrasts Optional list of additional hypotheses `L beta = b`.
#'   Each entry may be a matrix/vector `L` (zero null) or a list containing
#'   `L`, optional `b`, and optional `name`. Columns follow `pred` and use the
#'   reported coefficient units. LN/PGMM coefficients have input-scale slopes
#'   and a sample-centered intercept; corrected HL uses original-design
#'   coordinates. Internal scaling is applied to `L`, leaving `b` unchanged.
#'   An intercept-containing LN/PGMM contrast does not automatically define a
#'   fixed population target.
#' @param lrt_include_scale Logical; include the H/G scale factor in the summary.
#' @param lrt_details Logical; return a tidy data frame with raw LRT diagnostics.
#' @param keep_diagnostics Logical; retain per-gene diagnostic objects. For HL,
#'   `FALSE` also suppresses optional solver histories and cross-gene retention
#'   of full fits; required scientific computations and public summaries remain
#'   unchanged (see [scrope_hl()]).
#' @param model Explicit fitting route. `"LN"` uses the LN NBGMM pipeline,
#'   `"HL"` uses the corrected first-order nested HL NBGMM pipeline for every
#'   analyzed gene, and `"PGMM"` uses the exact Poisson-gamma mixed-model path.
#'   PGMM supports subject- and cell-level fixed effects but omits separate
#'   cell-level overdispersion, so robust results are preferred for cell-level
#'   predictors. No route is silently changed to another route.
#' @param ln_estimator LN estimator. `"legacy_post_pml"` is the default and
#'   preserves the released post-PML behavior. `"joint_ln"` is opt-in and
#'   estimates coefficients and dispersions together. It returns model-based
#'   and donor-robust Wald inference from the certified selected root for genes
#'   that pass its numerical checks, including a bounded same-basin review in
#'   the 0.2.5.9000 calibration-candidate build. Boundary, distinct-basin,
#'   hard-stop, and unresolved cases retain missing uncertainty and test fields.
#'   This argument applies only when `model = "LN"`.
#' @param ln_operational_fallback Optional independent legacy fallback for
#'   joint-LN fits. `"none"` (the default) performs no fallback fitting.
#'   `"legacy_post_pml"` independently refits a gene when joint-LN inference is
#'   unavailable because a dispersion estimate reached its allowed boundary.
#'   The sidecar exposes only an operational legacy coefficient, dispersions,
#'   and diagnostics when the legacy fit passes its numerical checks; it never
#'   exposes covariance or inferential output. This option requires
#'   `model = "LN"`, `ln_estimator = "joint_ln"`, and
#'   `allow_per_gene_switch = FALSE`.
#'
#' @details
#' The public default remains `ln_estimator = "legacy_post_pml"`. The opt-in
#' `"joint_ln"` estimator jointly estimates fixed effects and dispersions from
#' the LN approximate marginal criterion. It reports a numerically verified
#' local solution; global uniqueness is not claimed.
#'
#' In the 0.2.5.9000 calibration-candidate build, an extra nearby uncertified
#' endpoint may undergo one bounded same-basin review. A resolved endpoint is
#' not certified and never replaces the selected root. It only removes that
#' endpoint's veto; estimates, model-based covariance, donor-robust covariance,
#' and Wald inference remain functions of the fully certified selected root.
#' The original strict decision is retained in per-gene provenance. This build
#' is intended for simulation calibration and is not a calibrated public
#' release.
#'
#' Because the three routes share one function signature, non-LN routes accept
#' explicit default values of LN-only compatibility arguments. A non-default
#' LN-only value is rejected rather than silently ignored, except where the
#' route-specific documentation explicitly promises a warning.
#'
#' For genes passing the joint-LN numerical checks, model-based and robust Wald
#' inference is returned. Otherwise, covariance, standard errors, Wald
#' statistics, and p-values are `NA`. Adjusted Score and likelihood-ratio tests
#' are unavailable for `"joint_ln"` pending calibration.
#'
#' PGMM LRT requires accepted unrestricted and constrained endpoints; Score
#' requires an accepted constrained null endpoint. Rejected endpoints cannot
#' obtain genuine test success from a positive raw likelihood difference.
#' Raw quantities may be retained as diagnostics when inference is unavailable.
#'
#' When the fallback is requested, eligible genes are independently refitted
#' with the default legacy estimator. A successful fallback exposes only the
#' operational coefficient and dispersions from that fit. It is judged by its
#' own numerical checks and never overwrites the primary joint-LN result.
#'
#' The separate `operational_fallback` component contains per-gene fallback
#' status, a long coefficient table, fitted dispersions, and technical
#' diagnostics. Use `fallback_status == "coefficient_available"` and
#' `coefficient_available` to identify an operational point estimate.
#' `inference_supported` is always `FALSE`; there are deliberately no fallback
#' covariance, standard-error, test-statistic, p-value, or q-value fields.
#'
#' For joint-LN, `availability` gives exactly one row per analyzed gene/request.
#' Its three states are `"inference"` for certified primary Wald inference,
#' `"estimate_only"` for an explicitly requested successful operational
#' coefficient with no uncertainty, and `"unavailable"`. A finite diagnostic
#' joint-LN endpoint is never promoted to `"estimate_only"`. One consolidated
#' availability warning is emitted when any row lacks primary inference.
#' Columns include `gene_id`, `gene`, `estimator`, `certified_root`,
#' `inference_status`, `strict_inference_status`,
#' `enabled_inference_status`, `resolution_status`, `inference_source`,
#' `newly_enabled`, `availability_transition`, `resolution_policy_version`,
#' `availability`, `reason_group`, `status_message`, `inference_supported`,
#' `estimate_source`, `attention_flag`, and `technical_code`. The compact
#' `provenance$per_gene` table also retains strict and enabled route codes,
#' failure reasons, policy identity, mode, and all veto classes. Earlier
#' `route`, `inference_available`, and `failure_reason` fields are retained for
#' compatibility.
#'
#' FDR control for joint-LN Score/Pace/LRT remains calibration-pending. The
#' fallback cannot enter an FDR procedure because it contains no p-values.
#' Internal route codes and numerical checks are retained only in provenance
#' and diagnostics.
#'
#' With `output_re = TRUE`, joint-LN donor modes are auxiliary fixed-beta
#' conditional PML predictions evaluated at the plug-in joint-LN parameters.
#' They do not alter coefficients, dispersions, covariance, tests, or
#' convergence.
#'
#' @param pgmm_control PGMM main unconstrained `nlminb` budget: a named list
#'   with positive integer `iter.max` and/or `eval.max`. Missing entries default
#'   to 500 each. Use `list(iter.max = 150L, eval.max = 200L)` to reproduce
#'   the earlier optimizer limits. This does not change tolerances, starts,
#'   constrained/profile fit controls, or retry policy, and is not a total
#'   fitting-work cap. It applies only to `model = "PGMM"`.
#'
#' @return A list with elements `summary`, `overdispersion`, `convergence`,
#'   `algorithm`, optional `random_effect`, `diagnostics`, optional
#'   route-specific `lrt_details` or `hl_adj_details`, structured `provenance`,
#'   and, for joint-LN, the one-row-per-gene `availability` table and, only when
#'   requested, `operational_fallback`. The availability table applies the
#'   precedence certified primary inference, requested successful operational
#'   coefficient, then unavailable. The fallback component contains `per_gene`
#'   status, long `coefficients`, `overdispersion`, and technical `diagnostics`,
#'   all explicitly marked as inference-unsupported. Fallback values never
#'   overwrite the primary result fields.
#'
#' @useDynLib scRoPE, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
scrope <- function(
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
    cutoff_cell = 20,
    kappa = 800,
    allow_per_gene_switch = FALSE,
    use_betas = "betae", # or "final_betas"
    output_re = FALSE,
    additional_tests = character(0),
    lrt = FALSE,
    lrt_contrasts = NULL,
    lrt_include_scale = FALSE,
    lrt_details = FALSE,
    keep_diagnostics = TRUE,
    model = c("LN", "HL", "PGMM"),
    ln_estimator = c("legacy_post_pml", "joint_ln"),
    ln_operational_fallback = c("none", "legacy_post_pml"),
    pgmm_control = NULL) {
  cutoff_cell_missing <- missing(cutoff_cell)
  kappa_missing <- missing(kappa)
  ln_estimator_missing <- missing(ln_estimator)
  ln_estimator <- match.arg(ln_estimator)
  ln_operational_fallback <- match.arg(ln_operational_fallback)
  model <- match.arg(model)
  if (!identical(model, "PGMM") && !is.null(pgmm_control)) {
    stop("pgmm_control applies only to model='PGMM'.")
  }
  if (identical(model, "PGMM")) pgmm_control <- pmg_resolve_control(pgmm_control)
  if (length(opt) != 1L || is.na(opt) || !opt %in% c("lbfgs", "trust")) {
    stop("opt must be 'lbfgs' or 'trust'.")
  }
  if (length(use_betas) != 1L || is.na(use_betas) ||
      !use_betas %in% c("betae", "final_betas")) {
    stop("use_betas must be 'betae' or 'final_betas'.")
  }

  if (identical(model, "HL")) {
    if (!identical(ln_operational_fallback, "none")) {
      stop("ln_operational_fallback applies only to model='LN'.")
    }
    if (!ln_estimator_missing && !identical(ln_estimator, "legacy_post_pml")) {
      stop("ln_estimator applies only to model='LN'; do not request joint_ln for HL.")
    }
    if (isTRUE(allow_per_gene_switch)) {
      stop("allow_per_gene_switch is unavailable for model='HL'; the requested HL route is used for every gene.")
    }
    if (!cutoff_cell_missing && !isTRUE(all.equal(cutoff_cell, 20))) {
      stop("cutoff_cell controls legacy LN refinement and cannot be set for model='HL'.")
    }
    if (!kappa_missing && !isTRUE(all.equal(kappa, 800))) {
      stop("kappa controls legacy LN refinement and cannot be set for model='HL'.")
    }
    return(hl_fo_public_scrope_hl(
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
      keep_diagnostics = keep_diagnostics
    ))
  }

  if (identical(model, "PGMM")) {
    if (!identical(ln_operational_fallback, "none")) {
      stop("ln_operational_fallback applies only to model='LN'.")
    }
    if (!ln_estimator_missing && !identical(ln_estimator, "legacy_post_pml")) {
      stop("ln_estimator applies only to model='LN'; do not request joint_ln for PGMM.")
    }
    if (isTRUE(allow_per_gene_switch)) {
      stop("allow_per_gene_switch is unavailable for model='PGMM'.")
    }
    if (!cutoff_cell_missing && !isTRUE(all.equal(cutoff_cell, 20))) {
      stop("cutoff_cell controls legacy LN refinement and cannot be set for model='PGMM'.")
    }
    if (!kappa_missing && !isTRUE(all.equal(kappa, 800))) {
      stop("kappa controls legacy LN refinement and cannot be set for model='PGMM'.")
    }
    if (identical(use_betas, "final_betas")) {
      stop("use_betas='final_betas' applies only to model='LN'.")
    }
    pgmm_result <- scrope_pmg_internal(
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
      cutoff_cell = cutoff_cell,
      kappa = kappa,
      allow_per_gene_switch = allow_per_gene_switch,
      use_betas = use_betas,
      output_re = output_re,
      additional_tests = additional_tests,
      lrt = lrt,
      lrt_contrasts = lrt_contrasts,
      lrt_include_scale = lrt_include_scale,
      lrt_details = lrt_details,
      keep_diagnostics = keep_diagnostics,
      pgmm_control = pgmm_control
    )
    pgmm_result$provenance <- list(
      schema_version = 2L,
      route = "PGMM",
      requested_model = "PGMM",
      effective_model = "PGMM",
      estimator = "pgmm_exact_marginal",
      optimizer = list(
        primary_backend = "nlminb",
        main_control = pgmm_control,
        scope = "main_unconstrained_fit",
        retry_requested = FALSE
      ),
      model_covariance = "pgmm_observed_information",
      robust_covariance = "pgmm_donor_sandwich",
      test_capabilities = data.frame(
        test = c("wald", "score", "lrt"),
        public_available = rep(TRUE, 3L),
        status = rep("available_when_requested", 3L),
        stringsAsFactors = FALSE
      )
    )
    return(pgmm_result)
  }

  if (!identical(ln_operational_fallback, "none") &&
      !identical(ln_estimator, "joint_ln")) {
    stop(paste0(
      "ln_operational_fallback='legacy_post_pml' requires ",
      "ln_estimator='joint_ln'."
    ))
  }
  if (!identical(ln_operational_fallback, "none") &&
      isTRUE(allow_per_gene_switch)) {
    stop(paste0(
      "ln_operational_fallback='legacy_post_pml' does not support ",
      "allow_per_gene_switch=TRUE."
    ))
  }
  if (identical(ln_estimator, "joint_ln") && isTRUE(allow_per_gene_switch)) {
    stop("ln_estimator='joint_ln' does not support allow_per_gene_switch=TRUE.")
  }
  eps <- 1e-06
  cell_init <- 1

  ######################
  # (A) Process Input  #
  ######################
  if (is.vector(count)) {
    count <- matrix(count, nrow = 1)
  }
  count <- as(count, "dgCMatrix")
  nind <- ncol(count)
  if (nind < 2) {
    stop("The count matrix must have >= 2 columns.")
  }
  if (is.null(rownames(count))) {
    rownames(count) <- paste0("Gene", seq_len(nrow(count)))
  }
  gname <- rownames(count)
  ngene <- nrow(count)

  if (is.null(pred)) {
    pred <- matrix(1, nrow = nind, ncol = 1)
    sds <- 1
    intcol <- 1
    predn <- "(Intercept)"
  } else {
    predn <- colnames(pred)
    if ((!is.null(predn)) && length(unique(predn)) < ncol(pred)) {
      stop("All columns of the design matrix must have unique names.")
    }
    predC <- center_m(as.matrix(pred))
    sds <- predC$sds
    if ((sum(sds == 0) > 1) || any(sds < 0)) {
      stop("Some predictors have zero variation or are zero vectors.")
    }
    if (sum(sds == 0) == 0) {
      stop("Design matrix must have an intercept column.")
    }
    intcol <- which(sds == 0)
    pred <- as.matrix(predC$pred)
    if (Matrix::rankMatrix(pred) < ncol(pred)) {
      warning("Predictors are collinear or rank-deficient.")
    }
    if (nrow(pred) != nind) {
      stop("Design matrix row != count matrix column.")
    }
  }
  nb <- ncol(pred)
  public_scale <- as.numeric(sds)
  public_scale[intcol] <- 1
  coef_internal_to_public <- diag(1 / public_scale, nrow = nb, ncol = nb)

  extra_tests <- additional_tests
  if (!is.null(extra_tests)) {
    extra_tests <- unique(tolower(extra_tests))
  }
  if (length(extra_tests) == 1 && identical(extra_tests, "")) {
    extra_tests <- character(0)
  }
  if (lrt) {
    extra_tests <- union(extra_tests, "lrt")
  }
  allowed_tests <- c("score", "lrt")
  invalid_tests <- setdiff(extra_tests, allowed_tests)
  if (length(invalid_tests) > 0) {
    stop(sprintf(
      "Unsupported entries in additional_tests: %s. Allowed values are %s.",
      paste(invalid_tests, collapse = ", "),
      paste(allowed_tests, collapse = ", ")
    ))
  }

  want_score <- "score" %in% extra_tests
  want_lrt <- "lrt" %in% extra_tests
  need_constrained <- want_score || want_lrt

  if ((lrt_details || !is.null(lrt_contrasts)) && !need_constrained) {
    extra_tests <- union(extra_tests, "lrt")
    want_lrt <- TRUE
    need_constrained <- TRUE
  }
  requested_score <- want_score
  requested_lrt <- want_lrt
  optional_test_skip_reason <- NA_character_

  if (need_constrained && opt != "lbfgs" &&
      !identical(ln_estimator, "joint_ln")) {
    optional_test_skip_reason <- "not_computed_opt_requires_lbfgs"
    warning(
      "Legacy LN adjusted Score/LRT require opt='lbfgs'; skipping those requested tests.",
      call. = FALSE
    )
    extra_tests <- setdiff(extra_tests, c("score", "lrt"))
    want_score <- FALSE
    want_lrt <- FALSE
    need_constrained <- FALSE
  }
  if (need_constrained && allow_per_gene_switch) {
    optional_test_skip_reason <- "not_computed_per_gene_switch"
    warning(
      "Adjusted Score/LRT are unavailable when allow_per_gene_switch=TRUE; skipping those requested tests.",
      call. = FALSE
    )
    extra_tests <- setdiff(extra_tests, c("score", "lrt"))
    want_score <- FALSE
    want_lrt <- FALSE
    need_constrained <- FALSE
  }

  if (!want_lrt) {
    lrt_details <- FALSE
  }

  lrt_defs <- list()
  lrt_names <- character(0)
  if (want_score || want_lrt || !is.null(lrt_contrasts)) {
    coef_names <- if (is.null(predn)) paste0("X", seq_len(nb)) else predn
    non_intercept <- setdiff(seq_len(nb), intcol)
    add_contrast <- function(name, L, b) {
      L <- as.matrix(L)
      if (ncol(L) != nb) {
        stop("Each L must have ncol equal to number of predictors.")
      }
      if (is.null(dim(b))) {
        b <- rep(b, length.out = nrow(L))
      }
      if (length(b) != nrow(L)) {
        stop("Length of b must match number of rows of L.")
      }
      nm <- as.character(name)
      if (is.na(nm) || !nzchar(nm)) {
        nm <- paste0("contrast", length(lrt_defs) + 1)
      }
      base <- nm
      idx <- 1
      while (nm %in% names(lrt_defs)) {
        idx <- idx + 1
        nm <- paste0(base, "_", idx)
      }
      lrt_defs[[nm]] <<- list(L = L, b = b, df = nrow(L), label = nm)
    }
    for (j in non_intercept) {
      L <- matrix(0, nrow = 1, ncol = nb)
      L[1, j] <- 1
      add_contrast(coef_names[j], L, 0)
    }
    if (!is.null(lrt_contrasts)) {
      if (!is.list(lrt_contrasts)) {
        stop("lrt_contrasts must be a list.")
      }
      for (cc in lrt_contrasts) {
        if (is.null(cc)) next
        if (is.matrix(cc)) {
          add_contrast(paste0("contrast", length(lrt_defs) + 1), cc, rep(0, nrow(cc)))
        } else if (is.numeric(cc)) {
          add_contrast(paste0("contrast", length(lrt_defs) + 1), matrix(cc, nrow = 1), 0)
        } else if (is.list(cc)) {
          L <- cc$L
          if (is.null(L)) {
            stop("Each contrast list must include L.")
          }
          b <- if (!is.null(cc$b)) cc$b else rep(0, if (is.matrix(L)) nrow(L) else 1)
          add_contrast(if (!is.null(cc$name)) cc$name else paste0("contrast", length(lrt_defs) + 1), L, b)
        } else {
          stop("Unsupported element in lrt_contrasts.")
        }
      }
    }
    lrt_names <- names(lrt_defs)
  }

  # Offset
  if (is.null(offset)) {
    of_re <- cv_offset(0, 0, nind)
  } else {
    if (length(offset) != nind) {
      stop("Offset length != #columns of count.")
    }
    of_re <- cv_offset(as.double(offset), 1, nind)
  }
  offset <- of_re$offset
  moffset <- log(of_re$mexpoffset)
  cv <- of_re$cv
  cv2 <- cv * cv

  # IDs
  if (length(id) != nind) {
    stop("ID length != #columns of count.")
  }
  id <- as.character(id)
  levels <- unique(id)
  id <- as.numeric(factor(id, levels = levels))
  if (is.unsorted(id)) {
    stop("Cells for the same subject must be contiguous.")
  }
  k <- length(levels)
  fid <- which(c(1, diff(id)) == 1)
  fid <- as.integer(c(fid, nind + 1))

  maxcore <- max(c(1, length(parallelly::availableWorkers()) - 1))
  if (ncore > maxcore && verbose) {
    cat("Specified ncore exceeds available cores. Using detected # of cores.\n")
    ncore <- maxcore
  }
  options(future.globals.maxSize = fmaxsize)

  mfs <- nind / k
  if (mfs < 30 && verbose) {
    cat("Average #cells/subject = ", round(mfs, 2), "< 30 => LN may be suboptimal.\n")
  }

  cumsumy <- call_cumsumy(count, fid - 1, k, ngene)
  if (ngene == 1) {
    cumsumy <- matrix(cumsumy, ncol = k)
  }
  posind <- lapply(seq_len(ngene), function(x) which(cumsumy[x, ] > 0))

  gid2 <- which(tabulate(count@i + 1L, ngene) >= mincp)
  gidc <- which((rowSums(cumsumy) / nind) > cpc)
  gid <- intersect(gid2, gidc)
  lgid <- length(gid)
  if (verbose) {
    cat("Remove", ngene - lgid, "genes with low expression.\n")
    cat("Analyzing", lgid, "genes with", k, "subjects and", nind, "cells.\n")
  }
  if (lgid == 0) {
    stop("No gene passed the filtering.")
  }
  count <- Matrix::t(count)

  cell_ind <- get_cell(pred, fid - 1, nb, k)
  cell_ind <- which(cell_ind > 0)
  ncell <- length(cell_ind)

  ######################
  # Parallel Setup     #
  ######################
  registerDoFuture()
  if (ncore == 1) {
    plan(sequential)
  } else {
    cls <- parallelly::makeClusterPSOCK(
      ncore,
      rscript_libs = scrope_psock_library_paths()
    )
    plan(cluster, workers = cls)
    cluster_open <- TRUE
    on.exit({
      if (isTRUE(cluster_open)) {
        try(parallel::stopCluster(cls), silent = TRUE)
        plan(sequential)
      }
    }, add = TRUE)
  }

  ###################################################################
  # (D) Main Loop: LN => optional HL => LN+HL
  ###################################################################
  ctx <- list(
    pred = pred,
    offset = offset,
    fid = fid,
    cumsumy = cumsumy,
    posind = posind,
    nb = nb,
    nind = nind,
    k = k,
    intcol = intcol,
    moffset = moffset,
    min = min,
    max = max,
    allow_per_gene_switch = allow_per_gene_switch,
    use_betas = use_betas,
    cutoff_cell = cutoff_cell,
    cell_ind = cell_ind,
    ncell = ncell,
    cv2 = cv2,
    kappa = kappa,
    eps = eps,
    opt = opt,
    cell_init = cell_init,
    output_re = output_re,
    id = id,
    ln_estimator = ln_estimator,
    ln_candidate_resolution_mode = if (identical(
        ln_estimator,
        "joint_ln"
      )) {
      "calibration_candidate_enabled"
    } else {
      "not_applicable"
    },
    coef_internal_to_public = coef_internal_to_public
  )

  gene_results <- foreach(i = gid) %dorng% {
    (function() {
    posv <- call_posindy(count, i - 1, nind)
    if (identical(ctx$ln_estimator, "joint_ln")) {
      on.exit(
        ln_joint_objective_instance_registry(
          "remove_matching",
          gene_index = i,
          posv = posv,
          ctx = ctx
        ),
        add = TRUE
      )
    }
    unconstrained <- fit_gene_unconstrained(i, posv, ctx)

    if ((want_score || want_lrt) && length(lrt_names) > 0) {
      lrt_neg_tol <- 1e-6

      build_lrt_cache <- function(current_unconstrained) {
        cache <- vector("list", length(lrt_names))
        names(cache) <- lrt_names
        neg_detected <- FALSE
        joint_ln_pending <- identical(ctx$ln_estimator, "joint_ln")
        empty_score_result <- function(L, msg = NULL, reason = NULL) {
          d0 <- nrow(L)
          if (!is.null(msg)) {
            warning(msg)
          }
          failure_reason <- if (is.null(reason)) {
            NA_character_
          } else {
            normalise_failure(reason, "score")
          }
          list(
            method = "adj-score",
            df = d0,
            stat = NA_real_,
            pval = NA_real_,
            ridge = NA_real_,
            success = NA,
            score = matrix(NA_real_, d0, 1),
            Spsi = matrix(NA_real_, d0, d0),
            Gpsi = matrix(NA_real_, d0, d0),
            Hpsi = matrix(NA_real_, d0, d0),
            scalar_check = NA_real_,
            failure_reason = failure_reason
          )
        }
        empty_ridge <- list(Hpsi = NA_real_, Gpsi = NA_real_, wald = NA_real_)
        joint_pending_lrt <- function(unconstrained_fit, constrained_fit) {
          diagnostics <- constrained_fit$diagnostics
          if (is.null(diagnostics)) diagnostics <- list()
          gate <- constrained_fit$lrt_gate
          if (is.null(gate)) gate <- list()
          raw_wP <- if (is.finite(unconstrained_fit$loglik_ln) &&
              is.finite(constrained_fit$loglik_ln)) {
            2 * (unconstrained_fit$loglik_ln - constrained_fit$loglik_ln)
          } else {
            NA_real_
          }
          reduced_gradient <- if (!is.null(diagnostics$reduced_gradient_max_abs)) {
            as.numeric(diagnostics$reduced_gradient_max_abs[[1L]])
          } else {
            NA_real_
          }
          newton_step <- if (!is.null(diagnostics$max_abs_newton_step_scaled)) {
            as.numeric(diagnostics$max_abs_newton_step_scaled[[1L]])
          } else {
            NA_real_
          }
          hessian_rcond <- if (!is.null(diagnostics$H_red_rcond)) {
            as.numeric(diagnostics$H_red_rcond[[1L]])
          } else {
            NA_real_
          }
          constraint_residual <- if (!is.null(diagnostics$constraint_residual_scaled)) {
            as.numeric(diagnostics$constraint_residual_scaled[[1L]])
          } else {
            NA_real_
          }
          hessian_ok <- isTRUE(diagnostics$H_red_positive_definite)
          diagnostics_finite <- all(is.finite(c(
            reduced_gradient,
            newton_step,
            hessian_rcond,
            constraint_residual
          )))
          list(
            wP = raw_wP,
            wU = NA_real_,
            qP = NA_real_,
            wP_minus_qP = NA_real_,
            wP_over_qP = NA_real_,
            scale = NA_real_,
            wP_adjusted = NA_real_,
            score = NA_real_,
            Hpsi = NA_real_,
            Gpsi = NA_real_,
            Gbeta = NA_real_,
            Jpsi = NA_real_,
            Vpsi_model = NA_real_,
            Vpsi_godambe = NA_real_,
            method = "unavailable_joint_ln",
            ridge = empty_ridge,
            fallback = FALSE,
            success = FALSE,
            adjusted_available = FALSE,
            adjusted_failure_reason = "joint_ln_calibration_pending",
            gate = list(
              converged = isTRUE(unconstrained_fit$inference_available) &&
                isTRUE(constrained_fit$inference_available),
              convergence_basis = "joint_ln_endpoint_pair",
              unrestricted_route = unconstrained_fit$route,
              unrestricted_failure_reason = unconstrained_fit$failure_reason,
              constrained_route = constrained_fit$route,
              constrained_failure_reason = constrained_fit$failure_reason,
              legacy_converged = FALSE,
              legacy_conv_code = NA_integer_,
              information_scaled_eligible = FALSE,
              information_scaled_convergence_ok = FALSE,
              information_scaled_accepted = FALSE,
              H_red_positive_definite = hessian_ok,
              diagnostics_finite = diagnostics_finite,
              H_red_rcond = hessian_rcond,
              reduced_gradient_max_abs = reduced_gradient,
              max_abs_newton_step_scaled = newton_step,
              gradient_quadratic_lr_bound = NA_real_,
              newton_step_scaled_ok = is.finite(newton_step) &&
                newton_step <= ln_joint_thresholds()$newton,
              gradient_quadratic_lr_bound_ok = FALSE,
              reduced_parameterization = "zeta, log(sigma2), log(gamma_nb_size)",
              constraint_satisfied = isTRUE(gate$constraint_satisfied),
              constraint_residual_max_abs = constraint_residual,
              nuisance_hessian_invertible = hessian_ok,
              interior = isTRUE(gate$interior),
              nuisance_interior = isTRUE(gate$nuisance_interior),
              profiled_curvature_positive = hessian_ok,
              profiled_variability_positive = NA,
              scale_positive = FALSE
            ),
            failure_reason = normalise_failure(
              "joint_ln_calibration_pending",
              "lrt"
            )
          )
        }

        for (idx_nm in seq_along(lrt_names)) {
          nm <- lrt_names[idx_nm]
          def <- lrt_defs[[nm]]
          L_fit <- if (joint_ln_pending) {
            def$L %*% ctx$coef_internal_to_public
          } else {
            def$L
          }
          constrained_fit <- fit_gene_constrained(i, posv, ctx, L_fit, def$b)
          lrt_res <- NULL
          if (want_lrt) {
            if (joint_ln_pending) {
              lrt_res <- joint_pending_lrt(current_unconstrained, constrained_fit)
            } else {
              ln_gate <- ln_lrt_gate_override(constrained_fit, def$L, def$b)
              lrt_res <- compute_profile_lrt(
                current_unconstrained,
                constrained_fit,
                def$L,
                def$b,
                lrt_gate_override = ln_gate
              )
              if (!is.na(lrt_res$wP) && lrt_res$wP < -lrt_neg_tol) {
                neg_detected <- TRUE
              }
            }
          }
          gene_label <- if (!is.null(gname) && length(gname) >= i) gname[i] else paste0("gene_", i)
          score_res <- empty_score_result(def$L)
          if (want_score) {
            if (joint_ln_pending) {
              score_res$method <- "unavailable_joint_ln"
              score_res$success <- FALSE
              score_res$failure_reason <- normalise_failure(
                "joint_ln_calibration_pending",
                "score"
              )
            } else {
              score_res <- tryCatch(
                compute_adj_score(constrained_fit, def$L),
                error = function(e) {
                  empty_score_result(
                    def$L,
                    msg = sprintf(
                      "Adjusted score test failed for gene %s (contrast %s): %s",
                      gene_label,
                      def$label,
                      conditionMessage(e)
                    ),
                    reason = "undefined_stat"
                  )
                }
              )
            }
          }
          stat_adj <- if (!is.null(lrt_res)) lrt_res$wP_adjusted else NA_real_
          df <- def$df
          pval <- if (!is.null(lrt_res) && !is.na(stat_adj)) {
            pchisq(max(stat_adj, 0), df = df, lower.tail = FALSE)
          } else {
            NA_real_
          }
          cache[[idx_nm]] <- list(
            constrained = constrained_fit,
            df = df,
            raw = lrt_res,
            stat = stat_adj,
            pval = pval,
            success = if (!is.null(lrt_res)) lrt_res$success else NA,
            failure_reason = if (!is.null(lrt_res)) lrt_res$failure_reason else NA_character_,
            scale = if (!is.null(lrt_res)) lrt_res$scale else NA_real_,
            wP = if (!is.null(lrt_res)) lrt_res$wP else NA_real_,
            wU = if (!is.null(lrt_res)) lrt_res$wU else NA_real_,
            qP = if (!is.null(lrt_res)) lrt_res$qP else NA_real_,
            method = if (!is.null(lrt_res)) lrt_res$method else NA_character_,
            ridge = if (!is.null(lrt_res)) lrt_res$ridge else empty_ridge,
            fallback = if (!is.null(lrt_res)) lrt_res$fallback else NA,
            score_stat = score_res$stat,
            score_pval = score_res$pval,
            score_ridge = score_res$ridge,
            score_success = score_res$success,
            score_failure = score_res$failure_reason,
            score_test = if (want_score) score_res else NULL
          )
        }
        list(cache = cache, neg = neg_detected)
      }

      lrt_build <- build_lrt_cache(unconstrained)
      lrt_cache <- lrt_build$cache

      if (want_lrt && !identical(ctx$ln_estimator, "joint_ln") && lrt_build$neg) {
        extra_starts <- lapply(lrt_cache, function(entry) {
          theta_ln <- entry$constrained$theta_ln
          if (is.null(theta_ln) || is.null(theta_ln$beta)) {
            return(NULL)
          }
          c(theta_ln$beta, theta_ln$subVar, theta_ln$cellVar)
        })
        extra_starts <- extra_starts[!vapply(extra_starts, is.null, logical(1))]

        if (length(extra_starts) > 0) {
          refined_fit <- fit_gene_unconstrained(i, posv, ctx, extra_starts = extra_starts)

          merge_refinement_meta <- function(orig, extra) {
            if (is.null(extra)) {
              return(orig)
            }
            if (is.null(orig)) {
              return(extra)
            }
            orig$refined <- isTRUE(orig$refined) || isTRUE(extra$refined)
            orig$grad_norm_start <- extra$grad_norm_start
            orig$grad_norm_final <- extra$grad_norm_final
            if (is.null(orig$extra_attempts)) orig$extra_attempts <- 0L
            if (is.null(orig$extra_success)) orig$extra_success <- 0L
            if (is.null(extra$extra_attempts)) extra$extra_attempts <- 0L
            if (is.null(extra$extra_success)) extra$extra_success <- 0L
            orig$extra_attempts <- orig$extra_attempts + extra$extra_attempts
            orig$extra_success <- orig$extra_success + extra$extra_success
            orig
          }

          update_full <- is.finite(refined_fit$loglik) && refined_fit$loglik > unconstrained$loglik + 1e-8
          update_ln <- is.finite(refined_fit$loglik_ln) && refined_fit$loglik_ln > unconstrained$loglik_ln + 1e-8

          if (update_full) {
            unconstrained <- refined_fit
            lrt_build <- build_lrt_cache(unconstrained)
            lrt_cache <- lrt_build$cache
          } else {
            unconstrained$ln_refinement <- merge_refinement_meta(unconstrained$ln_refinement, refined_fit$ln_refinement)
            if (update_ln) {
              unconstrained$loglik_ln <- refined_fit$loglik_ln
              unconstrained$theta_ln <- refined_fit$theta_ln
              lrt_build <- build_lrt_cache(unconstrained)
              lrt_cache <- lrt_build$cache
            }
          }
        }
      }

      lrt_list <- list()
      for (idx_nm in seq_along(lrt_names)) {
        nm <- lrt_names[idx_nm]
        entry <- lrt_cache[[idx_nm]]
        entry$success <- if (want_lrt) isTRUE(entry$success) else NA
        entry$failure_reason <- if (want_lrt) entry$failure_reason else NA_character_
        entry$score_success <- if (want_score) isTRUE(entry$score_success) else NA
        entry$score_failure <- if (want_score) entry$score_failure else NA_character_
        lrt_list[[nm]] <- entry
      }
      unconstrained$lrt <- lrt_list
    }

    unconstrained
    })()
  }

  fallback_eligibility <- NULL
  fallback_results <- NULL
  if (identical(ln_operational_fallback, "legacy_post_pml")) {
    fallback_eligibility <- lapply(
      gene_results,
      ln_operational_fallback_eligibility
    )
    fallback_ctx <- ctx
    fallback_ctx$ln_estimator <- "legacy_post_pml"
    fallback_ctx$allow_per_gene_switch <- FALSE
    fallback_ctx$output_re <- FALSE
    fallback_results <- foreach(
      fallback_position = seq_along(gid),
      .options.RNG = 20260721L
    ) %dorng% {
      if (!isTRUE(fallback_eligibility[[fallback_position]]$eligible)) {
        ln_operational_empty_assessment()
      } else {
        fallback_gene <- gid[[fallback_position]]
        tryCatch({
          fallback_posv <- call_posindy(count, fallback_gene - 1L, nind)
          ln_operational_attempt_gene(
            fallback_gene,
            fallback_posv,
            fallback_ctx
          )
        }, error = function(e) {
          fallback_error <- ln_operational_empty_assessment("legacy_fit_error")
          fallback_error$diagnostics <- list(
            error_message = conditionMessage(e)
          )
          fallback_error
        })
      }
    }
  }

  if (ncore > 1) {
    parallel::stopCluster(cls)
    cluster_open <- FALSE
    plan(sequential)
    gc()
  }

  re_mat <- do.call(cbind, lapply(gene_results, `[[`, "vec"))
  re_mat <- t(re_mat)
  godambe_list <- lapply(gene_results, `[[`, "godambe")
  re_all <- as.data.frame(re_mat)

  beta_names <- paste0(
    "logFC_",
    if (is.null(predn)) paste0("X", seq_len(nb)) else predn
  )
  naive_var_names <- paste0("naive_var_", beta_names)
  robust_var_names <- paste0("robust_var_", beta_names)

  col_offset <- 0
  colnames(re_all)[1:nb] <- beta_names
  col_offset <- nb
  colnames(re_all)[col_offset + 1] <- "Subject"
  colnames(re_all)[col_offset + 2] <- "Cell"
  col_offset <- col_offset + 2

  iv_naive <- seq(col_offset + 1, col_offset + nb)
  for (j in seq_len(nb)) {
    colnames(re_all)[iv_naive[j]] <- naive_var_names[j]
  }
  col_offset <- col_offset + nb

  iv_robust <- seq(col_offset + 1, col_offset + nb)
  for (j in seq_len(nb)) {
    colnames(re_all)[iv_robust[j]] <- robust_var_names[j]
  }
  col_offset <- col_offset + nb

  colnames(re_all)[col_offset + 1] <- "convergence"
  colnames(re_all)[col_offset + 2] <- "algorithm"
  colnames(re_all)[col_offset + 3] <- "loglik"
  col_offset <- col_offset + 3

  if (output_re) {
    re_names <- paste0("random_", seq_len(k))
    for (rr in seq_len(k)) {
      colnames(re_all)[col_offset + rr] <- re_names[rr]
    }
    col_offset <- col_offset + k
  }

  re_all$gene_id <- gid
  re_all$gene <- gname[gid]

  joint_inference_available <- if (identical(ln_estimator, "joint_ln")) {
    vapply(gene_results, function(gr) isTRUE(gr$inference_available), logical(1))
  } else {
    rep(TRUE, length(gene_results))
  }
  joint_wald_failure <- if (identical(ln_estimator, "joint_ln")) {
    vapply(gene_results, function(gr) {
      reason <- gr$failure_reason
      if (is.null(reason) || !length(reason) || is.na(reason[[1L]])) {
        "no_certified_candidate"
      } else {
        as.character(reason[[1L]])
      }
    }, character(1))
  } else {
    rep(NA_character_, length(gene_results))
  }

  for (j in seq_len(nb)) {
    coef_label <- if (is.null(predn)) paste0("X", j) else predn[j]

    ## naive Wald
    stat_naive <- (re_all[[beta_names[j]]]^2) / re_all[[naive_var_names[j]]]
    p_naive_j <- pchisq(stat_naive, 1, lower.tail = FALSE)
    re_all[[paste0("p_naive_", coef_label)]] <- p_naive_j

    ## robust Wald
    var_robust <- re_all[[robust_var_names[j]]]
    stat_robust <- (re_all[[beta_names[j]]]^2) / var_robust
    valid_var <- is.finite(var_robust) & var_robust > 0 &
      joint_inference_available
    stat_robust[!valid_var] <- NA_real_
    p_robust_j <- pchisq(stat_robust, 1, lower.tail = FALSE)

    success_vec <- rep(NA, length(stat_robust))
    success_vec[valid_var & is.finite(stat_robust)] <- TRUE
    success_vec[!(valid_var & is.finite(stat_robust))] <- FALSE

    failure_vec <- rep(NA_character_, length(stat_robust))
    failure_vec[!valid_var] <- normalise_failure("singular_cov", "wald")
    joint_unavailable <- !joint_inference_available
    failure_vec[joint_unavailable] <- normalise_failure(
      joint_wald_failure[joint_unavailable],
      "wald"
    )
    failure_vec[valid_var & !is.finite(stat_robust)] <- normalise_failure("undefined_stat", "wald")

    re_all[[paste0("p_robust_", coef_label)]] <- p_robust_j
    re_all[[paste0("wald_stat_", coef_label)]] <- stat_robust
    re_all[[paste0("wald_success_", coef_label)]] <- success_vec
    re_all[[paste0("wald_failure_", coef_label)]] <- failure_vec
  }

  ###########################################
  # Convert naive/robust var => SE, rename  #
  ###########################################
  re_all[iv_naive] <- sqrt(re_all[iv_naive])
  re_all[iv_robust] <- sqrt(re_all[iv_robust])

  naive_se_names <- sub("naive_var_logFC_", "se_", naive_var_names)
  naive_se_names <- sub("naive_var_", "se_", naive_se_names)

  robust_se_names <- sub("robust_var_logFC_", "se_robust_", robust_var_names)
  robust_se_names <- sub("robust_var_", "se_robust_", robust_se_names)

  colnames(re_all)[iv_naive] <- naive_se_names
  colnames(re_all)[iv_robust] <- robust_se_names

  ###########################################
  # Final scaling for betas + SE columns    #
  ###########################################
  sds[intcol] <- 1
  scale_cols <- c(seq_len(nb), iv_naive, iv_robust)
  scvec <- c(sds, sds, sds)
  re_all[, scale_cols] <- t(t(re_all[, scale_cols]) / scvec)

  ####################################
  # Build final output => summary    #
  ####################################
  overdisp_df <- re_all[, c("Subject", "Cell")]
  conv_vec <- re_all$convergence
  alg_vec <- re_all$algorithm
  loglik_vec <- re_all$loglik

  drop_cols <- c("Subject", "Cell", "convergence", "algorithm", "loglik")
  if (output_re) {
    re_names <- paste0("random_", seq_len(k))
    drop_cols <- c(drop_cols, re_names)
  }
  keep_cols <- setdiff(names(re_all), drop_cols)
  summary_df <- re_all[, keep_cols, drop = FALSE]

  collect_status_rows <- function(success_prefix, failure_prefix, test_label, norm_type) {
    success_cols <- grep(paste0("^", success_prefix), names(summary_df), value = TRUE)
    if (!length(success_cols)) {
      return(NULL)
    }
    do.call(rbind, lapply(success_cols, function(success_col) {
      lab <- sub(paste0("^", success_prefix), "", success_col)
      failure_col <- paste0(failure_prefix, lab)
      failure_vals <- if (failure_col %in% names(summary_df)) {
        normalise_failure(summary_df[[failure_col]], norm_type)
      } else {
        rep(NA_character_, nrow(summary_df))
      }
      data.frame(
        gene_id = summary_df$gene_id,
        gene = summary_df$gene,
        contrast = lab,
        test = test_label,
        success = as.logical(summary_df[[success_col]]),
        failure_reason = failure_vals,
        stringsAsFactors = FALSE
      )
    }))
  }

  if ((want_score || want_lrt) && length(lrt_names) > 0) {
    sanitize_name <- function(x) {
      gsub("[^A-Za-z0-9]+", "_", x)
    }
    for (nm in lrt_names) {
      def <- lrt_defs[[nm]]
      key <- sanitize_name(def$label)
      if (want_lrt) {
        stats_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_real_)
          }
          gr$lrt[[nm]]$stat
        }, numeric(1))
        p_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_real_)
          }
          gr$lrt[[nm]]$pval
        }, numeric(1))
        summary_df[[paste0("lrt_adj_", key)]] <- stats_vec
        summary_df[[paste0("p_lrt_adj_", key)]] <- p_vec
        if (lrt_include_scale) {
          scale_vec <- vapply(gene_results, function(gr) {
            if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
              return(NA_real_)
            }
            gr$lrt[[nm]]$scale
          }, numeric(1))
          summary_df[[paste0("lrt_scale_", key)]] <- scale_vec
        }
        method_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_character_)
          }
          gr$lrt[[nm]]$method
        }, character(1))
        summary_df[[paste0("lrt_method_", key)]] <- method_vec
        success_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA)
          }
          gr$lrt[[nm]]$success
        }, logical(1))
        failure_raw <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_character_)
          }
          gr$lrt[[nm]]$failure_reason
        }, character(1))
        failure_norm <- rep(NA_character_, length(failure_raw))
        idx_fail <- !is.na(failure_raw) & nzchar(failure_raw)
        if (any(idx_fail)) {
          failure_norm[idx_fail] <- normalise_failure(failure_raw[idx_fail], "lrt")
        }
        summary_df[[paste0("lrt_success_", key)]] <- success_vec
        summary_df[[paste0("lrt_failure_", key)]] <- failure_norm
      }
      if (want_score) {
        score_stats_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_real_)
          }
          gr$lrt[[nm]]$score_stat
        }, numeric(1))
        score_p_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_real_)
          }
          gr$lrt[[nm]]$score_pval
        }, numeric(1))
        score_ridge_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_real_)
          }
          gr$lrt[[nm]]$score_ridge
        }, numeric(1))
        score_success_vec <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA)
          }
          gr$lrt[[nm]]$score_success
        }, logical(1))
        summary_df[[paste0("score_adj_", key)]] <- score_stats_vec
        summary_df[[paste0("p_score_adj_", key)]] <- score_p_vec
        summary_df[[paste0("score_ridge_", key)]] <- score_ridge_vec
        summary_df[[paste0("score_success_", key)]] <- score_success_vec
        score_failure_raw <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) {
            return(NA_character_)
          }
          gr$lrt[[nm]]$score_failure
        }, character(1))
        score_failure_norm <- rep(NA_character_, length(score_failure_raw))
        idx_score <- !is.na(score_failure_raw) & nzchar(score_failure_raw)
        if (any(idx_score)) {
          score_failure_norm[idx_score] <- normalise_failure(score_failure_raw[idx_score], "score")
        }
        summary_df[[paste0("score_failure_", key)]] <- score_failure_norm
      }
    }
  }

  status_blocks <- list(
    collect_status_rows("wald_success_", "wald_failure_", "wald", "wald"),
    collect_status_rows("lrt_success_", "lrt_failure_", "lrt", "lrt"),
    collect_status_rows("score_success_", "score_failure_", "score", "score")
  )
  status_blocks <- Filter(Negate(is.null), status_blocks)
  test_status_df <- if (length(status_blocks) > 0) {
    do.call(rbind, status_blocks)
  } else {
    data.frame(
      gene_id = integer(0),
      gene = character(0),
      contrast = character(0),
      test = character(0),
      success = logical(0),
      failure_reason = character(0),
      stringsAsFactors = FALSE
    )
  }

  status_vec <- rep("ok", nrow(summary_df))
  detail_vec <- rep(NA_character_, nrow(summary_df))
  if (nrow(test_status_df) > 0) {
    for (i in seq_len(nrow(summary_df))) {
      subset_rows <- test_status_df[test_status_df$gene_id == summary_df$gene_id[i], , drop = FALSE]
      if (nrow(subset_rows) == 0) {
        next
      }
      fail_rows <- subset_rows[!is.na(subset_rows$success) & !subset_rows$success, , drop = FALSE]
      if (nrow(fail_rows) == 0) {
        next
      }
      test_types <- unique(fail_rows$test)
      status_vec[i] <- paste(test_types, collapse = ",")
      detail_entries <- tapply(fail_rows$failure_reason, fail_rows$test, function(vals) {
        clean_vals <- unique(na.omit(vals))
        if (length(clean_vals) == 0) {
          "unknown"
        } else {
          paste(clean_vals, collapse = ",")
        }
      })
      detail_vec[i] <- paste(paste0(names(detail_entries), ":", detail_entries), collapse = ";")
    }
  }
  summary_df$test_status <- status_vec
  summary_df$test_failure_detail <- detail_vec

  if (want_lrt && length(lrt_names) > 0) {
    method_cols <- grep("^lrt_method_", names(summary_df), value = TRUE)
    if (length(method_cols) > 0) {
      method_vals <- unlist(summary_df[method_cols], use.names = FALSE)
      method_vals <- method_vals[!is.na(method_vals)]
      if (length(method_vals) > 0) {
        count_fallback <- sum(method_vals %in% c("wald", "wald_ridge", "lrt_fail"))
        count_ridge <- sum(method_vals %in% c("lrt_ridge", "wald_ridge"))
        if (count_fallback > 0) {
          warning(sprintf(
            "Adjusted LRT fell back to Wald for %d gene-contrast(s); see lrt_method_* columns.",
            count_fallback
          ))
        }
        if (count_ridge > 0) {
          warning(sprintf(
            "Adjusted LRT applied ridge regularisation for %d gene-contrast(s).",
            count_ridge
          ))
        }
      }
    }
  }

  result_list <- list(
    summary        = summary_df,
    overdispersion = overdisp_df,
    convergence    = conv_vec,
    algorithm      = alg_vec
  )
  if (output_re) {
    result_list$random_effect <- re_all[, re_names, drop = FALSE]
  } else {
    result_list$random_effect <- NULL
  }

  if (identical(ln_estimator, "joint_ln")) {
    joint_status <- lapply(gene_results, ln_joint_status_record)
    per_gene_provenance <- data.frame(
      gene_id = gid,
      gene = gname[gid],
      estimator = vapply(joint_status, `[[`, character(1), "estimator"),
      certified_root = vapply(
        joint_status, `[[`, logical(1), "certified_root"
      ),
      inference_status = vapply(
        joint_status, `[[`, character(1), "inference_status"
      ),
      strict_inference_status = vapply(
        joint_status, `[[`, character(1), "strict_inference_status"
      ),
      enabled_inference_status = vapply(
        joint_status, `[[`, character(1), "enabled_inference_status"
      ),
      resolution_status = vapply(
        joint_status, `[[`, character(1), "resolution_status"
      ),
      inference_source = vapply(
        joint_status, `[[`, character(1), "inference_source"
      ),
      newly_enabled = vapply(
        joint_status, `[[`, logical(1), "newly_enabled"
      ),
      availability_transition = vapply(
        joint_status, `[[`, character(1), "availability_transition"
      ),
      resolution_veto_classes = vapply(
        joint_status, `[[`, character(1), "resolution_veto_classes"
      ),
      resolution_policy_version = vapply(
        joint_status, `[[`, character(1), "resolution_policy_version"
      ),
      resolution_mode = vapply(
        joint_status, `[[`, character(1), "resolution_mode"
      ),
      route = vapply(gene_results, function(gr) gr$route, character(1)),
      route_code = vapply(
        joint_status, `[[`, character(1), "route_code"
      ),
      strict_route_code = vapply(
        joint_status, `[[`, character(1), "strict_route_code"
      ),
      enabled_route_code = vapply(
        joint_status, `[[`, character(1), "enabled_route_code"
      ),
      start_provenance = vapply(gene_results, function(gr) {
        value <- gr$start_provenance
        if (is.null(value) || !length(value)) NA_character_ else as.character(value[[1L]])
      }, character(1)),
      inference_available = vapply(gene_results, function(gr) {
        isTRUE(gr$inference_available)
      }, logical(1)),
      failure_reason = vapply(gene_results, function(gr) {
        value <- gr$failure_reason
        if (is.null(value) || !length(value)) NA_character_ else as.character(value[[1L]])
      }, character(1)),
      strict_failure_reason = vapply(
        joint_status, `[[`, character(1), "strict_failure_reason"
      ),
      enabled_failure_reason = vapply(
        joint_status, `[[`, character(1), "enabled_failure_reason"
      ),
      mode_method = vapply(gene_results, function(gr) gr$mode_method, character(1)),
      mode_status = vapply(gene_results, function(gr) gr$mode_status, character(1)),
      mode_max_score_residual = vapply(gene_results, function(gr) {
        value <- gr$mode$max_score_residual
        if (is.null(value) || !length(value)) NA_real_ else as.numeric(value[[1L]])
      }, numeric(1)),
      mode_max_newton_displacement = vapply(gene_results, function(gr) {
        value <- gr$mode$max_newton_displacement
        if (is.null(value) || !length(value)) NA_real_ else as.numeric(value[[1L]])
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
    mode_diagnostics <- lapply(gene_results, function(gr) gr$mode)
    names(mode_diagnostics) <- gname[gid]
    result_list$provenance <- list(
      schema_version = 3L,
      estimator = "ln_joint_root",
      candidate_resolution_policy =
        ln_joint_resolution_thresholds()$policy_version,
      candidate_resolution_mode = ctx$ln_candidate_resolution_mode,
      model_covariance = "ln_joint_curvature_inverse",
      robust_covariance = "ln_joint_donor_sandwich",
      random_effect = if (output_re) {
        "pml_mode_given_ln_joint_root"
      } else {
        "none"
      },
      parameterization = list(
        beta_internal = "standardized_design",
        beta_public = "centered_public_design",
        subject = "sigma2_natural",
        cell_internal = "gamma_nb_size",
        cell_public = "Cell_equals_1_over_gamma",
        donor_mode = "w_equals_log_Gamma_effect"
      ),
      coefficient_map_internal_to_public = coef_internal_to_public,
      per_gene = per_gene_provenance,
      mode_diagnostics = mode_diagnostics
    )
  } else {
    legacy_start <- vapply(gene_results, function(gr) {
      rescue <- gr$debug$rescue
      if (!is.null(rescue) && length(rescue) && !is.na(rescue[[1L]]) &&
          !identical(as.character(rescue[[1L]]), "none")) {
        "pql_start"
      } else {
        "direct"
      }
    }, character(1))
    legacy_mode_status <- if (output_re) {
      vapply(gene_results, function(gr) {
        values <- tail(gr$vec, ctx$k)
        if (length(values) == ctx$k && all(is.finite(values))) "ok" else "unavailable"
      }, character(1))
    } else {
      rep("not_requested", length(gene_results))
    }
    result_list$provenance <- list(
      schema_version = 2L,
      estimator = "legacy_post_pml_output_reproduction",
      model_covariance = "legacy_post_pml_curvature",
      robust_covariance = "legacy_ln_donor_sandwich_corrected_posindy",
      random_effect = if (output_re) "legacy_pml_mode" else "none",
      parameterization = list(
        beta_internal = "standardized_design_post_pml",
        beta_public = "centered_public_design",
        subject = "sigma2_natural",
        cell_internal = "gamma_nb_size",
        cell_public = "Cell_equals_1_over_gamma",
        donor_mode = "legacy_pml_logw"
      ),
      coefficient_map_internal_to_public = coef_internal_to_public,
      per_gene = data.frame(
        gene_id = gid,
        gene = gname[gid],
        route = rep("R1", length(gene_results)),
        start_provenance = legacy_start,
        inference_available = rep(FALSE, length(gene_results)),
        failure_reason = rep(NA_character_, length(gene_results)),
        mode_method = if (output_re) rep("legacy_pml_mode", length(gene_results)) else
          rep("none", length(gene_results)),
        mode_status = legacy_mode_status,
        stringsAsFactors = FALSE
      )
    )
  }

  result_list$provenance$route <- "LN"
  result_list$provenance$requested_model <- "LN"
  result_list$provenance$effective_model <- if (isTRUE(allow_per_gene_switch)) {
    "LN_with_explicit_per_gene_switch"
  } else {
    "LN"
  }
  if (identical(ln_estimator, "joint_ln")) {
    result_list$provenance$inference_policy <-
      "joint_ln_same_basin_wald_calibration_candidate"
    result_list$provenance$pace_status <- "joint_ln_calibration_pending"
    result_list$provenance$test_capabilities <- data.frame(
      test = c("wald", "score", "lrt"),
      public_available = c(TRUE, FALSE, FALSE),
      status = c(
        "available_when_inference_status_available",
        "joint_ln_calibration_pending",
        "joint_ln_calibration_pending"
      ),
      stringsAsFactors = FALSE
    )
  } else {
    result_list$provenance$inference_policy <- "legacy_compatibility_output"
    optional_status <- function(requested, computed) {
      if (!isTRUE(requested)) {
        return("available_when_requested")
      }
      if (!isTRUE(computed) && !is.na(optional_test_skip_reason)) {
        return(optional_test_skip_reason)
      }
      "requested_see_test_status"
    }
    result_list$provenance$test_capabilities <- data.frame(
      test = c("wald", "score", "lrt"),
      public_available = rep(TRUE, 3L),
      status = c(
        "legacy_compatibility_output",
        optional_status(requested_score, want_score),
        optional_status(requested_lrt, want_lrt)
      ),
      stringsAsFactors = FALSE
    )
  }

  if (keep_diagnostics) {
    result_list$diagnostics <- list(
      loglik = loglik_vec,
      godambe = godambe_list,
      raw = gene_results,
      lrt = if ((want_score || want_lrt) && length(lrt_names) > 0) lapply(gene_results, function(gr) gr$lrt) else NULL,
      test_status = test_status_df
    )
  } else {
    result_list$diagnostics <- NULL
  }

  if (want_lrt && lrt_details && length(lrt_names) > 0) {
    gene_labels <- summary_df$gene
    detail_rows <- lapply(seq_along(gene_results), function(idx) {
      gr <- gene_results[[idx]]
      if (is.null(gr$lrt)) {
        return(NULL)
      }
      do.call(rbind, lapply(lrt_names, function(nm) {
        res <- gr$lrt[[nm]]
        if (is.null(res)) {
          return(NULL)
        }
        gate_value <- function(name, default) {
          value <- res$raw$gate[[name]]
          if (is.null(value) || length(value) == 0L) default else value[[1L]]
        }
        data.frame(
          gene = gene_labels[idx],
          contrast = lrt_defs[[nm]]$label,
          stat = res$stat,
          p_value = res$pval,
          score_stat = res$score_stat,
          score_p_value = res$score_pval,
          score_ridge = res$score_ridge,
          score_success = res$score_success,
          scale = res$scale,
          wP = res$wP,
          wU = res$wU,
          qP = res$qP,
          wP_minus_qP = res$raw$wP_minus_qP,
          wP_over_qP = res$raw$wP_over_qP,
          df = res$df,
          method = res$method,
          success = res$success,
          adjusted_available = res$raw$adjusted_available,
          adjusted_failure_reason = res$raw$adjusted_failure_reason,
          failure_reason = res$failure_reason,
          gate_converged = res$raw$gate$converged,
          convergence_basis = gate_value("convergence_basis", NA_character_),
          unrestricted_route = gate_value("unrestricted_route", NA_character_),
          unrestricted_failure_reason = gate_value(
            "unrestricted_failure_reason",
            NA_character_
          ),
          constrained_route = gate_value("constrained_route", NA_character_),
          constrained_failure_reason = gate_value(
            "constrained_failure_reason",
            NA_character_
          ),
          legacy_converged = gate_value("legacy_converged", NA),
          legacy_conv_code = gate_value("legacy_conv_code", NA_integer_),
          information_scaled_eligible = gate_value("information_scaled_eligible", NA),
          information_scaled_convergence_ok = gate_value("information_scaled_convergence_ok", NA),
          information_scaled_accepted = gate_value("information_scaled_accepted", NA),
          H_red_positive_definite = gate_value("H_red_positive_definite", NA),
          reduced_diagnostics_finite = gate_value("diagnostics_finite", NA),
          H_red_rcond = gate_value("H_red_rcond", NA_real_),
          reduced_gradient_max_abs = gate_value("reduced_gradient_max_abs", NA_real_),
          max_abs_newton_step_scaled = gate_value("max_abs_newton_step_scaled", NA_real_),
          gradient_quadratic_lr_bound = gate_value("gradient_quadratic_lr_bound", NA_real_),
          newton_step_scaled_ok = gate_value("newton_step_scaled_ok", NA),
          gradient_quadratic_lr_bound_ok = gate_value("gradient_quadratic_lr_bound_ok", NA),
          reduced_parameterization = gate_value("reduced_parameterization", NA_character_),
          gate_constraint_satisfied = gate_value("constraint_satisfied", NA),
          constraint_residual_max_abs = gate_value("constraint_residual_max_abs", NA_real_),
          gate_nuisance_hessian_invertible = res$raw$gate$nuisance_hessian_invertible,
          gate_interior = res$raw$gate$interior,
          gate_nuisance_interior = res$raw$gate$nuisance_interior,
          gate_profiled_curvature_positive = res$raw$gate$profiled_curvature_positive,
          gate_profiled_variability_positive = gate_value("profiled_variability_positive", NA),
          gate_scale_positive = res$raw$gate$scale_positive,
          ridge_H = res$ridge$Hpsi,
          ridge_G = res$ridge$Gpsi,
          ridge_wald = res$ridge$wald,
          fallback = res$fallback,
          stringsAsFactors = FALSE
        )
      }))
    })
    detail_rows <- Filter(Negate(is.null), detail_rows)
    result_list$lrt_details <- if (length(detail_rows) > 0) do.call(rbind, detail_rows) else data.frame(
      gene = character(0),
      contrast = character(0),
      stat = numeric(0),
      p_value = numeric(0),
      score_stat = numeric(0),
      score_p_value = numeric(0),
      score_ridge = numeric(0),
      score_success = logical(0),
      scale = numeric(0),
      wP = numeric(0),
      wU = numeric(0),
      qP = numeric(0),
      wP_minus_qP = numeric(0),
      wP_over_qP = numeric(0),
      df = numeric(0),
      method = character(0),
      success = logical(0),
      adjusted_available = logical(0),
      adjusted_failure_reason = character(0),
      failure_reason = character(0),
      gate_converged = logical(0),
      convergence_basis = character(0),
      unrestricted_route = character(0),
      unrestricted_failure_reason = character(0),
      constrained_route = character(0),
      constrained_failure_reason = character(0),
      legacy_converged = logical(0),
      legacy_conv_code = integer(0),
      information_scaled_eligible = logical(0),
      information_scaled_convergence_ok = logical(0),
      information_scaled_accepted = logical(0),
      H_red_positive_definite = logical(0),
      reduced_diagnostics_finite = logical(0),
      H_red_rcond = numeric(0),
      reduced_gradient_max_abs = numeric(0),
      max_abs_newton_step_scaled = numeric(0),
      gradient_quadratic_lr_bound = numeric(0),
      newton_step_scaled_ok = logical(0),
      gradient_quadratic_lr_bound_ok = logical(0),
      reduced_parameterization = character(0),
      gate_constraint_satisfied = logical(0),
      constraint_residual_max_abs = numeric(0),
      gate_nuisance_hessian_invertible = logical(0),
      gate_interior = logical(0),
      gate_nuisance_interior = logical(0),
      gate_profiled_curvature_positive = logical(0),
      gate_profiled_variability_positive = logical(0),
      gate_scale_positive = logical(0),
      ridge_H = numeric(0),
      ridge_G = numeric(0),
      ridge_wald = numeric(0),
      fallback = logical(0),
      stringsAsFactors = FALSE
    )
  } else {
    result_list$lrt_details <- NULL
  }

  if (identical(ln_operational_fallback, "legacy_post_pml")) {
    fallback_terms <- if (is.null(predn)) {
      paste0("X", seq_len(nb))
    } else {
      predn
    }
    result_list$operational_fallback <- ln_operational_build_sidecar(
      primary_results = gene_results,
      eligibility = fallback_eligibility,
      fallback_results = fallback_results,
      gid = gid,
      gene_names = gname[gid],
      term_names = fallback_terms,
      ctx = ctx
    )
  }

  if (identical(ln_estimator, "joint_ln")) {
    operational_fallback <- if (
        "operational_fallback" %in% names(result_list)) {
      result_list$operational_fallback
    } else {
      NULL
    }
    result_list$availability <- ln_availability_build(
      primary_results = gene_results,
      gid = gid,
      gene_names = gname[gid],
      operational_fallback = operational_fallback
    )
    ln_availability_warn(
      result_list$availability,
      operational_fallback_requested = !is.null(operational_fallback)
    )
  }

  return(result_list)
}
