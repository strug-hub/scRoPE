#' sc-RoPE main entry point
#'
#' Fits the sc-RoPE negative binomial mixed model for multi-subject
#' single-cell data and, optionally, computes invariantly adjusted profile
#' likelihood ratio tests for fixed-effect contrasts.
#'
#' @param count Gene-by-cell count matrix (dense matrix or `dgCMatrix`).
#' @param id Vector of subject identifiers aligned with the columns of `count`.
#' @param pred Cell-level design matrix (defaults to intercept only).
#' @param offset Optional positive vector of exposure offsets.
#' @param min Lower bounds for the overdispersion parameters `(sigma^2, phi)`.
#' @param max Upper bounds for the overdispersion parameters `(sigma^2, phi)`.
#' @param opt Optimiser for the LN fit (`"lbfgs"` or `"trust"`).
#' @param verbose Logical; print progress messages.
#' @param cpc Minimum counts-per-cell threshold for gene filtering.
#' @param mincp Minimum number of non-zero cells required to keep a gene.
#' @param ncore Number of cores for parallel processing.
#' @param fmaxsize Maximum future.globals size (in bytes) when using parallel execution.
#' @param cutoff_cell Threshold (cells per subject × phi) below which HL refinement is triggered.
#' @param kappa Kappa threshold for adaptive HL refinement.
#' @param allow_per_gene_switch Logical; allow per-gene LN/HL switching. Defaults to
#'   `FALSE` so the LN pipeline (required for adjusted LRTs) is used for every gene.
#' @param use_betas Which fixed-effect anchor to use in LN variance-refinement
#'   starts: `"betae"` matches legacy `nebula` intercept-only anchoring, while
#'   `"final_betas"` uses the LN optimizer coefficients.
#' @param output_re Logical; return subject-level random effects.
#' @param additional_tests Optional character vector selecting extra
#'   contrast-level diagnostics to compute, in addition to the default
#'   adjusted Wald columns (`se_robust_*`, `p_robust_*`). Supported values
#'   are `"score"` and `"lrt"`. Defaults to character(0) (Wald only).
#' @param lrt Logical; legacy flag for adjusted profile LRTs. Equivalent to
#'   including `"lrt"` in `additional_tests`.
#' @param lrt_contrasts Optional list of additional contrasts (see details).
#' @param lrt_include_scale Logical; include the H/G scale factor in the summary.
#' @param lrt_details Logical; return a tidy data frame with raw LRT diagnostics.
#' @param keep_diagnostics Logical; retain per-gene diagnostic objects.
#' @param model Fitting path. `"LN"` uses the existing LN pipeline with optional
#'   HL refinement. `"PGMM"` uses the exact Poisson-gamma mixed model path and
#'   currently supports only subject-level predictors that are constant within
#'   subject.
#'
#' @return A list with elements `summary`, `overdispersion`, `convergence`,
#'   `algorithm`, optional `random_effect`, `diagnostics`, and optional
#'   `lrt_details`, mirroring the original `nebula_sand` output with the new
#'   adjusted LRT fields.
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
    model = c("LN", "PGMM")) {
  model <- match.arg(model)
  if (identical(model, "PGMM")) {
    return(scrope_pmg_internal(
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
      keep_diagnostics = keep_diagnostics
    ))
  }

  # fail-safe for use_betas
  if (!use_betas %in% c("betae", "final_betas")) {
    stop("use_betas must be 'betae' or 'final_betas'.")
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

  if (need_constrained && opt != "lbfgs") {
    if (verbose) {
      warning("Adjusted LRT/score require opt='lbfgs'; skipping those tests.")
    }
    extra_tests <- setdiff(extra_tests, c("score", "lrt"))
    want_score <- FALSE
    want_lrt <- FALSE
    need_constrained <- FALSE
  }
  if (need_constrained && allow_per_gene_switch) {
    if (verbose) {
      warning("Adjusted LRT/score not available when allow_per_gene_switch=TRUE; skipping those tests.")
    }
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
    cls <- parallelly::makeClusterPSOCK(ncore)
    plan(cluster, workers = cls)
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
    id = id
  )

  gene_results <- foreach(i = gid) %dorng% {
    posv <- call_posindy(count, i - 1, nind)
    unconstrained <- fit_gene_unconstrained(i, posv, ctx)

    if ((want_score || want_lrt) && length(lrt_names) > 0) {
      lrt_neg_tol <- 1e-6

      build_lrt_cache <- function(current_unconstrained) {
        cache <- vector("list", length(lrt_names))
        names(cache) <- lrt_names
        neg_detected <- FALSE
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

        for (idx_nm in seq_along(lrt_names)) {
          nm <- lrt_names[idx_nm]
          def <- lrt_defs[[nm]]
          constrained_fit <- fit_gene_constrained(i, posv, ctx, def$L, def$b)
          lrt_res <- NULL
          if (want_lrt) {
            lrt_res <- compute_profile_lrt(current_unconstrained, constrained_fit, def$L, def$b)
            if (!is.na(lrt_res$wP) && lrt_res$wP < -lrt_neg_tol) {
              neg_detected <- TRUE
            }
          }
          gene_label <- if (!is.null(gname) && length(gname) >= i) gname[i] else paste0("gene_", i)
          score_res <- empty_score_result(def$L)
          if (want_score) {
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

      if (want_lrt && lrt_build$neg) {
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
  }

  if (ncore > 1) {
    parallel::stopCluster(cls)
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

  for (j in seq_len(nb)) {
    coef_label <- if (is.null(predn)) paste0("X", j) else predn[j]

    ## naive Wald
    stat_naive <- (re_all[[beta_names[j]]]^2) / re_all[[naive_var_names[j]]]
    p_naive_j <- pchisq(stat_naive, 1, lower.tail = FALSE)
    re_all[[paste0("p_naive_", coef_label)]] <- p_naive_j

    ## robust Wald
    var_robust <- re_all[[robust_var_names[j]]]
    stat_robust <- (re_all[[beta_names[j]]]^2) / var_robust
    valid_var <- is.finite(var_robust) & var_robust > 0
    stat_robust[!valid_var] <- NA_real_
    p_robust_j <- pchisq(stat_robust, 1, lower.tail = FALSE)

    success_vec <- rep(NA, length(stat_robust))
    success_vec[valid_var & is.finite(stat_robust)] <- TRUE
    success_vec[!(valid_var & is.finite(stat_robust))] <- FALSE

    failure_vec <- rep(NA_character_, length(stat_robust))
    failure_vec[!valid_var] <- normalise_failure("singular_cov", "wald")
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
          df = res$df,
          method = res$method,
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
      df = numeric(0),
      method = character(0),
      ridge_H = numeric(0),
      ridge_G = numeric(0),
      ridge_wald = numeric(0),
      fallback = logical(0),
      stringsAsFactors = FALSE
    )
  } else {
    result_list$lrt_details <- NULL
  }

  return(result_list)
}
