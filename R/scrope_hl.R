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

  if (!use_betas %in% c("betae", "final_betas")) {
    stop("use_betas must be 'betae' or 'final_betas'.")
  }
  if (use_betas == "final_betas") {
    warning(
      "use_betas is currently ignored in scrope_hl(); retained for API compatibility.",
      call. = FALSE
    )
  }

  reml <- 0

  prep <- hl_prepare_data(
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
    output_re = output_re,
    reml = reml
  )
  ctx <- prep$ctx
  count <- prep$count
  gid <- prep$gid
  gname <- prep$gene_names
  predn <- prep$predictor_names
  sds <- prep$design_sds
  nb <- ctx$nb
  nind <- ctx$nind
  k <- ctx$k
  intcol <- ctx$intcol

  coef_names <- if (is.null(predn)) paste0("X", seq_len(nb)) else predn

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

  if ((lrt_details || !is.null(lrt_contrasts)) && !want_lrt) {
    extra_tests <- union(extra_tests, "lrt")
    want_lrt <- TRUE
  }

  has_user_contrasts <- !is.null(lrt_contrasts) && length(lrt_contrasts) > 0
  default_wald_only <- !want_score && !want_lrt && !has_user_contrasts

  sanitize_name <- function(x) {
    nm <- gsub("[^A-Za-z0-9]+", "_", x)
    if (!nzchar(nm)) {
      nm <- "contrast"
    }
    nm
  }

  contrast_defs <- list()
  if (!default_wald_only || has_user_contrasts) {
    add_contrast <- function(name, L, b) {
      L <- as.matrix(L)
      if (ncol(L) != nb) {
        stop("Each contrast L must have ncol equal to number of predictors.")
      }
      if (is.null(dim(b))) {
        b <- rep(b, length.out = nrow(L))
      }
      if (length(b) != nrow(L)) {
        stop("Length of b must match the number of rows of L.")
      }
      label_display <- if (is.null(name) || !nzchar(name)) {
        paste0("contrast", length(contrast_defs) + 1)
      } else {
        as.character(name)
      }
      label <- sanitize_name(label_display)
      base <- label
      idx <- 1
      while (label %in% names(contrast_defs)) {
        idx <- idx + 1
        label <- paste0(base, "_", idx)
      }
      contrast_defs[[label]] <<- list(
        label = label,
        display = label_display,
        L = L,
        b = b,
        df = nrow(L)
      )
    }

    non_intercept <- setdiff(seq_len(nb), intcol)
    for (j in non_intercept) {
      L <- matrix(0, nrow = 1, ncol = nb)
      L[1, j] <- 1
      add_contrast(coef_names[j], L, 0)
    }

    if (has_user_contrasts) {
      if (!is.list(lrt_contrasts)) {
        stop("lrt_contrasts must be a list when provided.")
      }
      for (cc in lrt_contrasts) {
        if (is.null(cc)) next
        if (is.matrix(cc)) {
          add_contrast(paste0("contrast", length(contrast_defs) + 1), cc, rep(0, nrow(cc)))
        } else if (is.numeric(cc)) {
          add_contrast(paste0("contrast", length(contrast_defs) + 1), matrix(cc, nrow = 1), 0)
        } else if (is.list(cc)) {
          L <- cc$L
          if (is.null(L)) {
            stop("Each contrast list must include L.")
          }
          b_vec <- if (!is.null(cc$b)) cc$b else rep(0, if (is.matrix(L)) nrow(L) else 1)
          add_contrast(if (!is.null(cc$name)) cc$name else paste0("contrast", length(contrast_defs) + 1), L, b_vec)
        } else {
          stop("Unsupported element in lrt_contrasts.")
        }
      }
    }
  }

  need_adj <- (!default_wald_only && length(contrast_defs) > 0) || has_user_contrasts
  contrast_payload <- if (need_adj) {
    lapply(contrast_defs, function(def) list(L = def$L, b = def$b, label = def$label))
  } else {
    list()
  }

  if (ncore > 1 && verbose) {
    warning("scrope_hl currently uses sequential execution; ignoring ncore > 1.")
  }

  gene_results <- lapply(gid, function(i) {
    posv <- call_posindy(count, i - 1, nind)
    ctx$gene_index <- i
    fit <- fit_gene_hl(i, posv, ctx)
    if (need_adj) {
      adj_res <- hl_adjusted_contrasts(
        unconstrained = fit$hl_unconstrained,
        ctx = ctx,
        gene_index = i,
        count = count,
        contrasts = contrast_payload,
        posv = posv,
        ridge_factor = adj_ridge_factor,
        cond_limit = adj_cond_limit
      )
      fit$hl_adj <- adj_res
    }
    fit$hl_unconstrained <- NULL
    fit
  })

  re_mat <- do.call(cbind, lapply(gene_results, `[[`, "vec"))
  re_mat <- t(re_mat)
  godambe_list <- lapply(gene_results, `[[`, "godambe")
  hl_list <- lapply(gene_results, `[[`, "hl")

  re_all <- as.data.frame(re_mat)
  beta_names <- paste0(
    "logFC_",
    if (is.null(predn)) paste0("X", seq_len(nb)) else predn
  )
  var_names <- paste0("var_", beta_names)

  col_offset <- 0
  colnames(re_all)[1:nb] <- beta_names
  col_offset <- nb
  colnames(re_all)[col_offset + 1] <- "Subject"
  colnames(re_all)[col_offset + 2] <- "Cell"
  col_offset <- col_offset + 2

  iv_var <- seq(col_offset + 1, col_offset + nb)
  for (j in seq_len(nb)) {
    colnames(re_all)[iv_var[j]] <- var_names[j]
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
    stat_naive <- (re_all[[beta_names[j]]]^2) / re_all[[var_names[j]]]
    p_naive_j <- pchisq(stat_naive, 1, lower.tail = FALSE)
    re_all[[paste0(
      "p_naive_",
      if (is.null(predn)) paste0("X", j) else predn[j]
    )]] <- p_naive_j
  }

  re_all[iv_var] <- sqrt(re_all[iv_var])

  se_names <- sub("var_logFC_", "se_", var_names)
  se_names <- sub("var_", "se_", se_names)

  colnames(re_all)[iv_var] <- se_names

  sds[intcol] <- 1
  scale_cols <- c(seq_len(nb), iv_var)
  scvec <- c(sds, sds)
  re_all[, scale_cols] <- t(t(re_all[, scale_cols]) / scvec)

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

  compute_beta_cov <- function(g) {
    if (is.null(g)) {
      return(list(success = FALSE, cov = NULL, reason = "missing_godambe"))
    }
    H <- g$H
    U <- g$per_subject_gradients
    if (is.null(H) || is.null(U)) {
      return(list(success = FALSE, cov = NULL, reason = "missing_blocks"))
    }
    S <- if (!is.null(g$S)) g$S else U %*% t(U)
    s_info <- invert_with_ridge(S,
      cond_limit = adj_cond_limit,
      base_factor = adj_ridge_factor
    )
    if (!s_info$success || is.null(s_info$inverse) || s_info$ridge > 0) {
      return(list(success = FALSE, cov = NULL, reason = "singular_S"))
    }
    S_inv <- s_info$inverse
    G_full <- H %*% S_inv %*% H
    idx_beta <- seq_len(nb)
    idx_lambda <- nb + seq_len(nrow(G_full) - nb)
    G_bb <- G_full[idx_beta, idx_beta, drop = FALSE]
    if (length(idx_lambda) > 0) {
      G_bl <- G_full[idx_beta, idx_lambda, drop = FALSE]
      G_ll <- G_full[idx_lambda, idx_lambda, drop = FALSE]
      G_ll_info <- invert_with_ridge(G_ll,
        cond_limit = adj_cond_limit,
        base_factor = adj_ridge_factor
      )
      if (!G_ll_info$success || is.null(G_ll_info$inverse) || G_ll_info$ridge > 0) {
        return(list(success = FALSE, cov = NULL, reason = "singular_G_ll"))
      }
      G_beta <- G_bb - G_bl %*% G_ll_info$inverse %*% t(G_bl)
    } else {
      G_beta <- G_bb
    }
    G_beta <- 0.5 * (G_beta + t(G_beta))
    G_beta_info <- invert_with_ridge(G_beta,
      cond_limit = adj_cond_limit,
      base_factor = adj_ridge_factor
    )
    if (!G_beta_info$success || is.null(G_beta_info$inverse) || G_beta_info$ridge > 0) {
      return(list(success = FALSE, cov = NULL, reason = "singular_G_beta"))
    }
    list(success = TRUE, cov = G_beta_info$inverse, reason = NULL)
  }

  if (length(godambe_list) > 0) {
    coef_names <- if (is.null(predn)) paste0("X", seq_len(nb)) else predn
    beta_mat <- as.matrix(summary_df[, beta_names, drop = FALSE])
    robust_var_mat <- matrix(NA_real_, nrow = nrow(beta_mat), ncol = ncol(beta_mat))
    failure_tracker <- list()
    row_failure_reason <- rep(NA_character_, length(godambe_list))
    scale_vec <- sds
    scale_vec[intcol] <- 1

    for (idx in seq_along(godambe_list)) {
      g <- godambe_list[[idx]]
      if (is.null(g) || is.null(g$H) || is.null(g$per_subject_gradients)) {
        row_failure_reason[idx] <- normalise_failure("missing_godambe", "wald")
        failure_tracker[[length(failure_tracker) + 1]] <- normalise_failure("missing_godambe", "wald")
        next
      }
      cov_info <- compute_beta_cov(g)
      if (!cov_info$success || !is.matrix(cov_info$cov)) {
        row_failure_reason[idx] <- normalise_failure(cov_info$reason, "wald")
        failure_tracker[[length(failure_tracker) + 1]] <- normalise_failure(cov_info$reason, "wald")
        next
      }
      var_diag <- diag(cov_info$cov)
      var_scaled <- var_diag / (scale_vec ^ 2)
      robust_var_mat[idx, ] <- var_scaled
    }

    robust_se_mat <- sqrt(robust_var_mat)
    robust_se_mat[!is.finite(robust_se_mat)] <- NA_real_
    for (j in seq_len(nb)) {
      se_name <- paste0("se_robust_", coef_names[j])
      summary_df[[se_name]] <- robust_se_mat[, j]
    }

    for (j in seq_len(nb)) {
      lab <- coef_names[j]
      beta_vals <- beta_mat[, j]
      var_vals <- robust_var_mat[, j]
      invalid_var <- !is.finite(var_vals) | var_vals <= 0
      wald_stat <- (beta_vals ^ 2) / var_vals
      wald_stat[invalid_var] <- NA_real_
      p_vals <- pchisq(wald_stat, df = 1, lower.tail = FALSE)
      success_vec <- rep(NA, length(wald_stat))
      success_vec[!invalid_var & is.finite(wald_stat)] <- TRUE
      success_vec[invalid_var] <- FALSE
      failure_vec <- rep(NA_character_, length(wald_stat))
      fail_rows <- which(!is.na(row_failure_reason))
      if (length(fail_rows) > 0) {
      failure_vec[fail_rows] <- normalise_failure(row_failure_reason[fail_rows], "wald")
      }
      failure_vec[is.na(failure_vec) & invalid_var] <- normalise_failure("singular_godambe", "wald")
      summary_df[[paste0("wald_stat_", lab)]] <- wald_stat
      summary_df[[paste0("p_robust_", lab)]] <- p_vals
      summary_df[[paste0("wald_success_", lab)]] <- success_vec
      summary_df[[paste0("wald_failure_", lab)]] <- failure_vec
    }

    failure_reasons <- unlist(failure_tracker, use.names = FALSE)
    failure_reasons <- failure_reasons[!is.na(failure_reasons)]
    if (length(failure_reasons) > 0) {
      reason_tab <- table(failure_reasons)
      warn_msg <- paste0(
        "HL Wald: unable to compute robust variance for ",
        sum(reason_tab),
        " gene(s) [",
        paste(names(reason_tab), reason_tab, sep = ": ", collapse = ", "),
        "]; see wald_failure_* columns."
      )
      warning(warn_msg, call. = FALSE)
    }
  }

  test_status_df <- NULL
  hl_adj_list <- if (need_adj) lapply(gene_results, `[[`, "hl_adj") else NULL

  if (need_adj) {
    safe_num <- function(expr) {
      val <- tryCatch(expr, error = function(...) NA_real_)
      if (length(val) == 0L) NA_real_ else as.numeric(val)[1]
    }
    safe_chr <- function(expr) {
      val <- tryCatch(expr, error = function(...) NA_character_)
      if (length(val) == 0L) NA_character_ else as.character(val)[1]
    }
    safe_lgl <- function(expr) {
      val <- tryCatch(expr, error = function(...) NA)
      if (length(val) == 0L) NA else as.logical(val)[1]
    }

    for (def_idx in seq_along(contrast_defs)) {
      def <- contrast_defs[[def_idx]]
      lab <- def$label
      extract_numeric <- function(fun) {
        vapply(hl_adj_list, function(entry) {
          if (is.null(entry)) {
            return(NA_real_)
          }
          val <- entry[[lab]]
          if (is.null(val)) {
            NA_real_
          } else {
            res <- tryCatch(fun(val), error = function(...) NA_real_)
            if (length(res) == 0L) {
              return(NA_real_)
            }
            as.numeric(res)[1]
          }
        }, numeric(1))
      }
      extract_character <- function(fun) {
        vapply(hl_adj_list, function(entry) {
          if (is.null(entry)) {
            return(NA_character_)
          }
          val <- entry[[lab]]
          if (is.null(val)) {
            NA_character_
          } else {
            res <- tryCatch(fun(val), error = function(...) NA_character_)
            if (length(res) == 0L) {
              return(NA_character_)
            }
            as.character(res)[1]
          }
        }, character(1))
      }
      extract_logical <- function(fun) {
        vapply(hl_adj_list, function(entry) {
          if (is.null(entry)) {
            return(NA)
          }
          val <- entry[[lab]]
          if (is.null(val)) {
            NA
          } else {
            res <- tryCatch(fun(val), error = function(...) NA)
            if (length(res) == 0L) {
              return(NA)
            }
            as.logical(res)[1]
          }
        }, logical(1))
      }

      wald_stat <- extract_numeric(function(val) val$adj$wald$stat)
      wald_p <- extract_numeric(function(val) val$p_values$wald)
      wald_success <- extract_logical(function(val) val$adj$wald$success)
      wald_failure <- extract_character(function(val) val$adj$wald$failure_reason)
      target_wald <- paste0("wald_stat_", lab)
      target_p <- paste0("p_robust_", lab)
      target_success <- paste0("wald_success_", lab)
      target_failure <- paste0("wald_failure_", lab)
      existing_wald <- summary_df[[target_wald]]
      if (is.null(existing_wald) || all(is.na(existing_wald))) {
        summary_df[[target_wald]] <- wald_stat
        summary_df[[target_p]] <- wald_p
        summary_df[[target_success]] <- wald_success
        summary_df[[target_failure]] <- wald_failure
      } else if (any(!is.na(wald_stat) & is.na(existing_wald))) {
        idx <- which(is.na(existing_wald))
        summary_df[[target_wald]][idx] <- wald_stat[idx]
        summary_df[[target_p]][idx] <- wald_p[idx]
        summary_df[[target_success]][idx] <- wald_success[idx]
        summary_df[[target_failure]][idx] <- wald_failure[idx]
      }

      if (want_lrt) {
        lrt_stat <- extract_numeric(function(val) val$adj$lr$adj)
        lrt_p <- extract_numeric(function(val) val$p_values$lr)
        lrt_method <- extract_character(function(val) val$adj$lr$method)
        lrt_success <- extract_logical(function(val) val$adj$lr$success)
        lrt_failure <- extract_character(function(val) val$adj$lr$failure_reason)
        summary_df[[paste0("hl_lrt_adj_", lab)]] <- lrt_stat
        summary_df[[paste0("p_hl_lrt_", lab)]] <- lrt_p
        summary_df[[paste0("hl_lrt_method_", lab)]] <- lrt_method
        summary_df[[paste0("hl_lrt_success_", lab)]] <- lrt_success
        summary_df[[paste0("hl_lrt_failure_", lab)]] <- lrt_failure
        if (lrt_include_scale) {
          lrt_scale <- extract_numeric(function(val) val$adj$lr$scale)
          summary_df[[paste0("hl_lrt_scale_", lab)]] <- lrt_scale
        }
      }

      if (want_score) {
        score_stat <- extract_numeric(function(val) val$adj$score$stat)
        score_p <- extract_numeric(function(val) val$p_values$score)
        score_success <- extract_logical(function(val) val$adj$score$success)
        score_failure <- extract_character(function(val) val$adj$score$failure_reason)
        summary_df[[paste0("hl_score_adj_", lab)]] <- score_stat
        summary_df[[paste0("p_hl_score_", lab)]] <- score_p
        summary_df[[paste0("hl_score_success_", lab)]] <- score_success
        summary_df[[paste0("hl_score_failure_", lab)]] <- score_failure
      }
    }
  }

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

  status_blocks <- list(
    collect_status_rows("wald_success_", "wald_failure_", "wald", "wald"),
    collect_status_rows("hl_lrt_success_", "hl_lrt_failure_", "lrt", "lrt"),
    collect_status_rows("hl_score_success_", "hl_score_failure_", "score", "score")
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
      types <- unique(fail_rows$test)
      status_vec[i] <- paste(types, collapse = ",")
      detail_entries <- tapply(fail_rows$failure_reason, fail_rows$test, function(x) {
        vals <- unique(na.omit(x))
        if (length(vals) == 0) "unknown" else paste(vals, collapse = ",")
      })
      detail_vec[i] <- paste(paste0(names(detail_entries), ":", detail_entries), collapse = ";")
    }
  }
  summary_df$test_status <- status_vec
  summary_df$test_failure_detail <- detail_vec

  result_list <- list(
    summary = summary_df,
    overdispersion = overdisp_df,
    convergence = conv_vec,
    algorithm = alg_vec
  )

  if (output_re) {
    result_list$random_effect <- re_all[, re_names, drop = FALSE]
  } else {
    result_list$random_effect <- NULL
  }

  hl_adj_details <- NULL
  if (need_adj && lrt_details) {
    detail_rows <- list()
    row_id <- 1L
    for (idx in seq_along(gid)) {
      adj_map <- hl_adj_list[[idx]]
      if (is.null(adj_map)) next
      gene_id <- gid[idx]
      gene_label <- if (!is.null(gname) && length(gname) >= gene_id) gname[gene_id] else paste0("gene_", gene_id)
      for (def in contrast_defs) {
        lab <- def$label
        val <- adj_map[[lab]]
        if (is.null(val)) next
        detail_rows[[row_id]] <- data.frame(
          gene_id = gene_id,
          gene = gene_label,
          contrast = def$display,
          label = lab,
          df = def$df,
          wald_stat = safe_num(val$adj$wald$stat),
          p_wald = safe_num(val$p_values$wald),
          lr_stat = safe_num(val$adj$lr$adj),
          p_lr = safe_num(val$p_values$lr),
          lr_method = safe_chr(val$adj$lr$method),
          lr_fallback = safe_lgl(val$adj$lr$fallback),
          lr_scale = safe_num(val$adj$lr$scale),
          lr_success = safe_lgl(val$adj$lr$success),
          lr_failure = safe_chr(val$adj$lr$failure_reason),
          score_stat = safe_num(val$adj$score$stat),
          p_score = safe_num(val$p_values$score),
          score_success = safe_lgl(val$adj$score$success),
          score_failure = safe_chr(val$adj$score$failure_reason),
          ridge_wald_H = safe_num(val$adj$wald$ridge$H),
          ridge_wald_G = safe_num(val$adj$wald$ridge$G),
          ridge_score = safe_num(val$adj$score$ridge),
          ridge_lambda = safe_num(val$adj$nuisances$H_lambda_ridge),
          wald_success = safe_lgl(val$adj$wald$success),
          wald_failure = safe_chr(val$adj$wald$failure_reason),
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }
    }
    if (length(detail_rows) > 0) {
      hl_adj_details <- do.call(rbind, detail_rows)
    } else {
      hl_adj_details <- data.frame()
    }
  }

  result_list$hl_adj_details <- hl_adj_details

  if (need_adj) {
    if (want_lrt) {
      lr_methods <- unlist(lapply(hl_adj_list, function(entry) {
        if (is.null(entry)) return(character(0))
        vapply(entry, function(x) {
          val <- tryCatch(x$adj$lr$method, error = function(...) NA_character_)
          if (length(val) == 0L) NA_character_ else as.character(val)[1]
        }, character(1))
      }), use.names = FALSE)
      lr_success_vec <- unlist(lapply(hl_adj_list, function(entry) {
        if (is.null(entry)) return(logical(0))
        vapply(entry, function(x) {
          val <- tryCatch(x$adj$lr$success, error = function(...) NA)
          if (length(val) == 0L) NA else isTRUE(val)
        }, logical(1))
      }), use.names = FALSE)
      if (length(lr_methods) > 0) {
        fallback_count <- sum(lr_methods %in% c("wald", "wald_ridge", "lrt_fail", "lrt_raw"), na.rm = TRUE)
        ridge_count <- sum(lr_methods %in% c("lrt_ridge", "wald_ridge"), na.rm = TRUE)
        failure_count <- if (length(lr_success_vec) > 0) {
          sum(is.na(lr_success_vec) | !lr_success_vec)
        } else {
          0
        }
        if (fallback_count > 0) {
          warning(sprintf(
            "HL adjusted LRT fell back to Wald for %d contrast(s); see hl_adj_details$lr_method.",
            fallback_count
          ))
        }
        if (ridge_count > 0) {
          warning(sprintf(
            "HL adjusted LRT applied ridge regularisation for %d contrast(s).",
            ridge_count
          ))
        }
        if (failure_count > 0) {
          warning(sprintf(
            "HL adjusted LRT returned non-finite results for %d contrast(s); see hl_adj_details$lr_failure.",
            failure_count
          ))
        }
      }
    }
    if (want_score) {
      score_success <- unlist(lapply(hl_adj_list, function(entry) {
        if (is.null(entry)) return(logical(0))
        vapply(entry, function(x) {
          val <- tryCatch(x$adj$score$success, error = function(...) NA)
          if (length(val) == 0L) NA else isTRUE(val)
        }, logical(1))
      }), use.names = FALSE)
      if (length(score_success) > 0 && any(is.na(score_success) | !score_success)) {
        warning("Adjusted HL score test failed or returned non-finite results for some contrasts; check hl_adj_details$score_failure.")
      }
    }
  }

  if (keep_diagnostics) {
    diag_list <- list(
      loglik = loglik_vec,
      godambe = godambe_list,
      hl = hl_list,
      test_status = test_status_df
    )
    if (need_adj) {
      diag_list$hl_adj <- hl_adj_list
    }
    result_list$diagnostics <- diag_list
  } else {
    result_list$diagnostics <- NULL
  }

  result_list
}
