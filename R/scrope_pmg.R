pmg_find_within_subject_predictors <- function(pred,
                                               fid,
                                               intcol,
                                               predictor_names = NULL,
                                               tol = 1e-10) {
  nb <- ncol(pred)
  non_intercept <- setdiff(seq_len(nb), intcol)
  if (length(non_intercept) == 0L) {
    return(character(0))
  }

  bad_cols <- logical(length(non_intercept))
  for (jj in seq_along(non_intercept)) {
    col_idx <- non_intercept[[jj]]
    for (ss in seq_len(length(fid) - 1L)) {
      idx <- fid[ss]:(fid[ss + 1L] - 1L)
      if (length(idx) <= 1L) {
        next
      }
      vals <- pred[idx, col_idx]
      if ((max(vals) - min(vals)) > tol) {
        bad_cols[[jj]] <- TRUE
        break
      }
    }
  }

  if (!any(bad_cols)) {
    return(character(0))
  }

  labels <- if (is.null(predictor_names)) {
    paste0("X", seq_len(nb))
  } else {
    predictor_names
  }
  labels[non_intercept[bad_cols]]
}

scrope_pmg_internal <- function(
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
    use_betas = "betae",
    output_re = FALSE,
    additional_tests = character(0),
    lrt = FALSE,
    lrt_contrasts = NULL,
    lrt_include_scale = FALSE,
    lrt_details = FALSE,
    keep_diagnostics = TRUE) {

  if (output_re && verbose) {
    warning("output_re=TRUE is not available for model='PGMM'; returning NULL random_effect.")
  }
  if (allow_per_gene_switch && verbose) {
    warning("allow_per_gene_switch is ignored for model='PGMM'.")
  }
  if (!identical(opt, "lbfgs") && verbose) {
    warning("opt is ignored for model='PGMM'; the exact PMG path uses nlminb with bobyqa fallback.")
  }
  if (!identical(use_betas, "betae") && verbose) {
    warning("use_betas is ignored for model='PGMM'.")
  }

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
    output_re = FALSE,
    reml = 0
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

  within_subject_predictors <- pmg_find_within_subject_predictors(
    pred = ctx$pred,
    fid = ctx$fid,
    intcol = intcol,
    predictor_names = predn
  )
  if (length(within_subject_predictors) > 0L && verbose) {
    warning(sprintf(
      paste(
        "model='PGMM' permits cell-level fixed predictors, matching legacy",
        "nebula(model='PMM'), but it only models subject-level overdispersion.",
        "Inference for predictors varying within subject should be used with caution.",
        "Detected within-subject variation in: %s."
      ),
      paste(within_subject_predictors, collapse = ", ")
    ))
  }

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

  sanitize_name <- function(x) {
    nm <- gsub("[^A-Za-z0-9]+", "_", x)
    if (!nzchar(nm)) {
      nm <- "contrast"
    }
    nm
  }

  coef_names <- if (is.null(predn)) paste0("X", seq_len(nb)) else predn
  contrast_defs <- list()
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
      paste0("contrast", length(contrast_defs) + 1L)
    } else {
      as.character(name)
    }
    label <- sanitize_name(label_display)
    base <- label
    idx <- 1L
    while (label %in% names(contrast_defs)) {
      idx <- idx + 1L
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
  if (!is.null(lrt_contrasts)) {
    if (!is.list(lrt_contrasts)) {
      stop("lrt_contrasts must be a list.")
    }
    for (cc in lrt_contrasts) {
      if (is.null(cc)) next
      if (is.matrix(cc)) {
        add_contrast(paste0("contrast", length(contrast_defs) + 1L), cc, rep(0, nrow(cc)))
      } else if (is.numeric(cc)) {
        add_contrast(paste0("contrast", length(contrast_defs) + 1L), matrix(cc, nrow = 1), 0)
      } else if (is.list(cc)) {
        L <- cc$L
        if (is.null(L)) {
          stop("Each contrast list must include L.")
        }
        b_vec <- if (!is.null(cc$b)) cc$b else rep(0, if (is.matrix(L)) nrow(L) else 1)
        add_contrast(if (!is.null(cc$name)) cc$name else paste0("contrast", length(contrast_defs) + 1L), L, b_vec)
      } else {
        stop("Unsupported element in lrt_contrasts.")
      }
    }
  }
  contrast_payload <- lapply(contrast_defs, function(def) {
    list(L = def$L, b = def$b, label = def$label)
  })
  contrast_names <- names(contrast_defs)

  maxcore <- max(c(1, length(parallelly::availableWorkers()) - 1))
  if (ncore > maxcore && verbose) {
    cat("Specified ncore exceeds available cores. Using detected # of cores.\n")
    ncore <- maxcore
  }
  options(future.globals.maxSize = fmaxsize)

  registerDoFuture()
  if (ncore == 1) {
    plan(sequential)
  } else {
    cls <- parallelly::makeClusterPSOCK(ncore)
    plan(cluster, workers = cls)
  }

  gene_results <- foreach(i = gid) %dorng% {
    posv <- call_posindy(count, i - 1, nind)
    fit <- fit_gene_pmg_unconstrained(i, posv, ctx)

    if ((want_score || want_lrt) && length(contrast_payload) > 0L) {
      adj_res <- pmg_adjusted_contrasts(
        ctx = ctx,
        gene_index = i,
        posv = posv,
        contrasts = contrast_payload,
        unconstrained = fit
      )
      fit$lrt <- lapply(adj_res, function(res) {
        list(
          constrained = res$constrained,
          df = res$df,
          raw = res$adj$lr,
          stat = res$adj$lr$wP_adjusted,
          pval = res$p_values$lr,
          success = res$adj$lr$success,
          adjusted_available = res$adj$lr$adjusted_available,
          adjusted_failure_reason = res$adj$lr$adjusted_failure_reason,
          failure_reason = res$adj$lr$failure_reason,
          scale = res$adj$lr$scale,
          wP = res$adj$lr$wP,
          wU = res$adj$lr$wU,
          qP = res$adj$lr$qP,
          method = res$adj$lr$method,
          ridge = res$adj$lr$ridge,
          fallback = res$adj$lr$fallback,
          gate = res$adj$lr$gate,
          score_stat = res$adj$score$stat,
          score_pval = res$p_values$score,
          score_ridge = res$adj$score$ridge,
          score_success = res$adj$score$success,
          score_failure = res$adj$score$failure_reason,
          score_test = res$adj$score,
          wald_stat = res$adj$wald$stat,
          wald_pval = res$p_values$wald,
          wald_success = res$adj$wald$success,
          wald_failure = res$adj$wald$failure_reason,
          wald_test = res$adj$wald
        )
      })
    }

    fit
  }

  if (ncore > 1) {
    parallel::stopCluster(cls)
    plan(sequential)
    gc()
  }

  beta_mat <- t(vapply(gene_results, function(gr) gr$theta$beta, numeric(nb)))
  subvar_vec <- vapply(gene_results, function(gr) gr$theta$subVar, numeric(1))
  naive_var_mat <- t(vapply(gene_results, function(gr) diag(gr$naive_cov), numeric(nb)))
  robust_var_mat <- t(vapply(gene_results, function(gr) diag(gr$robust_cov), numeric(nb)))
  conv_vec <- vapply(gene_results, function(gr) gr$convergence, integer(1))
  optimizer_vec <- vapply(gene_results, function(gr) gr$optimizer, character(1))
  loglik_vec <- vapply(gene_results, function(gr) gr$loglik, numeric(1))
  algorithm_vec <- rep("PGMM", length(gene_results))

  re_all <- data.frame(beta_mat, stringsAsFactors = FALSE)
  beta_names <- paste0("logFC_", coef_names)
  naive_var_names <- paste0("naive_var_", beta_names)
  robust_var_names <- paste0("robust_var_", beta_names)
  colnames(re_all) <- beta_names
  re_all$Subject <- subvar_vec
  for (j in seq_len(nb)) {
    re_all[[naive_var_names[j]]] <- naive_var_mat[, j]
    re_all[[robust_var_names[j]]] <- robust_var_mat[, j]
  }
  re_all$convergence <- conv_vec
  re_all$algorithm <- algorithm_vec
  re_all$optimizer <- optimizer_vec
  re_all$loglik <- loglik_vec
  re_all$gene_id <- gid
  re_all$gene <- gname

  for (j in seq_len(nb)) {
    coef_label <- coef_names[j]

    stat_naive <- (re_all[[beta_names[j]]]^2) / re_all[[naive_var_names[j]]]
    p_naive_j <- pchisq(stat_naive, 1, lower.tail = FALSE)
    re_all[[paste0("p_naive_", coef_label)]] <- p_naive_j

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

  naive_var_cols <- naive_var_names
  robust_var_cols <- robust_var_names
  re_all[naive_var_cols] <- sqrt(re_all[naive_var_cols])
  re_all[robust_var_cols] <- sqrt(re_all[robust_var_cols])

  naive_se_names <- sub("naive_var_logFC_", "se_", naive_var_names)
  naive_se_names <- sub("naive_var_", "se_", naive_se_names)
  robust_se_names <- sub("robust_var_logFC_", "se_robust_", robust_var_names)
  robust_se_names <- sub("robust_var_", "se_robust_", robust_se_names)

  colnames(re_all)[match(naive_var_names, names(re_all))] <- naive_se_names
  colnames(re_all)[match(robust_var_names, names(re_all))] <- robust_se_names

  sds[intcol] <- 1
  scale_beta <- beta_names
  scale_naive <- naive_se_names
  scale_robust <- robust_se_names
  re_all[, c(scale_beta, scale_naive, scale_robust)] <- t(t(re_all[, c(scale_beta, scale_naive, scale_robust)]) / c(sds, sds, sds))

  if ((want_score || want_lrt) && length(contrast_names) > 0L) {
    for (nm in contrast_names) {
      def <- contrast_defs[[nm]]
      key <- def$label
      if (want_lrt) {
        summary_stat <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_real_)
          gr$lrt[[nm]]$stat
        }, numeric(1))
        summary_p <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_real_)
          gr$lrt[[nm]]$pval
        }, numeric(1))
        summary_method <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_character_)
          gr$lrt[[nm]]$method
        }, character(1))
        summary_success <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA)
          gr$lrt[[nm]]$success
        }, logical(1))
        summary_failure <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_character_)
          gr$lrt[[nm]]$failure_reason
        }, character(1))
        summary_dfail <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_character_)
          gr$lrt[[nm]]$adjusted_failure_reason
        }, character(1))

        re_all[[paste0("lrt_adj_", key)]] <- summary_stat
        re_all[[paste0("p_lrt_adj_", key)]] <- summary_p
        if (lrt_include_scale) {
          re_all[[paste0("lrt_scale_", key)]] <- vapply(gene_results, function(gr) {
            if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_real_)
            gr$lrt[[nm]]$scale
          }, numeric(1))
        }
        re_all[[paste0("lrt_method_", key)]] <- summary_method
        re_all[[paste0("lrt_success_", key)]] <- summary_success
        re_all[[paste0("lrt_failure_", key)]] <- normalise_failure(summary_failure, "lrt")
        re_all[[paste0("lrt_adjusted_failure_", key)]] <- normalise_failure(summary_dfail, "lrt")
      }
      if (want_score) {
        re_all[[paste0("score_adj_", key)]] <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_real_)
          gr$lrt[[nm]]$score_stat
        }, numeric(1))
        re_all[[paste0("p_score_adj_", key)]] <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_real_)
          gr$lrt[[nm]]$score_pval
        }, numeric(1))
        re_all[[paste0("score_ridge_", key)]] <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_real_)
          gr$lrt[[nm]]$score_ridge
        }, numeric(1))
        re_all[[paste0("score_success_", key)]] <- vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA)
          gr$lrt[[nm]]$score_success
        }, logical(1))
        re_all[[paste0("score_failure_", key)]] <- normalise_failure(vapply(gene_results, function(gr) {
          if (is.null(gr$lrt) || is.null(gr$lrt[[nm]])) return(NA_character_)
          gr$lrt[[nm]]$score_failure
        }, character(1)), "score")
      }
    }
  }

  collect_status_rows <- function(success_prefix, failure_prefix, test_label, norm_type) {
    success_cols <- grep(paste0("^", success_prefix), names(re_all), value = TRUE)
    if (!length(success_cols)) {
      return(NULL)
    }
    do.call(rbind, lapply(success_cols, function(success_col) {
      lab <- sub(paste0("^", success_prefix), "", success_col)
      failure_col <- paste0(failure_prefix, lab)
      failure_vals <- if (failure_col %in% names(re_all)) {
        normalise_failure(re_all[[failure_col]], norm_type)
      } else {
        rep(NA_character_, nrow(re_all))
      }
      data.frame(
        gene_id = re_all$gene_id,
        gene = re_all$gene,
        contrast = lab,
        test = test_label,
        success = as.logical(re_all[[success_col]]),
        failure_reason = failure_vals,
        stringsAsFactors = FALSE
      )
    }))
  }

  status_blocks <- list(
    collect_status_rows("wald_success_", "wald_failure_", "wald", "wald"),
    collect_status_rows("lrt_success_", "lrt_failure_", "lrt", "lrt"),
    collect_status_rows("score_success_", "score_failure_", "score", "score")
  )
  status_blocks <- Filter(Negate(is.null), status_blocks)
  test_status_df <- if (length(status_blocks) > 0L) {
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

  status_vec <- rep("ok", nrow(re_all))
  detail_vec <- rep(NA_character_, nrow(re_all))
  if (nrow(test_status_df) > 0L) {
    for (ii in seq_len(nrow(re_all))) {
      subset_rows <- test_status_df[test_status_df$gene_id == re_all$gene_id[ii], , drop = FALSE]
      if (nrow(subset_rows) == 0L) {
        next
      }
      fail_rows <- subset_rows[!is.na(subset_rows$success) & !subset_rows$success, , drop = FALSE]
      if (nrow(fail_rows) == 0L) {
        next
      }
      test_types <- unique(fail_rows$test)
      status_vec[ii] <- paste(test_types, collapse = ",")
      detail_entries <- tapply(fail_rows$failure_reason, fail_rows$test, function(vals) {
        clean_vals <- unique(na.omit(vals))
        if (length(clean_vals) == 0L) {
          "unknown"
        } else {
          paste(clean_vals, collapse = ",")
        }
      })
      detail_vec[ii] <- paste(paste0(names(detail_entries), ":", detail_entries), collapse = ";")
    }
  }
  re_all$test_status <- status_vec
  re_all$test_failure_detail <- detail_vec

  if (want_lrt && length(contrast_names) > 0L) {
    method_cols <- grep("^lrt_method_", names(re_all), value = TRUE)
    if (length(method_cols) > 0L) {
      method_vals <- unlist(re_all[method_cols], use.names = FALSE)
      method_vals <- method_vals[!is.na(method_vals)]
      if (length(method_vals) > 0L) {
        count_fallback <- sum(method_vals %in% c("wald", "wald_ridge", "lrt_fail", "lrt_raw"))
        count_ridge <- sum(method_vals %in% c("lrt_ridge", "wald_ridge"))
        if (count_fallback > 0L) {
          warning(sprintf(
            "Adjusted LRT fell back or was gated out for %d gene-contrast(s); see lrt_method_* columns.",
            count_fallback
          ))
        }
        if (count_ridge > 0L) {
          warning(sprintf(
            "Adjusted LRT applied ridge regularisation for %d gene-contrast(s).",
            count_ridge
          ))
        }
      }
    }
  }

  overdisp_df <- re_all[, "Subject", drop = FALSE]
  drop_cols <- c("Subject", "convergence", "algorithm", "optimizer", "loglik")
  summary_df <- re_all[, setdiff(names(re_all), drop_cols), drop = FALSE]

  result_list <- list(
    summary = summary_df,
    overdispersion = overdisp_df,
    convergence = conv_vec,
    algorithm = algorithm_vec,
    random_effect = NULL
  )

  if (keep_diagnostics) {
    result_list$diagnostics <- list(
      loglik = loglik_vec,
      godambe = lapply(gene_results, `[[`, "godambe"),
      raw = gene_results,
      lrt = if ((want_score || want_lrt) && length(contrast_names) > 0L) lapply(gene_results, function(gr) gr$lrt) else NULL,
      test_status = test_status_df
    )
  } else {
    result_list$diagnostics <- NULL
  }

  if (want_lrt && lrt_details && length(contrast_names) > 0L) {
    detail_rows <- lapply(seq_along(gene_results), function(idx) {
      gr <- gene_results[[idx]]
      if (is.null(gr$lrt)) {
        return(NULL)
      }
      do.call(rbind, lapply(contrast_names, function(nm) {
        res <- gr$lrt[[nm]]
        if (is.null(res)) {
          return(NULL)
        }
        data.frame(
          gene = gname[idx],
          contrast = contrast_defs[[nm]]$display,
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
          success = res$success,
          adjusted_available = res$adjusted_available,
          adjusted_failure_reason = res$adjusted_failure_reason,
          failure_reason = res$failure_reason,
          gate_converged = res$gate$converged,
          gate_nuisance_hessian_invertible = res$gate$nuisance_hessian_invertible,
          gate_interior = res$gate$interior,
          gate_nuisance_interior = res$gate$nuisance_interior,
          gate_profiled_curvature_positive = res$gate$profiled_curvature_positive,
          gate_scale_positive = res$gate$scale_positive,
          ridge_H = res$ridge$Hpsi,
          ridge_G = res$ridge$Gpsi,
          ridge_wald = res$ridge$wald,
          fallback = res$fallback,
          stringsAsFactors = FALSE
        )
      }))
    })
    detail_rows <- Filter(Negate(is.null), detail_rows)
    result_list$lrt_details <- if (length(detail_rows) > 0L) {
      do.call(rbind, detail_rows)
    } else {
      data.frame(
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
        df = integer(0),
        method = character(0),
        success = logical(0),
        adjusted_available = logical(0),
        adjusted_failure_reason = character(0),
        failure_reason = character(0),
        gate_converged = logical(0),
        gate_nuisance_hessian_invertible = logical(0),
        gate_interior = logical(0),
        gate_nuisance_interior = logical(0),
        gate_profiled_curvature_positive = logical(0),
        gate_scale_positive = logical(0),
        ridge_H = numeric(0),
        ridge_G = numeric(0),
        ridge_wald = numeric(0),
        fallback = logical(0),
        stringsAsFactors = FALSE
      )
    }
  }

  result_list
}
