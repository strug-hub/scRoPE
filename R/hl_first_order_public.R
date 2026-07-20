# Public scrope_hl() adapter for the estimator-aligned first-order engine.

hl_fo_public_contrasts <- function(p, intcol, predictor_names,
                                   user_contrasts = NULL) {
  sanitize <- function(x) {
    out <- gsub("[^A-Za-z0-9]+", "_", x)
    if (!nzchar(out)) "contrast" else out
  }
  out <- list()
  add <- function(name, R, rhs) {
    R <- if (is.null(dim(R))) matrix(R, nrow = 1L) else as.matrix(R)
    if (ncol(R) != p) stop("Each contrast must have one column per predictor.")
    rhs <- rep(as.numeric(rhs), length.out = nrow(R))
    display <- if (is.null(name) || !nzchar(name)) paste0("contrast", length(out) + 1L) else
      as.character(name)
    label <- sanitize(display)
    base <- label
    suffix <- 1L
    while (label %in% names(out)) {
      suffix <- suffix + 1L
      label <- paste0(base, "_", suffix)
    }
    out[[label]] <<- list(label = label, display = display, R = R, rhs = rhs)
  }
  for (jj in setdiff(seq_len(p), intcol)) {
    rr <- matrix(0, 1L, p)
    rr[1L, jj] <- 1
    add(predictor_names[jj], rr, 0)
  }
  if (!is.null(user_contrasts)) {
    if (!is.list(user_contrasts)) stop("lrt_contrasts must be a list.")
    for (item in user_contrasts) {
      if (is.null(item)) next
      if (is.matrix(item)) {
        add(NULL, item, numeric(nrow(item)))
      } else if (is.numeric(item)) {
        add(NULL, item, 0)
      } else if (is.list(item) && !is.null(item$L)) {
        nr <- if (is.null(dim(item$L))) 1L else nrow(item$L)
        add(item$name, item$L, if (is.null(item$b)) numeric(nr) else item$b)
      } else {
        stop("Unsupported lrt_contrasts entry.")
      }
    }
  }
  out
}

hl_fo_public_failure_row <- function(p, predictor_names, gene_id, gene,
                                     status, reason) {
  row <- list()
  for (nm in predictor_names) {
    row[[paste0("logFC_", nm)]] <- NA_real_
    row[[paste0("se_", nm)]] <- NA_real_
    row[[paste0("p_naive_", nm)]] <- NA_real_
    row[[paste0("se_robust_", nm)]] <- NA_real_
    row[[paste0("wald_stat_", nm)]] <- NA_real_
    row[[paste0("p_robust_", nm)]] <- NA_real_
    row[[paste0("wald_df_", nm)]] <- NA_integer_
    row[[paste0("wald_m_", nm)]] <- NA_integer_
    row[[paste0("wald_contrast_rank_", nm)]] <- NA_integer_
    row[[paste0("wald_covariance_rank_", nm)]] <- NA_integer_
    row[[paste0("wald_success_", nm)]] <- FALSE
    row[[paste0("wald_failure_", nm)]] <- reason
  }
  row$gene_id <- gene_id
  row$gene <- gene
  row$fit_success <- FALSE
  row$fit_status <- status
  row$fit_failure <- reason
  row$covariance_success <- FALSE
  row$covariance_status <- NA_character_
  as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
}

hl_fo_public_scrope_hl <- function(
    count, id, pred = NULL, offset = NULL,
    min = c(1e-4, 1e-4), max = c(10, 1000), opt = "lbfgs",
    verbose = TRUE, cpc = 0.005, mincp = 5, ncore = 1,
    fmaxsize = Inf, use_betas = "betae", output_re = FALSE,
    additional_tests = character(0), lrt = FALSE, lrt_contrasts = NULL,
    lrt_include_scale = FALSE, lrt_details = FALSE,
    adj_ridge_factor = 1e-8, adj_cond_limit = 1e8,
    keep_diagnostics = TRUE) {

  if (!use_betas %in% c("betae", "final_betas")) {
    stop("use_betas must be 'betae' or 'final_betas'.")
  }
  if (identical(use_betas, "final_betas")) {
    warning(
      "use_betas is ignored by the estimator-aligned first-order HL path; retained for API compatibility.",
      call. = FALSE
    )
  }
  extra <- unique(tolower(additional_tests))
  extra <- extra[nzchar(extra)]
  if (isTRUE(lrt)) extra <- union(extra, "lrt")
  if (isTRUE(lrt_details)) extra <- union(extra, "lrt")
  if (!is.null(lrt_contrasts) && length(lrt_contrasts)) extra <- union(extra, "lrt")
  invalid <- setdiff(extra, c("score", "lrt"))
  if (length(invalid)) stop("Unsupported entries in additional_tests: ",
                            paste(invalid, collapse = ", "), ".")
  want_score <- "score" %in% extra
  want_lrt <- "lrt" %in% extra
  if (ncore > 1L && verbose) {
    message("scrope_hl first-order multistarts are deterministic and currently run sequentially.")
  }

  prep <- hl_prepare_data(
    count = count, id = id, pred = pred, offset = offset,
    min = min, max = max, opt = opt, verbose = verbose,
    cpc = cpc, mincp = mincp, output_re = output_re, reml = 0
  )
  ctx <- prep$ctx
  count_t <- prep$count
  gids <- prep$gid
  p <- ctx$nb
  predictor_names <- prep$predictor_names
  if (is.null(predictor_names)) predictor_names <- paste0("X", seq_len(p))
  predictor_names[ctx$intcol] <- if (is.null(prep$predictor_names)) "X1" else
    predictor_names[ctx$intcol]
  contrasts <- if (want_score || want_lrt ||
                   (!is.null(lrt_contrasts) && length(lrt_contrasts))) {
    hl_fo_public_contrasts(p, ctx$intcol, predictor_names, lrt_contrasts)
  } else list()

  gene_fits <- vector("list", length(gids))
  summary_rows <- vector("list", length(gids))
  overdispersion <- matrix(NA_real_, length(gids), 2L,
                           dimnames = list(NULL, c("Subject", "Cell")))
  convergence <- rep(-25, length(gids))
  algorithm <- rep(2, length(gids))
  random_effect <- if (output_re) matrix(NA_real_, length(gids), ctx$k) else NULL
  test_objects <- vector("list", length(gids))
  coefficient_wald <- vector("list", length(gids))

  for (gg in seq_along(gids)) {
    gene_id <- gids[gg]
    gene_name <- prep$gene_names[gg]
    posv <- call_posindy(count_t, gene_id - 1L, ctx$nind)
    y <- hl_build_counts(posv, ctx$nind)
    data <- hl_fo_make_data(y, ctx$pred, ctx$offset, ctx$id, ctx$intcol)
    fit <- hl_fo_fit(
      data,
      active = c("tau", "xi"),
      guard_min = min,
      guard_max = max
    )
    fit$gene_id <- gene_id
    fit$gene <- gene_name
    fit$design_transform <- ctx$design_transform
    gene_fits[[gg]] <- fit
    test_objects[[gg]] <- list()
    coefficient_wald[[gg]] <- list()

    if (!isTRUE(fit$point_success)) {
      summary_rows[[gg]] <- hl_fo_public_failure_row(
        p, predictor_names, gene_id, gene_name,
        if (is.null(fit$fit_status)) fit$status else fit$fit_status,
        if (is.null(fit$fit_failure_reason)) fit$reason else fit$fit_failure_reason
      )
      next
    }

    beta_report <- as.numeric(ctx$M_X %*% fit$beta)
    names(beta_report) <- predictor_names
    overdispersion[gg, ] <- c(exp(fit$rho_full[1]), exp(-fit$rho_full[2]))
    convergence[gg] <- 1
    if (output_re) random_effect[gg, ] <- fit$eta

    V_report <- NULL
    if (isTRUE(fit$covariance_success)) {
      V_report <- ctx$M_X %*% fit$covariance$V_beta %*% t(ctx$M_X)
    }
    row <- list()
    for (jj in seq_len(p)) {
      nm <- predictor_names[jj]
      row[[paste0("logFC_", nm)]] <- beta_report[jj]
      row[[paste0("se_", nm)]] <- NA_real_
      row[[paste0("p_naive_", nm)]] <- NA_real_
      rr <- matrix(0, 1L, p)
      rr[1L, jj] <- 1
      wald <- if (is.null(V_report)) {
        list(success = FALSE, statistic = NA_real_, p_value = NA_real_,
             status = fit$covariance_status, df = 1L,
             n_subjects = data$m, contrast_rank = 1L,
             covariance_rank = NA_integer_)
      } else {
        hl_fo_robust_wald(beta_report, V_report, rr, 0,
                          n_subjects = data$m)
      }
      row[[paste0("se_robust_", nm)]] <- if (is.null(V_report)) NA_real_ else
        sqrt(max(0, V_report[jj, jj]))
      row[[paste0("wald_stat_", nm)]] <- wald$statistic
      row[[paste0("p_robust_", nm)]] <- wald$p_value
      row[[paste0("wald_df_", nm)]] <- wald$df
      row[[paste0("wald_m_", nm)]] <- wald$n_subjects
      row[[paste0("wald_contrast_rank_", nm)]] <- wald$contrast_rank
      row[[paste0("wald_covariance_rank_", nm)]] <- wald$covariance_rank
      row[[paste0("wald_success_", nm)]] <- isTRUE(wald$success)
      row[[paste0("wald_failure_", nm)]] <- if (isTRUE(wald$success)) NA_character_ else
        wald$status
      coefficient_wald[[gg]][[nm]] <- wald
    }

    for (label in names(contrasts)) {
      def <- contrasts[[label]]
      wald <- if (is.null(V_report)) {
        list(success = FALSE, statistic = NA_real_, p_value = NA_real_,
             status = fit$covariance_status, df = NA_integer_,
             contrast_rank = NA_integer_, covariance_rank = NA_integer_)
      } else {
        hl_fo_robust_wald(beta_report, V_report, def$R, def$rhs,
                          n_subjects = data$m)
      }
      row[[paste0("wald_stat_", label)]] <- wald$statistic
      row[[paste0("p_robust_", label)]] <- wald$p_value
      row[[paste0("wald_df_", label)]] <- wald$df
      row[[paste0("wald_m_", label)]] <- if (is.null(wald$n_subjects)) data$m else
        wald$n_subjects
      row[[paste0("wald_contrast_rank_", label)]] <- wald$contrast_rank
      row[[paste0("wald_covariance_rank_", label)]] <- wald$covariance_rank
      row[[paste0("wald_success_", label)]] <- isTRUE(wald$success)
      row[[paste0("wald_failure_", label)]] <- if (isTRUE(wald$success)) NA_character_ else
        wald$status

      score <- NULL
      if (want_score) {
        R_internal <- def$R %*% ctx$M_X
        if (exists("hl_fo_score_test", mode = "function")) {
          score <- hl_fo_score_test(
            data = data, R = R_internal, rhs = def$rhs,
            unrestricted_fit = fit, guard_min = min, guard_max = max
          )
        } else {
          score <- list(success = FALSE, statistic = NA_real_, p_value = NA_real_,
                        status = "score_null_no_convergence",
                        reason = "score_solver_unavailable")
        }
        row[[paste0("hl_score_adj_", label)]] <- score$statistic
        row[[paste0("p_hl_score_", label)]] <- score$p_value
        row[[paste0("hl_score_success_", label)]] <- isTRUE(score$success)
        row[[paste0("hl_score_failure_", label)]] <- if (isTRUE(score$success))
          NA_character_ else if (!is.null(score$status)) score$status else score$reason
      }

      objective_ratio <- if (want_lrt) hl_fo_objective_ratio_unsupported(
        raw_difference = NA_real_, failure_reason = "objective_ratio_unsupported",
        wald = wald
      ) else NULL
      if (want_lrt) {
        row[[paste0("hl_lrt_adj_", label)]] <- NA_real_
        row[[paste0("p_hl_lrt_", label)]] <- NA_real_
        row[[paste0("hl_lrt_method_", label)]] <- "unsupported"
        row[[paste0("hl_lrt_success_", label)]] <- FALSE
        row[[paste0("hl_lrt_failure_", label)]] <- "objective_ratio_unsupported"
        if (lrt_include_scale) row[[paste0("hl_lrt_scale_", label)]] <- NA_real_
      }
      test_objects[[gg]][[label]] <- list(
        wald = wald, score = score, objective_ratio = objective_ratio,
        contrast_report = def,
        contrast_internal = list(R = def$R %*% ctx$M_X, rhs = def$rhs)
      )
    }

    row$gene_id <- gene_id
    row$gene <- gene_name
    row$fit_success <- TRUE
    row$fit_status <- fit$fit_status
    row$fit_failure <- NA_character_
    row$covariance_success <- isTRUE(fit$covariance_success)
    row$covariance_status <- fit$covariance_status
    summary_rows[[gg]] <- as.data.frame(row, check.names = FALSE,
                                        stringsAsFactors = FALSE)
  }

  all_names <- unique(unlist(lapply(summary_rows, names), use.names = FALSE))
  summary_rows <- lapply(summary_rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (nm in missing) row[[nm]] <- NA
    row[, all_names, drop = FALSE]
  })
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  overdispersion <- as.data.frame(overdispersion)

  collect_status <- function(success_prefix, failure_prefix, test_name) {
    success_columns <- grep(paste0("^", success_prefix), names(summary_df), value = TRUE)
    if (!length(success_columns)) return(NULL)
    blocks <- lapply(success_columns, function(success_column) {
      label <- substring(success_column, nchar(success_prefix) + 1L)
      failure_column <- paste0(failure_prefix, label)
      success <- as.logical(summary_df[[success_column]])
      failure <- if (failure_column %in% names(summary_df))
        as.character(summary_df[[failure_column]]) else rep(NA_character_, nrow(summary_df))
      data.frame(
        gene_id = summary_df$gene_id,
        gene = summary_df$gene,
        contrast = label,
        test = test_name,
        success = success,
        failure_reason = failure,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, blocks)
  }
  status_blocks <- Filter(Negate(is.null), list(
    collect_status("wald_success_", "wald_failure_", "wald"),
    collect_status("hl_score_success_", "hl_score_failure_", "score"),
    collect_status("hl_lrt_success_", "hl_lrt_failure_", "objective_ratio")
  ))
  test_status <- if (length(status_blocks)) do.call(rbind, status_blocks) else
    data.frame(gene_id = integer(), gene = character(), contrast = character(),
               test = character(), success = logical(), failure_reason = character(),
               stringsAsFactors = FALSE)

  summary_df$test_status <- "ok"
  summary_df$test_failure_detail <- NA_character_
  for (gg in seq_len(nrow(summary_df))) {
    rows <- test_status[test_status$gene_id == summary_df$gene_id[gg], , drop = FALSE]
    failures <- rows[!is.na(rows$success) & !rows$success, , drop = FALSE]
    if (!nrow(failures)) next
    types <- unique(failures$test)
    summary_df$test_status[gg] <- paste(types, collapse = ",")
    detail <- vapply(types, function(type) {
      values <- unique(na.omit(failures$failure_reason[failures$test == type]))
      paste0(type, ":", if (length(values)) paste(values, collapse = ",") else "unknown")
    }, character(1))
    summary_df$test_failure_detail[gg] <- paste(detail, collapse = ";")
  }

  details <- NULL
  if (isTRUE(lrt_details)) {
    details <- test_status
  }
  if (want_lrt && length(contrasts) && verbose) {
    message("HL beta-constraint objective-ratio inference is unsupported; robust Wald remains labeled as Wald.")
  }
  if (want_score && any(test_status$test == "score" & !test_status$success) && verbose) {
    message("Some first-order HL Score null roots were unavailable; see score status fields.")
  }

  list(
    summary = summary_df,
    overdispersion = overdispersion,
    convergence = convergence,
    algorithm = algorithm,
    random_effect = if (output_re) as.data.frame(random_effect) else NULL,
    hl_adj_details = details,
    diagnostics = if (keep_diagnostics) list(
      hl = gene_fits,
      first_order = gene_fits,
      godambe = lapply(gene_fits, function(x) if (!is.null(x$augmented))
        list(A = x$augmented$A, J = x$augmented$J,
             covariance = x$covariance, coordinate_tag = x$coordinate_tag) else NULL),
      hl_adj = test_objects,
      wald = coefficient_wald,
      test_status = test_status,
      design_transform = ctx$design_transform
    ) else NULL
  )
}
