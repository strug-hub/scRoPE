#!/usr/bin/env Rscript

# Targeted common-random-number validation for first-order HL Wald/Score
# equivalence.  This is deliberately outside testthat: the locked 250-replicate
# gate is a resumable validation job, not an ordinary package unit test.

.hl_mc_abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

.hl_mc_scalar_character <- function(x, default = NA_character_) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) default else as.character(x[[1L]])
}

.hl_mc_scalar_numeric <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || !is.finite(x[[1L]])) default else as.numeric(x[[1L]])
}

.hl_mc_scalar_integer <- function(x, default = NA_integer_) {
  value <- .hl_mc_scalar_numeric(x, default = default)
  if (is.na(value)) default else as.integer(value)
}

.hl_mc_flag <- function(x) {
  isTRUE(x)
}

.hl_mc_usage <- function() {
  paste(
    "Usage:",
    "  Rscript tests/validation/run_hl_first_order_wald_score_mc.R \\",
    "    --mode=smoke --replicates=2 --workers=2",
    "",
    "The locked acceptance run is intentionally guarded:",
    "  Rscript tests/validation/run_hl_first_order_wald_score_mc.R \\",
    "    --mode=full --replicates=250 --workers=8 \\",
    "    --confirm-full=RUN_LOCKED_250",
    "",
    "Optional: --repo-root=PATH --output=PATH",
    sep = "\n"
  )
}

.hl_mc_parse_args <- function(args) {
  out <- list(
    mode = "smoke",
    replicates = 2L,
    workers = 1L,
    repo_root = getwd(),
    output = NULL,
    confirm_full = ""
  )
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat(.hl_mc_usage(), "\n")
      quit(save = "no", status = 0L)
    }
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pair) != 2L) .hl_mc_abort("Malformed argument: ", arg)
    key <- gsub("-", "_", pair[[1L]], fixed = TRUE)
    if (!key %in% names(out)) .hl_mc_abort("Unknown argument: ", arg)
    out[[key]] <- pair[[2L]]
  }
  out$replicates <- suppressWarnings(as.integer(out$replicates))
  out$workers <- suppressWarnings(as.integer(out$workers))
  if (!out$mode %in% c("smoke", "full")) {
    .hl_mc_abort("--mode must be smoke or full.")
  }
  if (is.na(out$replicates) || out$replicates < 1L) {
    .hl_mc_abort("--replicates must be a positive integer.")
  }
  if (is.na(out$workers) || out$workers < 1L) {
    .hl_mc_abort("--workers must be a positive integer.")
  }
  if (out$mode == "smoke" && out$replicates > 2L) {
    .hl_mc_abort("Smoke mode is capped at two replicates. Use the explicitly guarded full mode later.")
  }
  if (out$mode == "full" &&
      (out$replicates != 250L || out$confirm_full != "RUN_LOCKED_250")) {
    .hl_mc_abort(
      "The acceptance run requires --replicates=250 and ",
      "--confirm-full=RUN_LOCKED_250."
    )
  }
  out$workers <- min(out$workers, out$replicates)
  out$repo_root <- normalizePath(out$repo_root, mustWork = TRUE)
  if (!file.exists(file.path(out$repo_root, "DESCRIPTION"))) {
    .hl_mc_abort("--repo-root is not an R package root: ", out$repo_root)
  }
  if (is.null(out$output)) {
    out$output <- file.path(out$repo_root, "temp", "hl_first_order_validation")
  }
  out$output <- normalizePath(out$output, mustWork = FALSE)
  out
}

.hl_mc_atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3L)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      .hl_mc_abort("Could not write artifact: ", path)
    }
    unlink(tmp)
  }
  invisible(path)
}

.hl_mc_atomic_write_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(object, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      .hl_mc_abort("Could not write artifact: ", path)
    }
    unlink(tmp)
  }
  invisible(path)
}

.hl_mc_rng_streams <- function(n, master_seed) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(as.integer(master_seed))
  streams <- vector("list", n)
  streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
  if (n > 1L) {
    for (ii in 2:n) streams[[ii]] <- parallel::nextRNGStream(streams[[ii - 1L]])
  }
  streams
}

.hl_mc_with_stream <- function(stream, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  assign(".Random.seed", stream, envir = .GlobalEnv)
  force(code)
}

.hl_mc_locked_design <- function() {
  list(
    schema_version = 1L,
    master_seed = 8675309L,
    rng = c(
      kind = "L'Ecuyer-CMRG",
      normal.kind = "Inversion",
      sample.kind = "Rejection"
    ),
    m_grid = c(250L, 500L, 1000L, 2000L),
    n_i = 5L,
    beta = c(`(Intercept)` = -1.2, z = 0.35, w = -0.25),
    sigma2 = 0.4,
    phi = 8,
    contrast = matrix(c(0, 1, 0), nrow = 1L,
      dimnames = list("z_at_truth", c("(Intercept)", "z", "w"))),
    rhs = c(z_at_truth = 0.35),
    z_distribution = "subject-level Bernoulli(0.5)",
    w_distribution = "cell-level standard Normal",
    log_offset_distribution = "cell-level Uniform(-0.3, 0.3)",
    subject_effect = paste(
      "Gamma(shape=1/expm1(sigma2),",
      "rate=shape*exp(-sigma2/2))"
    ),
    count_distribution = "NB(mean=offset*exp(X beta)*omega, size=phi)",
    active_dispersion = character(0),
    fixed_natural = c(sigma2 = 0.4, phi = 8),
    dispersion_choice = paste(
      "q=0: sigma2 and phi are fixed at their generating truth.",
      "This is an explicitly supported edge case and isolates the first-order",
      "equivalence of the unrestricted augmented-sandwich Wald statistic and",
      "the separate full-chi Score null root without outer-dispersion noise."
    ),
    common_random_numbers = paste(
      "Each replicate is generated exactly once at m=2000; the first m",
      "subjects and their five cells are used for each smaller m."
    ),
    discrepancy = "abs(W-Q_EE)/(1+W+Q_EE)",
    acceptance = list(
      replicates = 250L,
      D_2000_max = 0.05,
      log_log_slope_max = -0.25
    )
  )
}

.hl_mc_generate_max <- function(replicate_id, stream, design) {
  .hl_mc_with_stream(stream, {
    m <- max(design$m_grid)
    n_i <- design$n_i
    n <- m * n_i
    subject <- rep(seq_len(m), each = n_i)
    z_subject <- stats::rbinom(m, size = 1L, prob = 0.5)
    z <- rep(z_subject, each = n_i)
    w <- stats::rnorm(n)
    X <- cbind(`(Intercept)` = 1, z = z, w = w)
    log_offset <- stats::runif(n, min = -0.3, max = 0.3)
    offset <- exp(log_offset)
    alpha <- 1 / expm1(design$sigma2)
    lambda <- alpha * exp(-design$sigma2 / 2)
    omega <- stats::rgamma(m, shape = alpha, rate = lambda)
    eta <- log(omega)
    mu <- offset * exp(drop(X %*% design$beta)) * omega[subject]
    y <- stats::rnbinom(n, size = design$phi, mu = mu)
    list(
      schema_version = design$schema_version,
      replicate = as.integer(replicate_id),
      rng_stream = stream,
      rng = design$rng,
      m_max = as.integer(m),
      n_i = as.integer(n_i),
      subject = subject,
      z_subject = z_subject,
      z = z,
      w = w,
      X = X,
      log_offset = log_offset,
      offset = offset,
      alpha = alpha,
      lambda = lambda,
      omega = omega,
      eta = eta,
      mu = mu,
      y = y,
      truth = list(
        beta = design$beta,
        sigma2 = design$sigma2,
        phi = design$phi,
        tau = log(design$sigma2),
        xi = log(design$phi)
      )
    )
  })
}

.hl_mc_validate_generated <- function(object, replicate_id, design) {
  expected_n <- max(design$m_grid) * design$n_i
  is.list(object) && identical(object$schema_version, design$schema_version) &&
    identical(object$replicate, as.integer(replicate_id)) &&
    identical(object$m_max, as.integer(max(design$m_grid))) &&
    identical(object$n_i, as.integer(design$n_i)) &&
    identical(dim(object$X), c(expected_n, length(design$beta))) &&
    length(object$y) == expected_n && length(object$subject) == expected_n &&
    length(object$omega) == max(design$m_grid) &&
    identical(as.numeric(object$truth$beta), as.numeric(design$beta)) &&
    identical(object$truth$sigma2, design$sigma2) &&
    identical(object$truth$phi, design$phi)
}

.hl_mc_error_rows <- function(replicate_id, design, status, message) {
  do.call(rbind, lapply(design$m_grid, function(m) {
    data.frame(
      replicate = as.integer(replicate_id), m = as.integer(m),
      fit_point_success = FALSE, fit_status = status, fit_reason = message,
      covariance_success = FALSE, covariance_status = status,
      wald_success = FALSE, wald_status = "not_run", wald_statistic = NA_real_,
      wald_p_value = NA_real_, wald_df = NA_integer_,
      wald_contrast_rank = NA_integer_, wald_covariance_rank = NA_integer_,
      score_success = FALSE, score_status = "not_run", score_statistic = NA_real_,
      score_p_value = NA_real_, score_df = NA_integer_,
      score_contrast_rank = NA_integer_, score_covariance_rank = NA_integer_,
      score_constrained_fit_status = NA_character_, score_null_fit_status = NA_character_,
      discrepancy = NA_real_, elapsed_fit_seconds = NA_real_,
      elapsed_wald_seconds = NA_real_, elapsed_score_seconds = NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

.hl_mc_fit_one_m <- function(generated, m, design, api) {
  cells <- seq_len(as.integer(m) * design$n_i)
  data <- api$make_data(
    y = generated$y[cells],
    X = generated$X[cells, , drop = FALSE],
    log_offset = generated$log_offset[cells],
    subject = generated$subject[cells],
    intcol = 1L
  )

  fit_error <- NULL
  fit_time <- system.time({
    fit <- tryCatch(
      api$fit(
        data = data,
        active = character(0),
        fixed_natural = design$fixed_natural
      ),
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
  })[["elapsed"]]

  if (is.null(fit)) {
    row <- .hl_mc_error_rows(generated$replicate, design, "fit_error", fit_error)
    row <- row[row$m == m, , drop = FALSE]
    row$elapsed_fit_seconds <- fit_time
    return(row)
  }

  wald_error <- NULL
  wald_time <- system.time({
    wald <- tryCatch(
      api$wald(fit, R = design$contrast, rhs = design$rhs),
      error = function(e) {
        wald_error <<- conditionMessage(e)
        NULL
      }
    )
  })[["elapsed"]]

  score_error <- NULL
  score_time <- system.time({
    score <- tryCatch(
      api$score(
        data = data,
        R = design$contrast,
        rhs = design$rhs,
        unrestricted_fit = fit,
        active = character(0),
        fixed_natural = design$fixed_natural
      ),
      error = function(e) {
        score_error <<- conditionMessage(e)
        NULL
      }
    )
  })[["elapsed"]]

  wald_stat <- if (is.null(wald)) NA_real_ else .hl_mc_scalar_numeric(wald$stat)
  score_stat <- if (is.null(score)) NA_real_ else .hl_mc_scalar_numeric(score$stat)
  discrepancy <- if (isTRUE(!is.null(wald) && wald$success) &&
      isTRUE(!is.null(score) && score$success) &&
      is.finite(wald_stat) && is.finite(score_stat)) {
    abs(wald_stat - score_stat) / (1 + wald_stat + score_stat)
  } else {
    NA_real_
  }

  data.frame(
    replicate = as.integer(generated$replicate),
    m = as.integer(m),
    fit_point_success = .hl_mc_flag(fit$point_success),
    fit_status = .hl_mc_scalar_character(fit$fit_status,
      default = .hl_mc_scalar_character(fit$status, "missing_status")),
    fit_reason = .hl_mc_scalar_character(fit$fit_failure_reason,
      default = .hl_mc_scalar_character(fit$reason, NA_character_)),
    covariance_success = .hl_mc_flag(fit$covariance_success),
    covariance_status = .hl_mc_scalar_character(fit$covariance_status, "missing_status"),
    wald_success = .hl_mc_flag(if (is.null(wald)) FALSE else wald$success),
    wald_status = if (is.null(wald)) paste0("wald_error: ", wald_error) else
      .hl_mc_scalar_character(wald$status, "missing_status"),
    wald_statistic = wald_stat,
    wald_p_value = if (is.null(wald)) NA_real_ else .hl_mc_scalar_numeric(wald$p),
    wald_df = if (is.null(wald)) NA_integer_ else .hl_mc_scalar_integer(wald$df),
    wald_contrast_rank = if (is.null(wald)) NA_integer_ else
      .hl_mc_scalar_integer(wald$contrast_rank),
    wald_covariance_rank = if (is.null(wald)) NA_integer_ else
      .hl_mc_scalar_integer(wald$covariance_rank),
    score_success = .hl_mc_flag(if (is.null(score)) FALSE else score$success),
    score_status = if (is.null(score)) paste0("score_error: ", score_error) else
      .hl_mc_scalar_character(score$status, "missing_status"),
    score_statistic = score_stat,
    score_p_value = if (is.null(score)) NA_real_ else .hl_mc_scalar_numeric(score$p),
    score_df = if (is.null(score)) NA_integer_ else .hl_mc_scalar_integer(score$df),
    score_contrast_rank = if (is.null(score)) NA_integer_ else
      .hl_mc_scalar_integer(score$contrast_rank),
    score_covariance_rank = if (is.null(score)) NA_integer_ else
      .hl_mc_scalar_integer(score$covariance_rank),
    score_constrained_fit_status = if (is.null(score)) NA_character_ else
      .hl_mc_scalar_character(score$constrained_fit_status),
    score_null_fit_status = if (is.null(score)) NA_character_ else
      .hl_mc_scalar_character(score$null_fit_status),
    discrepancy = discrepancy,
    elapsed_fit_seconds = fit_time,
    elapsed_wald_seconds = wald_time,
    elapsed_score_seconds = score_time,
    stringsAsFactors = FALSE
  )
}

.hl_mc_run_replicate <- function(replicate_id, stream, design, output, api) {
  generated_path <- file.path(
    output, "generated", sprintf("replicate_%03d_m2000.rds", replicate_id)
  )
  result_path <- file.path(
    output, "per_replicate", sprintf("replicate_%03d_results.rds", replicate_id)
  )

  if (file.exists(generated_path)) {
    generated <- readRDS(generated_path)
    if (!.hl_mc_validate_generated(generated, replicate_id, design)) {
      .hl_mc_abort("Existing generated artifact failed validation: ", generated_path)
    }
  } else {
    generated <- .hl_mc_generate_max(replicate_id, stream, design)
    .hl_mc_atomic_save_rds(generated, generated_path)
  }
  generated_md5 <- unname(tools::md5sum(generated_path))

  if (file.exists(result_path)) {
    existing <- readRDS(result_path)
    if (identical(existing$generated_md5, generated_md5) &&
        identical(sort(as.integer(existing$rows$m)), sort(as.integer(design$m_grid)))) {
      return(list(replicate = replicate_id, resumed = TRUE, path = result_path))
    }
  }

  rows <- tryCatch(
    do.call(rbind, lapply(design$m_grid, function(m) {
      .hl_mc_fit_one_m(generated, m, design, api)
    })),
    error = function(e) .hl_mc_error_rows(
      replicate_id, design, "replicate_error", conditionMessage(e)
    )
  )
  artifact <- list(
    schema_version = design$schema_version,
    replicate = as.integer(replicate_id),
    generated_path = generated_path,
    generated_md5 = generated_md5,
    rows = rows
  )
  .hl_mc_atomic_save_rds(artifact, result_path)
  list(replicate = replicate_id, resumed = FALSE, path = result_path)
}

.hl_mc_summarize <- function(rows, design, requested_replicates, mode) {
  by_m <- do.call(rbind, lapply(design$m_grid, function(m) {
    part <- rows[rows$m == m, , drop = FALSE]
    valid <- part$wald_success & part$score_success & is.finite(part$discrepancy)
    data.frame(
      m = as.integer(m),
      n_replicates = nrow(part),
      n_fit_point_success = sum(part$fit_point_success),
      n_covariance_success = sum(part$covariance_success),
      n_wald_success = sum(part$wald_success),
      n_score_success = sum(part$score_success),
      n_valid_pairs = sum(valid),
      fit_point_success_rate = mean(part$fit_point_success),
      covariance_success_rate = mean(part$covariance_success),
      wald_success_rate = mean(part$wald_success),
      score_success_rate = mean(part$score_success),
      D_m = if (any(valid)) stats::median(part$discrepancy[valid]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  slope <- NA_real_
  if (all(is.finite(by_m$D_m)) && all(by_m$D_m > 0)) {
    slope <- unname(stats::coef(stats::lm(log(D_m) ~ log(m), data = by_m))[[2L]])
  }
  required <- design$acceptance$replicates
  enough <- identical(as.integer(requested_replicates), as.integer(required)) &&
    all(by_m$n_replicates == required) && all(by_m$n_valid_pairs == required)
  d2000 <- by_m$D_m[by_m$m == 2000L]
  evaluated <- identical(mode, "full") && enough && is.finite(d2000) && is.finite(slope)
  passed <- evaluated && d2000 <= design$acceptance$D_2000_max &&
    slope <= design$acceptance$log_log_slope_max

  list(
    by_m = by_m,
    slope = slope,
    acceptance_evaluated = evaluated,
    acceptance_passed = passed,
    acceptance_status = if (mode != "full") {
      "smoke_only_not_acceptance"
    } else if (!enough) {
      "full_run_incomplete_or_contains_failed_pairs"
    } else if (!is.finite(slope)) {
      "slope_not_estimable"
    } else if (passed) {
      "passed"
    } else {
      "failed_threshold"
    },
    thresholds = design$acceptance
  )
}

.hl_mc_load_worktree <- function(repo_root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    load_failure <- NULL
    loaded <- tryCatch({
      pkgload::load_all(repo_root, recompile = FALSE, helpers = FALSE, quiet = TRUE)
      TRUE
    }, error = function(e) {
      load_failure <<- conditionMessage(e)
      FALSE
    })
    if (isTRUE(loaded)) {
      return(list(
        method = "pkgload::load_all",
        reason = NULL,
        api = list(
          make_data = getFromNamespace("hl_fo_make_data", "scRoPE"),
          fit = getFromNamespace("hl_fo_fit", "scRoPE"),
          wald = getFromNamespace("hl_fo_unrestricted_wald", "scRoPE"),
          score = getFromNamespace("hl_fo_score_test", "scRoPE")
        )
      ))
    }
    load_reason <- load_failure
  } else {
    load_reason <- "pkgload is unavailable"
  }

  # A compiler is not required for the q=0 validation: the solver contains an
  # exact R reference path and dispatches to C++ only when its registered
  # symbol is loaded.  This fallback sources only the corrected first-order
  # worktree files into one private environment, in dependency order.
  worktree <- new.env(parent = globalenv())
  files <- c(
    "hl_first_order_utils.R",
    "hl_first_order_primitives.R",
    "hl_first_order_solver.R",
    "hl_first_order_inference.R",
    "hl_first_order_score.R"
  )
  paths <- file.path(repo_root, "R", files)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    .hl_mc_abort("Missing corrected first-order source file: ", missing[[1L]])
  }
  for (path in paths) sys.source(path, envir = worktree, keep.source = TRUE)
  list(
    method = "private sourced exact-R reference path",
    reason = load_reason,
    api = list(
      make_data = get("hl_fo_make_data", envir = worktree, inherits = FALSE),
      fit = get("hl_fo_fit", envir = worktree, inherits = FALSE),
      wald = get("hl_fo_unrestricted_wald", envir = worktree, inherits = FALSE),
      score = get("hl_fo_score_test", envir = worktree, inherits = FALSE)
    )
  )
}

.hl_mc_main <- function() {
  options(warn = 1)
  args <- .hl_mc_parse_args(commandArgs(trailingOnly = TRUE))
  if (.Platform$OS.type == "windows" && args$workers > 1L) {
    .hl_mc_abort("Parallel mode uses deterministic fork workers; use --workers=1 on Windows.")
  }
  dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(args$output, "generated"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(args$output, "per_replicate"), recursive = TRUE, showWarnings = FALSE)

  design <- .hl_mc_locked_design()
  metadata_path <- file.path(args$output, "design_metadata.rds")
  metadata <- list(
    artifact = "first_order_hl_wald_score_common_random_number_validation",
    locked_design = design,
    q_zero_is_intentional = TRUE,
    q_zero_rationale = design$dispersion_choice
  )
  if (file.exists(metadata_path)) {
    previous <- readRDS(metadata_path)
    if (!identical(previous$locked_design, metadata$locked_design)) {
      .hl_mc_abort("Existing output has incompatible locked design metadata: ", metadata_path)
    }
  } else {
    .hl_mc_atomic_save_rds(metadata, metadata_path)
  }

  loader <- .hl_mc_load_worktree(args$repo_root)
  api <- loader$api

  streams <- .hl_mc_rng_streams(args$replicates, design$master_seed)
  worker <- function(ii) {
    .hl_mc_run_replicate(ii, streams[[ii]], design, args$output, api)
  }
  ids <- seq_len(args$replicates)
  started <- Sys.time()
  run_records <- if (args$workers == 1L) {
    lapply(ids, worker)
  } else {
    parallel::mclapply(
      ids, worker, mc.cores = args$workers, mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  }

  result_paths <- file.path(
    args$output, "per_replicate", sprintf("replicate_%03d_results.rds", ids)
  )
  if (!all(file.exists(result_paths))) {
    .hl_mc_abort("One or more per-replicate result artifacts are missing.")
  }
  per_replicate <- lapply(result_paths, readRDS)
  rows <- do.call(rbind, lapply(per_replicate, `[[`, "rows"))
  rows <- rows[order(rows$replicate, rows$m), , drop = FALSE]
  rownames(rows) <- NULL
  summary <- .hl_mc_summarize(rows, design, args$replicates, args$mode)
  summary_table <- summary$by_m
  summary_table$log_log_slope <- summary$slope
  summary_table$acceptance_evaluated <- summary$acceptance_evaluated
  summary_table$acceptance_passed <- summary$acceptance_passed
  summary_table$acceptance_status <- summary$acceptance_status

  .hl_mc_atomic_save_rds(
    list(metadata = metadata, results = rows),
    file.path(args$output, "fit_results.rds")
  )
  .hl_mc_atomic_write_csv(rows, file.path(args$output, "fit_results.csv"))
  .hl_mc_atomic_save_rds(
    list(metadata = metadata, summary = summary),
    file.path(args$output, "wald_score_summary.rds")
  )
  .hl_mc_atomic_write_csv(
    summary_table, file.path(args$output, "wald_score_summary.csv")
  )
  manifest <- list(
    mode = args$mode,
    requested_replicates = args$replicates,
    workers = args$workers,
    started = started,
    finished = Sys.time(),
    resumed_replicates = sum(vapply(run_records, `[[`, logical(1), "resumed")),
    new_replicates = sum(!vapply(run_records, `[[`, logical(1), "resumed")),
    output = args$output,
    package_version = unname(read.dcf(
      file.path(args$repo_root, "DESCRIPTION"), fields = "Version"
    )[[1L]]),
    worktree_loader = loader$method,
    worktree_loader_fallback_reason = loader$reason,
    session_info = utils::capture.output(sessionInfo())
  )
  .hl_mc_atomic_save_rds(
    manifest,
    file.path(args$output, paste0("run_manifest_", args$mode, ".rds"))
  )

  cat("Artifacts:", args$output, "\n")
  print(summary_table, row.names = FALSE)
  cat("Acceptance status:", summary$acceptance_status, "\n")
  invisible(list(results = rows, summary = summary, manifest = manifest))
}

.hl_mc_main()
