auto_pgmm_stream_function <- function(name) {
  get(name, envir = environment(getFromNamespace("scrope_pmg_internal", "scRoPE")))
}

auto_pgmm_stream_fixture <- function() {
  # Distinct genes and a filtered first row exercise retained gene IDs/order.
  y <- c(1, 4, 1, 5, 4, 2, 9, 12, 0, 3, 11, 14,
         1, 1, 7, 11, 3, 2, 11, 17, 2, 4, 6, 12)
  list(count = rbind(filtered = 0 * y, first = y, second = y + 1),
       id = rep(paste0("s", 1:6), each = 4),
       pred = cbind(Intercept = 1, x = rep(c(-1.5, -.5, .5, 1.5), 6)))
}

auto_pgmm_stream_observer <- function(native, prep, gene_id, candidates) {
  # This callback has no mutable parent-process state. Its return value is the
  # only source of the records checked below, including with PSOCK workers.
  list(gene_id = gene_id, candidates = length(candidates),
       donor_count = nrow(native$per_subject_gradients),
       optional_tests = names(native$lrt),
       candidates_removed = is.null(native$working_candidates),
       worker_pid = Sys.getpid(), random_draw = runif(1))
}

auto_pgmm_stream_run <- function(keep_diagnostics, ncore = 1L, observer = TRUE) {
  fx <- auto_pgmm_stream_fixture()
  args <- c(fx, list(verbose = FALSE, cpc = 0, mincp = 1,
    keep_diagnostics = keep_diagnostics, ncore = ncore,
    pgmm_control = list(iter.max = 150L, eval.max = 200L),
    additional_tests = c("score", "lrt"), lrt_details = TRUE,
    lrt_include_scale = TRUE))
  if (observer) args$.fit_observer <- auto_pgmm_stream_observer
  set.seed(11226)
  suppressWarnings(do.call(auto_pgmm_stream_function("scrope_pmg_internal"), args))
}

test_that("PGMM compact optional tests and pure records equal full native output", {
  full <- auto_pgmm_stream_run(TRUE)
  native_full <- auto_pgmm_stream_run(TRUE, observer = FALSE)
  compact <- auto_pgmm_stream_run(FALSE)
  ordinary <- auto_pgmm_stream_run(FALSE, observer = FALSE)
  expect_null(compact$diagnostics)
  expect_null(ordinary$diagnostics)
  expect_identical(full$diagnostics, native_full$diagnostics)
  expect_identical(compact$summary, full$summary)
  expect_identical(compact$lrt_details, full$lrt_details)
  expect_identical(compact$overdispersion, full$overdispersion)
  expect_identical(compact$fit_provenance, full$fit_provenance)
  expect_identical(ordinary$summary, full$summary)
  expect_identical(ordinary$lrt_details, full$lrt_details)
  expect_identical(ordinary$fit_provenance, full$fit_provenance)
  expect_null(attr(ordinary, ".working_records"))
  records <- attr(compact, ".working_records")
  expect_identical(records, attr(full, ".working_records"))
  expect_identical(names(records), as.character(compact$summary$gene_id))
  expect_identical(vapply(records, `[[`, integer(1), "gene_id"),
                   stats::setNames(compact$summary$gene_id, names(records)))
  expect_true(all(vapply(records, function(x) x$candidates >= 1L, logical(1))))
  expect_true(all(vapply(records, `[[`, logical(1), "candidates_removed")))
  expect_true(all(vapply(records, function(x) identical(x$optional_tests, "x"), logical(1))))
  for (fit in full$diagnostics$raw) {
    expect_null(fit$working_candidates)
    expect_null(fit$working_record)
    expect_true(length(fit$per_subject_gradients) > 0)
    expect_true(length(fit$lrt$x$constrained) > 0)
  }
  expect_false(is.null(attr(full$diagnostics$raw, "rng")))
  for (j in seq_len(nrow(full$fit_provenance))) {
    row <- full$fit_provenance[j, , drop = FALSE]
    native <- full$diagnostics$raw[[j]]
    expect_identical(row$raw_convergence, native$raw_fit$convergence)
    expect_identical(row$optimizer_control[[1]], native$diagnostics$optimizer_control)
    expect_identical(row$requested_control[[1]], list(iter.max = 150L, eval.max = 200L))
    expect_identical(row$optimizer_evaluations[[1]], native$diagnostics$optimizer_evaluations)
    expect_identical(row$gradient_quadratic_lr_bound,
                     native$diagnostics$gradient_quadratic_lr_bound)
  }
})

test_that("PGMM worker compacts ordinary and observed payloads before aggregation", {
  fn <- auto_pgmm_stream_function("scrope_pmg_internal")
  env <- new.env(parent = environment(fn))
  environment(fn) <- env
  pred <- cbind(Intercept = 1, x = c(-1, -1, 1, 1))
  fit <- list(theta = list(beta = c(.4, .6), subVar = .3),
    naive_cov = diag(c(.2, .4)), robust_cov = diag(c(.3, .5)),
    convergence = -50L, optimizer = "nlminb", raw_fit = list(convergence = 1L),
    loglik = -10, per_subject_gradients = matrix(1, 50000L, 3L),
    godambe = list(U = matrix(2, 50000L, 3L)),
    working_candidates = list(list(par = c(.4, .6, .3))))
  env$hl_prepare_data <- function(...) list(
    ctx = list(nb = 2L, nind = c(2L, 2L), k = 2L, intcol = 1L,
      pred = pred, fid = c(1L, 3L, 5L)), count = matrix(1L, 2, 4),
    gid = c(2L, 4L), gene_names = c("two", "four"),
    predictor_names = colnames(pred), design_sds = c(1, 1))
  env$call_posindy <- function(...) NULL
  env$fit_gene_pmg_unconstrained <- function(...) fit
  env$registerDoFuture <- function(...) invisible(NULL)
  env$plan <- function(...) invisible(NULL)
  env$foreach <- function(gene_idx) gene_idx
  returned <- list()
  env$`%dorng%` <- function(obj, expr) {
    block <- substitute(expr)
    eval_env <- new.env(parent = parent.frame())
    lapply(obj, function(idx) {
      eval_env$gene_idx <- idx
      payload <- eval(block, eval_env)
      returned[[length(returned) + 1L]] <<- payload
      payload
    })
  }
  for (observed in c(FALSE, TRUE)) {
    returned <- list()
    res <- fn(count = NULL, id = NULL, verbose = FALSE, keep_diagnostics = FALSE,
      .fit_observer = if (observed) auto_pgmm_stream_observer else NULL)
    expect_length(returned, 2L)
    for (payload in returned) {
      expect_lt(as.numeric(object.size(payload)), 15000)
      expect_null(payload$fit$per_subject_gradients)
      expect_null(payload$fit$godambe)
      expect_null(payload$fit$working_candidates)
      expect_identical(payload$fit$raw_fit, list(convergence = 1L))
      if (observed) {
        expect_identical(payload$working_record$donor_count, 50000L)
        expect_identical(payload$working_record$candidates, 1L)
      } else expect_null(payload$working_record)
    }
    expect_false(any(res$summary$optimizer_success))
    expect_true(all(is.na(res$summary$p_robust_x)))
    expect_identical(res$fit_provenance$gene_id, c(2L, 4L))
    expect_identical(res$fit_provenance$gene, c("two", "four"))
  }
})

test_that("PGMM compact size excludes donor and optional null-fit arrays", {
  full <- auto_pgmm_stream_run(TRUE)
  fit <- full$diagnostics$raw[[1]]
  compact <- auto_pgmm_stream_function("pmg_compact_gene_result")
  baseline <- compact(fit)
  donor_array <- matrix(seq_len(150000L), 50000L, 3L)
  fit$per_subject_gradients <- donor_array
  fit$godambe <- list(U = donor_array)
  fit$lrt$x$constrained <- list(donor = donor_array)
  fit$lrt$x$raw <- list(donor = donor_array)
  fit$lrt$x$score_test <- list(donor = donor_array)
  fit$lrt$x$wald_test <- list(donor = donor_array)
  expect_identical(compact(fit), baseline)
  expect_lt(as.numeric(object.size(baseline)), 15000)
  expect_true(all(c("gate", "ridge", "score_stat", "stat") %in% names(baseline$lrt$x)))
})

test_that("PGMM two-worker observation preserves optional output and gene RNG streams", {
  skip_on_cran()
  skip_if(Sys.getenv("SCROPE_RUN_AUTO_PGMM_PARALLEL") != "true",
          "Set SCROPE_RUN_AUTO_PGMM_PARALLEL=true for the installed PSOCK gate.")
  serial <- auto_pgmm_stream_run(FALSE)
  parallel <- auto_pgmm_stream_run(FALSE, ncore = 2L)
  expect_identical(parallel$summary, serial$summary)
  expect_identical(parallel$lrt_details, serial$lrt_details)
  expect_identical(parallel$fit_provenance, serial$fit_provenance)
  serial_records <- attr(serial, ".working_records")
  parallel_records <- attr(parallel, ".working_records")
  pids <- vapply(parallel_records, `[[`, integer(1), "worker_pid")
  expect_true(all(pids != Sys.getpid()))
  expect_length(unique(pids), 2L)
  without_pid <- function(records) lapply(records, function(x) {
    x$worker_pid <- NULL
    x
  })
  expect_identical(without_pid(parallel_records), without_pid(serial_records))
  expect_null(parallel$diagnostics)
})
