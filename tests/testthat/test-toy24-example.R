toy24_helpers <- function() {
  directory <- system.file("examples", "toy24", package = "scRoPE")
  if (!nzchar(directory)) directory <- testthat::test_path("..", "..", "inst", "examples", "toy24")
  e <- new.env(parent = globalenv())
  sys.source(file.path(directory, "generate.R"), e)
  sys.source(file.path(directory, "extract_results.R"), e)
  list(functions = e, directory = directory)
}

test_that("toy data regenerate exactly and preserve the declared design and RNG", {
  data(toy_data)
  e <- toy24_helpers()$functions
  original_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) original_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(original_kind))
    if (had_seed) assign(".Random.seed", original_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(31)
  caller_seed <- .Random.seed
  expect_identical(e$generate_toy24(), toy_data)
  expect_identical(.Random.seed, caller_seed)
  expect_identical(RNGkind(), c("L'Ecuyer-CMRG", "Inversion", "Rejection"))
  rm(".Random.seed", envir = .GlobalEnv)
  expect_identical(e$generate_toy24(), toy_data)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))

  expect_identical(dim(toy_data$count), c(24L, 8000L))
  expect_type(toy_data$count, "integer")
  expect_true(all(toy_data$count >= 0))
  expect_identical(rownames(toy_data$count), toy_data$gene_metadata$gene)
  expect_identical(colnames(toy_data$count), toy_data$cell_metadata$cell)
  expect_identical(toy_data$id, toy_data$cell_metadata$donor)
  expect_identical(toy_data$pred[, "x"], toy_data$cell_metadata$x)
  expect_identical(rle(toy_data$id)$lengths, rep(40L, 200L))
  expect_true(all(table(toy_data$id, toy_data$pred[, "x"]) == 20L))
  expect_true(all(toy_data$offset == 1))
  g <- toy_data$gene_metadata
  expect_true(all(table(g$scenario, g$truth == 0) == 6L))
  expect_equal(log(g$mean_plus / g$mean_minus), g$truth, tolerance = 1e-14)
  expect_equal(g$fold_change, exp(g$truth), tolerance = 1e-14)
  expect_identical(toy_data$seed, 20261021L)
  for (i in 2:24) expect_identical(toy_data$streams[[i]], parallel::nextRNGStream(toy_data$streams[[i - 1L]]))

  # Independent numerical integration checks the mean-preserving Normal factor.
  for (x in c(-0.5, 0.5)) {
    mean_factor <- integrate(function(b) exp(b*x - x*x/2) * dnorm(b), -10, 10)$value
    expect_equal(mean_factor, 1, tolerance = 1e-10)
  }
})

test_that("the installed quick start reproduces the complete reviewed inference", {
  data(toy_data)
  h <- toy24_helpers()
  e <- new.env(parent = globalenv())
  expect_output(expect_no_warning(sys.source(file.path(h$directory, "quick_start.R"), e)))
  r <- e$results
  expected <- read.csv(file.path(h$directory, "reference_results.csv"), check.names = FALSE)
  expect_identical(r$gene, toy_data$gene_metadata$gene)
  expect_identical(r$availability, rep("inference", 24L))
  expect_identical(r$reason, expected$reason)
  # Numerical solvers/BLAS can differ across platforms. 1e-6 in the coefficient
  # or SE is tiny compared with the smallest SE (~0.02); never test p<.05 as a
  # proxy for numerical agreement. Data themselves must be exactly identical.
  for (nm in c("estimate", "model_se", "robust_se", "model_p_zero", "robust_p_zero"))
    expect_lte(max(abs(r[[nm]] - expected[[nm]])), 1e-6)
  for (method in c("model", "robust")) {
    se <- r[[paste0(method, "_se")]]
    expect_equal(r[[paste0(method, "_p_zero")]], 2 * pnorm(-abs(r$estimate / se)), tolerance = 1e-12)
    expect_equal(r[[paste0(method, "_lower")]], r$estimate - qnorm(.975) * se)
    expect_equal(r[[paste0(method, "_upper")]], r$estimate + qnorm(.975) * se)
  }
  # Both uncertainty calculations use precisely the one reported coefficient.
  expect_identical(r$estimate, e$fit$summary$logFC_x[match(r$gene, e$fit$summary$gene)])
})

test_that("toy export retains unavailable genes and withholds unsupported uncertainty", {
  data(toy_data)
  e <- toy24_helpers()$functions
  fake <- list(summary = data.frame(gene = toy_data$gene_metadata$gene,
    logFC_x = rep(0, 24), se_x = rep(1, 24), se_robust_x = rep(2, 24),
    p_naive_x = rep(1, 24), p_robust_x = rep(1, 24)),
    availability = data.frame(gene = toy_data$gene_metadata$gene,
      availability = rep("inference", 24), reason_group = rep("certified", 24)))
  fake$availability$availability[1] <- "unavailable"
  fake$availability$reason_group[1] <- "example_rejection"
  fake$summary$se_x[2] <- NA_real_
  fake$summary <- fake$summary[-3, ]
  r <- e$extract_toy_results(toy_data, fake)
  expect_identical(r$gene, toy_data$gene_metadata$gene)
  expect_identical(r$reason[1:3], c("example_rejection", "invalid_or_missing_inference_fields", "missing_summary"))
  for (nm in c("model_se", "robust_se", "model_p_zero", "robust_p_zero", "model_lower", "robust_upper"))
    expect_true(all(is.na(r[[nm]][1:3])))
  expect_true(all(r$estimate[1:2] == 0))
  failed <- e$extract_toy_results(toy_data, NULL, "simulated failure")
  expect_identical(failed$gene, r$gene)
  expect_true(all(failed$availability == "fit_error"))
  expect_true(all(is.na(failed$robust_lower)))
  expect_true(all(failed$reason == "fit_error: simulated failure"))
})
