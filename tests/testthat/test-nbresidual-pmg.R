test_that("nbresidual works for top-level PGMM fits with one-column overdispersion", {
  fx <- pmg_interior_test_fixture()

  fit <- scrope(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    model = "PGMM",
    verbose = FALSE,
    cpc = 0,
    mincp = 1
  )

  res <- nbresidual(
    nebulaSand = fit,
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    conditional = FALSE
  )

  expect_equal(dim(res$residuals), dim(fx$count))
  expect_equal(rownames(res$residuals), as.character(fit$summary$gene_id))
  expect_equal(colnames(res$residuals), colnames(fx$count))
  expect_true(all(is.finite(res$residuals)))
  expect_equal(res$gene, rownames(fx$count))
})

test_that("nbresidual PGMM conditional request warns and falls back to marginal residuals", {
  fx <- pmg_interior_test_fixture()

  fit <- scrope(
    count = fx$count,
    id = fx$id,
    pred = fx$pred,
    model = "PGMM",
    verbose = FALSE,
    cpc = 0,
    mincp = 1
  )

  res <- NULL
  expect_warning(
    res <- nbresidual(
      nebulaSand = fit,
      count = fx$count,
      id = fx$id,
      pred = fx$pred,
      conditional = TRUE
    ),
    "does not support PGMM yet"
  )

  expect_equal(dim(res$residuals), dim(fx$count))
  expect_true(all(is.finite(res$residuals)))
})
