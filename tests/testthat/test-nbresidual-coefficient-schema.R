nbresidual_schema_fixture <- function() {
  count <- matrix(
    c(1, 3, 2, 5),
    nrow = 1L,
    dimnames = list("gene_one", paste0("cell", seq_len(4L)))
  )
  pred <- cbind(
    "(Intercept)" = 1,
    x = c(-1, 0, 1, 2)
  )
  summary <- data.frame(
    "logFC_(Intercept)" = log(2),
    "se_(Intercept)" = NA_real_,
    "p_naive_(Intercept)" = NA_real_,
    "logFC_x" = 0.2,
    gene_id = 1L,
    check.names = FALSE
  )
  fit <- list(
    summary = summary,
    overdispersion = data.frame(Subject = 0.2, Cell = 0.3),
    convergence = 1L,
    algorithm = "HL",
    random_effect = NULL
  )
  list(count = count, pred = pred, id = c(1, 1, 2, 2), fit = fit)
}

test_that("nbresidual selects interleaved coefficients by exact design names", {
  fx <- nbresidual_schema_fixture()
  observed <- nbresidual(
    fx$fit,
    count = fx$count,
    id = fx$id,
    pred = fx$pred
  )$residuals

  centered_pred <- fx$pred
  centered_pred[, "x"] <- centered_pred[, "x"] - mean(centered_pred[, "x"])
  beta <- c(log(2), 0.2)
  subject_variance <- 0.2
  cell_variance <- 0.3
  expected_mean <- exp(
    centered_pred %*% beta + subject_variance / 2
  )
  expected_variance <- expected_mean + expected_mean^2 * (
    exp(subject_variance) * (1 + cell_variance) - 1
  )
  expected <- matrix(
    (as.numeric(fx$count) - expected_mean) / sqrt(expected_variance),
    nrow = 1L,
    dimnames = list("1", colnames(fx$count))
  )

  expect_equal(observed, expected, tolerance = 1e-12)
  expect_true(all(is.finite(observed)))
})

test_that("nbresidual preserves the unnamed legacy coefficient block", {
  fx <- nbresidual_schema_fixture()
  names(fx$fit$summary)[c(1L, 4L)] <- c("beta_intercept", "beta_x")
  fx$fit$summary <- fx$fit$summary[, c(
    "beta_intercept", "beta_x", "se_(Intercept)",
    "p_naive_(Intercept)", "gene_id"
  )]

  legacy <- nbresidual(
    fx$fit,
    count = fx$count,
    id = fx$id,
    pred = fx$pred
  )$residuals

  modern <- nbresidual(
    nbresidual_schema_fixture()$fit,
    count = fx$count,
    id = fx$id,
    pred = fx$pred
  )$residuals
  expect_equal(legacy, modern, tolerance = 0)
})

test_that("nbresidual rejects missing and duplicate named coefficients", {
  fx <- nbresidual_schema_fixture()
  missing <- fx$fit
  missing$summary$logFC_x <- NULL
  expect_error(
    nbresidual(missing, fx$count, fx$id, pred = fx$pred),
    "missing: logFC_x",
    fixed = TRUE
  )

  duplicate <- fx$fit
  duplicate$summary[["duplicate_intercept"]] <- duplicate$summary[[
    "logFC_(Intercept)"
  ]]
  names(duplicate$summary)[ncol(duplicate$summary)] <- "logFC_(Intercept)"
  expect_error(
    nbresidual(duplicate, fx$count, fx$id, pred = fx$pred),
    "duplicate coefficient columns: logFC_(Intercept)",
    fixed = TRUE
  )
})
