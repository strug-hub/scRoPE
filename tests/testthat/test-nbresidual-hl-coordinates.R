nbresidual_hl_fixture <- function() {
  count <- matrix(c(1, 4, 2, 7, 5, 2, 8, 3, 6, 1, 4, 9), nrow = 3L,
    dimnames = list(c("gene_a", "filtered", "gene_c"), paste0("cell", 1:4)))
  pred <- cbind(x = c(2, 3, 5, 6), constant = 2, z = c(-2, 1, -1, 4))
  beta <- rbind(c(0.12, 0.2, -0.07), c(-0.08, 0.4, 0.05))
  fit <- list(
    summary = data.frame(gene_id = c(3L, 1L), gene = c("gene_c", "gene_a"),
      logFC_z = beta[, 3L], se_robust_x = c(0.1, 0.2),
      logFC_constant = beta[, 2L], logFC_x = beta[, 1L],
      check.names = FALSE),
    overdispersion = data.frame(Subject = c(0.15, 0.35), Cell = c(0.2, 0.4)),
    algorithm = c(2, 2),
    random_effect = data.frame(V1 = c(-0.2, 0.1), V2 = c(0.3, -0.15)),
    diagnostics = NULL,
    provenance = list(schema_version = 2L, route = "HL", requested_model = "HL",
      effective_model = "HL", estimator = "first_order_nested_hl")
  )
  list(fit = fit, count = count, pred = pred, id = c("donor_b", "donor_b", "donor_a", "donor_a"),
    offset = c(0.5, 2, 1.5, 3), beta = beta)
}

nbresidual_hl_reference <- function(fx, conditional = FALSE) {
  # Compute each gene separately from the original design and NB-lognormal
  # moments; no package prediction/centering helper is used for this oracle.
  result <- t(vapply(seq_len(nrow(fx$beta)), function(g) {
    eta <- fx$pred[, "x"] * fx$beta[g, 1L] +
      fx$pred[, "constant"] * fx$beta[g, 2L] +
      fx$pred[, "z"] * fx$beta[g, 3L] + log(fx$offset)
    subject <- fx$fit$overdispersion$Subject[g]
    cell <- fx$fit$overdispersion$Cell[g]
    if (conditional) {
      mode <- as.numeric(fx$fit$random_effect[g, ])
      mean <- exp(eta + mode[c(1L, 1L, 2L, 2L)])
      variance <- mean + cell * mean^2
    } else {
      mean <- exp(eta + subject / 2)
      variance <- mean + mean^2 * (exp(subject) * (1 + cell) - 1)
    }
    observed <- as.numeric(fx$count[fx$fit$summary$gene_id[g], ])
    (observed - mean) / sqrt(variance)
  }, numeric(ncol(fx$count))))
  dimnames(result) <- list(as.character(fx$fit$summary$gene_id), colnames(fx$count))
  result
}

test_that("first-order HL marginal residuals use original coordinates and offsets", {
  fx <- nbresidual_hl_fixture()
  actual <- nbresidual(fx$fit, fx$count, fx$id, fx$pred, fx$offset)
  expect_equal(actual$residuals, nbresidual_hl_reference(fx), tolerance = 1e-12)
  expect_identical(actual$gene, c("gene_c", "gene_a"))
  # Predictor name matching must remain correct when input columns move.
  permuted <- nbresidual(fx$fit, fx$count, fx$id,
    fx$pred[, c("z", "x", "constant")], fx$offset)
  expect_equal(permuted, actual, tolerance = 1e-12)
  fx$fit$overdispersion <- fx$fit$overdispersion[, c("Cell", "Subject")]
  expect_equal(nbresidual(fx$fit, fx$count, fx$id, fx$pred, fx$offset), actual,
    tolerance = 0)
  # Compact saved fits contain the required provenance without diagnostics.
  expect_null(fx$fit$diagnostics)
})

test_that("first-order HL conditional residuals use stored modes in donor order", {
  fx <- nbresidual_hl_fixture()
  actual <- nbresidual(fx$fit, fx$count, fx$id, fx$pred, fx$offset, conditional = TRUE)
  expect_equal(actual$residuals, nbresidual_hl_reference(fx, TRUE), tolerance = 1e-12)
  expect_identical(actual$gene, c("gene_c", "gene_a"))
  fx$fit$random_effect <- NULL
  expect_error(nbresidual(fx$fit, fx$count, fx$id, fx$pred, fx$offset,
    conditional = TRUE), "output_re=TRUE", fixed = TRUE)
})

test_that("HL filtered and reordered gene rows align by retained public gene names", {
  fx <- nbresidual_hl_fixture()
  expected <- nbresidual(fx$fit, fx$count, fx$id, fx$pred, fx$offset)
  reordered <- fx$count[c(1L, 3L), , drop = FALSE]
  actual <- nbresidual(fx$fit, reordered, fx$id, fx$pred, fx$offset)
  expect_equal(actual, expected, tolerance = 0)
  sparse <- methods::as(reordered, "dgCMatrix")
  expect_equal(nbresidual(fx$fit, sparse, fx$id, fx$pred, fx$offset), expected,
    tolerance = 0)
  conditional <- nbresidual(fx$fit, reordered, fx$id, fx$pred, fx$offset,
    conditional = TRUE)
  expect_equal(conditional$residuals, nbresidual_hl_reference(fx, TRUE),
    tolerance = 1e-12)

  # Named rows may move only when the mapping is unambiguous.
  expect_error(nbresidual(fx$fit, fx$count[1L, , drop = FALSE], fx$id,
    fx$pred, fx$offset), "missing genes")
  rownames(reordered) <- c("gene_a", "gene_a")
  expect_error(nbresidual(fx$fit, reordered, fx$id, fx$pred, fx$offset),
    "duplicated gene names")
  unnamed <- unname(fx$count[c(1L, 3L), , drop = FALSE])
  expect_error(nbresidual(fx$fit, unnamed, fx$id, fx$pred, fx$offset),
    "gene_id values exceed")
})

test_that("HL intercept-only residuals preserve constant-column magnitude", {
  fx <- nbresidual_hl_fixture()
  fx$fit$summary <- fx$fit$summary[, c("gene_id", "gene", "logFC_constant")]
  pred <- matrix(-2, nrow = 4L, ncol = 1L, dimnames = list(NULL, "constant"))
  actual <- nbresidual(fx$fit, fx$count, fx$id, pred, fx$offset)$residuals
  eta <- outer(fx$fit$summary$logFC_constant, rep(-2, 4L))
  mean <- exp(sweep(eta, 2L, log(fx$offset), "+") +
    fx$fit$overdispersion$Subject / 2)
  variance <- mean + mean^2 * (exp(fx$fit$overdispersion$Subject) *
    (1 + fx$fit$overdispersion$Cell) - 1)
  expected <- (fx$count[fx$fit$summary$gene_id, , drop = FALSE] - mean) / sqrt(variance)
  rownames(expected) <- as.character(fx$fit$summary$gene_id)
  expect_equal(actual, expected, tolerance = 1e-12)

  names(fx$fit$summary)[3L] <- "logFC_(Intercept)"
  unit <- matrix(1, nrow = 4L, ncol = 1L,
    dimnames = list(NULL, "(Intercept)"))
  expect_equal(nbresidual(fx$fit, fx$count, fx$id, offset = fx$offset),
    nbresidual(fx$fit, fx$count, fx$id, unit, fx$offset), tolerance = 0)
})

test_that("HL residuals reject inconsistent metadata and recycling-prone inputs", {
  fx <- nbresidual_hl_fixture()
  bad <- fx$fit
  bad$provenance$effective_model <- "LN"
  expect_error(nbresidual(bad, fx$count, fx$id, fx$pred, fx$offset),
    "HL route and effective model")
  bad <- fx$fit
  bad$provenance$estimator <- "unknown_hl_coordinates"
  expect_error(nbresidual(bad, fx$count, fx$id, fx$pred, fx$offset),
    "supported prediction coordinates")
  bad$provenance$estimator <- NULL
  expect_error(nbresidual(bad, fx$count, fx$id, fx$pred, fx$offset),
    "supported prediction coordinates")
  for (offset in list(1, c(1, 1, 0, 1), c(1, 1, Inf, 1))) {
    expect_error(nbresidual(fx$fit, fx$count, fx$id, fx$pred, offset),
      "one per cell")
  }
  expect_error(nbresidual(fx$fit, fx$count, fx$id, fx$pred[-1L, ], fx$offset),
    "one row per count column")
  expect_error(nbresidual(fx$fit, fx$count, fx$id[-1L], fx$pred, fx$offset),
    "one nonmissing subject ID")
  bad <- fx$fit
  bad$overdispersion <- bad$overdispersion[1L, , drop = FALSE]
  expect_error(nbresidual(bad, fx$count, fx$id, fx$pred, fx$offset),
    "aligned with each fitted gene")
  bad <- fx$fit
  bad$random_effect <- bad$random_effect[, 1L, drop = FALSE]
  expect_error(nbresidual(bad, fx$count, fx$id, fx$pred, fx$offset,
    conditional = TRUE), "original donor order")
})

test_that("failed HL gene rows stay unavailable without hiding other genes", {
  fx <- nbresidual_hl_fixture()
  expected <- nbresidual_hl_reference(fx)
  fx$fit$summary[2L, c("logFC_x", "logFC_constant", "logFC_z")] <- NA_real_
  fx$fit$overdispersion[2L, ] <- NA_real_
  actual <- nbresidual(fx$fit, fx$count, fx$id, fx$pred, fx$offset)$residuals
  expect_equal(actual[1L, ], expected[1L, ], tolerance = 1e-12)
  expect_true(all(is.na(actual[2L, ])))
})

test_that("JLN PGMM and metadata-free legacy residual coordinates stay centered", {
  fx <- nbresidual_hl_fixture()
  centered <- fx$pred
  centered[, c("x", "z")] <- scale(centered[, c("x", "z")], scale = FALSE)
  centered[, "constant"] <- 1
  for (route in c("legacy", "JLN", "PGMM")) {
    fit <- fx$fit
    fit$provenance <- NULL
    reference <- fx
    reference$pred <- centered
    if (route == "JLN") {
      fit$provenance <- list(estimator = "ln_joint_root",
        parameterization = list(beta_public = "centered_public_design",
          subject = "sigma2_natural", cell_public = "Cell_equals_1_over_gamma"),
        random_effect = "pml_mode_given_ln_joint_root")
    }
    if (route == "PGMM") {
      fit$algorithm <- rep("PGMM", 2L)
      fit$overdispersion <- fit$overdispersion[, "Subject", drop = FALSE]
      reference$fit$overdispersion$Cell <- 0
    }
    actual <- nbresidual(fit, fx$count, fx$id, fx$pred, fx$offset)$residuals
    expect_equal(actual, nbresidual_hl_reference(reference), tolerance = 1e-12,
      info = route)
    if (route == "PGMM") {
      expect_warning(fallback <- nbresidual(fit, fx$count, fx$id, fx$pred,
        fx$offset, conditional = TRUE), "does not support PGMM yet")
      expect_equal(fallback$residuals, actual, tolerance = 0)
    } else {
      expect_equal(nbresidual(fit, fx$count, fx$id, fx$pred, fx$offset,
        conditional = TRUE)$residuals, nbresidual_hl_reference(reference, TRUE),
        tolerance = 1e-12, info = route)
    }
  }
})
