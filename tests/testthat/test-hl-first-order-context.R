hl_context_fixture <- function() {
  count <- matrix(
    c(2, 0, 1, 3, 1, 2),
    nrow = 1,
    dimnames = list("gene1", paste0("cell", seq_len(6)))
  )
  id <- rep(paste0("subject", seq_len(3)), each = 2)
  pred <- cbind(
    z = c(-1, 0, 1, 2, 3, 4),
    intercept_raw = rep(2, 6),
    w = c(0.5, -0.5, 1.5, -1.5, 2.5, -2.5)
  )

  getFromNamespace("hl_prepare_data", "scRoPE")(
    count = count,
    id = id,
    pred = pred,
    cpc = 0,
    mincp = 1,
    verbose = FALSE
  )
}

test_that("HL context retains an arbitrary raw intercept coordinate", {
  prep <- hl_context_fixture()
  tr <- prep$design_transform

  expect_equal(tr$intcol, 2L)
  expect_equal(tr$c_int, 2)
  expect_equal(tr$pred_original, prep$pred_original)
  expect_equal(tr$pred_original, prep$ctx$pred_original)
  expect_equal(tr$means, colMeans(tr$pred_original))
  expect_equal(tr$means, prep$ctx$design_means)
  expect_equal(tr$pred_internal[, tr$intcol], rep(1, nrow(tr$pred_internal)))
  expect_equal(tr$sds, prep$design_sds)
  expect_equal(tr$M_X, prep$ctx$M_X)
  expect_equal(tr$M_X_inv, prep$ctx$M_X_inv)
})

test_that("HL original and internal beta maps preserve fitted linear predictors", {
  prep <- hl_context_fixture()
  tr <- prep$design_transform
  beta_internal <- c(z = 0.7, intercept_raw = -1.1, w = 0.25)

  beta_report <- as.vector(tr$M_X %*% beta_internal)
  beta_roundtrip <- as.vector(tr$M_X_inv %*% beta_report)

  expect_equal(beta_roundtrip, as.vector(beta_internal), tolerance = 1e-14)
  expect_equal(
    as.vector(tr$pred_original %*% beta_report),
    as.vector(tr$pred_internal %*% beta_internal),
    tolerance = 1e-13
  )
})

test_that("HL covariance maps use the full intercept-recentering congruence", {
  prep <- hl_context_fixture()
  tr <- prep$design_transform
  V_internal <- matrix(
    c(
      0.40, 0.08, -0.03,
      0.08, 0.30, 0.04,
      -0.03, 0.04, 0.20
    ),
    nrow = 3,
    byrow = TRUE
  )

  V_report <- tr$M_X %*% V_internal %*% t(tr$M_X)
  V_roundtrip <- tr$M_X_inv %*% V_report %*% t(tr$M_X_inv)

  expect_equal(unname(V_roundtrip), V_internal, tolerance = 1e-13)
  expect_false(isTRUE(all.equal(
    V_report,
    diag(1 / c(tr$sds[1], 1, tr$sds[3])) %*%
      V_internal %*%
      diag(1 / c(tr$sds[1], 1, tr$sds[3]))
  )))
})

test_that("HL reporting-scale constraints map through C_report M_X", {
  prep <- hl_context_fixture()
  tr <- prep$design_transform
  beta_internal <- c(0.7, -1.1, 0.25)
  beta_report <- as.vector(tr$M_X %*% beta_internal)
  C_report <- rbind(
    c(1, 2, -1),
    c(0, 1, 0)
  )
  C_internal <- C_report %*% tr$M_X

  expect_equal(
    as.vector(C_internal %*% beta_internal),
    as.vector(C_report %*% beta_report),
    tolerance = 1e-14
  )
})

test_that("public HL defaults retain locked boundary and replacement arrays", {
  prepare <- getFromNamespace("hl_prepare_data", "scRoPE")
  d1 <- hl_fo_fixture("D1", fresh = TRUE)
  d2 <- hl_fo_fixture("D2", fresh = TRUE)
  fixtures <- list(
    Z1 = hl_fo_fixture("Z1", fresh = TRUE),
    D1_low_v_high_phi = d1$cases$low_v_high_phi,
    D1_low_v_low_phi = d1$cases$low_v_low_phi,
    Z2 = hl_fo_fixture("Z2", fresh = TRUE),
    D2_low_v_high_phi = d2$cases$low_v_high_phi,
    D2_low_v_low_phi = d2$cases$low_v_low_phi
  )

  for (fixture_name in names(fixtures)) {
    fixture <- fixtures[[fixture_name]]
    prepared <- prepare(
      count = fixture$count,
      id = fixture$subject_id,
      pred = fixture$X_raw,
      offset = fixture$pi,
      cpc = 0.005,
      mincp = 5,
      verbose = FALSE
    )
    expect_identical(prepared$gid, 1L, info = fixture_name)
    expect_identical(prepared$gene_names, fixture$id, info = fixture_name)
  }
})
