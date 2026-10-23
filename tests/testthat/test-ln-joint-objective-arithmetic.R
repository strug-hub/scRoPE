ln_joint_arithmetic_fixture <- function(y = c(3, 0, 1, 0, 0, 0, 0, 2, 0, 4)) {
  fx <- list(
    count = Matrix::Matrix(matrix(y, nrow = 1L), sparse = TRUE),
    id = rep(c("third", "first", "second"), c(3L, 4L, 3L)),
    pred = cbind(
      x = seq(-1, 1, length.out = 10),
      "(Intercept)" = 1,
      z = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 1)
    ),
    offset = exp(seq(-0.3, 0.2, length.out = 10))
  )
  ln_joint_prepare_fixture(fx)
}

ln_joint_arithmetic_legacy <- function(theta, prep) {
  ctx <- prep$ctx
  posv <- prep$posv
  ptmg_ll(
    theta, ctx$pred, ctx$offset, posv$Y, posv$n_onetwo, posv$ytwo,
    ctx$id, ctx$fid, ctx$cumsumy[1L, ], ctx$posind[[1L]], posv$posindy,
    ctx$nb, ctx$nind, ctx$k
  )
}

test_that("joint LN scalar keeps constants and sparse donor boundaries", {
  prep <- ln_joint_arithmetic_fixture()
  ctx <- prep$ctx
  posv <- prep$posv
  theta <- c(-0.2, 0.1, 0.3, 0.4, 1.3)
  map <- ln_joint_map(ctx$nb)
  theta_log <- c(theta[seq_len(ctx$nb)], log(tail(theta, 2L)))
  objective <- ln_joint_objective(theta, 1L, posv, ctx, map)
  donors <- ln_joint_donor_objective_values(theta_log, 1L, posv, ctx, map)

  expect_identical(posv$posindy[c(1L, length(posv$posindy))], c(0L, 9L))
  expect_equal(ctx$cumsumy[1L, ], c(4, 0, 6), tolerance = 0)
  expect_equal(objective, ln_joint_arithmetic_legacy(theta, prep),
               tolerance = 1e-12)
  expect_equal(sum(donors), objective, tolerance = 1e-13)
  expect_true(all(is.finite(donors)))
  for (j in seq_len(ctx$k)) {
    cells <- seq.int(ctx$fid[j], ctx$fid[j + 1L] - 1L)
    one <- ctx
    one$pred <- ctx$pred[cells, , drop = FALSE]
    one$offset <- ctx$offset[cells]
    one$id <- rep(1, length(cells))
    one$fid <- as.integer(c(1, length(cells) + 1L))
    one$cumsumy <- ctx$cumsumy[, j, drop = FALSE]
    one$posind <- list(which(one$cumsumy[1L, ] > 0))
    one$k <- 1L
    one$nind <- length(cells)
    one_posv <- call_posindy(Matrix::t(prep$count[, cells, drop = FALSE]),
                            0L, one$nind)
    expect_equal(
      donors[j], ln_joint_objective(theta, 1L, one_posv, one, map),
      tolerance = 1e-13, ignore_attr = TRUE
    )
  }

  # The requested gene's totals, not the first row, define the objective.
  ctx$cumsumy <- rbind(rep(0, ctx$k), ctx$cumsumy)
  ctx$posind <- c(list(integer(0)), ctx$posind)
  expect_equal(ln_joint_objective(theta, 2L, posv, ctx, map), objective,
               tolerance = 0)
})

test_that("joint LN scalar arithmetic resolves small genuine objective changes", {
  fx <- list(
    count = Matrix::Matrix(matrix(c(0, 1, 2, 4, 1, 0, 3, 2), nrow = 1L),
                           sparse = TRUE),
    id = rep(c("b", "a"), each = 4L),
    pred = matrix(1, nrow = 8L, ncol = 1L),
    offset = rep(1, 8L)
  )
  prep <- ln_joint_prepare_fixture(fx)
  map <- ln_joint_map(1L)
  q <- function(beta) {
    ln_joint_objective(c(beta, 1e-4, 10000), 1L, prep$posv, prep$ctx, map)
  }
  # Decimal evaluation at 80 digits, starting from exact binary doubles and
  # integer log1p gamma-ratio sums; these include the full working constants.
  expected <- c(7.02241581564605289983, 7.02241581537349957682)
  actual <- c(q(0.25), q(0.2500000001))
  expect_lte(max(abs(actual - expected)), 2e-12)
  expect_lte(abs(diff(actual) - (-2.7255332300771624e-10)), 1e-13)
  expect_lt(diff(actual), 0)
  increase <- q(0.24999) - actual[1L]
  expect_lte(abs(increase - 2.7255843442804422e-5), 1e-13)
  expect_gt(increase, 0)
})

test_that("joint LN scalar shares the unchanged derivative and coordinate paths", {
  prep <- ln_joint_arithmetic_fixture()
  ctx <- prep$ctx
  posv <- prep$posv
  beta <- c(-0.2, 0.1, 0.3)
  dispersion <- c(0.4, 1.3)
  maps <- list(
    ln_joint_map(ctx$nb),
    ln_joint_map(ctx$nb, matrix(c(1, -0.5, 0.25), nrow = 1L), 0.15),
    ln_joint_map(ctx$nb, diag(ctx$nb), beta)
  )
  for (map in maps) {
    zeta <- map$project(beta)
    theta <- c(zeta, dispersion)
    theta_log <- c(zeta, log(dispersion))
    full <- ln_joint_full_natural(theta, map)
    full_log <- c(full[seq_len(ctx$nb)], log(dispersion))
    state <- ln_joint_derivative_state(theta_log, 1L, posv, ctx, map)
    raw <- ptmg_ll_der_hes4(
      full_log, ctx$pred, ctx$offset, posv$Y, posv$n_onetwo, posv$ytwo,
      ctx$fid, ctx$cumsumy[1L, ], ctx$posind[[1L]], posv$posindy,
      ctx$nb, ctx$nind, ctx$k, correct_posindy = TRUE
    )
    expect_equal(state$raw$value,
                 ln_joint_objective(theta, 1L, posv, ctx, map),
                 tolerance = 1e-14)
    derivative_fields <- setdiff(names(raw), "value")
    expect_identical(state$raw[derivative_fields], raw[derivative_fields])
    expect_equal(ln_joint_objective(theta, 1L, posv, ctx, map),
                 ln_joint_objective(full, 1L, posv, ctx, maps[[1L]]),
                 tolerance = 0)
    expect_equal(sum(ln_joint_donor_objective_values(
      theta_log, 1L, posv, ctx, map
    )), state$raw$value, tolerance = 1e-14)
    expect_equal(state$q_gradient_reduced,
                 as.numeric(crossprod(map$derivative_map, raw$gradient)),
                 tolerance = 0)
  }

  permutation <- c(2L, 3L, 1L)
  permuted <- ctx
  permuted$pred <- ctx$pred[, permutation, drop = FALSE]
  expect_equal(
    ln_joint_objective(c(beta, dispersion), 1L, posv, ctx, maps[[1L]]),
    ln_joint_objective(c(beta[permutation], dispersion), 1L, posv,
                       permuted, maps[[1L]]),
    tolerance = 1e-14
  )
})

test_that("joint LN scalar remains finite across the default dispersion bounds", {
  for (y in list(rep(0, 10L), c(3, 0, 1, 0, 0, 0, 0, 2, 0, 4))) {
    prep <- ln_joint_arithmetic_fixture(y)
    map <- ln_joint_map(prep$ctx$nb)
    bounds <- expand.grid(sigma2 = c(1e-4, 10), gamma = c(1e-3, 10000))
    for (i in seq_len(nrow(bounds))) {
      dispersion <- unlist(bounds[i, ], use.names = FALSE)
      theta <- c(-0.2, 0.1, 0.3, dispersion)
      theta_log <- c(head(theta, -2L), log(dispersion))
      objective <- ln_joint_objective(theta, 1L, prep$posv, prep$ctx, map)
      state <- ln_joint_derivative_state(theta_log, 1L, prep$posv,
                                         prep$ctx, map)
      expect_true(is.finite(objective))
      expect_lte(abs(objective - state$raw$value),
                 2e-12 * (1 + abs(objective)))
      expect_lte(abs(sum(ln_joint_donor_objective_values(
        theta_log, 1L, prep$posv, prep$ctx, map
      )) - objective), 2e-12 * (1 + abs(objective)))
    }
  }
})

test_that("joint LN preserves existing prepared fractional-input arithmetic", {
  # Fractional counts are not promoted to a newly supported model. The current
  # extractor truncates cells but retains numeric donor totals; the repair must
  # neither silently change that preparation nor fail on its fractional totals.
  prep <- ln_joint_arithmetic_fixture(c(3.25, 0, 1.5, 0, 0, 0, 0, 2.75, 0, 4.25))
  expect_equal(prep$posv$Y, c(3L, 1L, 2L, 4L), tolerance = 0)
  expect_equal(prep$ctx$cumsumy[1L, ], c(4.75, 0, 7), tolerance = 0)
  theta <- c(-0.2, 0.1, 0.3, 0.4, 1.3)
  value <- ln_joint_objective(theta, 1L, prep$posv, prep$ctx,
                              ln_joint_map(prep$ctx$nb))
  expect_true(is.finite(value))
  expect_equal(value, ln_joint_arithmetic_legacy(theta, prep), tolerance = 1e-12)
})

test_that("joint LN handles donor totals exceeding the integer index range", {
  # Each cell remains within the existing extractor's integer representation.
  # A donor total is numeric and must not allocate one entry per count or wrap.
  prep <- ln_joint_arithmetic_fixture(c(1200000000, 1200000000, 0, rep(0, 7L)))
  expect_gt(prep$ctx$cumsumy[1L, 1L], .Machine$integer.max)
  theta <- c(-0.2, 0.1, 0.3, 0.4, 1.3)
  value <- ln_joint_objective(theta, 1L, prep$posv, prep$ctx,
                              ln_joint_map(prep$ctx$nb))
  legacy <- ln_joint_arithmetic_legacy(theta, prep)
  expect_true(is.finite(value))
  expect_lte(abs(value - legacy), 1e-3)
})

test_that("normalized gamma ratios agree across bounded arithmetic branches", {
  counts <- c(0, 1, 2, 31, 32, 33, 255, 256, 257, 1000)
  for (a in c(1e-3, 31.999999, 32, 32.000001, 10000)) {
    reference <- vapply(counts, function(n) {
      if (n < 2) 0 else sum(log1p(seq_len(n - 1L) / a))
    }, numeric(1))
    value <- ln_joint_rising_normalized_cpp(a, counts)
    expect_lte(max(abs(value - reference) / (1 + abs(reference))), 1e-13)
  }

  # 256-bit MPFR log-gamma differences from exact input doubles. The cases
  # straddle the argument shift and the small-ratio series boundary, and retain
  # fractional donor totals without rounding or a count-proportional vector.
  a <- c(1e-3, 31.999999, 32, 32.000001, 10000, 10000, 10000)
  n <- c(0.25, 4.75, 4.75, 4.75, 0.25, 1250, 1250.0001)
  reference <- c(
    -3.8964364172997525481, 0.26674764746158599638,
    0.26674763946564057671, 0.26674763146969566522,
    -9.3750781220701904377e-6, 75.000258690559746735,
    75.000270464419243594
  )
  value <- vapply(seq_along(a), function(i) {
    ln_joint_rising_normalized_cpp(a[i], n[i])
  }, numeric(1))
  expect_lte(max(abs(value - reference) / (1 + abs(reference))), 1e-13)
  expect_identical(ln_joint_rising_normalized_cpp(1, numeric(0)), numeric(0))
})
