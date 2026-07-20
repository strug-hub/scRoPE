test_that("first-order fixture registry matches the locked specification", {
  definitions <- hl_fo_fixture_definitions()

  expect_identical(names(definitions), hl_fo_fixture_ids())
  expect_identical(
    hl_fo_fixture_ids(),
    c("R1", "C1", "C2", "K0", "Q0", "Q1", "Z1", "D1", "W1", "X0", "N0", "B1", "S1")
  )
  expect_identical(
    vapply(definitions[c("R1", "C1", "C2", "K0", "Q0", "Q1", "Z1", "W1", "X0", "N0", "B1", "S1")],
      `[[`, integer(1), "seed"),
    c(R1 = 104729L, C1 = 104729L, C2 = 104759L, K0 = 104761L,
      Q0 = 104729L, Q1 = 104729L, Z1 = 104773L, W1 = 104789L,
      X0 = 104801L, N0 = 104803L, B1 = 104827L, S1 = 104831L)
  )
  expect_identical(as.integer(definitions$D1$seeds), 104779:104782)
  expect_identical(
    names(definitions$D1$cases),
    c("low_v_high_phi", "low_v_low_phi", "high_v_high_phi", "high_v_low_phi")
  )
})

test_that("fixture generation is deterministic and preserves caller RNG state", {
  original_kind <- RNGkind()
  had_original_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  original_seed <- if (had_original_seed) .Random.seed else NULL
  on.exit({
    do.call(RNGkind, as.list(original_kind))
    if (had_original_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(90210)
  seed_before <- .Random.seed
  kind_before <- RNGkind()

  first <- hl_fo_fixture("R1", fresh = TRUE)

  expect_identical(.Random.seed, seed_before)
  expect_identical(RNGkind(), kind_before)

  second <- hl_fo_fixture("R1", fresh = TRUE)
  expect_identical(first$y, second$y)
  expect_identical(first$X, second$X)
  expect_identical(first$log_pi, second$log_pi)
  expect_identical(first$omega, second$omega)
  expect_identical(first$mu, second$mu)
  expect_identical(
    first$rng,
    c(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  )
})

test_that("regular, constrained, fixed-dispersion, and edge fixtures share the declared arrays", {
  r1 <- hl_fo_fixture("R1", fresh = TRUE)
  c1 <- hl_fo_fixture("C1", fresh = TRUE)
  q0 <- hl_fo_fixture("Q0", fresh = TRUE)
  q1 <- hl_fo_fixture("Q1", fresh = TRUE)
  c2 <- hl_fo_fixture("C2", fresh = TRUE)
  k0 <- hl_fo_fixture("K0", fresh = TRUE)

  expect_identical(dim(r1$X), c(1200L, 3L))
  expect_identical(c(r1$m, r1$n_i, r1$n, r1$p), c(240L, 5L, 1200L, 3L))
  expect_identical(c2 |> with(c(m, n_i, n, p)), c(300L, 5L, 1500L, 3L))

  for (derived in list(c1, q0, q1)) {
    expect_identical(derived$y, r1$y)
    expect_identical(derived$X, r1$X)
    expect_identical(derived$pi, r1$pi)
    expect_identical(derived$omega, r1$omega)
  }
  expect_identical(k0$y, c2$y)
  expect_identical(k0$X, c2$X)
  expect_identical(k0$pi, c2$pi)
  expect_identical(k0$data_seed, 104759L)
  expect_identical(k0$seed, 104761L)

  expect_identical(length(r1$active_dispersion), 2L)
  expect_identical(length(q0$active_dispersion), 0L)
  expect_identical(length(q1$active_dispersion), 1L)
  expect_identical(q1$active_dispersion, "xi")
  expect_equal(q0$fixed_dispersion, c(sigma2 = 0.4, phi = 8))
  expect_equal(q1$fixed_dispersion, c(sigma2 = 0.4))

  expect_identical(hl_fo_num_rank(c1$constraint$C)$rank, 1L)
  expect_identical(hl_fo_num_rank(c2$constraint$C)$rank, 2L)
  expect_identical(hl_fo_num_rank(k0$constraint$C)$rank, 3L)
  expect_equal(unname(drop(c1$constraint$C %*% r1$truth$beta)), unname(c1$constraint$rhs))
  expect_equal(unname(drop(c2$constraint$C %*% c2$truth$beta)), unname(c2$constraint$rhs))
  expect_equal(unname(drop(k0$constraint$C %*% k0$truth$beta)), unname(k0$constraint$rhs))
})

test_that("locked generated arrays have reproducible dimensions, ranks, and count summaries", {
  ordinary_ids <- c("R1", "C1", "C2", "K0", "Q0", "Q1", "Z1", "W1", "X0", "N0", "B1", "S1")
  fixtures <- setNames(lapply(ordinary_ids, hl_fo_fixture, fresh = TRUE), ordinary_ids)

  expected_n <- c(
    R1 = 1200L, C1 = 1200L, C2 = 1500L, K0 = 1500L,
    Q0 = 1200L, Q1 = 1200L, Z1 = 1500L, W1 = 400L,
    X0 = 400L, N0 = 90L, B1 = 1500L, S1 = 30L
  )
  expected_p <- c(
    R1 = 3L, C1 = 3L, C2 = 3L, K0 = 3L, Q0 = 3L, Q1 = 3L,
    Z1 = 3L, W1 = 3L, X0 = 3L, N0 = 1L, B1 = 3L, S1 = 3L
  )
  expected_rank <- c(
    R1 = 3L, C1 = 3L, C2 = 3L, K0 = 3L, Q0 = 3L, Q1 = 3L,
    Z1 = 3L, W1 = 3L, X0 = 2L, N0 = 1L, B1 = 3L, S1 = 3L
  )
  expected_total <- c(
    R1 = 597, C1 = 597, C2 = 856, K0 = 856, Q0 = 597, Q1 = 597,
    Z1 = 57, W1 = 169, X0 = 134, N0 = 0, B1 = 558, S1 = 4
  )
  expected_nonzero <- c(
    R1 = 412L, C1 = 412L, C2 = 528L, K0 = 528L, Q0 = 412L, Q1 = 412L,
    Z1 = 56L, W1 = 115L, X0 = 108L, N0 = 0L, B1 = 458L, S1 = 4L
  )

  expect_identical(vapply(fixtures, `[[`, integer(1), "n"), expected_n)
  expect_identical(vapply(fixtures, `[[`, integer(1), "p"), expected_p)
  expect_identical(
    vapply(fixtures, function(fixture) hl_fo_num_rank(fixture$X)$rank, integer(1)),
    expected_rank
  )
  expect_equal(vapply(fixtures, function(fixture) sum(fixture$y), numeric(1)), expected_total)
  expect_identical(vapply(fixtures, function(fixture) sum(fixture$y > 0), integer(1)), expected_nonzero)

  expect_true(any(fixtures$Z1$y > 0))
  expect_identical(fixtures$Z1$requested_seed, 104773L)
  expect_identical(fixtures$Z1$actual_seed, 104773L)
  expect_true(all(fixtures$N0$y == 0))
  expect_true(all(fixtures$N0$pi == 1))
  expect_true(all(fixtures$N0$X == 1))
})

test_that("D1 uses the four ordered seed and dispersion combinations", {
  d1 <- hl_fo_fixture("D1", fresh = TRUE)

  expect_identical(length(d1$cases), 4L)
  expect_identical(unname(vapply(d1$cases, `[[`, integer(1), "seed")), 104779:104782)
  expect_equal(
    vapply(d1$cases, function(fixture) fixture$truth$sigma2, numeric(1)),
    c(low_v_high_phi = 0.02, low_v_low_phi = 0.02, high_v_high_phi = 2, high_v_low_phi = 2)
  )
  expect_equal(
    vapply(d1$cases, function(fixture) fixture$truth$phi, numeric(1)),
    c(low_v_high_phi = 100, low_v_low_phi = 0.3, high_v_high_phi = 100, high_v_low_phi = 0.3)
  )
  expect_identical(
    vapply(d1$cases, function(fixture) hl_fo_num_rank(fixture$X)$rank, integer(1)),
    setNames(rep(3L, 4L), names(d1$cases))
  )
  expect_equal(
    vapply(d1$cases, function(fixture) sum(fixture$y), numeric(1)),
    c(low_v_high_phi = 520, low_v_low_phi = 517, high_v_high_phi = 1028, high_v_low_phi = 2454)
  )
})

test_that("W1 records raw-display and package-analysis rank conventions explicitly", {
  w1 <- hl_fo_fixture("W1", fresh = TRUE)

  # The fixture-table columns are interpreted as raw package input for W1.
  # Direct-engine tests must use `X`, the package-centered/scaled analysis
  # matrix.  Under the global SVD rule the displayed raw matrix has rank two,
  # while the analysis matrix retains all three columns and can exercise the
  # intended weak-curvature gate.
  raw_rank <- hl_fo_num_rank(w1$X_raw)
  analysis_rank <- hl_fo_num_rank(w1$X)

  expect_identical(w1$display_coordinate, "raw_package_input")
  expect_identical(w1$design_coordinate, "package_centered_scaled_analysis")
  expect_identical(raw_rank$rank, 2L)
  expect_identical(analysis_rank$rank, 3L)
  expect_identical(w1$design_rank_raw$rank, 2L)
  expect_identical(w1$design_rank_analysis$rank, 3L)
  expect_lt(raw_rank$singular_values[[3L]], raw_rank$tolerance)
  expect_gt(analysis_rank$singular_values[[3L]], analysis_rank$tolerance)
  expect_equal(colMeans(w1$X[, -w1$intcol, drop = FALSE]), c(z = 0, w = 0), tolerance = 1e-14)
  expect_equal(
    sqrt(colMeans(w1$X[, -w1$intcol, drop = FALSE]^2)),
    c(z = 1, w = 1),
    tolerance = 1e-12
  )
  expect_equal(
    drop(w1$X_raw %*% w1$truth$beta_raw),
    drop(w1$X %*% w1$truth$beta),
    tolerance = 1e-12
  )
})

test_that("global scaled-discrepancy and SVD rank rules follow the specification", {
  a <- matrix(c(1, -2, 0, 4), nrow = 2L)
  b <- matrix(c(1.1, -1.5, 0.2, 4), nrow = 2L)
  expected <- max(abs(a - b) / (1 + pmax(abs(a), abs(b))))
  expect_equal(hl_fo_sdisc(a, b), expected)
  expect_equal(hl_fo_sdisc(numeric(0), numeric(0)), 0)
  expect_error(hl_fo_sdisc(1:4, matrix(1:4, nrow = 2L)), "equally shaped")

  X <- diag(c(4, 2, 1e-20))
  rank_info <- hl_fo_num_rank(X)
  expect_equal(rank_info$tolerance, 3 * .Machine$double.eps * 4)
  expect_identical(rank_info$rank, 2L)
  expect_identical(hl_fo_num_rank(matrix(0, 2L, 3L))$rank, 0L)
})

test_that("step-halving and Richardson traces retain the complete numerical evidence", {
  point <- c(0.4, -0.7)
  lower_order <- function(x) {
    c(first = exp(x[[1L]]) + x[[2L]]^3, second = x[[1L]] * x[[2L]])
  }
  analytic <- c(first = exp(point[[1L]]), second = point[[2L]])

  with_ad <- hl_fo_fd_trace(
    lower_order,
    point = point,
    coordinate = 1L,
    analytic = analytic,
    ad_reference = analytic
  )
  expected_h0 <- .Machine$double.eps^(1 / 3) * (1 + abs(point[[1L]]))
  expect_equal(with_ad$steps, expected_h0 / (2^(0:5)))
  expect_identical(length(with_ad$central), 6L)
  expect_identical(length(with_ad$richardson), 5L)
  expect_identical(with_ad$candidate_ells, 1:3)
  expect_true(with_ad$selected_ell %in% 1:3)
  expect_identical(with_ad$selector, "ad")
  expect_equal(
    with_ad$richardson[[1L]],
    (4 * with_ad$central[[2L]] - with_ad$central[[1L]]) / 3
  )
  expect_lt(with_ad$analytic_discrepancy, 1e-9)

  without_ad <- hl_fo_fd_trace(
    lower_order,
    point = point,
    coordinate = 1L,
    analytic = analytic
  )
  expect_identical(without_ad$candidate_ells, 1:2)
  expect_true(without_ad$selected_ell %in% 1:2)
  expect_identical(without_ad$selector, "finite_difference")
  expect_identical(
    without_ad$reference,
    without_ad$richardson[[without_ad$selected_ell + 2L]]
  )
  expect_lt(without_ad$analytic_discrepancy, 1e-8)
  expect_identical(length(without_ad$selection_discrepancies), 4L)
})
