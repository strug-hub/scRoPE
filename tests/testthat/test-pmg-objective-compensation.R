test_that("PGMM scalar retains accuracy with many small cell contributions", {
  # A deterministic, large-cell fixture exposes the serial addition error.
  # MPFR independently evaluates the original marginal likelihood in real
  # arithmetic; the binary64 inputs are converted without decimal rounding.
  skip_if_not_installed("Rmpfr")
  id <- rep(seq_len(200),each=400)
  x <- rep(rep(c(-.5,.5),each=200),200)
  y <- rep(c(0,1,0,2,0,1,0,0),10000)
  X <- cbind(Intercept=1,x=x)
  p <- scRoPE:::working_pmg_problem(y,X,id,rep(1,length(y)))
  M <- function(x)Rmpfr::mpfr(x,256)
  vals <- sort(unique(p$ctx$pred[,2]))
  for(v in c(1e-4,.02,.258,3)) {
    z <- c(-.7914,.009067,log(v)); st <- p$evaluate(z)
    beta <- M(st$par[1:2]); variance <- M(st$par[3])
    alpha <- 1/expm1(variance); rate <- alpha*exp(-variance/2)
    eta <- beta[1]+beta[2]*M(vals)
    total <- M(200)
    mu <- M(200)*sum(exp(eta))
    # The fixed pattern has 100 counts in each half of every donor.
    q <- M(200)*(-lgamma(total+alpha)+lgamma(alpha)-
      M(100)*sum(eta)-alpha*log(rate)+(total+alpha)*log(mu+rate))
    expect_lte(abs(as.numeric(M(st$native_q)-q)),st$scalar_tolerance)
    expect_lte(abs(st$q-st$native_q),st$scalar_tolerance)
  }
})

test_that("scalar identity failure remains a veto and has a specific diagnostic", {
  id <- rep(seq_len(30),each=5)
  totals <- rep(8:12,6)
  y <- unlist(lapply(totals,function(s)rep(s%/%5,5)+(seq_len(5)<=s%%5)))
  problem <- scRoPE:::working_pmg_problem(y,matrix(1,length(y),1),id,
    rep(1,length(y)),c(sigma2=1e-4))
  st <- problem$evaluate(log(mean(y))-1e-4/2)
  ok <- scRoPE:::working_pmg_certificate(st,problem)
  expect_true(ok$point && ok$covariance)
  st$native_q <- st$q+100*st$scalar_tolerance
  bad <- scRoPE:::working_pmg_certificate(st,problem)
  expect_false(bad$point)
  expect_false(bad$covariance)
  expect_identical(bad$reason,"free_scalar_identity")
  expect_identical(bad$metrics$failed_checks,"scalar_identity")
  st$native_q <- NA_real_
  nonfinite <- scRoPE:::working_pmg_certificate(st,problem)
  expect_false(nonfinite$point)
  expect_identical(nonfinite$reason,"free_scalar_identity")
})
