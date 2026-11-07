.pdc <- function(...) getFromNamespace("pmg_sandwich_covariance", "scRoPE")(...)
.pdc_scaled_error <- function(actual, expected) {
  scale <- outer(sqrt(diag(expected)), sqrt(diag(expected)))
  max(abs(actual - expected) / scale)
}

test_that("PGMM singular meat is distinct from singular sensitivity", {
  H <- diag(c(2, 3, 4))
  U <- rbind(c(1,-1,2,-2), c(2,-2,4,-4), c(0,1,0,-1))
  obj <- list(observed_info=H, per_subject_gradients=U)
  expected <- tcrossprod(diag(c(.5,1/3,.25)) %*% U)
  got <- .pdc(obj, 2)
  expect_lt(.pdc_scaled_error(got$Var_full, expected), 1e-10)
  expect_equal(got$Var_beta_adjusted, expected[1:2,1:2], tolerance=1e-12)
  expect_true(all(diag(got$Var_beta_adjusted)>0))
  expect_equal(qr(got$Var_beta_adjusted)$rank, 1L)
  expect_error(.pdc(list(observed_info=matrix(0,3,3),per_subject_gradients=U),2))
  expect_error(.pdc(list(observed_info=H,per_subject_gradients=U[1:2,]),2), "incompatible")
  expect_error(.pdc(obj,4), "coefficient block")
  expect_error(.pdc(obj,0), "coefficient block")
  obj$per_subject_gradients[1,1] <- Inf
  expect_error(.pdc(obj,2), "Nonfinite")
})

test_that("agreeing native PGMM covariance objects are exactly preserved", {
  obj <- list(observed_info=matrix(c(3,1,1,2),2),
              per_subject_gradients=matrix(c(1,2,-1,3,2,-2),2))
  legacy <- getFromNamespace("compute_sandwich_variance2", "scRoPE")
  for (full in c(TRUE,FALSE))
    expect_identical(.pdc(obj,1,compute_full=full),legacy(obj,1,compute_full=full))
})

test_that("tiny variances use relative agreement independently of BLAS error", {
  H <- diag(c(1e8,2)); u <- c(1,-1,1,-1)
  U <- rbind(u,u+1e-7*c(1,1,-1,-1))
  obj <- list(observed_info=H,per_subject_gradients=U)
  expected <- tcrossprod(diag(c(1e-8,.5)) %*% U)
  actual <- .pdc(obj,1)
  expect_lt(abs(actual$Var_beta_adjusted[1,1]/4e-16-1),1.01e-10)
  expect_lt(.pdc_scaled_error(actual$Var_full,expected),1.01e-10)
  # Do not require the platform's legacy calculation to be inaccurate.
  # Inject a deterministic one-percent error to verify the acceptance policy.
  f <- getFromNamespace("pmg_sandwich_covariance", "scRoPE")
  env <- new.env(parent=environment(f)); environment(f) <- env
  env$compute_sandwich_variance2 <- function(...) list(
    Var_beta_adjusted=matrix(4e-16*1.01),Var_full=expected)
  expect_equal(f(obj,1)$Var_beta_adjusted[1,1]/4e-16,1,tolerance=1e-14)
})

test_that("native preservation requires every requested covariance block to agree", {
  obj <- list(observed_info=diag(c(2,3)),per_subject_gradients=diag(2))
  expected <- diag(c(1/4,1/9))
  f <- getFromNamespace("pmg_sandwich_covariance", "scRoPE")
  env <- new.env(parent=environment(f)); environment(f) <- env
  old <- list(Var_beta_adjusted=matrix(1/4),Var_full=expected,sentinel="legacy")
  env$compute_sandwich_variance2 <- function(...) old
  expect_identical(f(obj,1),old)
  old$Var_full[2,2] <- 10
  expect_identical(f(obj,1)$Var_full,expected)
  expect_identical(f(obj,1,compute_full=FALSE),old)
  old$Var_beta_adjusted <- matrix(1/4,2,2)
  expect_equal(f(obj,1)$Var_beta_adjusted,matrix(1/4))
  obj$per_subject_gradients <- matrix(1e200,2,2)
  expect_error(f(obj,1),"Nonfinite PGMM influence covariance")
})

test_that("PGMM direct covariance obeys parameter transformations and equation scaling", {
  H <- matrix(c(3,.2,.1,.4,2,.3,.1,.2,4),3)
  U <- rbind(c(1,-1,2,-2),c(2,-2,4,-4),c(0,1,0,-1))
  baseline <- tcrossprod(solve(H,U))
  maps <- list(diag(c(-1e3,.01,2)), diag(3)[,c(2,1,3)],
               matrix(c(1,.2,0,-.3,2,0,0,0,3),3))
  for (T in maps) {
    obj <- list(observed_info=t(T)%*%H%*%T,per_subject_gradients=t(T)%*%U)
    got <- .pdc(obj,2)$Var_full
    expect_lt(.pdc_scaled_error(T%*%got%*%t(T),baseline),1.01e-10)
  }
  for (scale in c(1e-30,1e30)) {
    got <- .pdc(list(observed_info=H*scale,per_subject_gradients=U*scale),2)$Var_full
    expect_lt(.pdc_scaled_error(got,baseline),1.01e-10)
  }
})

test_that("public PGMM covariance maps predictor permutations and nonunit intercepts", {
  n <- 16L; id <- rep(seq_len(n),each=4L)
  x <- 4+2.7*rep(seq(-1.5,1.5,length.out=n),each=4L)
  z <- -1+.8*rep(c(-1.2,-.4,.4,1.2),n)
  omega <- rep(qgamma((((seq_len(n)*7)%%17)+.5)/18,shape=2,rate=1.7),each=4L)
  mu <- exp(1+.23*(x-4)-.2*(z+1))*omega
  y <- matrix(qpois((((seq_along(id)*29)%%67)+.5)/68,mu),1L)
  args <- list(count=y,id=id,pred=cbind(x=x,Intercept=2,z=z),model="PGMM",
               cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE,recovery="none")
  base <- do.call(scrope,args)
  args$pred <- cbind(z=-2.3*z,Intercept=-2.5,x=x/4)
  mapped <- do.call(scrope,args)
  # PGMM reports a unit-intercept, centered-predictor coordinate system.
  for (nm in c("x","z","Intercept")) {
    factor <- switch(nm,x=4,z=-1/2.3,Intercept=1)
    expect_equal(mapped$summary[[paste0("logFC_",nm)]]/factor,
                 base$summary[[paste0("logFC_",nm)]],tolerance=2e-7)
    expect_equal(mapped$summary[[paste0("se_robust_",nm)]]/abs(factor),
                 base$summary[[paste0("se_robust_",nm)]],tolerance=2e-7)
    expect_equal(mapped$summary[[paste0("p_robust_",nm)]],
                 base$summary[[paste0("p_robust_",nm)]],tolerance=2e-7)
    expect_true(base$summary[[paste0("wald_success_",nm)]])
    expect_true(mapped$summary[[paste0("wald_success_",nm)]])
  }
})

test_that("singular PGMM joint contrasts are withheld while scalar marginals survive", {
  H <- diag(c(2,3,4)); U <- rbind(c(1,-1,2,-2),c(2,-2,4,-4),c(0,1,0,-1))
  V <- .pdc(list(observed_info=H,per_subject_gradients=U),2)$Var_beta_adjusted
  fit <- list(theta=list(beta=c(1,0)),robust_cov=V,convergence=1L,
              optimizer="nlminb",raw_fit=list(convergence=0L))
  wald <- getFromNamespace("pmg_robust_wald","scRoPE")
  for (scale in c(1e-12,1,1e12)) {
    joint <- wald(fit,diag(c(scale,1)),c(0,0))
    expect_false(joint$success)
    expect_identical(joint$failure_reason,"singular_cov")
    expect_true(is.na(joint$stat));expect_true(is.na(joint$ridge))
  }
  for (j in 1:2) {
    scalar <- wald(fit,diag(2)[j,,drop=FALSE],0)
    expect_true(scalar$success)
    expect_equal(scalar$stat,fit$theta$beta[j]^2/V[j,j],tolerance=1e-12)
  }
  fit$robust_cov <- diag(c(4e-16,1))
  expect_true(wald(fit,matrix(c(1,0),1),0)$success)
  # Do not apply a raw absolute rank tolerance to a very small valid variance.
  out <- wald(fit,diag(c(1e8,1)),c(0,0))
  expect_true(out$success)
  expect_equal(out$stat,1/4e-16,tolerance=1e-12)
})

test_that("PGMM adjusted tests cannot re-expose a rejected joint Wald as LRT fallback", {
  f <- getFromNamespace("pmg_adjusted_contrasts","scRoPE")
  env <- new.env(parent=environment(f));environment(f)<-env
  fit <- list(theta=list(beta=c(1,0)),robust_cov=matrix(1,2,2),convergence=1L,
    optimizer="nlminb",raw_fit=list(convergence=0L),diagnostics=list(),
    lrt_gate=list(interior=TRUE,nuisance_interior=TRUE))
  env$pmg_fit_constrained_pair <- function(...) list(unconstrained=fit,constrained=fit,
    diagnostics=list(raw_lr_final=0))
  env$pmg_assess_constrained_convergence <- function(...) list(accepted=TRUE)
  env$compute_adj_score <- function(...) list(stat=NA_real_,success=FALSE)
  for (method in c("wald", "wald_ridge", "lrt_raw")) {
    env$compute_profile_lrt <- function(...) list(method=method,wP=2,wP_adjusted=999,
      success=FALSE,adjusted_available=FALSE,fallback=TRUE,gate=list())
    out <- f(list(k=4),1L,NULL,list(list(L=diag(2),b=c(1,0))),
             unconstrained=fit,retry_constrained=FALSE)[[1]]
    expect_false(out$adj$wald$success);expect_false(out$adj$lr$success)
    expect_false(out$adj$lr$fallback);expect_identical(out$adj$lr$method,"lrt_fail")
    expect_identical(out$adj$lr$failure_reason,"other_singular_cov")
    expect_true(is.na(out$p_values$lr));expect_true(is.na(out$p_values$wald))
    expect_equal(out$adj$lr$wP,2)
  }
  env$compute_profile_lrt <- function(...) list(method="lrt",wP_adjusted=2,
    success=TRUE,adjusted_available=TRUE,fallback=FALSE,gate=list())
  out <- f(list(k=4),1L,NULL,list(list(L=diag(2),b=c(1,0))),
           unconstrained=fit,retry_constrained=FALSE)[[1]]
  expect_true(out$adj$lr$success)
  expect_equal(out$p_values$lr,pchisq(2,2,lower.tail=FALSE))
})
