.nor_native <- function() {
  x <- pmg_interior_test_fixture()
  x$native <- scRoPE:::fit_gene_pmg_unconstrained(x$gene_index,x$posv,x$prep$ctx,
    keep_candidates=TRUE)
  x
}
.nor_assess <- function(x, full=FALSE) scRoPE:::working_assess_native("PGMM",
  x$native,x$prep,x$gene_index,x$native$working_candidates,full)
.nor_call <- function(x, ...) suppressWarnings(scrope(x$count,x$id,x$pred,
  model="PGMM",cpc=0,mincp=1,verbose=FALSE,...))

test_that("complete native PGMM uses no recovery certificate, refit or precision trigger", {
  x <- .nor_native()
  # Neither Newton displacement nor a nuisance robust variance is a new
  # availability gate for already qualified coefficient inference.
  x$native$diagnostics$newton_step[] <- 1e20
  x$native$robust_cov_full[nrow(x$native$robust_cov_full),ncol(x$native$robust_cov_full)] <- 0
  testthat::local_mocked_bindings(
    working_pmg_certificate=function(...)stop("unexpected certificate"),
    working_recover=function(...)stop("unexpected recovery"),
    working_pmg_problem=function(...)stop("unexpected reference evaluator"),.package="scRoPE")
  for(full in c(FALSE,TRUE)) {
    a <- .nor_assess(x,full)
    expect_identical(a$route,"native")
    expect_true(a$native_regular)
    expect_true(a$model_inference_available)
    expect_false(a$recovery_attempted)
    expect_null(a$native_certificate)
    expect_null(a$audit$native_certificate)
  }
  a <- .nor_call(x)
  b <- .nor_call(x,recovery="none")
  expect_identical(a$summary,b$summary)
  expect_identical(a$overdispersion,b$overdispersion)
  expect_null(a$recovery$primary)
})

test_that("native qualification retains optimizer, geometry and curvature rules", {
  x <- .nor_native()
  q <- function(n=x$native,unresolved=FALSE)
    scRoPE:::working_native_qualification("PGMM",n,x$prep$ctx,unresolved)
  expect_true(q()$point)
  n <- x$native;n$raw_fit$convergence <- 1L
  expect_false(q(n)$point)
  n <- x$native;n$diagnostics$boundary_any <- TRUE
  expect_false(q(n)$point)
  n <- x$native;n$diagnostics$information_positive_definite <- FALSE
  expect_false(q(n)$point)
  n <- x$native;n$diagnostics$information_rcond <- 0
  expect_false(q(n)$point)
  expect_false(q(unresolved=TRUE)$point)
})

test_that("native partial survives failed and incomplete recovery without splicing", {
  x <- .nor_native();x$native$robust_cov[2,2] <- NA_real_
  for(reply in list(list(success=FALSE,point_success=FALSE,reason="trial_failed"),
    list(success=TRUE,point_success=TRUE,inference_available=FALSE,reason="covariance_not_certified",
      beta_report=c(90,91),V_report=diag(2),estimable=c(TRUE,TRUE)),
    list(success=TRUE,point_success=TRUE,inference_available=TRUE,model_inference_available=FALSE,
      reason="model_covariance_unavailable",beta_report=c(90,91),V_report=diag(2),estimable=c(TRUE,TRUE)))) {
    local({
      testthat::local_mocked_bindings(working_recover=function(...)reply,.package="scRoPE")
      a <- .nor_assess(x);full <- .nor_assess(x,TRUE)
      expect_identical(a$route,"native_partial")
      expect_true(a$native_selected);expect_false(a$native_regular)
      expect_true(a$recovery_attempted)
      expect_identical(a$recovery_reason,reply$reason)
      expect_null(a[["recovery_attempt",exact=TRUE]])
      expect_identical(full$recovery_attempt,reply)
      expect_null(a$beta_report)
    })
  }
  fitter <- scRoPE:::fit_gene_pmg_unconstrained
  testthat::local_mocked_bindings(fit_gene_pmg_unconstrained=function(...) {
    n<-fitter(...);n$robust_cov[2,2]<-NA_real_;n
  },working_recover=function(...)stop("controlled recovery failure"),.package="scRoPE")
  a <- .nor_call(x);n <- .nor_call(x,recovery="none")
  expect_identical(a$summary,n$summary)
  expect_identical(a$overdispersion,n$overdispersion)
  expect_identical(a$availability$selected_source,"native")
  expect_identical(a$availability$status,"partial")
  expect_true(a$summary$model_wald_success_x)
  expect_false(a$summary$wald_success_x)
  expect_null(a$recovery$primary)
})

test_that("complete recovered model and robust inference replaces a partial fit as a whole", {
  x <- .nor_native();fitter<-scRoPE:::fit_gene_pmg_unconstrained
  recovered <- list(success=TRUE,point_success=TRUE,inference_available=TRUE,
    route="certified_free_refit",reason="certified",beta_report=c(8,9),
    model_inference_available=TRUE,model_V_report=diag(c(.16,.64)),
    V_report=diag(c(.25,1)),estimable=c(TRUE,TRUE),par=c(8,9,.2),
    fixed=numeric(0),sampling_target="regular_working_equation_target")
  testthat::local_mocked_bindings(fit_gene_pmg_unconstrained=function(...) {
    n<-fitter(...);n$robust_cov[2,2]<-NA_real_;n
  },working_recover=function(...)recovered,.package="scRoPE")
  a <- .nor_call(x);b <- .nor_call(x,keep_diagnostics=TRUE)
  expect_identical(a$availability$selected_source,"recovery")
  expect_equal(a$summary$logFC_x,9)
  expect_equal(a$summary$se_robust_x,1)
  expect_equal(a$overdispersion$Subject,.2)
  expect_equal(a$summary$se_x,.8)
  expect_equal(a$summary$p_naive_x,pchisq((9/.8)^2,1,lower.tail=FALSE))
  expect_true(a$summary$model_wald_success_x)
  expect_null(a$recovery$primary)
  expect_true(is.list(b$recovery$primary))
  expect_false(b$recovery$primary$summary$wald_success_x)
  expect_false(identical(b$recovery$primary$summary$logFC_x,b$summary$logFC_x))
  expect_identical(a$summary,b$summary)
  expect_lt(as.numeric(object.size(a)),as.numeric(object.size(b)))
  expect_null(a$diagnostics);expect_null(a$recovery$diagnostics)
})

test_that("model-only missing inference retains native when recovery fails", {
  x <- .nor_native();x$native$naive_cov[2,2]<-NA_real_
  testthat::local_mocked_bindings(working_recover=function(...)list(success=FALSE,reason="controlled_failure"),.package="scRoPE")
  a <- .nor_assess(x)
  expect_identical(a$route,"native_partial")
  expect_identical(a$reason,"native_partial_recovery_incomplete")
  expect_true(a$native_selected);expect_true(a$recovery_attempted)
  expect_true(a$native_qualification$robust_complete)
  expect_false(a$native_qualification$model_complete)
})

test_that("HL reuses its resolved native geometry and recovers missing covariance", {
  x <- .nor_native();ctx<-x$prep$ctx
  data<-scRoPE:::hl_fo_make_data(as.numeric(x$prep$count),ctx$pred,ctx$offset,ctx$id,ctx$intcol)
  # Use the established regular fixture rather than assert this small fixture
  # necessarily has a fully estimable two-dispersion HL covariance.
  n<-list(point_success=TRUE,working_inference=list(model=c(TRUE,TRUE),robust=c(TRUE,TRUE)))
  testthat::local_mocked_bindings(hl_fo_inner_recession=function(...)stop("duplicate geometry"),
    working_recover=function(...)list(success=FALSE,reason="controlled_covariance_failure"),.package="scRoPE")
  a<-scRoPE:::working_assess_native("HL",n,x$prep,1L)
  expect_identical(a$route,"native");expect_false(a$recovery_attempted)
  expect_identical(a$geometry_source,"native_hl_qualification")
  n$working_inference$robust[]<-FALSE
  a<-scRoPE:::working_assess_native("HL",n,x$prep,1L)
  expect_identical(a$route,"native_partial");expect_true(a$recovery_attempted)
})

test_that("both interfaces distinguish model-only inference from point-only counts", {
  x<-.nor_native();fitter<-scRoPE:::fit_gene_pmg_unconstrained
  testthat::local_mocked_bindings(fit_gene_pmg_unconstrained=function(...) {
    n<-fitter(...);n$robust_cov[,]<-NA_real_;n
  },working_recover=function(...)list(success=FALSE,reason="failed"),.package="scRoPE")
  capture<-function(fun) {
    notices<-list()
    value<-withCallingHandlers(fun(),warning=function(w) {
      if(inherits(w,"scrope_recovery_warning")) {
        notices[[length(notices)+1L]]<<-w;invokeRestart("muffleWarning")
      }
    });list(value=value,notices=notices)
  }
  a<-capture(function()scrope(x$count,x$id,x$pred,model="PGMM",cpc=0,mincp=1,verbose=FALSE))
  b<-capture(function()scrope_working(x$count,x$id,x$pred,model="PGMM",cpc=0,mincp=1,verbose=FALSE))
  expect_length(a$notices,1);expect_length(b$notices,1)
  expect_identical(a$notices[[1]]$counts,b$notices[[1]]$counts)
  expect_equal(unname(a$notices[[1]]$counts["point_only"]),0)
  expect_equal(unname(a$notices[[1]]$counts["all_model_wald"]),1)
  expect_null(b$value$primary)
})

.nor_hl_data <- function() {
  kind<-RNGkind();on.exit(do.call(RNGkind,as.list(kind)))
  withr::with_preserve_seed({
    RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(90211)
    m<-60L;k<-8L;id<-rep(seq_len(m),each=k)
    x<-rnorm(length(id));z<-rep(rnorm(m),each=k)
    y<-rpois(length(id),exp(.3+.2*x+.1*z+rep(rnorm(m,0,.5),each=k)))
    list(y=y,id=id,X=cbind(Intercept=1,x=x,z=z))
  })
}

test_that("public native HL reuses recession qualification and preserves its result", {
  d<-.nor_hl_data();calls<-c(recession=0L,face=0L,recovery=0L,certificate=0L)
  recession<-scRoPE:::hl_fo_inner_recession;face<-scRoPE:::working_face_detect
  testthat::local_mocked_bindings(hl_fo_inner_recession=function(...) {
    calls['recession']<<-calls['recession']+1L;recession(...)
  },working_face_detect=function(...) {
    calls['face']<<-calls['face']+1L;face(...)
  },working_recover=function(...) {
    calls['recovery']<<-calls['recovery']+1L;stop('unexpected recovery')
  },working_hl_audit=function(...) {
    calls['certificate']<<-calls['certificate']+1L;stop('unexpected certificate')
  },.package='scRoPE')
  f<-scrope(d$y,d$id,d$X,model='HL',cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE)
  expect_identical(f$availability$selected_route,'native')
  expect_equal(unname(calls),c(1L,1L,0L,0L))
  expect_identical(f$summary,f$recovery$primary$summary)
  expect_identical(f$overdispersion,f$recovery$primary$overdispersion)
})

test_that("real HL recovery refusal preserves a native covariance-only partial fit", {
  d<-.nor_hl_data();fitter<-scRoPE:::hl_fo_fit;calls<-0L
  testthat::local_mocked_bindings(hl_fo_fit=function(...) {
    n<-fitter(...);calls<<-calls+1L
    # Perturb only the first native inference result. Recovery, derivatives,
    # history qualification and fit selection remain real.
    if(calls==1L) {
      expect_true(n$point_success);expect_true(n$covariance_success)
      n$covariance_success<-FALSE;n$covariance_status<-'controlled_native_covariance_failure'
    }
    n
  },.package='scRoPE')
  f<-suppressWarnings(scrope(d$y,d$id,d$X,model='HL',cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE))
  a<-f$recovery$diagnostics[[1]]
  expect_identical(f$availability$selected_route,'native_partial')
  expect_identical(f$availability$selected_source,'native')
  expect_true(a$recovery_attempted)
  expect_false(scRoPE:::working_recovery_complete(a$recovery_attempt))
  expect_identical(a$recovery_reason,a$recovery_attempt$reason)
  expect_identical(f$summary,f$recovery$primary$summary)
  expect_false(any(f$coefficient_availability$available))
  expect_true(all(f$coefficient_availability$point_available))
})


test_that("real PGMM recovery certifies a complete replacement for partial native inference", {
  x<-.nor_native();x$native$robust_cov[2,2]<-NA_real_
  a<-.nor_assess(x,TRUE)
  expect_identical(a$route,"certified_free_refit")
  expect_false(a$native_selected)
  expect_true(scRoPE:::working_recovery_complete(a))
  expect_true(all(vapply(a$certificates,function(z)isTRUE(z$point)&&isTRUE(z$covariance),logical(1))))
  expect_equal(a$V,tcrossprod(solve(a$A,a$U)),tolerance=1e-12)
  scale<-x$prep$design_sds;scale[x$prep$ctx$intcol]<-1
  expect_equal(a$beta_report,a$beta/scale,tolerance=1e-12)
  expect_equal(a$V_report,a$Vbeta/outer(scale,scale),tolerance=1e-12)
})

test_that("model-only native incompleteness reports consistently across recovery modes", {
  x<-.nor_native();fitter<-scRoPE:::fit_gene_pmg_unconstrained
  testthat::local_mocked_bindings(fit_gene_pmg_unconstrained=function(...) {
    n<-fitter(...);n$naive_cov[2,2]<-NA_real_;n
  },working_recover=function(...)list(success=FALSE,reason="controlled_model_recovery_failure"),
    working_pmg_certificate=function(...)stop("model-only reporting must not certify"),.package="scRoPE")
  capture<-function(mode) {
    notices<-list()
    fit<-withCallingHandlers(scrope(x$count,x$id,x$pred,model="PGMM",cpc=0,mincp=1,
      verbose=FALSE,recovery=mode),warning=function(w) {
      if(inherits(w,"scrope_recovery_warning")) {
        notices[[length(notices)+1L]]<<-w;invokeRestart("muffleWarning")
      }
    });list(fit=fit,notices=notices)
  }
  a<-capture("auto");n<-capture("none")
  expect_identical(a$fit$summary,n$fit$summary)
  expect_identical(a$fit$overdispersion,n$fit$overdispersion)
  expect_identical(a$fit$availability$status,"partial")
  expect_identical(n$fit$availability$status,"partial")
  expect_identical(n$fit$availability$failure_reason,"native_model_inference_unavailable")
  expect_true(all(n$fit$coefficient_availability$available))
  expect_true(n$fit$summary$wald_success_x)
  expect_false(n$fit$summary$model_wald_success_x)
  expect_length(a$notices,1);expect_length(n$notices,1)
  legacy_counts<-setdiff(names(n$notices[[1]]$counts),c("recovery_attempted","recovery_failed","recovery_incomplete"))
  expect_identical(a$notices[[1]]$counts[legacy_counts],n$notices[[1]]$counts[legacy_counts])
  expect_identical(unname(a$notices[[1]]$counts["recovery_attempted"]),1L)
  expect_identical(unname(n$notices[[1]]$counts["recovery_attempted"]),0L)
  expect_identical(a$fit$availability$recovery_outcome,"failed")
  expect_identical(n$fit$availability$recovery_outcome,"disabled")
  expect_identical(unname(n$notices[[1]]$counts['native_incomplete']),1L)
  expect_null(n$fit$recovery$primary)
  # Absent legacy flags and optional-test failures are not a new requirement.
  b<-n$fit;b$summary<-b$summary[,!grepl('^model_wald_',names(b$summary)),drop=FALSE]
  b$summary$score_success_x<-FALSE;b$summary$score_failure_x<-"controlled_optional_failure"
  b<-scRoPE:::scrope_finish_recovery(b,"PGMM","none",notify=FALSE)
  expect_identical(b$availability$status,"available")
  # The remaining model family prevents a false point-only label.
  b<-n$fit
  for(nm in grep('^wald_success_',names(b$summary),value=TRUE))b$summary[[nm]]<-FALSE
  b<-scRoPE:::scrope_finish_recovery(b,"PGMM","none",notify=FALSE)
  expect_identical(b$availability$status,"partial")
})


test_that("genuine model-only failure can recover both families at one selected fit", {
  x<-.nor_native();x$native$naive_cov[2,2]<-NA_real_
  a<-.nor_assess(x,TRUE)
  expect_identical(a$route,"certified_free_refit")
  expect_false(a$native_selected);expect_true(a$recovery_attempted)
  expect_true(scRoPE:::working_recovery_complete(a))
  expect_equal(a$model_covariance,solve(a$A),tolerance=1e-10)
  expect_equal(a$V,tcrossprod(solve(a$A,a$U)),tolerance=1e-10)
})
