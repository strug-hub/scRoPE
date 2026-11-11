.apr_fixture <- function(combined=FALSE,face=TRUE) {
  old_kind<-RNGkind();old<-if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)).Random.seed else NULL
  on.exit({do.call(RNGkind,as.list(old_kind));if(is.null(old))rm(".Random.seed",envir=.GlobalEnv) else assign(".Random.seed",old,.GlobalEnv)})
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(if(combined)4003 else 90211)
  m<-if(combined)16L else 60L;k<-if(combined)6L else 8L
  id<-rep(seq_len(m),each=k);x<-rnorm(length(id));z<-rep(rnorm(m),each=k)
  y<-if(combined)rbinom(length(id),4,plogis(.1+.3*x+.15*z)) else
    rpois(length(id),exp(.3+.2*x+.1*z+rep(rnorm(m,0,.5),each=k)))
  X<-cbind(Intercept=1,x=x,region=0,z=z)
  if(face)list(y=c(y,0,0),id=c(id,m+1:2),X=rbind(X,c(1,0,1,0),c(1,0,1,0))) else
    list(y=y,id=id,X=X[,c(1,2,4),drop=FALSE])
}
.apr_call <- function(d,model,...)suppressWarnings(scrope(d$y,d$id,d$X,model=model,cpc=0,mincp=1,verbose=FALSE,...))
.apr_expect_native_refit <- function(actual,expected) {
  expect_identical(names(actual),names(expected))
  exact<-setdiff(names(expected),c("summary","overdispersion","coefficient_availability"))
  expect_identical(actual[exact],expected[exact])
  expect_recovery_refit_table(actual$summary,expected$summary)
  expect_recovery_refit_table(actual$overdispersion,expected$overdispersion,
    numerical=intersect(c("Subject","Cell"),names(expected$overdispersion)))
  expect_recovery_refit_table(actual$coefficient_availability,expected$coefficient_availability,
    numerical=c("beta","robust_se","statistic","p_value","model_se","model_statistic","model_p_value"))
}

test_that("automatic recovery and compact defaults have an explicit compatible opt-out", {
  expect_identical(formals(scrope)$keep_diagnostics,FALSE)
  expect_identical(formals(scrope_hl)$keep_diagnostics,FALSE)
  expect_identical(eval(formals(scrope)$recovery),c("auto","none"))
  d<-.apr_fixture(face=FALSE)
  expect_error(.apr_call(d,"LN",recovery="auto",ln_recovery="none"),"must agree")
  a<-.apr_call(d,"PGMM",recovery="none",pgmm_control=list(iter.max=600L))
  expect_identical(a$provenance$optimizer$main_control$iter.max,600L)
  expect_identical(a$recovery$mode,"none")
  expect_null(a$diagnostics);expect_null(a$recovery$primary)
  expect_true(all(a$availability$status=="available"))
  expect_true(all(a$fit_provenance$optimizer_success))
})

test_that("public native results and compact/detail inference are preserved", {
  d<-.apr_fixture(face=FALSE)
  for(model in c("HL","PGMM")) {
    native<-.apr_call(d,model,recovery="none")
    compact<-.apr_call(d,model)
    detailed<-.apr_call(d,model,keep_diagnostics=TRUE)
    # Native output copied within this call is still preserved bit for bit.
    expect_null(compact$recovery$primary)
    expect_identical(detailed$summary,detailed$recovery$primary$summary)
    expect_identical(detailed$overdispersion,detailed$recovery$primary$overdispersion)
    expect_recovery_refit_table(compact$summary,native$summary)
    expect_recovery_refit_table(compact$summary,detailed$summary)
    expect_identical(compact$availability,detailed$availability)
    expect_recovery_refit_table(compact$coefficient_availability,detailed$coefficient_availability,
      numerical=c("beta","robust_se","statistic","p_value","model_se","model_statistic","model_p_value"))
    expect_null(compact$diagnostics);expect_null(compact$recovery$diagnostics)
    expect_true(all(compact$coefficient_availability$available))
    expect_true(all(compact$availability$selected_source=="native"))
  }
})

test_that("public nuisance limits retain points and withhold unsupported inference", {
  d<-.apr_fixture()
  for(model in c("HL","PGMM")) {
    compact<-.apr_call(d,model)
    detailed<-.apr_call(d,model,keep_diagnostics=TRUE)
    # Separate recovery fits may differ at roundoff scale across BLAS builds.
    # Every other column, including status and missingness, remains exact.
    expect_recovery_refit_table(compact$summary,detailed$summary,
      numerical=as.vector(outer(c("logFC_","se_robust_","p_robust_","wald_stat_","se_","p_naive_","model_wald_stat_"),
        c("Intercept","x","region","z"),paste0)))
    expect_recovery_refit_table(compact$coefficient_availability,detailed$coefficient_availability,
      numerical=c("beta","robust_se","statistic","p_value","model_se","model_statistic","model_p_value"))
    expect_true(all(compact$availability$selected_route=="identifiable_nuisance_limit"))
    expect_false(any(compact$coefficient_availability$selection_adjusted))
    if(model=="HL") {
      expect_true(all(is.na(compact$summary[,grep("^se_(?!robust)",names(compact$summary),perl=TRUE),drop=FALSE])))
      expect_identical(compact$availability$status,"point_only")
      expect_equal(sum(compact$coefficient_availability$point_available),3L)
      expect_false(any(compact$coefficient_availability$available))
      expect_true(all(is.na(compact$summary[,grep("^se_robust_|^p_robust_",names(compact$summary)),drop=FALSE])))
      h<-compact$recovery$audit[[1]]$history_audit
      expect_null(h$observations);expect_null(h$history_evidence)
      expect_equal(h$processed,h$retained)
      expect_true(length(detailed$recovery$audit[[1]]$history_audit$observations)>0L)
    } else {
      expect_equal(sum(compact$coefficient_availability$available),2L)
      expect_equal(sum(compact$coefficient_availability$model_available),2L)
      expect_true(all(is.na(detailed$diagnostics$loglik)))
      unavailable<-!detailed$coefficient_availability$covariance_available
      V<-detailed$diagnostics$godambe[[1]]$robust_cov
      expect_true(all(is.na(V[unavailable,,drop=FALSE])))
      expect_true(all(is.na(V[,unavailable,drop=FALSE])))
      expect_identical(detailed$diagnostics$godambe[[1]]$inferential_scope,"coefficient_wald_only")
    }
    expect_lt(as.numeric(object.size(compact)),250000)
  }
})

test_that("all-failed native HL requests retain explicit optional test statuses after recovery", {
  d<-.apr_fixture(combined=TRUE)
  f<-.apr_call(d,"HL",additional_tests=c("score","lrt"),lrt_details=TRUE,output_re=TRUE,keep_diagnostics=TRUE)
  expect_true(any(f$availability$selected_source=="recovery"))
  extra<-f$test_status[f$test_status$test%in%c("score","lrt"),,drop=FALSE]
  expect_true(nrow(extra)>0)
  expect_true(all(!extra$success))
  expect_true(all(extra$failure_reason=="not_supported_for_recovery"))
  expect_true(all(extra$availability_status=="unsupported"))
  expect_true(all(is.na(f$random_effect)))
  expect_identical(f$convergence,f$recovery$primary$convergence)
  original<-f$recovery$primary$summary
  coefficient_failures<-unlist(original[grep("^wald_failure_",names(original))],use.names=FALSE)
  expect_true(all(coefficient_failures==original$fit_failure))
  expect_identical(f$provenance$per_gene$fit_success,f$summary$fit_success)
  expect_true(all(f$coefficient_availability$available[f$coefficient_availability$estimable]))
})

test_that("custom boxes preserve native controls without an implicit working mask", {
  d<-.apr_fixture(face=FALSE)
  for(model in c("HL","PGMM")) {
    automatic<-.apr_call(d,model,min=c(2e-4,1e-4),max=c(9,999))
    native<-.apr_call(d,model,min=c(2e-4,1e-4),max=c(9,999),recovery="none")
    expect_recovery_refit_table(automatic$summary,native$summary)
    expect_false(automatic$recovery$eligible)
    expect_identical(automatic$recovery$skip_reason,"recovery_requires_default_box")
    expect_null(automatic$recovery$primary)
    expect_identical(automatic$provenance$recovery$controls$max,c(9,999))
  }
})

test_that("the aggregate recovery condition has a stable code and truthful counts", {
  a<-data.frame(requested_model="HL",selected_source=c("native","recovery","unavailable"),
    status=c("available","point_only","unavailable"),
    selected_wald_coefficients=c(3,0,0),total_coefficients=3)
  condition<-NULL
  withCallingHandlers(scRoPE:::scrope_recovery_notice(a),warning=function(w){condition<<-w;invokeRestart("muffleWarning")})
  expect_s3_class(condition,"scrope_recovery_warning")
  expect_identical(condition$code,"scrope_recovery_summary")
  expect_equal(unname(condition$counts[c("recovered","point_only","unavailable")]),c(1,1,1))
})

test_that("a missing reported SE cannot advertise an available interval", {
  x<-list(summary=data.frame(gene_id=1L,gene="g",logFC_x=1e200,
    se_robust_x=NA_real_,p_robust_x=NA_real_,wald_stat_x=NA_real_,
    wald_success_x=FALSE,wald_failure_x="nonfinite_wald_statistic"),
    convergence=-25L,provenance=list(),recovery=list(
      routes=data.frame(gene_id=1L,route="fixed_dispersion",inferential_scope="data_selected_working_target"),
      coefficients=data.frame(gene_id=1L,coefficient="x",point_available=TRUE,
        covariance_available=TRUE,estimable=TRUE,sampling_target="working_target")))
  x<-scRoPE:::scrope_finish_recovery(x,"HL","auto",notify=FALSE)
  expect_true(x$coefficient_availability$covariance_available)
  expect_false(x$coefficient_availability$interval_available)
  expect_false(x$coefficient_availability$available)
})
