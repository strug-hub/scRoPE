.rsc_primary <- function(model, partial=FALSE) {
  s<-data.frame(gene_id=1L,gene="controlled",logFC_x=.2,se_x=.1,p_naive_x=.04,
    model_wald_stat_x=4,model_wald_success_x=TRUE,model_wald_failure_x=NA_character_,
    se_robust_x=if(partial)NA_real_ else .12,p_robust_x=if(partial)NA_real_ else .1,
    wald_stat_x=if(partial)NA_real_ else (.2/.12)^2,wald_success_x=!partial,
    wald_failure_x=if(partial)paste0(model,"_controlled_covariance_failure") else NA_character_,
    optimizer_success=TRUE,fit_success=TRUE,fit_failure=NA_character_)
  list(summary=s,overdispersion=data.frame(Subject=.2),convergence=-10L,
    provenance=list(estimator=model,per_gene=data.frame(gene_id=1L,inference_available=TRUE,
      failure_reason=NA_character_,enabled_failure_reason=NA_character_)))
}

test_that("native status uses the same additive contract without new acceptance checks", {
  testthat::local_mocked_bindings(working_recover=function(...)stop("unexpected recovery"),
    working_pmg_certificate=function(...)stop("unexpected certificate"),.package="scRoPE")
  for(model in c("PGMM","HL","LN"))for(partial in c(FALSE,TRUE)) {
    p<-.rsc_primary(model,partial)
    f<-scRoPE:::scrope_finish_recovery(p,model,"none",notify=FALSE)
    a<-f$availability;c<-f$coefficient_availability
    expect_identical(f$convergence,p$convergence)
    expect_identical(f$summary[names(p$summary)],p$summary)
    expect_identical(a$native_qualification,if(partial)"partial" else "complete")
    expect_identical(a$recovery_outcome,if(partial)"disabled" else "not_needed")
    expect_false(a$recovery_attempted);expect_identical(a$selected_source,"native")
    expect_identical(c$native_available,!partial);expect_true(c$native_model_available)
    expect_identical(c$native_failure_reason,p$summary$wald_failure_x)
    expect_identical(c$native_available,c$available)
    expect_identical(c$native_model_available,c$model_available)
    expect_match(a$availability_message,"Native qualification:.*recovery:.*selected source:")
    expect_null(f$recovery$primary);expect_null(attr(f,".native_reporting"))
    expect_true(all(f$test_status$native_success==f$test_status$success))
    if(partial)expect_identical(c$failure_stage,"covariance")
  }
})

test_that("failed recovery preserves compact native partial evidence and the exact error", {
  for(model in c("PGMM","HL")) {
    p<-.rsc_primary(model,TRUE)
    a<-list(native_selected=TRUE,native_regular=FALSE,point_success=TRUE,
      route="native_partial",reason="native_partial_recovery_incomplete",
      native_qualification=list(point=TRUE,robust_complete=FALSE,model_complete=TRUE,complete=FALSE),
      recovery_attempted=TRUE,recovery_outcome="failed",recovery_reason="recovery_evaluation_error",
      recovery_error="controlled optimizer exception",elapsed_seconds=0)
    f<-scRoPE:::scrope_apply_working(p,list('1'=a),model,FALSE)
    f<-scRoPE:::scrope_finish_recovery(f,model,"auto",notify=FALSE)
    expect_identical(f$summary[names(p$summary)],p$summary)
    expect_identical(f$availability$native_qualification,"partial")
    expect_true(f$availability$recovery_attempted)
    expect_identical(f$availability$recovery_outcome,"failed")
    expect_identical(f$availability$recovery_error,"controlled optimizer exception")
    expect_identical(f$availability$selected_source,"native")
    expect_identical(f$coefficient_availability$native_failure_reason,p$summary$wald_failure_x)
    expect_null(f$recovery$primary);expect_null(f$recovery$diagnostics)
    expect_null(attr(f,".native_reporting"))
    again<-scRoPE:::scrope_finish_recovery(f,model,"auto",notify=FALSE)
    expect_identical(again,f)
    condition<-NULL
    withCallingHandlers(scRoPE:::scrope_recovery_notice(f$availability),warning=function(w) {
      condition<<-w;invokeRestart("muffleWarning")
    })
    expect_equal(unname(condition$counts[c("recovered","recovery_attempted","recovery_failed","native_incomplete_before")]),c(0,1,1,1))
  }
})

test_that("normalizing a failure repeatedly preserves its code", {
  for(type in c("wald","lrt","score")) {
    x<-c(NA_character_,"",unlist(scRoPE:::failure_codes[[type]],use.names=FALSE),
      "controlled_unknown","other_controlled_unknown")
    once<-scRoPE:::normalise_failure(x,type)
    expect_identical(scRoPE:::normalise_failure(once,type),once)
    expect_identical(tail(once,2),rep("other_controlled_unknown",2))
    expect_true(all(is.na(head(once,2))))
  }
})

test_that("JLN recovery separates covariance from overflowing or otherwise unavailable Wald", {
  x<-readRDS(test_path("fixtures","jln-research.rds"))$zero$input
  base<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE))
  native<-base$recovery$primary;original<-base$recovery$diagnostics[[1]]$result
  for(full in c(FALSE,TRUE))for(case in c("overflow","missing_point","zero_variance","negative_variance")) {
    f<-original
    if(case=="overflow"){f$beta[2]<-1e200;f$Vbeta[2,2]<-1}
    if(case=="missing_point")f$beta[2]<-NA_real_
    if(case=="zero_variance")f$Vbeta[2,2]<-0
    if(case=="negative_variance")f$Vbeta[2,2]<- -1
    answer<-list(route=base$recovery$routes$route,result=f)
    if(!full)answer<-scRoPE:::ln_compact_recovery(answer)
    out<-scRoPE:::ln_recovery_apply(native,native$diagnostics$raw,list(pred=x$X),TRUE,full,
      recovery_results=list(answer),.notify=FALSE)
    out<-scRoPE:::scrope_finish_recovery(out,"LN","auto",keep_diagnostics=full,notify=FALSE)
    a<-out$coefficient_availability[2,]
    reason<-switch(case,overflow="nonfinite_wald_statistic",missing_point="nonfinite_coefficient_estimate",
      "mapped_variance_not_positive_finite")
    expect_false(a$available);expect_false(a$reference_available);expect_false(a$statistic_available)
    expect_true(is.na(out$summary$p_robust_x));expect_true(is.na(out$summary$wald_stat_x))
    expect_identical(a$failure_reason,reason)
    expect_identical(out$summary$wald_failure_x,reason)
    row<-out$test_status[out$test_status$test=="wald"&out$test_status$contrast=="x",]
    expect_false(row$success);expect_identical(row$failure_reason,reason)
    expect_identical(row$availability_stage,a$failure_stage)
    expect_identical(a$native_failure_reason,native$coefficient_availability$failure_reason[2])
    expect_identical(out$availability$native_convergence,native$convergence)
    expect_true(out$availability$recovery_attempted)
    expect_identical(out$availability$recovery_outcome,"selected")
    expect_identical(out$recovery$routes$route,base$recovery$routes$route)
    if(case=="overflow") {
      expect_true(a$point_available);expect_true(a$covariance_available)
      expect_equal(a$robust_se,1);expect_true(a$interval_available)
      expect_equal(out$summary$se_robust_x,1)
      expect_identical(a$failure_stage,"statistic")
    } else expect_true(is.na(a$robust_se))
    if(!full)expect_null(out$recovery$primary)
  }
})


test_that("reported native flags do not erase an existing geometry refusal", {
  p<-.rsc_primary("PGMM")
  records<-list('1'=list(native_qualification=list(point=FALSE,robust_complete=FALSE,
    model_complete=FALSE,complete=FALSE,reason="native_nuisance_limit")))
  s<-scRoPE:::scrope_native_snapshot(p,"PGMM",records)
  expect_identical(s$genes$native_qualification,"unavailable")
  expect_identical(s$genes$native_failure_reason,"native_nuisance_limit")
  expect_false(s$genes$native_point_qualified)
  expect_true(s$coefficients$native_available)
})

test_that("native non-convergence retains method reasons within the common status", {
  reasons<-c(PGMM="optimizer_not_converged",HL="inner_armijo_failure",LN="competitive_hard_stop_candidate")
  for(model in names(reasons)) {
    p<-.rsc_primary(model);s<-p$summary
    s$optimizer_success<-s$fit_success<-s$wald_success_x<-s$model_wald_success_x<-FALSE
    s$wald_failure_x<-s$model_wald_failure_x<-unname(reasons[model])
    s$se_x<-s$se_robust_x<-s$p_robust_x<-s$p_naive_x<-s$wald_stat_x<-s$model_wald_stat_x<-NA_real_
    p$summary<-s;p$provenance$per_gene$inference_available<-FALSE
    p$convergence<- -50L
    f<-scRoPE:::scrope_finish_recovery(p,model,"none",notify=FALSE)
    expect_identical(f$availability$status,"unavailable")
    expect_identical(f$availability$native_qualification,"unavailable")
    expect_identical(f$availability$native_failure_reason,unname(reasons[model]))
    expect_identical(f$availability$recovery_outcome,"disabled")
    expect_false(f$availability$recovery_attempted)
    expect_identical(f$coefficient_availability$failure_stage,"point_qualification")
    expect_identical(f$convergence,p$convergence)
    expect_identical(f$summary[names(s)],s)
  }
})

test_that("a genuine failed JLN recovery preserves the original native failure", {
  x<-readRDS(test_path("fixtures","jln-face-competitor.rds"))
  f<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE))
  a<-f$availability
  expect_true(a$recovery_attempted);expect_identical(a$recovery_outcome,"failed")
  expect_identical(a$selected_source,"unavailable")
  expect_identical(a$recovery_reason,"reduced_root_worse_than_known_full_candidate")
  expect_false(is.na(a$native_failure_reason));expect_false(any(f$coefficient_availability$available))
  expect_null(f$recovery$primary);expect_null(f$recovery$diagnostics);expect_null(f$diagnostics)
})

test_that("Armijo and non-converged fit reasons share the point qualification stage", {
  reasons<-c("inner_armijo_failure","nonconverged_unconstrained",
    "nonconverged_constrained","score_null_no_convergence")
  expect_identical(vapply(reasons,scRoPE:::scrope_failure_stage,character(1)),
    setNames(rep("point_qualification",length(reasons)),reasons))
  expect_identical(scRoPE:::scrope_failure_stage("singular_cov"),"covariance")
  expect_identical(scRoPE:::scrope_failure_stage("nonfinite_wald_statistic"),"statistic")
  expect_identical(scRoPE:::scrope_failure_stage("unrecognized_method_reason"),"inference")
  expect_true(is.na(scRoPE:::scrope_failure_stage("inner_armijo_failure",TRUE)))
})
