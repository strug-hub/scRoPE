#!/usr/bin/env Rscript
# Validate the optional test observer in an isolated installed-package process.
# Synthetic failures below test capture, not the unresolved Linux failure cause.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Usage: run_hl_history_diagnostic_gate.R OUTPUT_RDS EXPECTED_PACKAGE_PATH')
output <- args[1]
pkg <- normalizePath(args[2],mustWork=TRUE)
.libPaths(unique(c(dirname(pkg),.libPaths())))
library(scRoPE,lib.loc=dirname(pkg))
stopifnot(identical(normalizePath(getNamespaceInfo(asNamespace('scRoPE'),'path')),pkg))
source(file.path(pkg,'tests','testthat','helper-hl-history-diagnostic.R'))
ns <- asNamespace('scRoPE')
getfun <- function(n)get(n,ns,inherits=FALSE)
bindings <- c('working_hl_history_audit','working_hl_profile','working_hl_refine')
original <- setNames(lapply(bindings,getfun),bindings)
seal <- getfun('.hl_score_implementation')()
backend <- getfun('.hl_score_execution_identity')()
strip_time <- function(x){x$elapsed_seconds<-NULL;x}
record_dir <- tempfile('hl-observer-validation-')
dir.create(record_dir)
on_exit_cleanup <- function() {
  for(n in bindings)assignInNamespace(n,original[[n]],ns='scRoPE')
  unlink(record_dir,recursive=TRUE)
}
run <- function() {
  on.exit(on_exit_cleanup(),add=TRUE)
  set.seed(90216);id<-rep(1:12,each=5);X<-cbind(1,rnorm(60))
  y<-rpois(60,exp(.3+.2*X[,2]+rep(rnorm(12,0,.4),each=5)))
  data<-getfun('hl_fo_make_data')(y,X,rep(0,60),id,1L)
  collector<-getfun('working_hl_history')()
  fit<-getfun('hl_fo_fit')(data,active=integer(0),fixed_natural=c(sigma2=.4,phi=7),
    keep_diagnostics=FALSE,.history=collector)
  history<-collector$snapshot()
  stopifnot(fit$point_success,fit$covariance_success)
  audit<-function()getfun('working_hl_history_audit')(history,fit,keep_details=FALSE)
  assert_restored<-function()stopifnot(all(vapply(bindings,function(n)
    identical(getfun(n),original[[n]]),logical(1))))
  expected<-audit()
  stopifnot(expected$success)
  local({
    restore<-.hl_history_diagnostic_install(file.path(record_dir,'real'))
    on.exit(restore())
    stopifnot(identical(seal,getfun('.hl_score_implementation')()),
      identical(backend,getfun('.hl_score_execution_identity')()),
      identical(strip_time(audit()),strip_time(expected)))
  })
  assert_restored()
  cat('PASS real audit and implementation/backend identities unchanged\n')
  failures<-list()
  for(kind in c('profile','refinement','resolution','write'))local({
    on.exit(for(n in bindings)assignInNamespace(n,original[[n]],ns='scRoPE'))
    if(kind%in%c('refinement','resolution')) {
      assignInNamespace('working_hl_profile',function(...) {
        ans<-original$working_hl_profile(...)
        ans$augmented$A<-ans$augmented$A*1.1
        ans
      },ns='scRoPE')
      assignInNamespace('working_hl_refine',function(...) {
        if(kind=='resolution')return(list(success=TRUE,state=fit,
          reason='coalesced_stalled_history_root',trace=list('observer validation only'),
          stationarity=getfun('hl_fo_outer_stationarity')(fit),
          coalescence=getfun('working_hl_history_compare')(fit,fit)))
        list(success=FALSE,reason='synthetic_refinement_failure',
          trace=list('observer validation only'))
      },ns='scRoPE')
    } else {
      assignInNamespace('working_hl_profile',function(fit,rho,active,...) {
        ans<-list(success=FALSE,status='diagnostic_synthetic_failure',
          reason='injected for observer validation only')
        stop('Synthetic cold profile failure')
      },ns='scRoPE')
    }
    expected<-audit()
    stopifnot(identical(expected$reason,if(kind=='resolution')'history_supported' else
      'history_required_replay_failed'))
    directory<-file.path(record_dir,kind)
    restore<-.hl_history_diagnostic_install(directory)
    on.exit(restore(),add=TRUE,after=FALSE)
    if(kind=='write')unlink(directory,recursive=TRUE)
    warnings<-character(0)
    printed<-capture.output(ans<-withCallingHandlers(audit(),warning=function(w) {
      warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
    }))
    stopifnot(identical(strip_time(ans),strip_time(expected)),!length(warnings))
    if(kind=='write') {
      stopifnot(any(grepl('HL_HISTORY_DIAGNOSTIC_WRITE_ERROR',printed)))
    } else {
      paths<-list.files(directory,pattern='\\.rds$',full.names=TRUE)
      stopifnot(length(paths)==1L)
      record<-readRDS(paths)
      stopifnot(identical(record$schema,1L),identical(record$history,history),
        identical(record$selected_fit,fit),identical(record$audit,ans),
        identical(record$environment$backend,backend),
        identical(record$audit$failure_witness,expected$failure_witness))
      if(kind=='profile')stopifnot(length(record$original_profile_errors)==1L,
        identical(record$original_profile_errors[[1]]$original_locals$ans$reason,
          'injected for observer validation only'))
      if(kind=='refinement')stopifnot(length(record$original_refinement_failures)==1L,
        identical(record$original_refinement_failures[[1]]$reason,'synthetic_refinement_failure'))
      if(kind=='resolution')stopifnot(length(record$original_refinement_resolutions)==1L,
        identical(record$original_refinement_resolutions[[1]]$reason,'coalesced_stalled_history_root'),
        record$original_refinement_resolutions[[1]]$coalescence$point)
      failures[[kind]]<<-list(original_reason=ans$reason,
        original_error=ans$failure_witness$observation$error,capture_verified=TRUE)
    }
    cat('PASS ',kind,' diagnostic preserves original result and warning state\n',sep='')
  })
  assert_restored()
  # The test entry point uses this same finally restoration when tests error.
  caught<-tryCatch(local({
    restore<-.hl_history_diagnostic_install(file.path(record_dir,'unwind'))
    tryCatch(stop('synthetic test runner error'),finally=restore())
  }),error=function(e)conditionMessage(e))
  stopifnot(identical(caught,'synthetic test runner error'))
  assert_restored()
  stopifnot(identical(seal,getfun('.hl_score_implementation')()),
    identical(backend,getfun('.hl_score_execution_identity')()))
  cat('PASS all bindings restored, including error unwinding; final identities unchanged\n')
  list(status='passed',purpose='diagnostic_behavior_gate_not_linux_cause_reproduction',
    R=R.version.string,platform=R.version$platform,failures=failures,
    real_audit_preserved=TRUE,write_failure_preserved=TRUE,
    error_unwind_restored=TRUE,implementation_and_backend_unchanged=TRUE)
}
result<-run()
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
saveRDS(result,output)
