# Synthetic derivative states isolate refinement decisions under a numerical
# score floor. Curvature, stationarity, full-state and covariance checks are real;
# only profile evaluation and scalar reconstruction are replaced below.
whr_get <- function(n)getFromNamespace(n,'scRoPE')
whr_target <- function(tiny_covariance=FALSE) {
  s<-list(data=list(n=20L,m=10L,p=1L),active=1L,beta_star=0,K=matrix(1,1,1),
    beta=0,gamma=0,eta=rep(0,10),rho=0,rho_full=c(0,0),F=0,g=0,
    H_F=matrix(-1,1,1),success=TRUE,status='ok',
    guard_bounds=list(lower=-5,upper=5),
    augmented=list(B=matrix(-1,1,1),A=diag(3),
      J=diag(if(tiny_covariance)rep(1e-12,3) else rep(1,3)),
      psi_subject=matrix(0,3,10),augmented_residual=rep(0,3),chi=0,L_gamma=0))
  whr_get('hl_fo_finalize_state')(s)
}
whr_profile <- function(target,covariance_mismatch=FALSE,better_probe=FALSE) {
  function(fit,rho,active=fit$active,history=NULL,source='synthetic_profile') {
    out<-target;out$rho_full<-rho;out$rho<-rho[active]
    out$g<-if(abs(out$rho)>1e-8)-out$rho else 2e-11
    out$F<-if(better_probe&&abs(out$rho)>0&&abs(out$rho)<1e-8)1 else -out$rho^2/2
    out$augmented$augmented_residual<-c(0,0,out$g)
    if(covariance_mismatch)out$augmented$J[1,1]<-2*out$augmented$J[1,1]
    if(!is.null(history)) {
      key<-history$register(out$data,active,out$beta_star,out$K,TRUE,rho)
      history$record(key,rho,out,TRUE,source)
    }
    out
  }
}

test_that('only historical coalescence can resolve a certified polishing stall', {
  target<-whr_target();profile<-whr_profile(target)
  state<-profile(target,c(0,0))
  local_mocked_bindings(working_hl_profile=profile,.package='scRoPE')
  refine<-whr_get('working_hl_refine')
  strict<-refine(state,target)
  expect_false(strict$success);expect_identical(strict$reason,'outer_refinement_stalled')
  history<-whr_get('working_hl_history')()
  recovered<-refine(state,target,history,coalesce_with_template=TRUE)
  expect_true(recovered$success)
  expect_identical(recovered$reason,'coalesced_stalled_history_root')
  expect_identical(recovered$state,state)
  expect_true(recovered$stationarity$success)
  expect_true(recovered$coalescence$point)
  expect_true(recovered$coalescence$covariance)
  expect_gt(max(abs(tail(recovered$trace,1)[[1]]$step)),1e-11)
  # All thirteen rejected line-search probes are retained, not dropped.
  expect_equal(history$snapshot()$observed,13L)
  expect_equal(history$snapshot()$retained,13L)
})

test_that('true stalls, invalid roots and untied targets remain failures', {
  target<-whr_target();base<-whr_profile(target)(target,c(0,0))
  cases<-list(
    different_beta=function(s){s$beta<-1e-4;s},
    different_modes=function(s){s$eta[1]<-1e-4;s},
    different_sensitivity=function(s){s$augmented$A[1,1]<-1.1;s},
    different_meat=function(s){s$augmented$J[1,1]<-1.1;s},
    nonstationary=function(s){s$g<-1e-3;s},
    not_tied=function(s){s$F<--1;s},
    better_scalar=function(s){s$F<-1;s})
  for(name in names(cases))local({
    state<-cases[[name]](base)
    local_mocked_bindings(working_hl_profile=function(...)state,.package='scRoPE')
    ans<-whr_get('working_hl_refine')(state,target,coalesce_with_template=TRUE)
    expect_false(ans$success,info=name)
    expect_identical(ans$reason,'outer_refinement_stalled',info=name)
  })
  local({
    bad<-base;bad$H_F[]<-1
    local_mocked_bindings(working_hl_profile=function(...)stop('must not evaluate'),.package='scRoPE')
    ans<-whr_get('working_hl_refine')(bad,target,coalesce_with_template=TRUE)
    expect_false(ans$success);expect_identical(ans$reason,'outer_refinement_curvature')
  })
  for(field in c('g','H_F'))local({
    bad_target<-target
    if(field=='g')bad_target$g<-1e-3 else bad_target$H_F[]<-1
    local_mocked_bindings(working_hl_profile=function(...)base,.package='scRoPE')
    ans<-whr_get('working_hl_refine')(base,bad_target,coalesce_with_template=TRUE)
    expect_false(ans$success);expect_identical(ans$reason,'outer_refinement_stalled')
  })
})

test_that('a failed required cold trial cannot be replaced by target agreement', {
  target<-whr_target();state<-whr_profile(target)(target,c(0,0))
  local_mocked_bindings(working_hl_profile=function(...)stop('cold trial failed'),.package='scRoPE')
  ans<-whr_get('working_hl_refine')(state,target,coalesce_with_template=TRUE)
  expect_false(ans$success);expect_identical(ans$reason,'outer_refinement_cold_failure')
  expect_match(ans$error,'cold trial failed')
})

test_that('resolved history stalls retain all probes and separate covariance support', {
  for(mismatch in c(FALSE,TRUE))local({
    target<-whr_target(tiny_covariance=TRUE)
    profile<-whr_profile(target,covariance_mismatch=mismatch)
    local_mocked_bindings(working_hl_profile=profile,
      working_hl_history_scalar=function(config,event)-event$rho[1]^2/2,.package='scRoPE')
    collector<-whr_get('working_hl_history')()
    profile(target,c(-1e-6,0),history=collector,source='initial_noncoalesced_root')
    before<-collector$snapshot()
    audit<-whr_get('working_hl_history_audit')(before,target)
    expect_true(audit$success)
    expect_identical(audit$covariance,!mismatch)
    expect_identical(audit$reason,if(mismatch)'history_tied_covariance_unresolved' else 'history_supported')
    expect_equal(audit$generated_evaluations,14L)
    expect_equal(audit$processed,15L);expect_equal(audit$retained,15L)
    expect_equal(length(audit$observations),15L)
    expect_length(audit$pending_event_ids,0L)
    expect_true(all(vapply(audit$observations[-1],function(e)
      identical(e$parent_event,1L)&&identical(e$source,'history_tied_refinement'),logical(1))))
    expect_identical(collector$snapshot(),before)
    if(mismatch) {
      expect_true(audit$observations[[1]]$replay$comparison$point)
      expect_false(audit$observations[[1]]$replay$comparison$covariance)
      expect_gt(audit$observations[[1]]$replay$comparison$covariance_agreement,1e-8)
    }
  })
})

test_that('a better rejected probe still vetoes an otherwise coalesced stalled root', {
  target<-whr_target();profile<-whr_profile(target,better_probe=TRUE)
  local_mocked_bindings(working_hl_profile=profile,
    working_hl_history_scalar=function(config,event)
      if(abs(event$rho[1])>0&&abs(event$rho[1])<1e-8)1 else -event$rho[1]^2/2,
    .package='scRoPE')
  collector<-whr_get('working_hl_history')()
  profile(target,c(-1e-6,0),history=collector,source='initial_noncoalesced_root')
  audit<-whr_get('working_hl_history_audit')(collector$snapshot(),target)
  expect_false(audit$success)
  expect_identical(audit$reason,'history_better_original_profile')
  expect_equal(audit$generated_evaluations,14L)
  expect_equal(audit$retained,15L)
  expect_identical(audit$failure_witness$observation$source,'history_tied_refinement')
  expect_gt(audit$pending_events,0L)
})
