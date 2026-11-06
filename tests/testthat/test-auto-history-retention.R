# Synthetic history records exercise retention without fitting any model.
ahr_get <- function(name)getFromNamespace(name,"scRoPE")
ahr_fixture <- function(n=1L) {
  data<-list(n=4L,m=2L,p=2L)
  fit<-list(data=data,active=integer(0),beta_star=c(0,0),K=diag(2),
    beta=c(0,0),rho_full=c(0,0),F=0,point_success=TRUE)
  collector<-ahr_get("working_hl_history")()
  key<-collector$register(data,fit$active,fit$beta_star,fit$K,TRUE,fit$rho_full)
  for(i in seq_len(n))collector$record(key,fit$rho_full,
    list(success=TRUE,beta=fit$beta,F=fit$F,status="synthetic"),TRUE,"synthetic")
  list(fit=fit,history=collector$snapshot(),collector=collector)
}
ahr_mocks <- function(.env=parent.frame()) {
  local_mocked_bindings(
    working_hl_history_scalar=function(config,event)sum(event$beta),
    working_hl_profile=function(fit,...)fit,
    hl_fo_outer_stationarity=function(...)list(maximum=TRUE,scaled_score=0),
    working_hl_history_compare=function(state,fit)list(point=TRUE,
      covariance=TRUE,agreement=0,covariance_agreement=0,F=state$F),
    .package="scRoPE",.env=.env)
}
ahr_pair <- function(history,fit,...) {
  audit<-ahr_get("working_hl_history_audit")
  full<-audit(history,fit,...,keep_details=TRUE)
  compact<-audit(history,fit,...,keep_details=FALSE)
  fields<-setdiff(names(full),c("elapsed_seconds","observations","history_evidence",
    "pending_event_ids","details_retained"))
  expect_identical(compact[fields],full[fields])
  expect_null(compact$observations)
  expect_null(compact$pending_event_ids)
  expect_null(compact$history_evidence$events)
  expect_false(compact$details_retained)
  expect_equal(sum(compact$classifications),compact$processed)
  expect_equal(compact$processed+compact$pending_events,length(full$observations)+
    length(full$pending_event_ids))
  list(full=full,compact=compact)
}

test_that("compact history audit preserves gates without retaining processed rows", {
  ahr_mocks();d<-ahr_fixture(80L);before<-d$history
  last<-d$history$events[[80L]]
  last$F<-NA_real_;last$beta<-NULL;last$success<-FALSE
  d$history$events[[80L]]<-last
  a<-ahr_pair(d$history,d$fit)
  expect_true(a$compact$success)
  expect_identical(a$compact$reason,"history_supported")
  expect_identical(a$compact$processed,80L)
  expect_identical(a$compact$unique_replays,1L)
  expect_identical(unname(a$compact$classifications),c(1L,79L,0L,0L))
  expect_identical(a$compact$maxima[["scalar_discrepancy"]],0)
  expect_identical(a$compact$maxima[["cold_profile_discrepancy"]],0)
  expect_length(a$full$observations,80L)
  expect_identical(a$full$observations[[1L]]$id,1L)
  expect_identical(a$full$observations[[80L]]$classification,"failed_probe_without_profile")
  expect_lt(as.numeric(object.size(a$compact)),as.numeric(object.size(a$full))/3)
  # The audit does not mutate the caller's input or discard mandatory events.
  expect_identical(d$history$events[1:79],before$events[1:79])
  expect_identical(d$history$retained,80L)
  default<-ahr_get("working_hl_history_audit")(d$history,d$fit)
  expect_identical(default$observations,a$full$observations)
})

test_that("compact early failures keep truthful counts and one failure witness", {
  ahr_mocks();d<-ahr_fixture(3L)
  mismatch<-d$history;mismatch$events[[2L]]$F<-1
  a<-ahr_pair(mismatch,d$fit)
  expect_identical(a$compact$reason,"history_scalar_mismatch")
  expect_identical(a$compact$processed,2L)
  expect_identical(a$compact$pending_events,1L)
  expect_identical(a$compact$failure_witness$observation$id,2L)
  expect_identical(a$compact$maxima[["scalar_discrepancy"]],1)
  expect_length(a$full$history_evidence$events,3L)
  expect_identical(a$compact$history_evidence$event_id_bounds,list(first=1L,last=3L))
  better<-d$history;better$events[[2L]]$F<-1;better$events[[2L]]$beta<-c(1,0)
  a<-ahr_pair(better,d$fit)
  expect_identical(a$compact$reason,"history_better_original_profile")
  expect_identical(a$compact$maxima[["comparable_original_F_gap"]],1)
  overflow<-d$history;overflow$overflow<-TRUE;overflow$observed<-4L
  a<-ahr_pair(overflow,d$fit)
  expect_identical(a$compact$reason,"history_collection_incomplete")
  expect_identical(a$compact$processed,0L)
  expect_identical(a$compact$pending_events,3L)
  expect_identical(a$compact$dropped,1L)
  expect_null(a$compact$failure_witness$observation)
  a<-ahr_pair(d$history,d$fit,max_replays=0L)
  expect_identical(a$compact$reason,"history_replay_budget")
  expect_identical(a$compact$processed,1L)
  expect_identical(a$compact$pending_events,2L)
  local_mocked_bindings(working_hl_profile=function(...)stop("synthetic replay failure"),
    .package="scRoPE")
  a<-ahr_pair(d$history,d$fit)
  expect_identical(a$compact$reason,"history_required_replay_failed")
  expect_match(a$compact$failure_witness$observation$error,"synthetic replay failure")
})

test_that("compact projected and tied-root ledgers retain the same decisions", {
  ahr_mocks();d<-ahr_fixture()
  h<-d$history;config<-h$configurations[[1L]]
  config$held_log<-c(1,0);key<-ahr_get(".hl_score_fingerprint")(config)
  h$configurations<-setNames(list(config),key)
  h$events[[1L]]$configuration<-key;h$events[[1L]]$rho<-c(1,0)
  local_mocked_bindings(working_hl_profile=function(fit,...) {
    fit$F<-1;fit
  },.package="scRoPE")
  a<-ahr_pair(h,d$fit)
  expect_identical(a$compact$reason,"history_better_profile")
  expect_identical(a$compact$classifications[["projected_dispersion_probe"]],1L)
  expect_identical(a$compact$maxima[["replay_F_gap"]],1)
  local_mocked_bindings(working_hl_profile=function(fit,...)fit,
    working_hl_history_compare=function(state,fit)list(point=TRUE,covariance=FALSE,
      agreement=0,covariance_agreement=.25,F=state$F),.package="scRoPE")
  a<-ahr_pair(d$history,d$fit)
  expect_true(a$compact$success)
  expect_false(a$compact$covariance)
  expect_identical(a$compact$reason,"history_tied_covariance_unresolved")
  expect_identical(a$compact$maxima[["covariance_agreement"]],.25)
})

test_that("compact audit accounts for generated probes and overflow without dropping checks", {
  ahr_mocks();d<-ahr_fixture()
  local_mocked_bindings(
    working_hl_profile=function(fit,...) {fit$needs_refinement<-TRUE;fit},
    working_hl_history_compare=function(state,fit)list(
      point=!isTRUE(state$needs_refinement),covariance=TRUE,
      agreement=if(isTRUE(state$needs_refinement)).1 else 0,covariance_agreement=0,F=state$F),
    working_hl_refine=function(state,template,history,source) {
      key<-history$register(template$data,template$active,template$beta_star,
        template$K,TRUE,template$rho_full)
      for(i in 1:2)history$record(key,template$rho_full,
        list(success=TRUE,beta=template$beta,F=template$F,status="generated"),TRUE,source)
      list(success=TRUE,state=template)
    },.package="scRoPE")
  a<-ahr_pair(d$history,d$fit)
  expect_true(a$compact$success)
  expect_identical(a$compact$initial_retained,1L)
  expect_identical(a$compact$generated_evaluations,2L)
  expect_identical(a$compact$processed,3L)
  expect_identical(a$compact$retained,3L)
  expect_identical(a$compact$pending_events,0L)
  expect_identical(a$full$observations[[2L]]$parent_event,1L)
  h<-d$history;h$limits$events<-2L
  a<-ahr_pair(h,d$fit)
  expect_false(a$compact$success)
  expect_identical(a$compact$reason,"history_collection_incomplete")
  expect_identical(a$compact$generated_evaluations,2L)
  expect_identical(a$compact$observed,3L)
  expect_identical(a$compact$retained,2L)
  expect_identical(a$compact$processed,1L)
  expect_identical(a$compact$pending_events,1L)
  expect_identical(a$compact$dropped,1L)
  expect_identical(a$compact$overflow_source,"audit_tied_refinement")
  expect_identical(a$compact$history_evidence$overflow_source,"audit_tied_refinement")
  expect_identical(a$compact$failure_witness$observation$id,1L)
})

test_that("working refit requests compact audit and failure history by default", {
  ahr_mocks();d<-ahr_fixture();refit<-ahr_get("working_hl_refit")
  native<-list(point_success=FALSE,fit_status="synthetic_native_failure")
  full<-refit(d$fit$data,native,TRUE,d$collector)
  compact<-refit(d$fit$data,native,FALSE,d$collector)
  expect_identical(compact$reason,full$reason)
  expect_length(full$history_evidence$events,1L)
  expect_null(compact$history_evidence$events)
  expect_identical(compact$history_evidence$observed,full$history_evidence$observed)
  local_mocked_bindings(working_hl_refine=function(state,...)list(success=TRUE,state=state,trace=list()),
    hl_fo_finalize_state=function(state)state,
    working_hl_audit=function(...)list(point=TRUE,covariance=TRUE,reason="certified",Vbeta=diag(2)),
    .package="scRoPE")
  full<-refit(d$fit$data,d$fit,TRUE,d$collector)
  compact<-refit(d$fit$data,d$fit,FALSE,d$collector)
  expect_true(full$success);expect_true(compact$success)
  expect_identical(compact$beta,full$beta)
  expect_identical(compact$Vbeta,full$Vbeta)
  expect_identical(compact$inference_available,full$inference_available)
  expect_length(full$history_audit$observations,1L)
  expect_null(compact$history_audit$observations)
  expect_false(compact$history_audit$details_retained)
})
