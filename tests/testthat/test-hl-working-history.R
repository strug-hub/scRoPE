wh <- function(name)getFromNamespace(name,'scRoPE')
wh_fixture <- local({
 saved<-NULL
 function() {
  if(!is.null(saved))return(saved)
  set.seed(90216);id<-rep(1:12,each=5);X<-cbind(1,rnorm(60))
  y<-rpois(60,exp(.3+.2*X[,2]+rep(rnorm(12,0,.4),each=5)))
  data<-wh('hl_fo_make_data')(y,X,rep(0,60),id,1L)
  history<-wh('working_hl_history')()
  fit<-wh('hl_fo_fit')(data,active=integer(0),fixed_natural=c(sigma2=.4,phi=7),
    keep_diagnostics=FALSE,.history=history)
  stopifnot(fit$point_success,fit$covariance_success)
  saved<<-list(fit=fit,history=history$snapshot(),data=data)
  saved
 }
})
wh_audit <- function(h,f)wh('working_hl_history_audit')(h,f,keep_details=TRUE)

test_that('installed history capture and audit work without attached utils', {
 d<-wh_fixture();input<-tempfile(fileext='.rds');script<-tempfile(fileext='.R')
 on.exit(unlink(c(input,script)),add=TRUE);saveRDS(d$data,input)
 quote_r<-function(x)paste(deparse(x),collapse='')
 lib<-dirname(getNamespaceInfo(asNamespace('scRoPE'),'path'))
 writeLines(c(sprintf('library(scRoPE,lib.loc=%s)',quote_r(lib)),
   'if("package:utils"%in%search())detach("package:utils")',
   'stopifnot(!"package:utils"%in%search())',
   sprintf('data<-readRDS(%s)',quote_r(input)),
   'h<-scRoPE:::working_hl_history()',
   'fit<-scRoPE:::hl_fo_fit(data,active=integer(0),fixed_natural=c(sigma2=.4,phi=7),keep_diagnostics=FALSE,.history=h)',
   'stopifnot(fit$point_success,fit$covariance_success)',
   's<-h$snapshot();a<-scRoPE:::working_hl_history_audit(s,fit,keep_details=TRUE)',
   'expected<-sum(vapply(c(s$configurations,s$events),function(x)as.numeric(utils::object.size(x)),numeric(1)))',
   'stopifnot(a$success,identical(s$bytes,expected),!"package:utils"%in%search())',
   'cat("HISTORY_WITHOUT_ATTACHED_UTILS_VERIFIED\\n")'),script)
 out<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(script)),stdout=TRUE,stderr=TRUE)
 expect_null(attr(out,'status'),info=paste(out,collapse='\n'))
 expect_true(any(grepl('HISTORY_WITHOUT_ATTACHED_UTILS_VERIFIED',out)),info=paste(out,collapse='\n'))
})

test_that('captured historical scalars reproduce actual HL profile values', {
 d<-wh_fixture();h<-d$history;e<-h$events[[1]];cfg<-h$configurations[[e$configuration]]
 expect_equal(wh('working_hl_history_scalar')(cfg,e),d$fit$F,tolerance=1e-12)
 a<-wh_audit(h,d$fit)
 expect_true(a$success);expect_identical(a$reason,'history_supported')
 expect_true(a$covariance);expect_equal(a$unique_replays,1L)
 expect_true(a$observations[[1]]$scalar_comparable)
 expect_length(a$configurations,1L)
 expect_null(a$configurations[[1]][['data',exact=TRUE]])
 expect_identical(a$observations[[1]]$rho,d$fit$rho_full)
})

test_that('incomplete and inconsistent finite history cannot certify recovery', {
 d<-wh_fixture()
 cases<-list(
  missing_coordinates=function(h){h$events[[1]]$rho<-NULL;h},
  scalar_disagreement=function(h){h$events[[1]]$F<-h$events[[1]]$F+1;h},
  configuration=function(h){h$configurations[[1]]$constant<-FALSE;h},
  backend=function(h){h$execution$inner<-'altered';h},
  missing_configuration=function(h){h$events[[1]]$configuration<-'missing';h},
  overflow=function(h){h$overflow<-TRUE;h},
  missing_event=function(h){h$observed<-h$observed+1L;h},
  missing_identity=function(h){h$events[[1]]$id<-2L;h})
 reasons<-c('history_profile_coordinates_missing','history_scalar_mismatch',
  'history_configuration_mismatch','history_execution_mismatch','history_configuration_missing',
  'history_collection_incomplete','history_collection_incomplete','history_event_identity_mismatch')
 for(i in seq_along(cases)) {
  a<-wh_audit(cases[[i]](d$history),d$fit)
  expect_false(a$success,info=names(cases)[i]);expect_identical(a$reason,reasons[i])
 }
 duplicate<-d$history;e<-duplicate$events[[1]];e$id<-2L;e$F<-e$F+1
 duplicate$events[[2]]<-e;duplicate$retained<-duplicate$observed<-2L
 expect_identical(wh_audit(duplicate,d$fit)$reason,'history_scalar_mismatch')
})

test_that('same-criterion higher observations are not hidden behind cold replay', {
 d<-wh_fixture();f<-d$fit;f$F<-f$F-1
 expect_identical(wh_audit(d$history,f)$reason,'history_better_original_profile')
 # An active-mask change at identical full rho does not change scalar comparability.
 h<-d$history;config<-h$configurations[[1]];config$active<-1:2;config$held_log<-numeric(0)
 key<-wh('.hl_score_fingerprint')(config);h$configurations<-setNames(list(config),key)
 h$events[[1]]$configuration<-key
 expect_identical(wh_audit(h,f)$reason,'history_better_original_profile')
 a<-wh_audit(h,d$fit)
 expect_true(a$success);expect_true(a$observations[[1]]$scalar_comparable)
 expect_identical(a$observations[[1]]$classification,'projected_dispersion_probe')
})

test_that('failed no-profile trials are explicit but required replay failures refuse', {
 d<-wh_fixture();h<-d$history;e<-h$events[[1]]
 e$id<-2L;e$F<-NA_real_;e$beta<-NULL;e$success<-FALSE;e$status<-'deliberate_failure'
 h$events[[2]]<-e;h$retained<-h$observed<-2L
 a<-wh_audit(h,d$fit)
 expect_true(a$success)
 expect_identical(a$observations[[2]]$classification,'failed_probe_without_profile')
 expect_identical(a$observations[[2]]$status,'deliberate_failure')
 local_mocked_bindings(working_hl_profile=function(...)stop('required replay failed'),.package='scRoPE')
 a<-wh_audit(d$history,d$fit)
 expect_false(a$success);expect_identical(a$reason,'history_required_replay_failed')
 expect_match(a$observations[[1]]$error,'required replay failed')
})

test_that('projected competitors use the selected configuration and retain strict rank gates', {
 d<-wh_fixture();profile<-wh('working_hl_profile')
 h<-d$history;config<-h$configurations[[1]]
 # Changed held coordinate: original F is evaluated at its original coordinates.
 e<-h$events[[1]];e$rho[1]<-e$rho[1]+.01;config$held_log<-e$rho
 key<-wh('.hl_score_fingerprint')(config);h$configurations<-setNames(list(config),key)
 e$configuration<-key;e$F<-wh('working_hl_history_scalar')(config,e);h$events[[1]]<-e
 a<-wh_audit(h,d$fit)
 expect_true(a$success);expect_false(a$observations[[1]]$scalar_comparable)
 expect_identical(a$observations[[1]]$selected_rho,d$fit$rho_full)
 local_mocked_bindings(working_hl_profile=function(...) {
  out<-profile(...);out$F<-out$F+1;out
 },.package='scRoPE')
 expect_identical(wh_audit(h,d$fit)$reason,'history_better_profile')
})

test_that('tied roots require complete point and covariance agreement', {
 d<-wh_fixture();profile<-wh('working_hl_profile')
 run<-function(field) {
  local_mocked_bindings(working_hl_profile=function(...) {
   out<-profile(...)
   if(field=='beta')out$beta[1]<-out$beta[1]+.1
   if(field=='meat')out$augmented$J<-out$augmented$J*1.1
   out
  },.package='scRoPE')
  wh_audit(d$history,d$fit)
 }
 expect_identical(run('beta')$reason,'history_cold_profile_mismatch')
 expect_identical(run('meat')$reason,'history_unresolved_tied_root')
 finalize<-wh('hl_fo_finalize_state')
 local_mocked_bindings(hl_fo_finalize_state=function(...) {
  out<-finalize(...);out$covariance$V_beta_rho<-out$covariance$V_beta_rho*1.1;out
 },.package='scRoPE')
 a<-wh_audit(d$history,d$fit)
 expect_true(a$success);expect_false(a$covariance)
 expect_identical(a$reason,'history_tied_covariance_unresolved')
})

test_that('unadopted audit-refinement probes join the same bounded audit queue', {
 d<-wh_fixture();refine<-wh('working_hl_refine');profile<-wh('working_hl_profile')
 local_mocked_bindings(working_hl_profile=function(...) {
  out<-profile(...);out$augmented$A<-out$augmented$A*1.1;out
 },working_hl_refine=function(state,template,history=NULL,source='unused',...) {
  out<-refine(d$fit,template)
  config<-history$register(template$data,template$active,template$beta_star,template$K,TRUE,template$rho_full)
  bad<-state;bad$F<-bad$F+1
  history$record(config,state$rho_full,bad,TRUE,'unadopted_injected_polish')
  out
 },.package='scRoPE')
 a<-wh_audit(d$history,d$fit)
 expect_false(a$success);expect_identical(a$reason,'history_scalar_mismatch')
 expect_equal(a$generated_evaluations,1L)
 expect_identical(tail(a$observations,1)[[1]]$source,'unadopted_injected_polish')
})

test_that('ledger and replay budgets cannot silently discard evidence', {
 d<-wh_fixture();h<-wh('working_hl_history')(max_events=0L)
 invisible(wh('hl_fo_fit')(d$data,active=integer(0),fixed_natural=c(sigma2=.4,phi=7),.history=h))
 expect_true(h$snapshot()$overflow)
 expect_identical(wh_audit(h$snapshot(),d$fit)$reason,'history_collection_incomplete')
 expect_identical(wh('working_hl_history_audit')(d$history,d$fit,max_replays=0L)$reason,
   'history_replay_budget')
})

test_that('ordinary native results bypass new history vetoes even on ledger overflow', {
 set.seed(90211);m<-60L;k<-8L;id<-rep(seq_len(m),each=k)
 x<-rnorm(length(id));z<-rep(rnorm(m),each=k)
 y<-rpois(length(id),exp(.3+.2*x+.1*z+rep(rnorm(m,0,.5),each=k)))
 X<-cbind(Intercept=1,x=x,z=z)
 native<-scrope(y,id,X,model='HL',cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE,recovery='none')
 collector<-wh('working_hl_history')
 local_mocked_bindings(working_hl_history=function(...)collector(max_events=0L),.package='scRoPE')
 a<-scrope_working(y,id,X,model='HL',cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE)
 core<-c('summary','overdispersion','convergence','algorithm','random_effect')
 expect_identical(a$primary[core],native[core])
 expect_identical(a$primary$diagnostics$hl,native$diagnostics$hl)
 expect_identical(a$routes$route,'native')
 expect_true(a$routes$native_regular)
 expect_identical(a$summary$beta,as.numeric(native$summary[1,grep('^logFC_',names(native$summary))]))
 expect_identical(a$summary$robust_se,as.numeric(native$summary[1,grep('^se_robust_',names(native$summary))]))
 expect_null(a$audit[[1]]$history_audit)
})

test_that('profile capture is independent of optional traces and labels every phase', {
 d<-wh_fixture();h<-wh('working_hl_history')()
 e<-wh('hl_fo_make_outer_evaluator')(d$data,d$fit$rho_full,integer(0),d$fit$beta_star,
   d$fit$K,d$fit$beta,keep_diagnostics=FALSE,.history=h,.history_source='round1')
 sources<-c('start','trust_proposal','total_gradient_polish','stationary_replay','cold_polish','final_cold')
 for(source in sources) {
  st<-wh('hl_fo_history_evaluate')(e,numeric(0),cold=TRUE,source=source)
  expect_true(st$success)
 }
 expect_length(e$trace(),0L)
 captured<-h$snapshot()
 expect_equal(captured$retained,length(sources))
 expect_true(all(vapply(captured$events,function(x)is.finite(x$F)&&length(x$beta)==2L,logical(1))))
 expect_identical(vapply(captured$events,`[[`,character(1),'source'),
                  paste('round1',sources,seq_along(sources),sep=':'))
 expect_true(wh_audit(captured,d$fit)$success)
})

test_that('exceptional profile attempts are recorded before their error propagates', {
 d<-wh_fixture();h<-wh('working_hl_history')()
 e<-wh('hl_fo_make_outer_evaluator')(d$data,d$fit$rho_full,integer(0),d$fit$beta_star,
   d$fit$K,d$fit$beta,keep_diagnostics=FALSE,.history=h)
 local_mocked_bindings(hl_fo_profile_evaluate=function(...)stop('injected evaluation exception'),.package='scRoPE')
 expect_error(e$evaluate(numeric(0),cold=TRUE),'injected evaluation exception')
 record<-h$snapshot()$events[[1]]
 expect_identical(record$status,'evaluation_error')
 expect_identical(record$rho,d$fit$rho_full)
 expect_false(record$success)
 expect_error(wh('working_hl_profile')(d$fit,d$fit$rho_full,d$fit$active,h),'injected evaluation exception')
 expect_equal(h$snapshot()$retained,2L)
})

test_that('already certified complete-state ties avoid redundant polishing', {
 d<-wh_fixture()
 local_mocked_bindings(working_hl_refine=function(...)stop('unnecessary polishing'),.package='scRoPE')
 a<-wh_audit(d$history,d$fit)
 expect_true(a$success);expect_equal(a$generated_evaluations,0L)
})

test_that('byte and generated-event exhaustion retain complete refusal evidence', {
 d<-wh_fixture();h<-wh('working_hl_history')(max_bytes=1)
 invisible(wh('hl_fo_fit')(d$data,active=integer(0),fixed_natural=c(sigma2=.4,phi=7),.history=h))
 a<-wh_audit(h$snapshot(),d$fit)
 expect_identical(a$reason,'history_collection_incomplete')
 expect_equal(a$history_evidence$observed,1L)
 expect_equal(a$history_evidence$retained,0L)
 # Force a noncoalesced first replay, then generate a trial beyond the event cap.
 profile<-wh('working_hl_profile');refine<-wh('working_hl_refine')
 local_mocked_bindings(working_hl_profile=function(...) {
  out<-profile(...);out$augmented$A<-out$augmented$A*1.1;out
 },working_hl_refine=function(state,template,history=NULL,source='unused',...) {
  config<-history$register(template$data,template$active,template$beta_star,template$K,TRUE,template$rho_full)
  history$record(config,state$rho_full,state,TRUE,'generated_before_failure')
  history$record(config,state$rho_full,state,TRUE,'generated_over_budget')
  list(success=FALSE,reason='injected_refinement_failure')
 },.package='scRoPE')
 h<-d$history;h$limits$events<-2L
 a<-wh_audit(h,d$fit)
 expect_identical(a$reason,'history_collection_incomplete')
 expect_equal(a$history_evidence$observed,3L)
 expect_equal(a$history_evidence$retained,2L)
 expect_identical(a$history_evidence$overflow_source,'audit_tied_refinement')
 expect_identical(a$pending_event_ids,2L)
 expect_identical(a$history_evidence$events[[2]]$source,'generated_before_failure')
 # With adequate budget, both observations survive an early refinement failure.
 h$limits$events<-10L;a<-wh_audit(h,d$fit)
 expect_identical(a$reason,'history_required_replay_failed')
 expect_equal(a$history_evidence$retained,3L)
 expect_identical(a$pending_event_ids,2:3)
})

test_that('criterion identity is checked independently of its consistency digest', {
 d<-wh_fixture();h<-d$history;config<-h$configurations[[1]]
 config$criterion$direction<-'minimize'
 key<-wh('.hl_score_fingerprint')(config)
 h$configurations<-setNames(list(config),key);h$events[[1]]$configuration<-key
 expect_identical(wh_audit(h,d$fit)$reason,'history_criterion_mismatch')
})

test_that('overflow diagnostics report retained event identities rather than positions', {
 d<-wh_fixture();h<-d$history;last<-h$events[[1]];last$id<-3L
 h$events[[2]]<-last;h$observed<-3L;h$retained<-2L;h$overflow<-TRUE
 a<-wh_audit(h,d$fit)
 expect_identical(a$reason,'history_collection_incomplete')
 expect_identical(a$pending_event_ids,c(1L,3L));expect_equal(a$dropped,1L)
})
