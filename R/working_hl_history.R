working_hl_history_criterion <- function()list(
  id="HL_profile_Laplace_p_root",version=1L,direction="maximize",
  beta_equation="p_beta=0",mode_equation="h_eta=0",
  coordinates="conditional_beta_log_sigma2_log_NB_precision",
  empty_donors="retain_Gamma_Laplace_r(v)")

# Private, bounded capture. It observes profile calls without changing their
# starts, caches or return objects; ordinary fits do not create a collector.
working_hl_history <- function(max_events=10000L,max_bytes=64*1024^2,known_configurations=character(0),parent_event=NULL) {
  configs<-events<-list();bytes<-0;seen<-0L;overflow<-FALSE
  execution<-.hl_score_execution_identity()
  register<-function(data,active,beta_star,K,constant=TRUE,rho_template=c(NA_real_,NA_real_)) {
    config<-list(data=data,active=as.integer(active),beta_star=as.numeric(beta_star),
                 K=K,criterion=working_hl_history_criterion(),constant=isTRUE(constant),held_log=rho_template[setdiff(1:2,active)],execution=execution)
    key<-.hl_score_fingerprint(config)
    if(is.null(configs[[key]])) {
      cost<-if(key%in%known_configurations)0 else as.numeric(utils::object.size(config))
      if(bytes+cost>max_bytes)overflow<<-TRUE else {
        configs[[key]]<<-config;bytes<<-bytes+cost
      }
    }
    key
  }
  record<-function(config,rho,ans,cold,source) {
    seen<<-seen+1L
    event<-list(id=seen,origin_id=seen,parent_event=parent_event,configuration=config,rho=as.numeric(rho),
      beta=if(isTRUE(ans$success)||length(ans$F))as.numeric(ans$beta) else NULL,
      F=if(length(ans$F)==1L)as.numeric(ans$F) else NA_real_,
      success=isTRUE(ans$success),status=ans$status,reason=ans$reason,
      cold=isTRUE(cold),source=source)
    cost<-as.numeric(utils::object.size(event))
    if(length(events)>=max_events||bytes+cost>max_bytes)overflow<<-TRUE else {
      events[[length(events)+1L]]<<-event;bytes<<-bytes+cost
    }
    invisible(NULL)
  }
  snapshot<-function()list(schema=1L,configurations=configs,events=events,
    observed=seen,retained=length(events),bytes=bytes,overflow=overflow,
    limits=list(events=max_events,bytes=max_bytes),execution=execution,
    byte_accounting="sum_object_size_unique_configurations_and_events")
  list(register=register,record=record,snapshot=snapshot)
}

# Optional phase labels leave synthetic evaluators and native results intact.
hl_fo_history_evaluate <- function(evaluator,rho,cold=FALSE,source="outer") {
  if(is.function(evaluator$mark))evaluator$mark(source)
  evaluator$evaluate(rho,cold=cold)
}

# Independently reconstruct the recorded scalar at its original beta/rho;
# recompute donor modes but do not refit beta or compare across working masks.
working_hl_history_scalar <- function(config,event) {
  data<-config$data;rho<-event$rho;beta<-event$beta
  modes<-hl_fo_solve_modes(data,beta,rho,eta_start=NULL)
  if(!isTRUE(modes$success))stop("Historical scalar donor modes failed.")
  base<-as.numeric(data$log_offset+data$X%*%beta)
  zeros<-if(is.null(data$zero_donors))0L else data$zero_donors
  terms<-numeric(data$m+zeros)
  for(i in seq_along(terms)) {
    keep<-if(i<=data$m)data$subject_rows[[i]] else integer(0)
    eta<-if(i<=data$m)modes$eta[i] else exp(rho[1])/2
    s<-hl_fo_mode_components(data$y[keep],base[keep],eta,rho[1],rho[2])
    if(!isTRUE(s$ok)||!is.finite(s$h_etaeta)||s$h_etaeta>=0)
      stop("Historical scalar component failed.")
    terms[i]<-sum(s$nb$loglik)+s$gamma$g-.5*log(-s$h_etaeta)+
      if(config$constant).5*log(2*pi) else 0
  }
  F<-sum(terms)
  if(!is.finite(F))stop("Historical scalar is nonfinite.")
  F
}

# Record identities are consistency checks, not authentication against a caller
# forging both a record and its digest. No rounded-coordinate cache is used.
working_hl_history_audit <- function(history,fit,max_replays=5000L,keep_details=TRUE) {
  began<-proc.time()[["elapsed"]];rows<-list();profiles<-scalars<-list()
  covariance<-TRUE;unique_replays<-0L;initial_retained<-history$retained
  generated_evaluations<-0L;processed<-0L;current<-NULL;failure_configuration<-NULL
  classifications<-setNames(integer(4L),c("failed_probe_without_profile",
    "same_target_probe","projected_dispersion_probe","unclassified"))
  maxima<-setNames(rep(NA_real_,6L),c("scalar_discrepancy","cold_profile_discrepancy",
    "comparable_original_F_gap","replay_F_gap","point_agreement","covariance_agreement"))
  # The current observation is sufficient for all existing gates. Retaining
  # prior processed rows is optional; the mandatory input queue is unchanged.
  retain_current<-function() {
    if(is.null(current))return(invisible(NULL))
    processed<<-processed+1L
    if(isTRUE(keep_details))rows[[processed]]<<-current
    classification<-current$classification
    if(is.null(classification))classification<-"unclassified"
    classifications[classification]<<-classifications[classification]+1L
    values<-list(scalar_discrepancy=current$scalar_discrepancy,
      cold_profile_discrepancy=current$cold_profile_discrepancy,
      comparable_original_F_gap=if(isTRUE(current$scalar_comparable))current$original_F-fit$F else NULL,
      replay_F_gap=if(!is.null(current$replay$F))current$replay$F-fit$F else NULL,
      point_agreement=current$replay$comparison$agreement,
      covariance_agreement=current$replay$comparison$covariance_agreement)
    for(name in names(values)) {
      value<-values[[name]]
      if(is.numeric(value)&&length(value)==1L&&!is.na(value))
        maxima[name]<<-if(is.na(maxima[name]))value else max(maxima[name],value)
    }
    current<<-NULL
    invisible(NULL)
  }
  finish<-function(success,reason) {
    witness<-if(!success)list(reason=reason,configuration=failure_configuration,
      observation=current) else NULL
    retain_current()
    list(success=success,reason=reason,
      configurations=working_hl_history_config_summary(history),
      covariance=covariance,observations=if(isTRUE(keep_details))rows else NULL,
      unique_replays=unique_replays,
      elapsed_seconds=proc.time()[["elapsed"]]-began,retained=history$retained,
      observed=history$observed,dropped=history$observed-history$retained,
      bytes=history$bytes,byte_accounting=history$byte_accounting,limits=history$limits,
      initial_retained=initial_retained,generated_evaluations=generated_evaluations,
      processed=processed,pending_event_ids=if(isTRUE(keep_details))vapply(
        history$events[setdiff(seq_along(history$events),seq_len(processed))],
        function(x)if(length(x$id)==1L)as.integer(x$id) else NA_integer_,integer(1)) else NULL,
      history_evidence=if(!success) {
        if(isTRUE(keep_details))working_hl_history_compact(history) else
          working_hl_history_summary(history)
      } else NULL,
      details_retained=isTRUE(keep_details),pending_events=length(history$events)-processed,
      overflow=isTRUE(history$overflow),overflow_source=history$overflow_source,
      classifications=classifications,maxima=maxima,failure_witness=witness)
  }

  if(!identical(history$schema,1L)||isTRUE(history$overflow)||
     !identical(history$observed,history$retained)||history$retained!=length(history$events))
    return(finish(FALSE,"history_collection_incomplete"))
  if(!length(history$events))return(finish(FALSE,"history_evidence_missing"))
  execution<-.hl_score_execution_identity()
  if(!identical(history$execution,execution))return(finish(FALSE,"history_execution_mismatch"))
  validate_configuration<-function(key,config) {
    if(!identical(config$criterion,working_hl_history_criterion()))return("history_criterion_mismatch")
    if(!identical(config$execution,execution)||!identical(.hl_score_fingerprint(config),key))
      return("history_configuration_mismatch")
    NULL
  }
  for(key in names(history$configurations)) {
    reason<-validate_configuration(key,history$configurations[[key]])
    if(!is.null(reason)) {
      failure_configuration<-key
      return(finish(FALSE,reason))
    }
  }
  target_data<-fit$data;target_data$inference_m<-NULL
  target_identity<-.hl_score_fingerprint(target_data)
  ids<-vapply(history$events,function(x)x$id,integer(1))
  if(!identical(ids,seq_along(history$events)))return(finish(FALSE,"history_event_identity_mismatch"))
  i<-0L
  while(i<length(history$events)) {
    retain_current()
    i<-i+1L
    event<-history$events[[i]];config<-history$configurations[[event$configuration]]
    current<-event
    if(is.null(config))return(finish(FALSE,"history_configuration_missing"))
    valid_rho<-is.numeric(event$rho)&&length(event$rho)==2L&&all(is.finite(event$rho))
    valid_beta<-is.numeric(event$beta)&&length(event$beta)==config$data$p&&all(is.finite(event$beta))
    finite_F<-is.numeric(event$F)&&length(event$F)==1L&&is.finite(event$F)
    required<-isTRUE(event$success)||finite_F||valid_beta
    if(!required) {
      current$classification<-"failed_probe_without_profile"
      next
    }
    if(!valid_rho||!valid_beta)return(finish(FALSE,"history_profile_coordinates_missing"))
    if(!identical(event$rho[setdiff(1:2,config$active)],as.numeric(config$held_log)))
      return(finish(FALSE,"history_fixed_configuration_mismatch"))
    if(length(config$beta_star)!=length(event$beta)||
       hl_fo_sdisc(event$beta,as.numeric(config$beta_star+config$K%*%
         crossprod(config$K,event$beta-config$beta_star)))>1e-8)
      return(finish(FALSE,"history_affine_coordinates_mismatch"))
    scalar_key<-.hl_score_fingerprint(list(event$configuration,event$rho,event$beta))
    scalar<-scalars[[scalar_key]]
    if(is.null(scalar)) {
      scalar<-tryCatch(working_hl_history_scalar(config,event),error=function(e)e)
      scalars[[scalar_key]]<-scalar
    }
    if(inherits(scalar,"error")) {
      current$error<-conditionMessage(scalar)
      return(finish(FALSE,"history_scalar_evaluation_failed"))
    }
    discrepancy<-if(finite_F)abs(event$F-scalar) else 0
    current$original_F<-scalar;current$scalar_discrepancy<-discrepancy
    if(discrepancy>1e-8*(1+abs(scalar)))return(finish(FALSE,"history_scalar_mismatch"))
    rho<-event$rho;held<-setdiff(1:2,fit$active);rho[held]<-fit$rho_full[held]
    source_data<-config$data;source_data$inference_m<-NULL
    comparable<-identical(.hl_score_fingerprint(source_data),target_identity)&&
      identical(event$rho,rho)&&
      identical(config$K,fit$K)&&identical(config$beta_star,fit$beta_star)&&isTRUE(config$constant)
    projected<-!comparable||!identical(config$active,fit$active)
    current$classification<-if(projected)"projected_dispersion_probe" else "same_target_probe"
    current$scalar_comparable<-comparable
    current$selected_rho<-rho
    current$selected_margin<-1e-10*(1+abs(fit$F))
    if(comparable&&scalar>fit$F+current$selected_margin)
      return(finish(FALSE,"history_better_original_profile"))
    key<-paste(sprintf("%a",rho),collapse=":")
    replay<-profiles[[key]]
    if(is.null(replay)) {
      if(unique_replays>=max_replays)return(finish(FALSE,"history_replay_budget"))
      unique_replays<-unique_replays+1L
      replay<-tryCatch(working_hl_profile(fit,rho,fit$active),error=function(e)e)
      # Cache compact comparisons, never the repeated donor derivative tensors.
      if(!inherits(replay,"error")) {
        tied<-abs(replay$F-fit$F)<=1e-10*(1+max(abs(replay$F),abs(fit$F)))
        station<-hl_fo_outer_stationarity(replay)
        root_like<-isTRUE(station$maximum)&&is.finite(station$scaled_score)&&station$scaled_score<=1e-4
        comparison<-NULL
        if(tied&&root_like) {
          comparison<-working_hl_history_compare(replay,fit)
          if(!isTRUE(comparison$point)) {
            refinement_history<-working_hl_history(
              max_events=history$limits$events-length(history$events),
              max_bytes=history$limits$bytes-history$bytes,
              known_configurations=names(history$configurations),parent_event=event$id)
            refined<-working_hl_refine(replay,fit,refinement_history,"history_tied_refinement")
            extra<-refinement_history$snapshot()
            generated_evaluations<-generated_evaluations+extra$observed
            for(name in names(extra$configurations))
              if(is.null(history$configurations[[name]]))
                history$configurations[[name]]<-extra$configurations[[name]]
            for(event_new in extra$events) {
              event_new$id<-length(history$events)+1L
              history$events[[event_new$id]]<-event_new
            }
            history$observed<-history$observed+extra$observed
            history$retained<-length(history$events)
            history$bytes<-history$bytes+extra$bytes
            if(isTRUE(extra$overflow)) {
              history$overflow<-TRUE
              history$overflow_source<-"audit_tied_refinement"
              return(finish(FALSE,"history_collection_incomplete"))
            }
            for(name in names(extra$configurations)) {
              reason<-validate_configuration(name,extra$configurations[[name]])
              if(!is.null(reason)) {
                failure_configuration<-name
                return(finish(FALSE,reason))
              }
            }
            if(!isTRUE(refined$success))replay<-simpleError("Historical tied root refinement failed.") else {
              comparison<-working_hl_history_compare(refined$state,fit)
            }
          }
        }
        if(!inherits(replay,"error"))replay<-list(F=replay$F,beta=replay$beta,rho=replay$rho_full,
                                                 tied=tied,root_like=root_like,
                                                 comparison=comparison)
      }
      profiles[[key]]<-replay
    }
    if(inherits(replay,"error")) {
      current$error<-conditionMessage(replay)
      return(finish(FALSE,"history_required_replay_failed"))
    }
    current$replay<-replay
    if(comparable) {
      discrepancy<-max(hl_fo_sdisc(event$beta,replay$beta),hl_fo_sdisc(scalar,replay$F))
      current$cold_profile_discrepancy<-discrepancy
      if(!is.finite(discrepancy)||discrepancy>1e-8)
        return(finish(FALSE,"history_cold_profile_mismatch"))
    }
    if(!is.finite(replay$F)||replay$F>fit$F+1e-10*(1+abs(fit$F)))
      return(finish(FALSE,"history_better_profile"))
    comparison<-replay$comparison
    if(!is.null(comparison)) {
      if(!comparison$point||comparison$F>fit$F+1e-10*(1+abs(fit$F)))
        return(finish(FALSE,"history_unresolved_tied_root"))
      covariance<-covariance&&comparison$covariance
    }
  }
  finish(TRUE,if(covariance)"history_supported" else "history_tied_covariance_unresolved")
}

working_hl_history_config_summary <- function(history) {
  lapply(history$configurations,function(config) {
    data<-config$data
    out<-config;out$data<-NULL
    out$data_identity<-.hl_score_fingerprint(data)
    out$dimensions<-list(cells=data$n,donors=data$m,p=data$p,
      zero_donors=if(is.null(data$zero_donors))0L else data$zero_donors)
    out
  })
}

working_hl_history_compact <- function(history) {
  history$configurations<-working_hl_history_config_summary(history)
  history
}

# Resource and identity evidence for compact failures. Complete event histories
# remain internal to the computation or are returned only in detailed mode.
working_hl_history_summary <- function(history) {
  fields<-c("schema","observed","retained","bytes","overflow","limits",
    "execution","byte_accounting","overflow_source")
  out<-history[intersect(fields,names(history))]
  out$configurations<-working_hl_history_config_summary(history)
  out$dropped<-history$observed-history$retained
  out$event_count<-length(history$events)
  out$event_id_bounds<-if(length(history$events))list(
    first=history$events[[1L]]$id,last=history$events[[length(history$events)]]$id) else NULL
  out$details_retained<-FALSE
  out
}

working_hl_history_compare <- function(state,fit) {
  cert<-hl_fo_finalize_state(state)
  agreement<-max(hl_fo_sdisc(state$beta,fit$beta),hl_fo_sdisc(state$rho_full,fit$rho_full),
    hl_fo_sdisc(state$eta,fit$eta),hl_fo_sdisc(state$augmented$A,fit$augmented$A),
    hl_fo_sdisc(state$augmented$J,fit$augmented$J))
  cv<-if(isTRUE(cert$covariance_success)&&isTRUE(fit$covariance_success))
    working_covariance_discrepancy(cert$covariance$V_beta_rho,fit$covariance$V_beta_rho) else Inf
  list(point=isTRUE(cert$point_success)&&agreement<=1e-8,
    agreement=agreement,covariance=cv<=1e-8,covariance_agreement=cv,F=state$F)
}
