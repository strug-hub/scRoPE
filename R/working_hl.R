working_covariance_discrepancy <- function(a,b) {
  if(!identical(dim(a),dim(b))||any(!is.finite(c(a,b)))||any(diag(b)<0))return(Inf)
  scale<-outer(sqrt(diag(b)),sqrt(diag(b)));difference<-abs(a-b)
  max(ifelse(scale>0,difference/scale,ifelse(difference==0,0,Inf)))
}

# Independent augmented-equation audit: re-solve only donor modes under
# perturbations, leaving gamma and chi as independent coordinates.
working_hl_psi <- function(z, state) {
  data<-state$data;K<-state$K;k<-ncol(K);active<-state$active;q<-length(active)
  beta<-as.numeric(state$beta_star+K%*%z[seq_len(k)])
  chi<-z[k+seq_len(k)];rho<-state$rho_full
  if(q)rho[active]<-z[2L*k+seq_len(q)]
  modes<-hl_fo_solve_modes(data,beta,rho,eta_start=NULL)
  if(!isTRUE(modes$success))stop("Cold perturbed donor-mode solve failed.")
  transform<-hl_fo_scientific_transform(data$p,K,active)
  zeros<-if(is.null(data$zero_donors))0L else data$zero_donors
  out<-matrix(0,2L*k+q,data$m+zeros)
  for(ii in seq_len(ncol(out))) {
    keep<-if(ii<=data$m)data$subject_rows[[ii]] else integer(0)
    eta<-if(ii<=data$m)modes$eta[ii] else exp(rho[1L])/2
    raw<-hl_fo_subject_raw(data$y[keep],data$X[keep,,drop=FALSE],data$log_offset[keep],
      beta,eta,rho,include_laplace_constant=TRUE)
    if(!isTRUE(raw$ok))stop("Perturbed donor derivative failed.")
    s<-hl_fo_contract_subject(raw,transform)
    B<-s$p_hess[seq_len(k),seq_len(k),drop=FALSE]
    out[seq_len(k),ii]<-s$p_grad[seq_len(k)]
    out[k+seq_len(k),ii]<-as.numeric(t(B)%*%chi)-s$ell_grad[seq_len(k)]
    if(q) {
      E<-s$p_hess[seq_len(k),k+seq_len(q),drop=FALSE]
      out[2L*k+seq_len(q),ii]<-s$ell_grad[k+seq_len(q)]-as.numeric(crossprod(E,chi))
    }
  }
  out
}

working_hl_audit <- function(fit) {
  fail<-function(reason,metrics=list())list(point=FALSE,covariance=FALSE,reason=reason,metrics=metrics)
  if(!isTRUE(fit$point_success))return(fail(fit$fit_status))
  z<-c(fit$gamma,fit$augmented$chi,fit$rho)
  U<-fit$augmented$psi_subject;A<-fit$augmented$A
  base<-working_hl_psi(z,fit)
  metrics<-list(donor_equations=working_scaled_error(base,U),
    additivity=working_scaled_error(rowSums(U),fit$augmented$augmented_residual),
    meat=working_scaled_error(tcrossprod(U),fit$augmented$J),
    root=max(abs(rowSums(U))/(1+sqrt(rowSums(U^2)))),
    jacobian_fd=0,jacobian_stability=0,summed_jacobian_fd=0)
  if(any(!is.finite(c(A,U)))||rcond(A)<sqrt(.Machine$double.eps))return(fail("singular_augmented_sensitivity",metrics))
  for(j in seq_along(z)) {
    h<-.Machine$double.eps^(1/3)*(1+abs(z[j]))
    derivative<-function(step) {
      a<-b<-z;a[j]<-a[j]+step;b[j]<-b[j]-step
      (working_hl_psi(a,fit)-working_hl_psi(b,fit))/(2*step)
    }
    first<-derivative(h);half<-derivative(h/2)
    analytic<-vapply(fit$augmented$jacobian_subject,function(x)x[,j],numeric(length(z)))
    metrics$jacobian_fd<-max(metrics$jacobian_fd,working_scaled_error(first,analytic))
    metrics$jacobian_stability<-max(metrics$jacobian_stability,working_scaled_error(first,half))
    metrics$summed_jacobian_fd<-max(metrics$summed_jacobian_fd,
      working_scaled_error(-rowSums(first),A[,j]))
  }
  if(metrics$donor_equations>1e-8||metrics$additivity>1e-10||metrics$meat>1e-10||
     metrics$root>1e-6||metrics$jacobian_fd>1e-4||metrics$jacobian_stability>1e-4||
     metrics$summed_jacobian_fd>1e-4)return(fail("augmented_derivative_or_root_audit",metrics))
  delta<-solve(A,rowSums(U));V<-tcrossprod(solve(A,U))
  metrics$newton<-max(abs(delta)/(1+abs(z)))
  metrics$newton_in_se<-max(ifelse(diag(V)>0,abs(delta)/sqrt(diag(V)),ifelse(delta==0,0,Inf)))
  if(metrics$newton>1e-6)return(fail("augmented_newton_displacement",metrics))
  covariance<-isTRUE(fit$covariance_success)&&all(is.finite(V))&&metrics$newton_in_se<=1e-4
  k<-ncol(fit$K)
  Vbeta<-fit$K%*%V[seq_len(k),seq_len(k),drop=FALSE]%*%t(fit$K)
  if(covariance) {
    metrics$covariance_agreement<-working_covariance_discrepancy(fit$covariance$V_beta,Vbeta)
    covariance<-metrics$covariance_agreement<=1e-10
  }
  list(point=TRUE,covariance=covariance,reason=if(covariance)"certified" else "covariance_not_certified",
       metrics=metrics,Vbeta=Vbeta,V=V)
}

working_hl_profile <- function(fit,rho,active=seq_len(2L),history=NULL,source="working_profile") {
  ans<-tryCatch(hl_fo_profile_evaluate(fit$data,rho[active],rho,active,
    fit$beta_star,fit$K,beta_start=fit$beta,cold=TRUE,
    inner_acceptance_tier="derivative",keep_diagnostics=FALSE),error=function(e) {
      if(!is.null(history)) {
        configuration<-history$register(fit$data,active,fit$beta_star,fit$K,TRUE,rho)
        history$record(configuration,rho,list(success=FALSE,status="evaluation_error",
          reason=conditionMessage(e)),TRUE,source)
      }
      stop(e)
    })
  if(!is.null(history)) {
    configuration<-history$register(fit$data,active,fit$beta_star,fit$K,TRUE,rho)
    history$record(configuration,rho,ans,TRUE,source)
  }
  if(!isTRUE(ans$success))stop("Cold HL profile evaluation failed: ",ans$status)
  ans
}

# Working-only refinement goes beyond native outer stopping so that complete
# candidate roots can be compared at the prespecified coalescence tolerance.
working_hl_refine <- function(state,template,history=NULL,source="working_refinement",
                              coalesce_with_template=FALSE) {
  trace<-list()
  if(!length(template$active))return(list(success=TRUE,state=state,trace=trace))
  for(iteration in seq_len(12L)) {
    if(!isTRUE(hl_fo_finalize_state(state)$point_success))
      return(list(success=FALSE,reason="outer_refinement_curvature",trace=trace))
    step<-tryCatch(solve(state$H_F,state$g),error=function(e)NULL)
    if(is.null(step)||any(!is.finite(step)))return(list(success=FALSE,reason="outer_refinement_singular",trace=trace))
    trace[[iteration]]<-list(rho=state$rho_full,F=state$F,g=state$g,step=step)
    if(max(abs(step)/(1+abs(state$rho)))<=1e-11)
      return(list(success=TRUE,state=state,trace=trace))
    accepted<-FALSE
    for(rate in 2^-(0:12)) {
      rho<-state$rho_full;rho[template$active]<-state$rho-rate*step
      if(any(rho[template$active]<=template$guard_bounds$lower|
             rho[template$active]>=template$guard_bounds$upper))next
      trial<-tryCatch(working_hl_profile(template,rho,template$active,history,source),error=function(e)e)
      if(inherits(trial,"error"))return(list(success=FALSE,
        reason="outer_refinement_cold_failure",error=conditionMessage(trial),trace=trace))
      if(trial$F>=state$F-1e-12*(1+abs(state$F))&&max(abs(trial$g))<max(abs(state$g))) {
        state<-trial;accepted<-TRUE;break
      }
    }
    if(!accepted) {
      # A historical replay can already satisfy the complete target agreement
      # required by its caller even when further strict score decrease stalls.
      # This is not a general small-gradient exit: both roots must be stationary
      # and certified, and the unchanged full-state/F agreement must hold.
      # All attempted probes remain in history; covariance is checked separately
      # by the caller. Selected/candidate-root polishing keeps its strict exit.
      if(isTRUE(coalesce_with_template)) {
        station<-hl_fo_outer_stationarity(state)
        target_station<-hl_fo_outer_stationarity(template)
        target_cert<-hl_fo_finalize_state(template)
        comparison<-working_hl_history_compare(state,template)
        tied<-is.finite(state$F)&&is.finite(template$F)&&
          abs(state$F-template$F)<=1e-10*(1+max(abs(state$F),abs(template$F)))
        if(isTRUE(station$success)&&isTRUE(target_station$success)&&
           isTRUE(target_cert$point_success)&&isTRUE(comparison$point)&&tied)
          return(list(success=TRUE,reason="coalesced_stalled_history_root",
            state=state,trace=trace,stationarity=station,coalescence=comparison))
      }
      return(list(success=FALSE,reason="outer_refinement_stalled",trace=trace))
    }
  }
  list(success=FALSE,reason="outer_refinement_budget",trace=trace)
}

working_hl_refit <- function(data,native,keep_diagnostics=FALSE,history=NULL) {
  fixed<-numeric(0);fit<-native;face_evidence<-NULL
  if(is.null(history))return(list(success=FALSE,point_success=FALSE,reason="history_evidence_missing"))
  fail<-function(reason,...)list(success=FALSE,point_success=FALSE,reason=reason,fixed=fixed,
    history_evidence=if(isTRUE(keep_diagnostics))working_hl_history_compact(history$snapshot()) else
      working_hl_history_summary(history$snapshot()),...)
  if(!isTRUE(native$point_success)) {
    if(!identical(native$fit_status,"outer_boundary"))return(fail(native$fit_status))
    faces<-native$boundary_faces
    if(!length(faces)||any(!faces%in%c("1:lower","2:upper")))return(fail("unsupported_dispersion_face"))
    if("1:lower"%in%faces)fixed<-c(fixed,sigma2=1e-4)
    if("2:upper"%in%faces)fixed<-c(fixed,phi=1000)
    face_evidence<-list(faces=faces,guard_bounds=native$guard_bounds,
      guard_expansions=native$guard_expansions,
      terminal_runs=lapply(native$outer_runs,hl_fo_compact_outer_run))
    active<-c("tau","xi")[!c("sigma2","phi")%in%names(fixed)]
    fit<-hl_fo_fit(data,active=active,fixed_natural=fixed,keep_diagnostics=keep_diagnostics,
      .history=history,.history_source="fixed_working")
  }
  if(!isTRUE(fit$point_success))return(fail(fit$fit_status,fit=fit,face_evidence=face_evidence))
  selected_refinement<-working_hl_refine(working_hl_profile(fit,fit$rho_full,fit$active,
    history,"selected_initial"),fit,history,"selected_refinement")
  if(!isTRUE(selected_refinement$success))return(fail(selected_refinement$reason,
    fit=fit,selected_refinement=selected_refinement))
  fit[names(selected_refinement$state)]<-selected_refinement$state
  fit<-hl_fo_finalize_state(fit)
  if(!isTRUE(fit$point_success))return(fail("refined_point_not_certified",fit=fit))
  # Every retained certified candidate is cold-replayed. A duplicate label is
  # insufficient: complete roots must coalesce, including nuisance and covariance.
  competitors<-list();roots<-fit$candidates;competitor_covariance<-TRUE
  if(length(roots))for(j in seq_along(roots)) {
    candidate<-roots[[j]]
    st<-working_hl_profile(fit,candidate$rho_full,fit$active,history,paste0("candidate",j))
    polished<-working_hl_refine(st,fit,history,paste0("candidate_refinement",j))
    if(!isTRUE(polished$success))return(fail("competitor_cold_refinement_failed",
      fit=fit,competitor_refinement=polished))
    st<-polished$state
    cert<-hl_fo_finalize_state(st)
    agreement<-max(hl_fo_sdisc(st$beta,fit$beta),hl_fo_sdisc(st$rho_full,fit$rho_full),
      hl_fo_sdisc(st$eta,fit$eta),hl_fo_sdisc(st$augmented$A,fit$augmented$A),
      hl_fo_sdisc(st$augmented$J,fit$augmented$J))
    tied<-abs(st$F-fit$F)<=1e-10*(1+max(abs(st$F),abs(fit$F)))
    supported<-isTRUE(cert$point_success)&&st$F<=fit$F+1e-10*(1+abs(fit$F))&&(!tied||agreement<=1e-8)
    covariance_agreement<-if(tied&&isTRUE(cert$covariance_success)&&isTRUE(fit$covariance_success))
      working_covariance_discrepancy(cert$covariance$V_beta_rho,fit$covariance$V_beta_rho) else Inf
    covariance_supported<-!tied||covariance_agreement<=1e-8
    competitor_covariance<-competitor_covariance&&covariance_supported
    competitors[[j]]<-list(covariance_supported=covariance_supported,
      covariance_agreement=covariance_agreement,refinement=polished$trace,original=hl_fo_compact_state(candidate),
      replay=hl_fo_compact_state(st),point_certified=cert$point_success,
      agreement=agreement,tied=tied,success=supported)
    if(!supported)return(fail("unresolved_recovered_competitor",fit=fit,competitors=competitors))
  }
  # Compare every retained terminal endpoint after cold inner profiling under
  # the selected working mask. Raw native F values from another mask are not
  # comparable. Keep a ledger entry for every attempt, including failed ones.
  runs<-c(native$outer_runs,unlist(lapply(native$boundary_rounds,`[[`,"runs"),recursive=FALSE),
          fit$outer_runs)
  terminal_audit<-list();cache<-list()
  for(j in seq_along(runs)) {
    run<-runs[[j]];endpoint<-run$state
    if(is.null(endpoint))endpoint<-run$inner_failure
    rho<-endpoint$rho_full
    if(length(rho)!=2L||any(!is.finite(rho)))
      return(fail("retained_endpoint_missing",fit=fit,terminal_audit=terminal_audit))
    if(length(fixed))rho[match(names(fixed),c("sigma2","phi"))]<-log(fixed)
    key<-paste(format(rho,digits=17),collapse=":")
    replay<-cache[[key]]
    if(is.null(replay)) {
      replay<-tryCatch(working_hl_profile(fit,rho,fit$active),error=function(e)e)
      cache[[key]]<-replay
    }
    if(inherits(replay,"error")) {
      terminal_audit[[j]]<-list(original=hl_fo_compact_outer_run(run),
        rho=rho,success=FALSE,error=conditionMessage(replay))
      return(fail("retained_endpoint_cold_failure",fit=fit,terminal_audit=terminal_audit))
    }
    supported<-is.finite(replay$F)&&replay$F<=fit$F+1e-10*(1+abs(fit$F))
    terminal_audit[[j]]<-list(original=hl_fo_compact_outer_run(run),rho=rho,
      F=replay$F,selected_F=fit$F,delta=replay$F-fit$F,success=supported)
    if(!supported)return(fail("retained_endpoint_better_profile",fit=fit,terminal_audit=terminal_audit))
  }
  history_audit<-working_hl_history_audit(history$snapshot(),fit,keep_details=keep_diagnostics)
  if(!isTRUE(history_audit$success))return(fail(history_audit$reason,fit=fit,history_audit=history_audit,
    terminal_audit=terminal_audit,competitors=competitors))
  competitor_covariance<-competitor_covariance&&history_audit$covariance
  audit<-working_hl_audit(fit)
  if(!audit$point)return(fail(audit$reason,fit=fit,certificate=audit,competitors=competitors))
  if(!competitor_covariance) {
    audit$covariance<-FALSE;audit$reason<-"competing_covariance_not_certified"
  }
  boundary<-NULL
  if(length(fixed)) {
    profile<-working_hl_profile(fit,fit$rho_full)
    if(max(hl_fo_sdisc(profile$beta,fit$beta),hl_fo_sdisc(profile$F,fit$F))>1e-8)
      return(fail("held_coordinate_profile_mismatch",fit=fit))
    held<-match(names(fixed),c("sigma2","phi"));boundary<-vector("list",length(held))
    for(h in seq_along(held)) {
      j<-held[h];direction<-if(j==1L)1 else -1
      rho<-fit$rho_full;rho[j]<-rho[j]+direction*log(2)
      inward<-working_hl_profile(fit,rho)
      fd<-ln_research_fd(function(r)working_hl_profile(fit,r)$F,fit$rho_full)[j]
      k<-ncol(fit$K);U<-profile$augmented$psi_subject[2L*k+j,]
      slope<-direction*profile$g[j]/(1+sqrt(sum(U^2)))
      boundary[[h]]<-list(coordinate=j,outward=slope,inward_gain=inward$F-profile$F,
        derivative_fd=fd,derivative_error=abs(fd-profile$g[j])/(1+abs(profile$g[j])))
      if(slope>1e-6||inward$F-profile$F>1e-10*(1+abs(profile$F))||boundary[[h]]$derivative_error>1e-5)
        return(fail("working_mask_direction_not_supported",fit=fit,boundary=boundary))
    }
  }
  list(success=TRUE,point_success=TRUE,inference_available=audit$covariance,
       model_inference_available=isTRUE(fit$model_covariance_success),
       model_inference_reason=if(isTRUE(fit$model_covariance_success))NA_character_ else fit$model_covariance_status,
       model_Vbeta=fit$model_covariance$V_beta,
       reason=audit$reason,fit=fit,fixed=fixed,beta=fit$beta,Vbeta=audit$Vbeta,
       certificate=audit,competitors=competitors,terminal_audit=terminal_audit,history_audit=history_audit,
       selected_refinement=selected_refinement$trace,
       boundary=boundary,face_evidence=face_evidence)
}
