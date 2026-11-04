# Supported JLN-first recovery. Each promoted result retains its native evidence.
ln_recovery_route <- function(primary,y,X,id,offset,target=2L,allow_face=TRUE) {
  if(isTRUE(primary$success))return(list(route="ordinary_jln",result=primary))
  fail<-function(reason)list(route="withheld",result=list(success=FALSE,reason=reason))
  dat<-tryCatch(ln_research_context(y,X,id,offset),error=function(e)NULL)
  if(is.null(dat))return(fail("recovery_input_or_design_failure"))
  refinement<-tryCatch(ln_recovery_refine(primary,dat,target),
    error=function(e)list(success=FALSE,reason="refinement_error",error=conditionMessage(e)))
  if(isTRUE(refinement$success))return(list(route="joint_refined",result=refinement,refinement=refinement))
  ans<-ln_research_route(primary,y,X,id,offset,target,TRUE,TRUE)
  ans$refinement<-refinement
  if(isTRUE(ans$result$success)||!allow_face)return(ans)
  face<-tryCatch(ln_recovery_face_fit(y,X,id,offset),
    error=function(e)list(success=FALSE,reason="face_evaluation_error",error=conditionMessage(e)))
  # Include native recovery history and probe endpoints, reevaluated with the
  # same scalar criterion, before trusting a local reduced solution.
  if(isTRUE(face$success)){
    audit<-ln_recovery_competitor_audit(primary,dat,face$q,face$objective_error)
    face$competitor_audit<-audit
    face$known_objective_min<-audit$known_objective_min
    if(!isTRUE(audit$success)){
      face$success<-FALSE
      face$reason<-if(identical(audit$reason,"known_better_full_point"))
        "reduced_root_worse_than_known_full_candidate" else audit$reason
    }
  }
  ans$face_audit<-face
  if(!isTRUE(face$success) && isTRUE(face$face$success))ans$result$reason<-face$reason
  if(isTRUE(face$success)){ans$route<-"identifiable_nuisance_limit";ans$result<-face}
  ans
}

ln_recovery_apply <- function(result,raw,input,enabled,keep_diagnostics) {
  primary<-result;nr<-nrow(result$summary);p<-ncol(input$pred)
  coefnames<-sub("^logFC_","",grep("^logFC_",names(result$summary),value=TRUE))
  colnames(input$pred)<-coefnames
  rows<-coefrows<-audits<-vector("list",nr)
  for(j in seq_len(nr)){
    gid<-result$summary$gene_id[j];gene<-result$summary$gene[j]
    native<-list(success=isTRUE(result$provenance$per_gene$inference_available[j]),joint=raw[[j]]$ln_joint)
    ans<-if(native$success)list(route="ordinary_jln",result=native) else if(!enabled)
      list(route="withheld",result=list(success=FALSE,reason="recovery_requires_default_box_and_lbfgs")) else
      tryCatch(ln_recovery_route(native,as.numeric(input$count[gid,]),input$pred,input$id,input$offset,
                        target=if(p>1L)2L else 1L),error=function(e)
        list(route="withheld",result=list(success=FALSE,reason="recovery_evaluation_error",error=conditionMessage(e))))
    f<-ans$result;added<-isTRUE(f$success)&&!native$success;route<-ans$route
    estimable<-rep(isTRUE(f$success),p);fixed<-list();q<-NA_real_
    scope<-if(route%in%c("ordinary_jln","joint_refined"))"regular_joint_JLN_target" else
      if(route=="ordinary_jln_with_search_warning")"intended_interior_branch_unresolved_search" else
        if(route=="identifiable_nuisance_limit")"identifiable_limiting_face_working_target" else
          if(startsWith(route,"fixed_"))"fixed_mask_epsilon_working_target" else "unavailable"
    caution<-switch(route,
      ordinary_jln_with_search_warning="Dominated boundary-start trajectories remain unresolved; unvisited competing branches are not excluded.",
      identifiable_nuisance_limit="Only identifiable contrasts are reported; a data-selected limiting face has no selective coverage guarantee.",
      fixed_v="Fixed-dispersion working target; scientific-target equality and selective coverage are not established.",
      fixed_phi="Fixed-dispersion working target; scientific-target equality and selective coverage are not established.",
      fixed_v_phi="Fixed-dispersion working target; scientific-target equality and selective coverage are not established.","")
    if(added){
      beta<-f$beta;V<-f$Vbeta
      if(route=="identifiable_nuisance_limit"){
        estimable<-f$face$estimable;full<-f$reduced_full;q<-f$q
        fixed<-f$reduced$result$fixed
      }else{full<-f$state$full;q<-f$state$q;fixed<-f$fixed}
      se<-sqrt(diag(V));ok<-estimable&is.finite(beta)&is.finite(se)&se>0
      beta[!ok]<-NA_real_;se[!ok]<-NA_real_
      for(k in seq_len(p)){
        nm<-coefnames[k]
        result$summary[j,paste0("logFC_",nm)]<-beta[k]
        result$summary[j,paste0("se_robust_",nm)]<-se[k]
        result$summary[j,paste0("p_robust_",nm)]<-if(ok[k])2*pnorm(-abs(beta[k]/se[k]))else NA_real_
        result$summary[j,paste0("wald_stat_",nm)]<-if(ok[k])(beta[k]/se[k])^2 else NA_real_
        result$summary[j,paste0("wald_success_",nm)]<-ok[k]
        result$summary[j,paste0("wald_failure_",nm)]<-if(ok[k])NA_character_ else "nonestimable_limiting_coefficient"
        result$summary[j,paste0("se_",nm)]<-NA_real_
        result$summary[j,paste0("p_naive_",nm)]<-NA_real_
      }
      # Optional tests were computed for the native estimator, not the recovery.
      optional<-grep("^(p_)?(score|lrt)_",names(result$summary),value=TRUE)
      for(nm in optional){
        result$summary[j,nm]<-if(grepl("success",nm))FALSE else if(is.numeric(result$summary[[nm]]))NA_real_ else
          if(grepl("failure|status|method",nm))"not_supported_for_recovery" else NA
      }
      result$summary$test_status[j]<-if(all(ok))"recovery_wald_available" else "recovery_wald_partially_available"
      result$summary$test_failure_detail[j]<-"Recovery supports coefficient Wald tests only; model-based SE, Score/LRT and random effects are unavailable."
      ts<-which(result$test_status$gene_id==gid)
      for(ti in ts){
        iswald<-result$test_status$test[ti]=="wald";k<-match(result$test_status$contrast[ti],coefnames)
        success<-iswald&&!is.na(k)&&ok[k]
        result$test_status$success[ti]<-success
        result$test_status$failure_reason[ti]<-if(success)NA_character_ else if(iswald)"nonestimable_limiting_coefficient" else "not_supported_for_recovery"
        result$test_status$status[ti]<-if(success)"available_with_route_limitations" else "unavailable"
        result$test_status$inference[ti]<-scope
        result$test_status$method[ti]<-if(iswald)"matched_donor_robust_wald" else "not_supported_for_recovery"
        result$test_status$failure_stage[ti]<-if(success)NA_character_ else
          if(iswald)"coefficient_estimability" else "recovery_test_scope"
      }
      result$overdispersion$Subject[j]<-exp(tail(full,2)[1]);result$overdispersion$Cell[j]<-exp(-tail(full,2)[2])
      result$convergence[j]<-1 # Certified free root; see route for full-coordinate estimability.
      if(!is.null(result$random_effect)&&nrow(result$random_effect)==nr)result$random_effect[j,]<-NA_real_
      if(!is.null(result$lrt_details) && "gene_id" %in% names(result$lrt_details)){
        di<-which(result$lrt_details$gene_id==gid)
        for(nm in setdiff(names(result$lrt_details),c("gene_id","gene","contrast"))){
          col<-result$lrt_details[[nm]]
          result$lrt_details[di,nm]<-if(is.logical(col))FALSE else if(is.numeric(col))NA_real_ else "not_supported_for_recovery"
        }
      }
      fullroot<-route%in%c("joint_refined","ordinary_jln_with_search_warning")
      for(nm in c("provenance","availability")){
        tab<-if(nm=="provenance")result$provenance$per_gene else result$availability
        tab$estimator[j]<-if(fullroot)"joint_ln" else if(route=="identifiable_nuisance_limit")"identifiable_limiting_joint_ln" else "fixed_dispersion_joint_ln"
        tab$certified_root[j]<-fullroot
        tab$inference_status[j]<-if(all(ok))"available_with_route_limitations" else "partially_available"
        tab$enabled_inference_status[j]<-tab$inference_status[j]
        tab$resolution_status[j]<-route;tab$inference_source[j]<-route
        tab$newly_enabled[j]<-TRUE;tab$availability_transition[j]<-"recovered"
        tab$resolution_policy_version[j]<-"jln_public_recovery_20261103"
        if(nm=="provenance"){
          tab$mode_method[j]<-"none";tab$mode_status[j]<-"not_supported_for_recovery"
          tab$mode_max_score_residual[j]<-tab$mode_max_newton_displacement[j]<-NA_real_
          tab$inference_available[j]<-any(ok);tab$failure_reason[j]<-if(all(ok))NA_character_ else "some_coefficients_nonestimable"
          tab$enabled_failure_reason[j]<-tab$failure_reason[j]
          tab$resolution_mode[j]<-"automatic_recovery"
          tab$resolution_veto_classes[j]<-NA_character_
          tab$route[j]<-tab$route_code[j]<-tab$enabled_route_code[j]<-route
          result$provenance$per_gene<-tab
        }else{
          tab$availability[j]<-if(all(ok))"inference" else "partial_inference"
          tab$reason_group[j]<-"supported_recovery";tab$status_message[j]<-scope
          tab$inference_supported[j]<-any(ok);tab$estimate_source[j]<-route;tab$attention_flag[j]<-TRUE
          tab$technical_code[j]<-route;result$availability<-tab
        }
      }
      result$provenance$mode_diagnostics[[j]]<-list(status="not_supported_for_recovery")
      if(keep_diagnostics){
        # The native raw covariance uses standardized coefficient coordinates.
        # Elementwise scaling preserves the estimable block when a face has NA
        # rows/columns (matrix products would propagate those NAs throughout).
        scale<-diag(result$provenance$coefficient_map_internal_to_public)
        Vinternal<-V/outer(scale,scale)
        result$diagnostics$loglik[j]<- -q
        result$diagnostics$raw[[j]]<-list(route=route,inference_available=any(ok),
          certified_full_root=fullroot,recovery=f,robust_cov=Vinternal,loglik= -q,
          failure_reason=if(all(ok))NA_character_ else "some_coefficients_nonestimable")
        result$diagnostics$godambe[[j]]<-list(robust_cov=Vinternal,
          robust_cov_public=V,robust_cov_coordinates="standardized_design",
          robust_cov_public_coordinates="centered_public_design",route=route)
      }
    }else{
      ok<-vapply(coefnames,function(nm)isTRUE(result$summary[[paste0("wald_success_",nm)]][j]),logical(1))
    }
    if(keep_diagnostics){
      scale<-diag(result$provenance$coefficient_map_internal_to_public)
      result$diagnostics$raw[[j]]$robust_cov_public<-
        result$diagnostics$raw[[j]]$robust_cov*outer(scale,scale)
      result$diagnostics$raw[[j]]$robust_cov_coordinates<-"standardized_design"
      result$diagnostics$raw[[j]]$robust_cov_public_coordinates<-"centered_public_design"
    }
    rows[[j]]<-data.frame(gene_id=gid,gene=gene,route=route,native_available=native$success,
      added=added,available_coefficients=sum(ok),total_coefficients=p,
      inferential_scope=scope,selection_adjusted=FALSE,finite_cluster_multiplier=1,
      fixed_donor_variance=if(is.null(fixed$v))NA_real_ else fixed$v,
      fixed_cell_overdispersion=if(is.null(fixed$phi))NA_real_ else fixed$phi,
      reason=if(isTRUE(f$success))scope else f$reason,warning=caution,stringsAsFactors=FALSE)
    coefrows[[j]]<-data.frame(gene_id=gid,gene=gene,coefficient=coefnames,
      available=ok,estimable=if(isTRUE(f$success))estimable else NA,route=route,stringsAsFactors=FALSE)
    if(keep_diagnostics)audits[[j]]<-ans
  }
  routes<-do.call(rbind,rows)
  if(any(routes$added)){
    result$provenance$estimator<-"jln_first_with_recovery"
    result$provenance$model_covariance<-"native_JLN_only"
    result$provenance$robust_covariance<-"route_matched_free_parameter_donor_sandwich"
    result$provenance$random_effect<-"native_JLN_only_when_requested"
    result$provenance$inference_policy<-"route_specific_coefficient_Wald_with_recorded_limitations"
  }
  result$provenance$recovery_policy<-"jln_public_recovery_20261103"
  result$recovery<-list(policy="jln_public_recovery_20261103",mode="auto",routes=routes,
    coefficients=do.call(rbind,coefrows),primary=primary,
    diagnostics=if(keep_diagnostics)audits else NULL,
    limitations=paste("Fixed-mask and limiting-face sandwiches concern regular working targets.",
      "Route/face selection, scientific-target equality and FDR calibration are not established.",
      "Finite-donor undercoverage has been observed; no finite-cluster correction is applied."))
  if(keep_diagnostics)result$diagnostics$test_status<-result$test_status
  if(any(routes$added))warning(sprintf(
    "JLN recovery: %d/%d genes have all coefficient Wald results, %d partial, %d unavailable; %d recovered. Inspect recovery$routes and recovery$coefficients for target and selection limitations.",
    sum(routes$available_coefficients==p),nr,sum(routes$available_coefficients>0&routes$available_coefficients<p),
    sum(routes$available_coefficients==0),sum(routes$added)),call.=FALSE) else
      ln_availability_warn(result$availability, operational_fallback_requested=FALSE)
  result
}
