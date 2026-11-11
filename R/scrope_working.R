# Contrast estimability is invariant to changing predictor/reporting units.
working_estimable <- function(face, report_map) {
  d <- face$direction
  residual <- abs(as.numeric(report_map %*% d))
  bound <- 128*.Machine$double.eps*length(d)*sqrt(rowSums(report_map^2))*sqrt(sum(d^2))
  residual <= bound
}

# The common face detector uses a first-column unit intercept. Canonicalize
# only its input, then map the detected internal coordinates back. Native
# fits and their public coordinate conventions remain untouched.
working_face_detect <- function(y, X, id) {
  p <- ncol(X)
  constant <- vapply(seq_len(p),function(j)all(X[,j]==X[1,j]),logical(1))
  intercept <- unname(which(constant & X[1,]!=0))
  if(length(intercept)!=1L || any(constant & X[1,]==0))
    return(list(success=FALSE,reason="invalid_intercept_design"))
  permutation <- c(intercept,setdiff(seq_len(p),intercept))
  canonical <- X[,permutation,drop=FALSE];canonical[,1L] <- 1
  face <- ln_recovery_face_detect(y,canonical,id)
  face$canonicalization <- list(permutation=permutation,
    intercept_column=intercept,intercept_value=unname(X[1,intercept]))
  if(!isTRUE(face$success))return(face)
  P <- diag(p)[,permutation,drop=FALSE]
  face$columns <- permutation[face$columns]
  face$direction <- as.numeric(P %*% face$direction)
  face$beta_map <- P %*% face$beta_map
  face$public_map <- P %*% face$public_map %*% t(P)
  face$estimable <- as.logical(P %*% face$estimable)
  face
}

working_compact_audit <- function(ans) {
  # No cell arrays, derivative tensors or repeated whole native histories.
  compact_trial<-function(t) {
    fields<-c("success","reason","error","start","optimizer","trace")
    out<-t[intersect(fields,names(t))]
    if(!is.null(t$state))out$endpoint<-t$state[c("x","par","q","g")]
    if(!is.null(t$optimizer_endpoint))out$optimizer_endpoint<-t$optimizer_endpoint[c("x","par","q","g")]
    out
  }
  fields<-c("success","point_success","inference_available","reason","route","fixed","free",
    "model_inference_available","model_inference_reason","model_covariance_convention",
    "sampling_target","inference_reason","error","full_competitors","mapping","estimable","beta_report","V_report",
    "working_competitors","objective","start_agreement","boundary","face_evidence","competitors","native_candidates",
    "history_capture","history_audit","history_evidence","native_preservation","terminal_audit","selected_refinement","competitor_refinement","identity","native_regular","native_selected","native_qualification","recovery_attempted","recovery_outcome","recovery_reason","recovery_error","geometry_source","elapsed_seconds","input_identity","point_failure")
  out<-ans[intersect(fields,names(ans))]
  if(!is.null(ans$face))out$face<-ans$face[intersect(c("success","reason","direction","columns",
    "estimable","original_donors","retained_donors","retained_donor_indices","singular_values",
    "positive_rank","canonicalization"),names(ans$face))]
  if(!is.null(ans$trials))out$trials<-lapply(ans$trials,compact_trial)
  if(!is.null(ans$certificates))out$certificates<-lapply(ans$certificates,function(x)
    x[intersect(c("point","covariance","reason","metrics"),names(x))])
  if(!is.null(ans$native_certificate))out$native_certificate<-ans$native_certificate[
    intersect(c("point","covariance","reason","metrics","error"),names(ans$native_certificate))]
  if(!is.null(ans$certificate))out$certificate<-ans$certificate[
    intersect(c("point","covariance","reason","metrics"),names(ans$certificate))]
  if(!is.null(ans$fit)) {
    f<-ans$fit
    if(!is.null(f$theta))out$fit<-list(theta=f$theta,loglik=f$loglik,
      diagnostics=f$diagnostics,optimizer=f$optimizer)
    else out$fit<-c(hl_fo_compact_state(f),list(certification=f$certification,
      active=f$active,point_success=f$point_success,covariance_success=f$covariance_success,
      covariance_status=f$covariance_status,guard_bounds=f$guard_bounds,
      guard_expansions=f$guard_expansions,outer_runs=lapply(f$outer_runs,hl_fo_compact_outer_run),
      candidates=lapply(f$candidates,hl_fo_compact_state)))
  }
  out
}

working_recover <- function(model,y,X,id,offset,native,face,candidates=NULL,
                            keep_diagnostics=FALSE) {
  keep<-if(isTRUE(face$success))face$retained else rep(TRUE,length(y))
  J<-if(isTRUE(face$success))face$columns else seq_len(ncol(X))
  prep<-hl_prepare_data(y[keep],id[keep],X[keep,J,drop=FALSE],offset[keep],
                       cpc=0,mincp=1,verbose=FALSE)
  original<-hl_prepare_data(y,id,X,offset,cpc=0,mincp=1,verbose=FALSE)$ctx
  if(model=="HL") {
    ctx<-prep$ctx
    data<-hl_fo_make_data(as.numeric(prep$count),ctx$pred,ctx$offset,ctx$id,ctx$intcol)
    if(isTRUE(face$success)) {
      data$zero_donors<-face$original_donors-face$retained_donors
      native<-hl_fo_fit(data,keep_diagnostics=keep_diagnostics,
        .history=candidates,.history_source="nuisance_face")
    }
    ans<-working_hl_refit(data,native,keep_diagnostics,history=candidates)
    if(!isTRUE(ans$point_success))return(ans)
    mapping<-if(isTRUE(face$success))original$M_X%*%face$beta_map else original$M_X
    estimability_map<-original$M_X
  } else {
    ans<-working_pmg_refit(y[keep],X[keep,J,drop=FALSE],id[keep],offset[keep],
      native=if(isTRUE(face$success))NULL else native,
      candidates=if(isTRUE(face$success))NULL else candidates)
    if(!isTRUE(ans$point_success))return(ans)
    if(isTRUE(face$success)) {
      # The exact limiting criterion cannot lose to any retained full-data
      # point under the same dispersion mask. Audit original native endpoints.
      full<-working_pmg_problem(y,X,id,offset,ans$fixed)
      if(!length(candidates))return(list(success=FALSE,reason="full_candidate_evidence_missing"))
      audit<-lapply(candidates,function(c) {
        if(!is.null(c$error)||is.null(c$par))return(list(success=FALSE,reason="native_attempt_failed"))
        par<-c$par;p<-ncol(X)
        if(length(par)!=p+1L||any(!is.finite(par))||par[p+1L]<=0)return(list(success=FALSE,reason="invalid_native_candidate"))
        z<-c(par[seq_len(p)],log(par[p+1L]))[full$free]
        q<-full$evaluate(z)$q
        list(success=is.finite(q)&&q>=ans$objective-1e-8*(1+abs(q)),q=q,
             limiting_q=ans$objective,delta=q-ans$objective,par=par)
      })
      ans$full_competitors<-audit
      if(any(!vapply(audit,function(x)isTRUE(x$success),logical(1)))) {
        ans$success<-ans$point_success<-FALSE;ans$reason<-"limiting_root_worse_than_full_candidate"
        return(ans)
      }
      mapping<-face$public_map%*%face$beta_map;estimability_map<-face$public_map
    } else {
      scale<-prep$design_sds;scale[prep$ctx$intcol]<-1
      mapping<-diag(1/scale);estimability_map<-mapping
    }
  }
  ans$beta_report<-as.numeric(mapping%*%ans$beta)
  ans$V_report<-mapping%*%ans$Vbeta%*%t(mapping)
  ans$model_V_report<-if(!is.null(ans$model_Vbeta))mapping%*%ans$model_Vbeta%*%t(mapping) else NULL
  ans$model_covariance_convention<-scrope_model_convention(model)
  ans$estimable<-if(isTRUE(face$success))working_estimable(face,estimability_map) else rep(TRUE,ncol(X))
  ans$mapping<-mapping;ans$face<-face
  ans$sampling_target<-"regular_working_equation_target"
  if(model=="HL"&&isTRUE(face$success)&&!"sigma2"%in%names(ans$fixed)) {
    ans$inference_available<-FALSE
    ans$inference_reason<-"hl_limit_sampling_target_unresolved"
    ans$model_inference_available<-FALSE
    ans$model_inference_reason<-"hl_limit_sampling_target_unresolved"
    ans$sampling_target<-"point_only_hl_limit"
  }
  ans$route<-if(isTRUE(face$success))"identifiable_nuisance_limit" else
    if(length(ans$fixed))"fixed_positive_dispersion" else "certified_free_refit"
  if(isTRUE(face$success)&&length(ans$fixed))ans$route<-"identifiable_limit_fixed_dispersion"
  ans
}

working_assess_native <- function(model,native,prep,gene_id,candidates=NULL,
                                  keep_diagnostics=FALSE,execution=NULL) {
  began<-proc.time()[["elapsed"]];ctx<-prep$ctx;X<-ctx$pred_original
  y<-as.numeric(prep$count[,gene_id]);exposure<-exp(ctx$offset)
  # A point-qualified native HL fit already passed the unrestricted recession,
  # geometry, curvature and cold-replay checks on these same prepared data.
  # Failed HL fits and PGMM require the recession check here. Always retain
  # the separate canonical face detector: its exact structural checks can
  # detect faces missed by numerical rank/LP checks in standardized coordinates.
  reuse_geometry <- model=="HL" && isTRUE(native$point_success)
  recession <- if(reuse_geometry) list(recession=FALSE,status="native_qualified") else
    hl_fo_inner_recession(hl_fo_make_data(y,ctx$pred,ctx$offset,ctx$id,ctx$intcol),diag(ncol(X)))
  face <- working_face_detect(y,X,ctx$id)
  unresolved <- isTRUE(face$success) || identical(face$reason,"partial_donor_face_outside_scope") ||
    isTRUE(recession$recession) || identical(recession$status,"nonattainment_check_failure")
  qualification <- working_native_qualification(model,native,ctx,unresolved)
  # Retain the already observed geometry refusal, without rerunning a check.
  if(unresolved)qualification$reason<-if(isTRUE(face$success))"native_nuisance_limit" else
    if(identical(face$reason,"partial_donor_face_outside_scope"))face$reason else
      if(isTRUE(recession$recession))"native_recession_detected" else recession$status
  regular <- isTRUE(qualification$complete)
  attempted <- !regular && (!unresolved || isTRUE(face$success))
  ans <- if(regular) list(success=TRUE,point_success=TRUE,inference_available=TRUE,
      route="native",reason="preserved_native") else
    if(!attempted) list(success=FALSE,point_success=FALSE,reason="unsupported_or_unverified_nuisance_limit") else
    tryCatch(working_recover(model,y,X,ctx$id,exposure,native,face,candidates,keep_diagnostics),
      error=function(e) {
        out<-list(success=FALSE,point_success=FALSE,reason="recovery_evaluation_error",error=conditionMessage(e))
        if(model=="HL"&&!is.null(candidates))out$history_evidence<-working_hl_history_compact(candidates$snapshot())
        out
      })
  recovery_reason <- if(regular) NA_character_ else ans$reason
  recovery_point <- isTRUE(ans$point_success)
  recovery_error <- if(regular)NULL else ans$error
  selected_native <- regular
  if(!regular && isTRUE(qualification$point) && !working_recovery_complete(ans)) {
    # Recovery has not supplied complete inference: retain the qualified native
    # point and whichever native inference is supported, as one coherent fit.
    attempt <- ans
    ans <- list(success=TRUE,point_success=TRUE,
      inference_available=qualification$robust_complete,route="native_partial",
      reason="native_partial_recovery_incomplete",
      recovery_attempt=if(keep_diagnostics)attempt else NULL)
    selected_native <- TRUE
  }
  if(selected_native) {
    ans$model_inference_available<-isTRUE(qualification$model_complete)
    ans$model_inference_reason<-if(isTRUE(qualification$model_complete))NA_character_ else "native_model_inference_unavailable"
    ans$model_covariance_convention<-scrope_model_convention(model)
  }
  ans$native_qualification <- qualification
  ans$native_selected <- selected_native
  ans$recovery_attempted <- attempted
  ans$recovery_outcome <- scrope_recovery_outcome(attempted,selected_native,recovery_point,regular)
  ans$recovery_reason <- recovery_reason
  ans$recovery_error <- recovery_error
  ans$geometry_source <- if(reuse_geometry)"native_hl_qualification" else "unrestricted_recession_check"
  if(model=="HL"&&!is.null(candidates)) {
    snapshot<-candidates$snapshot()
    ans$history_capture<-c(snapshot[c("schema","observed","retained","bytes","overflow","limits")],
      list(policy=if(regular)"native_preservation_bypass" else "recovery_history_audit"))
  }
  ans$native_regular<-regular;ans$face<-face
  if(!isTRUE(ans$success))ans$route<-"withheld"
  ans$elapsed_seconds<-proc.time()[["elapsed"]]-began
  ans$input_identity<-if(attempted || keep_diagnostics)
    .hl_score_fingerprint(list(y=y,X=X,id=ctx$id,offset=exposure,gene_id=gene_id,model=model)) else NULL
  ans$identity<-execution
  ans$dimensions<-c(cells=length(y),donors=length(unique(ctx$id)),coefficients=ncol(X))
  if(keep_diagnostics)ans else working_minimal_record(ans)
}

working_observer <- function(model,keep_diagnostics=FALSE) {
  force(model);force(keep_diagnostics)
  execution<-list(version=as.character(getNamespaceVersion("scRoPE")),
    R=R.version.string,platform=R.version$platform,backend=.hl_score_execution_identity(),
    policy="native_once_recovery_20261110")
  function(native,prep,gene_id,candidates) {
    working_assess_native(model,native,prep,gene_id,candidates,keep_diagnostics,execution)
  }
}

#' Explicit working-target HL and PGMM recovery
#'
#' Runs the requested native estimator first. Its result is retained in `primary`
#' only with `keep_diagnostics=TRUE`.
#' An explicit call opts into separately labelled coefficient estimates from
#' certified positive-dispersion working fits or identifiable limiting equations.
#' Ordinary `scrope()` and `scrope_hl()` now use these eligible working routes
#' automatically; `recovery="none"` selects their native-only behavior. This
#' interface remains available for its explicit long-format working output.
#' @inheritParams scrope
#' @param model Either `"HL"` or `"PGMM"`; no switching between estimators.
#' @param keep_diagnostics Retain full native and recovery objects when `TRUE`.
#'   The default `FALSE` retains compact acceptance/failure summaries and a
#'   bounded failure witness, without full event ledgers, donor arrays or traces. Computation
#'   and acceptance are unchanged. Native fitting is sequential in this interface.
#' @details Recovery uses the default dispersion guards and fixed donor parameter
#'   `sigma2=1e-4` and, for eligible HL cell boundaries, NB precision `phi=1000`.
#'   These positive values define working models, not zero-variance boundary laws.
#'   A nuisance limit must be one-dimensional and remove only whole all-zero
#'   donors. Nonestimable coefficients remain missing. Free-donor-dispersion HL
#'   limits return points only, without SEs or p-values. Recovered inference is
#'   model and donor-robust coefficient Wald from the same selected fit. PGMM
#'   uses inverse free marginal information; HL uses inverse profile h-likelihood
#'   coefficient information at fitted dispersions. Optional tests, arbitrary
#'   contrasts and random effects are not supplied for recovery.
#'
#'   Outcomes select the route. Working-target equality to the scientific target,
#'   selective coverage, finite-donor calibration and FDR control are not
#'   established by these computational checks. No finite-cluster multiplier or
#'   pseudoinverse is used. Successful regular native fits are retained exactly.
#' @return A list with long coefficient `summary`, per-gene `routes`, untouched
#'   native `primary` only with diagnostics (otherwise `NULL`), a compact per-gene
#'   `audit`, optional full `diagnostics`,
#'   and `provenance`. Point, covariance and Wald availability are separate.
#'   `native_available` is the original reported Wald flag; `native_regular`
#'   records complete native qualification; `routes$native_selected` also
#'   includes selected partial native fits. Use `available` together with `scope` and
#'   `sampling_target` for the selected output. A reported coefficient variance
#'   does not establish joint covariance rank or calibration.
#' @seealso [scrope()], [scrope_hl()]. The installed guide is available at
#'   `system.file("doc", "hl-pgmm-working-recovery.md", package="scRoPE")`.
#' @export
scrope_working <- function(count,id,pred=NULL,offset=NULL,model=c("HL","PGMM"),
                           cpc=.005,mincp=5,keep_diagnostics=FALSE,verbose=TRUE) {
  model<-match.arg(model)
  if(!is.logical(keep_diagnostics)||length(keep_diagnostics)!=1L||is.na(keep_diagnostics))
    stop("keep_diagnostics must be TRUE or FALSE.")
  args<-list(count=count,id=id,pred=pred,offset=offset,cpc=cpc,mincp=mincp,
             keep_diagnostics=keep_diagnostics,verbose=verbose,
             .fit_observer=working_observer(model,keep_diagnostics))
  primary<-if(model=="HL")do.call(hl_fo_public_scrope_hl,args) else do.call(scrope_pmg_internal,args)
  records<-attr(primary,".working_records");attr(primary,".working_records")<-NULL
  if(model=="PGMM")primary$provenance<-pmg_public_provenance(pmg_resolve_control())
  result<-working_selected_output(primary,records,model,keep_diagnostics)
  scrope_working_notice(result,model)
  result
}

working_selected_output <- function(primary,records,model,keep_diagnostics=FALSE) {
  audits<-lapply(records,function(a)if(keep_diagnostics)working_compact_audit(a) else a$audit)
  full<-if(keep_diagnostics)records else NULL
  labels<-sub("^logFC_","",grep("^logFC_",names(primary$summary),value=TRUE))
  rows<-routes<-list()
  for(j in seq_len(nrow(primary$summary))) {
    s<-primary$summary[j,,drop=FALSE];a<-records[[as.character(s$gene_id)]]
    regular<-working_native_selected(a);recovered<-!regular&&isTRUE(a$point_success)
    scope<-if(regular)"native_estimator" else if(recovered)"data_selected_working_target" else "unavailable"
    if(identical(a$sampling_target,"point_only_hl_limit"))scope<-"point_only_hl_limit_sampling_unresolved"
    npoint<-nwald<-nmodel<-0L
    for(k in seq_along(labels)) {
      nm<-labels[k];native_available<-isTRUE(s[[paste0("wald_success_",nm)]])
      beta<-if(regular)s[[paste0("logFC_",nm)]] else if(recovered)a$beta_report[k] else NA_real_
      se<-if(regular)s[[paste0("se_robust_",nm)]] else if(recovered)sqrt(max(0,if(!is.null(a$variance_report))a$variance_report[k] else a$V_report[k,k])) else NA_real_
      estimable<-if(regular)TRUE else recovered&&isTRUE(a$estimable[k])
      point<-estimable&&is.finite(beta)
      cov<-if(regular)is.finite(se)&&se>0 else point&&isTRUE(a$inference_available)&&is.finite(se)&&se>0
      ok<-point&&cov&&if(regular)native_available else TRUE
      pv<-if(regular)s[[paste0("p_robust_",nm)]] else if(ok)2*pnorm(-abs(beta/se)) else NA_real_
      if(ok&&!is.finite((beta/se)^2))ok<-FALSE
      reason<-if(ok)NA_character_ else if(recovered&&!estimable)"nonestimable_coefficient" else
        if(recovered&&!is.null(a$inference_reason))a$inference_reason else
        if(regular)s[[paste0("wald_failure_",nm)]] else if(recovered) {
          # A fit-level certificate does not certify every mapped coefficient.
          # Keep numerical gates unchanged and report the actual refusal stage.
          if(!point)"nonfinite_coefficient_estimate" else
            if(!isTRUE(a$inference_available)) {
              if(is.null(a$reason)||identical(a$reason,"certified"))"working_covariance_unavailable" else a$reason
            } else if(!is.finite(se)||se<=0)"mapped_variance_not_positive_finite" else
              if(!is.finite((beta/se)^2))"nonfinite_wald_statistic" else "working_wald_unavailable"
        } else a$reason
      model_se<-if(regular)s[[paste0("se_",nm)]] else if(recovered) {
        variance<-if(!is.null(a$model_variance_report))a$model_variance_report[k] else
          if(!is.null(a$model_V_report))a$model_V_report[k,k] else NA_real_
        sqrt(max(0,variance))
      } else NA_real_
      model_cov<-point&&is.finite(model_se)&&model_se>0&&
        (regular||isTRUE(a$model_inference_available))
      model_ok<-model_cov&&is.finite((beta/model_se)^2)&&
        if(regular)isTRUE(s[[paste0("model_wald_success_",nm)]]) else TRUE
      model_reason<-if(model_ok)NA_character_ else if(recovered&&!estimable)"nonestimable_coefficient" else
        if(regular)s[[paste0("model_wald_failure_",nm)]] else
          if(!point)"nonfinite_coefficient_estimate" else
            if(!isTRUE(a$model_inference_available)) {
              if(length(a$model_inference_reason)==1L&&!is.na(a$model_inference_reason))a$model_inference_reason else "model_covariance_unavailable"
            } else if(!is.finite(model_se)||model_se<=0)"mapped_model_variance_not_positive_finite" else
              if(!is.finite((beta/model_se)^2))"nonfinite_model_wald_statistic" else "model_wald_unavailable"
      rows[[length(rows)+1L]]<-data.frame(gene_id=s$gene_id,gene=s$gene,coefficient=nm,
        beta=if(point)beta else NA_real_,robust_se=if(cov)se else NA_real_,p_value=if(ok)pv else NA_real_,
        native_available=native_available,native_regular=isTRUE(a$native_regular),point_available=point,
        covariance_available=cov,available=ok,estimable=estimable,route=a$route,scope=scope,
        model_covariance_available=model_cov,model_available=model_ok,
        model_se=if(model_cov)model_se else NA_real_,
        model_p_value=if(model_ok)if(regular)s[[paste0("p_naive_",nm)]] else
          pchisq((beta/model_se)^2,1,lower.tail=FALSE) else NA_real_,
        model_failure=model_reason,model_covariance_convention=scrope_model_convention(model),
        sampling_target=if(regular)"native_estimator" else if(is.null(a$sampling_target))"unavailable" else a$sampling_target,
        selection_adjusted=FALSE,failure=reason,stringsAsFactors=FALSE)
      npoint<-npoint+point;nwald<-nwald+ok;nmodel<-nmodel+model_ok
    }
    fixed<-a$fixed
    routes[[j]]<-data.frame(gene_id=s$gene_id,gene=s$gene,route=a$route,native_regular=isTRUE(a$native_regular),
      native_selected=regular,native_complete=isTRUE(a$native_qualification$complete),
      recovery_attempted=isTRUE(a$recovery_attempted),recovery_outcome=if(is.null(a$recovery_outcome))"not_recorded" else a$recovery_outcome,recovery_reason=if(is.null(a$recovery_reason))NA_character_ else a$recovery_reason,
      recovery_error=if(is.null(a$recovery_error))NA_character_ else a$recovery_error,
      point_coefficients=npoint,wald_coefficients=nwald,model_wald_coefficients=nmodel,total_coefficients=length(labels),
      fixed_donor_parameter=if("sigma2"%in%names(fixed))unname(fixed["sigma2"]) else NA_real_,
      fixed_nb_precision=if("phi"%in%names(fixed))unname(fixed["phi"]) else NA_real_,
      inferential_scope=scope,reason=a$reason,recovery_seconds=a$elapsed_seconds,
      selection_adjusted=FALSE,finite_cluster_multiplier=1,stringsAsFactors=FALSE)
  }
  native_status<-scrope_native_snapshot(primary,model,records)
  native_coefficients<-native_status$coefficients[,setdiff(names(native_status$coefficients),"native_available"),drop=FALSE]
  list(summary=scrope_reporting_join(do.call(rbind,rows),native_coefficients,c("gene_id","coefficient")),
       routes=scrope_reporting_join(do.call(rbind,routes),native_status$genes,"gene_id"),primary=if(keep_diagnostics)primary else NULL,audit=audits,
       diagnostics=if(keep_diagnostics)full else NULL,
       provenance=list(model=model,policy="native_once_recovery_20261110",execution=if(length(records))records[[1]]$identity else NULL,
         coefficient_coordinates=if(model=="HL")"original_design" else "centered_unit_intercept_original_slope_units",
         model_covariance_convention=scrope_model_convention(model),
         recovery_capabilities=c(coefficient_robust_wald=TRUE,coefficient_model_wald=TRUE,score=FALSE,lrt=FALSE),
         limitations=c("working_targets","outcome_selected_routes","no_selective_coverage_guarantee",
                       "no_general_finite_sample_calibration","no_scientific_target_equality_claim")))
}

