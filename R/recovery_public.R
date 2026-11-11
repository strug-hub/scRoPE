# Public recovery reporting is separate from the method-specific acceptance
# rules. Compact records contain no evaluators, donor arrays or search ledgers.
working_minimal_record <- function(ans) {
  fields <- c("success", "point_success", "inference_available", "route", "reason",
    "fixed", "free", "estimable", "beta_report", "sampling_target", "inference_reason",
    "model_inference_available", "model_inference_reason", "model_covariance_convention",
    "native_regular", "native_selected", "native_qualification", "recovery_attempted", "recovery_outcome", "recovery_reason", "recovery_error", "geometry_source", "elapsed_seconds", "input_identity", "identity", "dimensions")
  out <- ans[intersect(fields, names(ans))]
  if (!is.null(ans$V_report)) out$variance_report <- diag(ans$V_report)
  if (!is.null(ans$model_V_report)) out$model_variance_report <- diag(ans$model_V_report)
  out$dispersion <- working_report_dispersion(ans)
  small_certificate <- function(x) {
    if (is.null(x)) return(NULL)
    x[intersect(c("point", "covariance", "reason", "metrics", "error"), names(x))]
  }
  h <- ans$history_audit
  if (!is.null(h)) {
    # Full rows may be present in direct callers; the public compact contract
    # always removes them, including failure evidence and pending-event lists.
    h$observations <- h$history_evidence <- NULL
    h$pending_events <- if (!is.null(h[["pending_events"]])) h[["pending_events"]] else length(h$pending_event_ids)
    h$pending_event_ids <- NULL
  }
  history <- ans$history_evidence
  if (!is.null(history) && !is.null(history$events)) {
    history <- history[intersect(c("schema", "observed", "retained", "bytes",
      "overflow", "limits", "execution", "byte_accounting"), names(history))]
  }
  out$audit <- list(success=isTRUE(ans$success), point_success=isTRUE(ans$point_success),
    inference_available=isTRUE(ans$inference_available), reason=ans$reason,
    model_inference_available=isTRUE(ans$model_inference_available),
    model_inference_reason=ans$model_inference_reason,
    native_regular=ans$native_regular, native_qualification=ans$native_qualification,
    native_selected=ans$native_selected,recovery_attempted=ans$recovery_attempted,recovery_outcome=ans$recovery_outcome,
    recovery_reason=ans$recovery_reason,recovery_error=ans$recovery_error,geometry_source=ans$geometry_source,
    native_certificate=small_certificate(ans$native_certificate),
    certificate=small_certificate(ans$certificate),
    certificates=lapply(ans$certificates, small_certificate),
    history_capture=ans$history_capture, history_audit=h, history_evidence=history,
    identity=ans$identity, input_identity=ans$input_identity, dimensions=ans$dimensions,
    retained_details=FALSE)
  out
}

working_report_dispersion <- function(ans) {
  if (!isTRUE(ans$point_success)) return(NULL)
  if (!is.null(ans$par)) return(c(Subject=tail(ans$par, 1L)))
  if (length(ans$fit$rho_full) == 2L)
    return(c(Subject=exp(ans$fit$rho_full[1L]), Cell=exp(-ans$fit$rho_full[2L])))
  NULL
}

scrope_recovery_policy <- function(model, recovery, ln_recovery,
                                   recovery_missing, ln_missing) {
  recovery <- match.arg(recovery, c("auto", "none"))
  ln_recovery <- match.arg(ln_recovery, c("auto", "none"))
  if (identical(model, "LN")) {
    if (!recovery_missing && !ln_missing && !identical(recovery, ln_recovery))
      stop("recovery and ln_recovery must agree when both are supplied for model='LN'.")
    if (recovery_missing && !ln_missing) recovery <- ln_recovery
  }
  recovery
}

scrope_recovery_default_box <- function(min, max) {
  identical(as.numeric(min), c(1e-4,1e-4)) &&
    identical(as.numeric(max), c(10,1000))
}

scrope_complete_working <- function(result,model,mode,eligible,keep_diagnostics,controls) {
  records<-attr(result,".working_records");attr(result,".working_records")<-NULL
  if(!is.null(records))result<-scrope_apply_working(result,records,model,keep_diagnostics)
  scrope_finish_recovery(result,model,mode,eligible,
    skip_reason=if(!eligible)"recovery_requires_default_box" else NA_character_,
    keep_diagnostics=keep_diagnostics,controls=controls)
}

scrope_clear_selected_tests <- function(result, row, reason) {
  s <- result$summary
  fields <- grep("^(se_|p_|wald_|model_wald_|score_|hl_score_|lrt_|hl_lrt_)",
                 names(s), value=TRUE)
  for (name in fields) {
    value <- s[[name]]
    s[row,name] <- if (grepl("success",name)) FALSE else
      if (is.numeric(value)) NA_real_ else if (is.logical(value)) FALSE else
        if (grepl("failure|status|method",name)) reason else NA_character_
  }
  result$summary <- s
  for (name in c("lrt_details","hl_adj_details")) {
    d <- result[[name]]
    if (!is.data.frame(d) || !nrow(d)) next
    ix <- if ("gene_id" %in% names(d)) d$gene_id == s$gene_id[row] else
      d$gene == s$gene[row]
    for (field in setdiff(names(d),c("gene_id","gene","contrast","test"))) {
      value <- d[[field]]
      d[ix,field] <- if (is.logical(value)) FALSE else if (is.numeric(value)) NA_real_ else
        if (grepl("failure|status|method|reason",field)) reason else NA_character_
    }
    result[[name]] <- d
  }
  if (!is.null(result$random_effect)) result$random_effect[row,] <- NA_real_
  result
}

scrope_apply_working <- function(result, records, model, keep_diagnostics) {
  selected <- working_selected_output(result, records, model, keep_diagnostics)
  primary <- if(keep_diagnostics)result else NULL
  attr(result,".native_reporting") <- scrope_native_snapshot(result,model,records)
  labels <- sub("^logFC_", "", grep("^logFC_",names(result$summary),value=TRUE))
  for (j in seq_len(nrow(result$summary))) {
    id <- result$summary$gene_id[j]
    a <- records[[as.character(id)]]
    if (working_native_selected(a)) next
    reason <- if (isTRUE(a$point_success)) "not_supported_for_recovery" else a$reason
    if (is.null(reason) || !length(reason)) reason <- "recovery_unavailable"
    result <- scrope_clear_selected_tests(result,j,reason)
    rows <- selected$summary[selected$summary$gene_id==id,,drop=FALSE]
    rows <- rows[match(labels,rows$coefficient),,drop=FALSE]
    for (k in seq_along(labels)) {
      nm <- labels[k];r <- rows[k,,drop=FALSE]
      result$summary[j,paste0("logFC_",nm)] <- r$beta
      result$summary[j,paste0("se_",nm)] <- r$model_se
      result$summary[j,paste0("p_naive_",nm)] <- r$model_p_value
      result$summary[j,paste0("model_wald_stat_",nm)] <- if(r$model_available)(r$beta/r$model_se)^2 else NA_real_
      result$summary[j,paste0("model_wald_success_",nm)] <- r$model_available
      result$summary[j,paste0("model_wald_failure_",nm)] <- r$model_failure
      result$summary[j,paste0("se_robust_",nm)] <- r$robust_se
      result$summary[j,paste0("p_robust_",nm)] <- r$p_value
      result$summary[j,paste0("wald_stat_",nm)] <- if (r$available) (r$beta/r$robust_se)^2 else NA_real_
      result$summary[j,paste0("wald_success_",nm)] <- r$available
      result$summary[j,paste0("wald_failure_",nm)] <- r$failure
      for (prefix in c("wald_df_","wald_contrast_rank_","wald_covariance_rank_")) {
        field <- paste0(prefix,nm)
        if (field %in% names(result$summary)) result$summary[j,field] <- if(r$available)1 else NA_real_
      }
      field <- paste0("wald_m_",nm)
      if(field %in% names(result$summary))result$summary[j,field] <- unname(a$dimensions["donors"])
    }
    if ("fit_success" %in% names(result$summary)) {
      result$summary$fit_success[j] <- any(rows$point_available)
      result$summary$fit_status[j] <- if(any(rows$point_available))a$route else "recovery_unavailable"
      result$summary$fit_failure[j] <- if(any(rows$point_available))NA_character_ else reason
      result$summary$covariance_success[j] <- any(rows$covariance_available)
      result$summary$covariance_status[j] <- if(any(rows$covariance_available))"available_with_route_limitations" else
        if(!is.null(a$inference_reason))a$inference_reason else reason
    }
    dispersion <- if(!is.null(a$dispersion))a$dispersion else working_report_dispersion(a)
    result$overdispersion[j,] <- NA_real_
    for(nm in intersect(names(dispersion),names(result$overdispersion)))
      result$overdispersion[j,nm] <- unname(dispersion[nm])
    # convergence is intentionally the unchanged native compatibility code.
    if (keep_diagnostics && !is.null(result$diagnostics)) {
      if (model=="HL") {
        result$diagnostics$hl[[j]] <- result$diagnostics$first_order[[j]] <-
          list(route=a$route,native_convergence=result$convergence[j],recovery=a)
        result$diagnostics$hl_adj[[j]] <- result$diagnostics$wald[[j]] <- list()
      } else {
        result$diagnostics$raw[[j]] <- list(route=a$route,native_convergence=result$convergence[j],recovery=a)
        if(length(result$diagnostics$loglik))result$diagnostics$loglik[j] <- NA_real_
        if(length(result$diagnostics$lrt))result$diagnostics$lrt[[j]] <- list()
      }
      reported_covariance <- if(any(rows$covariance_available))a$V_report else NULL
      if(!is.null(reported_covariance)) {
        reported_covariance[!rows$covariance_available,] <- NA_real_
        reported_covariance[,!rows$covariance_available] <- NA_real_
      }
      result$diagnostics$godambe[[j]] <- list(route=a$route,
        inference_available=any(rows$available), covariance_coordinates="reported_coefficients",
        inferential_scope="coefficient_wald_only", robust_cov=reported_covariance,
        model_cov=if(any(rows$model_covariance_available)) {
          V<-a$model_V_report;V[!rows$model_covariance_available,]<-NA_real_
          V[,!rows$model_covariance_available]<-NA_real_;V
        } else NULL,model_covariance_convention=scrope_model_convention(model))
    }
  }
  result$recovery <- list(policy="native_once_recovery_20261110",mode="auto",
    routes=selected$routes,coefficients=selected$summary,primary=primary,
    audit=selected$audit,diagnostics=selected$diagnostics,
    limitations=selected$provenance$limitations)
  if(model=="HL"&&!is.null(result$provenance$per_gene)) {
    fields<-intersect(names(result$provenance$per_gene),names(result$summary))
    result$provenance$per_gene[fields]<-result$summary[fields]
  }
  result$provenance$selected_routes<-selected$routes
  result$provenance$native_estimator<-result$provenance$estimator
  if(any(!selected$routes$native_selected)) {
    result$provenance$model_covariance<-scrope_model_convention(model)
    result$provenance$robust_covariance<-"route_specific_coefficient_donor_sandwich"
    result$provenance$inference_policy<-"native_tests_or_selected_coefficient_wald"
  }
  result
}

scrope_collect_test_status <- function(summary,model) {
  definitions <- list(c("wald_success_","wald_failure_","wald"),
    c("model_wald_success_","model_wald_failure_","wald_model"),
    c("score_success_","score_failure_","score"),c("lrt_success_","lrt_failure_","lrt"),
    c("hl_score_success_","hl_score_failure_","score"),c("hl_lrt_success_","hl_lrt_failure_","lrt"))
  blocks <- list()
  for(d in definitions)for(nm in grep(paste0("^",d[1]),names(summary),value=TRUE)) {
    label<-substring(nm,nchar(d[1])+1L);failure<-paste0(d[2],label)
    blocks[[length(blocks)+1L]]<-data.frame(gene_id=summary$gene_id,gene=summary$gene,
      contrast=label,test=d[3],success=as.logical(summary[[nm]]),
      failure_reason=if(failure %in% names(summary))as.character(summary[[failure]]) else NA_character_,
      stringsAsFactors=FALSE)
  }
  if(!length(blocks))return(data.frame(gene_id=integer(),gene=character(),contrast=character(),
    test=character(),success=logical(),failure_reason=character()))
  do.call(rbind,blocks)
}

scrope_finish_recovery <- function(result,model,mode,eligible=TRUE,skip_reason=NA_character_,
                                   keep_diagnostics=FALSE,controls=list(),notify=TRUE) {
  native_reporting<-attr(result,".native_reporting");attr(result,".native_reporting")<-NULL
  if(is.null(native_reporting)) {
    if(!is.null(result$recovery$routes)&&"native_qualification"%in%names(result$availability)) {
      native_reporting<-list(
        genes=result$availability[,c("gene_id",grep("^native_",names(result$availability),value=TRUE)),drop=FALSE],
        coefficients=result$coefficient_availability[,c("gene_id","coefficient",grep("^native_",names(result$coefficient_availability),value=TRUE)),drop=FALSE],
        tests=result$test_status[,c("gene_id","contrast","test","native_success","native_failure_reason"),drop=FALSE])
    } else native_reporting<-scrope_native_snapshot(result,model)
  }
  if(!is.null(result$recovery$primary))result$recovery$primary<-scrope_finish_recovery(
    result$recovery$primary,model,"none",eligible,skip_reason,keep_diagnostics,controls,notify=FALSE)
  s<-result$summary;labels<-sub("^logFC_","",grep("^logFC_",names(s),value=TRUE))
  model_convention<-if(model=="LN"&&identical(result$provenance$model_covariance,"legacy_post_pml_curvature"))
    "legacy_post_pml_curvature" else scrope_model_convention(model)
  if(!"gene_id"%in%names(s))stop("Recovery reporting requires original gene identifiers.")
  existing<-result$recovery
  selected<-if(!is.null(existing$coefficients)&&"point_available"%in%names(existing$coefficients)) existing$coefficients else NULL
  rows<-list();genes<-list()
  for(j in seq_len(nrow(s))) {
    id<-s$gene_id[j]
    route<-if(!is.null(existing$routes))as.character(existing$routes$route[match(id,existing$routes$gene_id)]) else "native"
    native<-route %in% c("native","native_partial","ordinary_jln")
    source<-if(native)"native" else if(route=="withheld")"unavailable" else "recovery"
    scope<-if(!is.null(existing$routes)&&"inferential_scope"%in%names(existing$routes))
      as.character(existing$routes$inferential_scope[match(id,existing$routes$gene_id)]) else "native_estimator"
    sampling<-if(native)"native_estimator" else scope
    gene_rows<-list()
    for(k in seq_along(labels)) {
      nm<-labels[k];beta<-s[[paste0("logFC_",nm)]][j]
      se<-s[[paste0("se_robust_",nm)]][j];pv<-s[[paste0("p_robust_",nm)]][j]
      ok<-isTRUE(s[[paste0("wald_success_",nm)]][j])&&is.finite(pv)
      if(!is.null(selected)) {
        r<-selected[selected$gene_id==id&selected$coefficient==nm,,drop=FALSE]
        point<-isTRUE(r$point_available);cov<-isTRUE(r$covariance_available)
        estimable<-isTRUE(r$estimable);sampling<-r$sampling_target
      } else {
        fit_ok<-if(model=="HL")isTRUE(s$fit_success[j]) else if(model=="PGMM")isTRUE(s$optimizer_success[j]) else
          isTRUE(result$provenance$per_gene$inference_available[j])
        point<-is.finite(beta)&&(fit_ok||source=="recovery")
        cov<-point&&is.finite(se)&&se>0;estimable<-point
      }
      stat<-s[[paste0("wald_stat_",nm)]][j]
      model_supported<-paste0("model_wald_success_",nm)%in%names(s)
      model_se<-s[[paste0("se_",nm)]][j]
      model_pv<-s[[paste0("p_naive_",nm)]][j]
      model_stat<-s[[paste0("model_wald_stat_",nm)]][j]
      model_cov<-point&&length(model_se)==1L&&is.finite(model_se)&&model_se>0
      model_ok<-model_supported&&isTRUE(s[[paste0("model_wald_success_",nm)]][j])&&is.finite(model_pv)
      model_reason<-s[[paste0("model_wald_failure_",nm)]][j]
      family_status<-function(available)if(available)if(native)"available" else "available_with_limitations" else
        if(point)"point_only" else "unavailable"
      status<-if(ok&&(!model_supported||model_ok))family_status(TRUE) else
        if(ok||model_ok)"partial" else family_status(FALSE)
      reason<-s[[paste0("wald_failure_",nm)]][j]
      gene_rows[[k]]<-data.frame(gene_id=id,gene=s$gene[j],coefficient=nm,model=model,
        route=route,source=source,status=status,point_available=point,covariance_available=cov,
        statistic_available=ok&&is.finite(stat),reference_available=ok,
        interval_available=point&&cov&&is.finite(se)&&se>0,available=ok,estimable=estimable,
        beta=if(point)beta else NA_real_,robust_se=if(cov)se else NA_real_,
        statistic=if(ok)stat else NA_real_,p_value=if(ok)pv else NA_real_,
        robust_status=family_status(ok),model_status=if(model_supported)family_status(model_ok) else "not_reported",
        model_supported=model_supported,model_covariance_available=model_cov,
        model_available=model_ok,model_se=if(model_cov)model_se else NA_real_,
        model_statistic=if(model_ok)model_stat else NA_real_,model_p_value=if(model_ok)model_pv else NA_real_,
        model_interval_available=model_cov,model_covariance_convention=model_convention,
        model_failure_reason=if(model_ok)NA_character_ else if(length(model_reason))model_reason else "not_reported",
        scope=scope,sampling_target=sampling,selection_adjusted=FALSE,
        failure_reason=if(ok)NA_character_ else if(length(reason))reason else "unavailable",
        failure_stage=scrope_failure_stage(reason,ok),
        model_failure_stage=scrope_failure_stage(model_reason,model_ok),
        stringsAsFactors=FALSE)
    }
    cr<-do.call(rbind,gene_rows);rows[[j]]<-cr
    np<-sum(cr$point_available);nc<-sum(cr$covariance_available);nw<-sum(cr$available);p<-nrow(cr)
    nmw<-sum(cr$model_available);nmc<-sum(cr$model_covariance_available)
    model_incomplete<-any(cr$model_supported & !cr$model_available)
    status<-if(nw==p&&!model_incomplete)if(native)"available" else "available_with_limitations" else
      if(nw>0||nmw>0)"partial" else if(np>0)"point_only" else "unavailable"
    if(identical(route,"native_partial")&&nw==p)status<-"partial"
    genes[[j]]<-data.frame(gene_id=id,gene=s$gene[j],requested_model=model,effective_model=model,
      selected_route=route,selected_source=source,status=status,recovery_mode=mode,
      recovery_eligible=eligible,recovery_skip_reason=skip_reason,
      native_convergence=if(!is.null(existing$native_convergence))existing$native_convergence[j] else
        if(!is.null(existing$primary))existing$primary$convergence[j] else result$convergence[j],
      selected_point_coefficients=np,selected_covariance_coefficients=nc,selected_wald_coefficients=nw,
      selected_model_covariance_coefficients=nmc,selected_model_wald_coefficients=nmw,
      model_covariance_convention=model_convention,
      total_coefficients=p,inferential_scope=scope,sampling_target=sampling,
      selection_adjusted=FALSE,failure_reason=if(identical(route,"native_partial") && nw==p)
        as.character(existing$routes$reason[match(id,existing$routes$gene_id)]) else
        if(native&&model_incomplete&&nw==p)"native_model_inference_unavailable" else if(nw==p&&!model_incomplete)NA_character_ else
        paste(unique(na.omit(c(cr$failure_reason,cr$model_failure_reason[cr$model_supported]))),collapse=";"),stringsAsFactors=FALSE)
  }
  availability<-scrope_reporting_join(do.call(rbind,genes),native_reporting$genes,"gene_id")
  coefficients<-scrope_reporting_join(do.call(rbind,rows),native_reporting$coefficients,c("gene_id","coefficient"))
  for(j in seq_len(nrow(availability))) {
    a<-availability[j,,drop=FALSE];r<-existing$routes
    ix<-if(!is.null(r))match(a$gene_id,r$gene_id) else NA_integer_
    value<-function(nm,default=NA_character_)if(!is.null(r)&&nm%in%names(r)&&!is.na(ix))r[[nm]][ix] else default
    attempted<-if(!is.null(r)&&"recovery_attempted"%in%names(r))isTRUE(value("recovery_attempted",FALSE)) else
      model=="LN"&&mode=="auto"&&isTRUE(eligible)&&!isTRUE(value("native_available",TRUE))
    outcome<-value("recovery_outcome")
    if(is.na(outcome)||outcome=="not_recorded"||!attempted)outcome<-if(!attempted) {
      if(a$native_qualification=="complete")"not_needed" else if(mode=="none")"disabled" else
        if(!isTRUE(eligible))"ineligible" else "not_attempted"
    } else if(a$selected_source=="recovery")"selected" else "failed"
    reason<-value("recovery_reason",value("reason"))
    if(outcome=="not_needed")reason<-NA_character_ else if(outcome=="disabled")reason<-"recovery_disabled" else
      if(outcome=="ineligible")reason<-skip_reason
    availability$recovery_attempted[j]<-attempted
    availability$recovery_outcome[j]<-outcome
    availability$recovery_reason[j]<-reason
    availability$recovery_error[j]<-value("recovery_error")
    availability$availability_message[j]<-paste0("Native qualification: ",a$native_qualification,
      "; recovery: ",outcome,"; selected source: ",a$selected_source,
      "; robust Wald: ",a$selected_wald_coefficients,"/",a$total_coefficients,
      "; model Wald: ",a$selected_model_wald_coefficients,"/",a$total_coefficients,".",
      if(!is.na(a$native_failure_reason))paste0(" Native reason: ",a$native_failure_reason,".") else "",
      if(!is.na(a$failure_reason)&&nzchar(a$failure_reason))paste0(" Selected reason: ",a$failure_reason,".") else "",
      if(!is.na(reason))paste0(" Recovery reason: ",reason,".") else "")
  }
  # Preserve existing LN availability columns/enumerations; common fields are
  # additive and have identical names and meanings across the three models.
  if(is.data.frame(result$availability)) {
    old<-result$availability[match(availability$gene_id,result$availability$gene_id),,drop=FALSE]
    for(nm in names(availability))old[[nm]]<-availability[[nm]]
    availability<-old
  }
  result$availability<-availability
  result$coefficient_availability<-coefficients
  if(model!="LN")result$test_status<-scrope_collect_test_status(s,model) else {
    # Keep LN optional-test metadata, adding the selected model family in the
    # same schema without relabelling native Score/LRT results as recovered.
    tab<-result$test_status
    if(is.null(tab))tab<-scrope_collect_test_status(s,model)
    tab<-tab[tab$test!="wald_model",,drop=FALSE]
    model_tab<-scrope_collect_test_status(s,model)
    model_tab<-model_tab[model_tab$test=="wald_model",,drop=FALSE]
    for(nm in setdiff(names(tab),names(model_tab)))model_tab[[nm]]<-rep(NA,nrow(model_tab))
    for(nm in setdiff(names(model_tab),names(tab)))tab[[nm]]<-rep(NA,nrow(tab))
    if(nrow(model_tab)) {
      if("status"%in%names(model_tab))model_tab$status<-ifelse(model_tab$success,"available","unavailable")
      if("failure_stage"%in%names(model_tab))model_tab$failure_stage<-ifelse(model_tab$success,NA_character_,"model_coefficient_inference")
      if("method"%in%names(model_tab))model_tab$method<-"selected_fit_inverse_curvature_wald"
      if("inference"%in%names(model_tab))model_tab$inference<-availability$inferential_scope[match(model_tab$gene_id,availability$gene_id)]
    }
    result$test_status<-rbind(tab,model_tab[,names(tab),drop=FALSE])
  }
  if(!is.null(result$test_status)) {
    tab<-scrope_reporting_join(result$test_status,native_reporting$tests,c("gene_id","contrast","test"))
    idx<-match(tab$gene_id,availability$gene_id)
    tab$selected_route<-availability$selected_route[idx]
    tab$selected_source<-availability$selected_source[idx]
    unified_status<-ifelse(tab$success%in%TRUE,"available",
      ifelse(grepl("unsupported|not_supported",tab$failure_reason),"unsupported","unavailable"))
    unified_status[tab$success%in%TRUE&tab$selected_source=="recovery"]<-"available_with_limitations"
    tab$availability_status<-unified_status
    tab$availability_stage<-vapply(seq_len(nrow(tab)),function(i)
      scrope_failure_stage(tab$failure_reason[i],isTRUE(tab$success[i])),character(1))
    if(!"status"%in%names(tab))tab$status<-unified_status
    result$test_status<-tab
    for(j in seq_len(nrow(s))) {
      failures<-tab[tab$gene_id==s$gene_id[j]&!tab$success%in%TRUE,,drop=FALSE]
      s$test_status[j]<-if(!nrow(failures))"ok" else paste(unique(failures$test),collapse=",")
      s$test_failure_detail[j]<-if(!nrow(failures))NA_character_ else
        paste(unique(paste(failures$test,failures$failure_reason,sep=":")),collapse=";")
    }
    if(keep_diagnostics&&!is.null(result$diagnostics))result$diagnostics$test_status<-tab
  }
  result$summary<-s
  if(is.null(result$recovery))result$recovery<-list(policy="unified_recovery_20261106",
    mode=mode,primary=NULL,diagnostics=NULL)
  result$recovery$mode<-mode;result$recovery$eligible<-eligible;result$recovery$skip_reason<-skip_reason
  result$recovery$report_schema_version<-3L
  result$provenance$result_schema_version<-6L
  result$provenance$model_covariance_convention<-model_convention
  result$provenance$recovery<-list(mode=mode,eligible=eligible,skip_reason=skip_reason,
    controls=controls,details_retained=keep_diagnostics,
    convergence_semantics="legacy_method_specific_use_availability_for_selected_results")
  if(notify)scrope_recovery_notice(availability)
  result
}

scrope_recovery_notice <- function(availability, fields="availability and coefficient_availability") {
  if(!nrow(availability))return(invisible(NULL))
  recovered<-sum(availability$selected_source=="recovery")
  attention<-availability$status %in% c("partial","point_only","unavailable")
  if(!recovered&&!any(attention))return(invisible(NULL))
  count_field<-function(nm,value=TRUE)if(nm%in%names(availability))sum(availability[[nm]]%in%value) else 0L
  counts<-c(genes=nrow(availability),recovered=recovered,
    all_wald=sum(availability$selected_wald_coefficients==availability$total_coefficients),
    partial_wald=sum(availability$selected_wald_coefficients>0&availability$selected_wald_coefficients<availability$total_coefficients),
    all_model_wald=sum(availability$selected_model_wald_coefficients==availability$total_coefficients),
    partial_model_wald=sum(availability$selected_model_wald_coefficients>0&availability$selected_model_wald_coefficients<availability$total_coefficients),
    point_only=sum(availability$status=="point_only"),unavailable=sum(availability$status=="unavailable"),
    native_incomplete=sum(availability$selected_source=="native" & attention),
    recovery_attempted=count_field("recovery_attempted"),
    recovery_failed=count_field("recovery_outcome","failed"),
    recovery_incomplete=count_field("recovery_outcome","incomplete"),
    native_incomplete_before=count_field("native_qualification",c("partial","point_only","unavailable")))
  message<-sprintf("scRoPE %s recovery: %d/%d genes recovered; robust coefficient Wald available for all coefficients in %d genes, partially in %d; %d point-only and %d unavailable. See %s.",
    availability$requested_model[1],counts["recovered"],counts["genes"],counts["all_wald"],
    counts["partial_wald"],counts["point_only"],counts["unavailable"],fields)
  message<-paste(message,sprintf("Model coefficient Wald available for all coefficients in %d genes, partially in %d.",counts["all_model_wald"],counts["partial_model_wald"]))
  if("recovery_attempted"%in%names(availability))message<-paste(message,sprintf(
    "Native inference incomplete before recovery in %d genes; recovery attempted for %d (%d failed, %d incomplete).",
    counts["native_incomplete_before"],counts["recovery_attempted"],counts["recovery_failed"],counts["recovery_incomplete"]))
  if(counts["native_incomplete"]>0)message<-paste(message,"Some selected native fits have incomplete inference; supported native results are retained.")
  if(recovered)message<-paste(message,"Recovered targets are route-specific; selection adjustment and general finite-sample calibration are not established.")
  warning(structure(list(message=message,call=NULL,code="scrope_recovery_summary",counts=counts),
    class=c("scrope_recovery_warning","warning","condition")))
  invisible(NULL)
}

# Adapt the explicit interface to the same aggregate warning contract.
scrope_working_notice <- function(result, model) {
  r <- result$routes
  status <- ifelse(r$wald_coefficients==r$total_coefficients&r$model_wald_coefficients==r$total_coefficients,"available",
    ifelse(r$wald_coefficients>0|r$model_wald_coefficients>0,"partial",ifelse(r$point_coefficients>0,"point_only","unavailable")))
  status[r$route=="native_partial" & r$wald_coefficients==r$total_coefficients] <- "partial"
  scrope_recovery_notice(data.frame(requested_model=model,
    selected_source=ifelse(r$native_selected,"native",ifelse(r$route=="withheld","unavailable","recovery")),
    status=status,selected_wald_coefficients=r$wald_coefficients,
    selected_model_wald_coefficients=r$model_wald_coefficients,total_coefficients=r$total_coefficients,
    native_qualification=r$native_qualification,recovery_attempted=r$recovery_attempted,recovery_outcome=r$recovery_outcome),
    fields="routes and summary")
}
