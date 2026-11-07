# Public recovery reporting is separate from the method-specific acceptance
# rules. Compact records contain no evaluators, donor arrays or search ledgers.
working_minimal_record <- function(ans) {
  fields <- c("success", "point_success", "inference_available", "route", "reason",
    "fixed", "free", "estimable", "beta_report", "sampling_target", "inference_reason",
    "native_regular", "elapsed_seconds", "input_identity", "identity", "dimensions")
  out <- ans[intersect(fields, names(ans))]
  if (!is.null(ans$V_report)) out$variance_report <- diag(ans$V_report)
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
    native_regular=ans$native_regular, native_preservation=ans$native_preservation,
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
  primary <- result
  labels <- sub("^logFC_", "", grep("^logFC_",names(result$summary),value=TRUE))
  for (j in seq_len(nrow(result$summary))) {
    id <- result$summary$gene_id[j]
    a <- records[[as.character(id)]]
    if (isTRUE(a$native_regular)) next
    reason <- if (isTRUE(a$point_success)) "not_supported_for_recovery" else a$reason
    if (is.null(reason) || !length(reason)) reason <- "recovery_unavailable"
    result <- scrope_clear_selected_tests(result,j,reason)
    rows <- selected$summary[selected$summary$gene_id==id,,drop=FALSE]
    rows <- rows[match(labels,rows$coefficient),,drop=FALSE]
    for (k in seq_along(labels)) {
      nm <- labels[k];r <- rows[k,,drop=FALSE]
      result$summary[j,paste0("logFC_",nm)] <- r$beta
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
      reported_covariance <- if(any(rows$available))a$V_report else NULL
      if(!is.null(reported_covariance)) {
        reported_covariance[!rows$covariance_available,] <- NA_real_
        reported_covariance[,!rows$covariance_available] <- NA_real_
      }
      result$diagnostics$godambe[[j]] <- list(route=a$route,
        inference_available=any(rows$available), covariance_coordinates="reported_coefficients",
        inferential_scope="coefficient_wald_only", robust_cov=reported_covariance)
    }
  }
  result$recovery <- list(policy="hl_pgmm_working_20261105",mode="auto",
    routes=selected$routes,coefficients=selected$summary,primary=primary,
    audit=selected$audit,diagnostics=selected$diagnostics,
    limitations=selected$provenance$limitations)
  if(model=="HL"&&!is.null(result$provenance$per_gene)) {
    fields<-intersect(names(result$provenance$per_gene),names(result$summary))
    result$provenance$per_gene[fields]<-result$summary[fields]
  }
  result$provenance$selected_routes<-selected$routes
  result$provenance$native_estimator<-result$provenance$estimator
  if(any(!selected$routes$native_regular)) {
    result$provenance$model_covariance<-"native_rows_only"
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
  if(!is.null(result$recovery$primary))result$recovery$primary<-scrope_finish_recovery(
    result$recovery$primary,model,"none",eligible,skip_reason,keep_diagnostics,controls,notify=FALSE)
  s<-result$summary;labels<-sub("^logFC_","",grep("^logFC_",names(s),value=TRUE))
  if(!"gene_id"%in%names(s))stop("Recovery reporting requires original gene identifiers.")
  existing<-result$recovery
  selected<-if(!is.null(existing$coefficients)&&"point_available"%in%names(existing$coefficients)) existing$coefficients else NULL
  rows<-list();genes<-list()
  for(j in seq_len(nrow(s))) {
    id<-s$gene_id[j]
    route<-if(!is.null(existing$routes))as.character(existing$routes$route[match(id,existing$routes$gene_id)]) else "native"
    native<-route %in% c("native","ordinary_jln")
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
      status<-if(ok)if(native)"available" else "available_with_limitations" else if(point)"point_only" else "unavailable"
      reason<-s[[paste0("wald_failure_",nm)]][j]
      gene_rows[[k]]<-data.frame(gene_id=id,gene=s$gene[j],coefficient=nm,model=model,
        route=route,source=source,status=status,point_available=point,covariance_available=cov,
        statistic_available=ok&&is.finite(stat),reference_available=ok,
        interval_available=point&&cov&&is.finite(se)&&se>0,available=ok,estimable=estimable,
        beta=if(point)beta else NA_real_,robust_se=if(cov)se else NA_real_,
        statistic=if(ok)stat else NA_real_,p_value=if(ok)pv else NA_real_,
        scope=scope,sampling_target=sampling,selection_adjusted=FALSE,
        failure_reason=if(ok)NA_character_ else if(length(reason))reason else "unavailable",
        stringsAsFactors=FALSE)
    }
    cr<-do.call(rbind,gene_rows);rows[[j]]<-cr
    np<-sum(cr$point_available);nc<-sum(cr$covariance_available);nw<-sum(cr$available);p<-nrow(cr)
    status<-if(nw==p)if(native)"available" else "available_with_limitations" else
      if(nw>0)"partial" else if(np>0)"point_only" else "unavailable"
    genes[[j]]<-data.frame(gene_id=id,gene=s$gene[j],requested_model=model,effective_model=model,
      selected_route=route,selected_source=source,status=status,recovery_mode=mode,
      recovery_eligible=eligible,recovery_skip_reason=skip_reason,
      native_convergence=if(!is.null(existing$primary))existing$primary$convergence[j] else result$convergence[j],
      selected_point_coefficients=np,selected_covariance_coefficients=nc,selected_wald_coefficients=nw,
      total_coefficients=p,inferential_scope=scope,sampling_target=sampling,
      selection_adjusted=FALSE,failure_reason=if(nw==p)NA_character_ else
        paste(unique(na.omit(cr$failure_reason)),collapse=";"),stringsAsFactors=FALSE)
  }
  availability<-do.call(rbind,genes);coefficients<-do.call(rbind,rows)
  # Preserve existing LN availability columns/enumerations; common fields are
  # additive and have identical names and meanings across the three models.
  if(is.data.frame(result$availability)) {
    old<-result$availability[match(availability$gene_id,result$availability$gene_id),,drop=FALSE]
    for(nm in names(availability))old[[nm]]<-availability[[nm]]
    availability<-old
  }
  result$availability<-availability
  result$coefficient_availability<-coefficients
  if(model!="LN")result$test_status<-scrope_collect_test_status(s,model)
  if(!is.null(result$test_status)) {
    tab<-result$test_status;idx<-match(tab$gene_id,availability$gene_id)
    tab$selected_route<-availability$selected_route[idx]
    tab$selected_source<-availability$selected_source[idx]
    unified_status<-ifelse(tab$success%in%TRUE,"available",
      ifelse(grepl("unsupported|not_supported",tab$failure_reason),"unsupported","unavailable"))
    unified_status[tab$success%in%TRUE&tab$selected_source=="recovery"]<-"available_with_limitations"
    tab$availability_status<-unified_status
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
  result$recovery$report_schema_version<-1L
  result$provenance$result_schema_version<-4L
  result$provenance$recovery<-list(mode=mode,eligible=eligible,skip_reason=skip_reason,
    controls=controls,details_retained=keep_diagnostics,
    convergence_semantics="legacy_method_specific_use_availability_for_selected_results")
  if(notify)scrope_recovery_notice(availability)
  result
}

scrope_recovery_notice <- function(availability) {
  if(!nrow(availability))return(invisible(NULL))
  recovered<-sum(availability$selected_source=="recovery")
  attention<-availability$status %in% c("partial","point_only","unavailable")
  if(!recovered&&!any(attention))return(invisible(NULL))
  counts<-c(genes=nrow(availability),recovered=recovered,
    all_wald=sum(availability$selected_wald_coefficients==availability$total_coefficients),
    partial_wald=sum(availability$selected_wald_coefficients>0&availability$selected_wald_coefficients<availability$total_coefficients),
    point_only=sum(availability$status=="point_only"),unavailable=sum(availability$status=="unavailable"))
  message<-sprintf("scRoPE %s recovery: %d/%d genes recovered; coefficient Wald available for all coefficients in %d genes, partially in %d; %d point-only and %d unavailable. See availability and coefficient_availability.",
    availability$requested_model[1],counts["recovered"],counts["genes"],counts["all_wald"],
    counts["partial_wald"],counts["point_only"],counts["unavailable"])
  if(recovered)message<-paste(message,"Recovered targets are route-specific; selection adjustment and general finite-sample calibration are not established.")
  warning(structure(list(message=message,call=NULL,code="scrope_recovery_summary",counts=counts),
    class=c("scrope_recovery_warning","warning","condition")))
  invisible(NULL)
}
