# Reporting only: summarize existing native flags before selection overwrites
# public columns. No estimator, derivative, qualification or recovery is run.
scrope_status_reason <- function(x) {
  x<-unique(as.character(x));x<-x[!is.na(x)&nzchar(x)]
  if(length(x))paste(x,collapse=";") else NA_character_
}

scrope_failure_stage <- function(reason, available=FALSE) {
  if(isTRUE(available)||!length(reason)||is.na(reason)||!nzchar(reason))return(NA_character_)
  if(grepl("unsupported|not_supported",reason))return("scope")
  if(grepl("nonestimable",reason))return("estimability")
  if(grepl("nonfinite_.*statistic|undefined_stat",reason))return("statistic")
  if(grepl("covariance|variance|singular_cov|godambe",reason))return("covariance")
  if(identical(reason,"inner_armijo_failure")||
     grepl("converg|optimizer|root|candidate|curvature|boundary|geometry|point|coefficient_estimate",reason))return("point_qualification")
  "inference"
}

scrope_native_snapshot <- function(result, model, records=NULL) {
  s<-result$summary;labels<-sub("^logFC_","",grep("^logFC_",names(s),value=TRUE))
  genes<-coefficients<-vector("list",nrow(s))
  for(j in seq_len(nrow(s))) {
    val<-function(nm,default=NA)if(nm%in%names(s))s[[nm]][j] else default
    pg<-result$provenance$per_gene
    ix<-if(is.data.frame(pg)&&"gene_id"%in%names(pg))match(s$gene_id[j],pg$gene_id) else j
    pgval<-function(nm)if(is.data.frame(pg)&&nm%in%names(pg))pg[[nm]][ix] else NA_character_
    point_ok<-if(model=="HL")isTRUE(val("fit_success")) else if(model=="PGMM")
      isTRUE(val("optimizer_success")) else isTRUE(pgval("inference_available"))
    rows<-lapply(labels,function(nm) {
      point<-point_ok&&is.finite(val(paste0("logFC_",nm),NA_real_))
      se<-val(paste0("se_robust_",nm),NA_real_);ms<-val(paste0("se_",nm),NA_real_)
      robust<-isTRUE(val(paste0("wald_success_",nm)))&&is.finite(val(paste0("p_robust_",nm),NA_real_))
      supported<-paste0("model_wald_success_",nm)%in%names(s)
      mod<-supported&&isTRUE(val(paste0("model_wald_success_",nm)))&&is.finite(val(paste0("p_naive_",nm),NA_real_))
      data.frame(gene_id=s$gene_id[j],coefficient=nm,native_point_available=point,
        native_covariance_available=point&&is.finite(se)&&se>0,
        native_model_covariance_available=point&&is.finite(ms)&&ms>0,
        native_available=robust,native_model_available=mod,native_model_supported=supported,
        native_failure_reason=if(robust)NA_character_ else as.character(val(paste0("wald_failure_",nm),NA_character_)),
        native_model_failure_reason=if(mod||!supported)NA_character_ else as.character(val(paste0("model_wald_failure_",nm),NA_character_)),
        stringsAsFactors=FALSE)
    })
    cc<-do.call(rbind,rows);coefficients[[j]]<-cc
    q<-records[[as.character(s$gene_id[j])]]$native_qualification
    point<-if(!is.null(q))isTRUE(q$point) else any(cc$native_point_available)
    robust<-if(!is.null(q))isTRUE(q$robust_complete) else all(cc$native_available)
    mod<-if(!is.null(q))isTRUE(q$model_complete) else all(!cc$native_model_supported|cc$native_model_available)
    qualification<-if(point&&robust&&mod)"complete" else if(point&&any(cc$native_available|cc$native_model_available))"partial" else
      if(point)"point_only" else "unavailable"
    reason<-scrope_status_reason(c(q$reason,val("fit_failure",NA_character_),pgval("enabled_failure_reason"),
      pgval("failure_reason"),cc$native_failure_reason,cc$native_model_failure_reason))
    if(qualification=="complete")reason<-NA_character_ else if(is.na(reason))reason<-if(!point)"native_point_unqualified" else "native_inference_incomplete"
    genes[[j]]<-data.frame(gene_id=s$gene_id[j],native_qualification=qualification,
      native_point_qualified=point,native_robust_complete=robust,native_model_complete=mod,
      native_point_coefficients=sum(cc$native_point_available),
      native_covariance_coefficients=sum(cc$native_covariance_available),
      native_wald_coefficients=sum(cc$native_available),
      native_model_covariance_coefficients=sum(cc$native_model_covariance_available),
      native_model_wald_coefficients=sum(cc$native_model_available),
      native_failure_reason=reason,stringsAsFactors=FALSE)
  }
  tests<-scrope_collect_test_status(s,model)
  tests<-tests[,c("gene_id","contrast","test","success","failure_reason"),drop=FALSE]
  names(tests)[4:5]<-c("native_success","native_failure_reason")
  list(genes=do.call(rbind,genes),coefficients=do.call(rbind,coefficients),tests=tests)
}

scrope_reporting_join <- function(tab, native, keys) {
  if(is.null(native)||!nrow(native))return(tab)
  key<-function(x)do.call(paste,c(x[keys],sep="\r"))
  idx<-match(key(tab),key(native))
  for(nm in setdiff(names(native),keys))tab[[nm]]<-native[[nm]][idx]
  tab
}

scrope_recovery_outcome <- function(attempted, selected_native, point, complete_native=FALSE) {
  if(isTRUE(complete_native))"not_needed" else if(!isTRUE(attempted))"not_attempted" else
    if(!isTRUE(point))"failed" else if(isTRUE(selected_native))"incomplete" else "selected"
}
