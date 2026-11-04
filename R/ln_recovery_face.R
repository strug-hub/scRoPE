# General one-dimensional whole-zero-donor face. No gene identifiers or data-set
# names enter detection. Uses the installed 0.3.0.9003 criterion/derivatives.
ln_recovery_refine <- function(primary,dat,target=2L) {
  fail<-function(reason)list(success=FALSE,reason=reason)
  rec<-primary$joint$candidate_records
  cert<-Filter(function(r)isTRUE(r$certified)&&all(r$checks),rec)
  cores<-Filter(function(r)isTRUE(r$core)&&isTRUE(r$required),cert)
  if(length(unique(lapply(cores,`[[`,'start')))<2L)return(fail('no_corroborated_regular_anchor'))
  anchor<-cert[[which.min(vapply(cert,`[[`,numeric(1),'objective'))]]
  if(any(vapply(cert,function(r)max(abs(r$theta_reduced_log-anchor$theta_reduced_log)/(1+abs(anchor$theta_reduced_log)))>1e-6,logical(1))))
    return(fail('distinct_certified_roots'))
  bad<-Filter(function(r)!isTRUE(r$certified),rec)
  allowed<-function(r)isTRUE(r$finite_endpoint)&&isTRUE(r$interior)&&
    length(names(r$checks)[!r$checks])>0&&all(names(r$checks)[!r$checks]%in%c('root','newton'))
  if(!length(bad)||!all(vapply(bad,allowed,logical(1))))return(fail('failure_outside_local_refinement'))
  fits<-lapply(bad,function(r)ln_research_fixed_fit(dat,list(),list(r$theta_reduced_log,anchor$theta_reduced_log),target))
  if(!all(vapply(fits,function(f)isTRUE(f$success)&&
    max(abs(f$state$full-anchor$theta_reduced_log)/(1+abs(anchor$theta_reduced_log)))<1e-6,logical(1))))
    return(list(success=FALSE,reason='local_refinement_not_certified',attempts=fits))
  f<-fits[[which.min(vapply(fits,function(f)f$state$q,numeric(1)))]]
  natural<-f$state$full;natural[dat$ctx$nb+seq_len(2L)]<-exp(tail(natural,2L))
  objective_error<-ln_recovery_scalar_evidence(natural,dat)$error
  audit<-ln_recovery_competitor_audit(primary,dat,f$state$q,objective_error)
  if(!isTRUE(audit$success))return(list(success=FALSE,
    reason='refinement_retained_competitor_not_resolved',competitor_audit=audit))
  f$competitor_audit<-audit
  f$reason<-'certified_joint_refinement';f$refined_labels<-vapply(bad,`[[`,'','label');f
}

ln_recovery_face_detect <- function(y, X, id) {
  fail<-function(reason,extra=list())c(list(success=FALSE,reason=reason),extra)
  p<-ncol(X);positive<-y>0
  if(p<2L || any(X[,1L]!=1))return(fail('face_requires_first_unit_intercept'))
  if(sum(positive)<p)return(fail('insufficient_positive_rows'))
  pc<-center_m(X);S<-as.matrix(pc$pred)
  sv<-svd(S[positive,,drop=FALSE],nu=0,nv=p)
  rank<-sum(sv$d>1e-10*max(sv$d))
  if(rank!=p-1L)return(fail('no_unique_one_dimensional_face',list(positive_rank=rank)))
  d<-sv$v[,p];eta<-as.numeric(S%*%d)
  tol<-128*.Machine$double.eps*p*pmax(1,rowSums(abs(S)))
  orient<-function(a)all(a<=tol)&&any(a< -100*tol)
  if(!orient(eta)){d<- -d;eta<- -eta}
  if(!orient(eta))return(fail('no_separating_direction'))
  if(any(abs(eta[positive])>tol[positive]))return(fail('positive_rows_not_exact_face'))
  retained<-abs(eta)<=tol;removed<-eta< -100*tol
  if(any(!retained&!removed))return(fail('ambiguous_face_membership'))
  if(any(y[removed]!=0)||!any(removed))return(fail('removed_counts_not_zero'))
  ids<-match(id,unique(id));bydonor<-split(seq_along(ids),ids)
  if(any(vapply(bydonor,function(i)any(removed[i])&&any(retained[i]),logical(1))))
    return(fail('partial_donor_face_outside_scope'))
  # Independent columns in original order, retaining the intercept. The final
  # linear identity is verified at machine-scale tolerance below.
  J<-integer();last_rank<-0L
  for(j in seq_len(p)){
    trial<-c(J,j);r<-qr(S[retained,trial,drop=FALSE],tol=1e-10)$rank
    if(r>last_rank){J<-trial;last_rank<-r}
  }
  if(length(J)!=p-1L||J[1]!=1L)return(fail('retained_design_not_rank_p_minus_one'))
  sr<-center_m(X[retained,J,drop=FALSE])
  Bsmall<-qr.solve(S[retained,J,drop=FALSE],as.matrix(sr$pred),tol=1e-10)
  B<-matrix(0,p,length(J));B[J,]<-Bsmall
  if(max(abs(S[retained,,drop=FALSE]%*%B-as.matrix(sr$pred)))>1e-11)
    return(fail('reduced_predictor_mapping_failed'))
  estimable<-abs(d)<=128*.Machine$double.eps*p
  scale<-as.numeric(pc$sds);scale[scale==0]<-1
  list(success=TRUE,reason='one_dimensional_whole_zero_donor_face',retained=retained,
       removed=removed,direction=d,direction_predictor=eta,tolerance=tol,
       columns=J,beta_map=B,public_map=diag(1/scale),estimable=estimable,
       original_donors=length(bydonor),retained_donors=length(unique(ids[retained])),
       retained_donor_indices=unique(ids[retained]),singular_values=sv$d,
       positive_rank=rank)
}

ln_recovery_face_fit <- function(y,X,id,offset) {
  face<-ln_recovery_face_detect(y,X,id)
  if(!face$success)return(list(success=FALSE,reason=face$reason,face=face))
  J<-face$columns;keep<-face$retained;p<-ncol(X);r<-length(J)
  targets<-which(face$estimable&seq_len(p)>1L)
  if(!length(targets))return(list(success=FALSE,reason='no_estimable_nonintercept_coefficient',face=face))
  target<-targets[1];target_reduced<-match(target,J)
  if(is.na(target_reduced))stop('Estimable coefficient absent from reduced basis.')
  native <- suppressWarnings(scrope(y[keep], id[keep], X[keep,J,drop=FALSE],
    offset=offset[keep], model="LN", ln_recovery="none", cpc=0, mincp=1,
    keep_diagnostics=TRUE, verbose=FALSE))
  primary <- list(success=isTRUE(native$provenance$per_gene$inference_available[1]),
                  joint=native$diagnostics$raw[[1]]$ln_joint)
  ans <- ln_recovery_route(primary,y[keep],X[keep,J,drop=FALSE],id[keep],offset[keep],
                            target_reduced,allow_face=FALSE)
  if(!isTRUE(ans$result$success))return(list(success=FALSE,reason='reduced_fit_unavailable',face=face,reduced=ans))
  route<-ans$route; f<-ans$result; refined<-ans$refinement
  if(route=='ordinary_jln'){
    state<-primary$joint$state
    full<-state$full_log;A<-state$A_full;U<-state$U_subject_full;free<-seq_len(r+2L)
    q<-state$raw$value
  }else{
    full<-f$state$full;A<-f$state$A;U<-f$state$U;q<-f$state$q
    free<-if(is.null(f$spec$free))seq_len(r+2L)else f$spec$free
  }
  IF<-solve(A,U);V<-tcrossprod(IF)
  dat<-ln_research_context(y[keep],X[keep,J,drop=FALSE],id[keep],offset[keep])
  natural<-full;natural[r+seq_len(2L)]<-exp(tail(natural,2L))
  objective_error<-ln_recovery_scalar_evidence(natural,dat)$error
  T<-matrix(0,p+2L,r+2L);T[seq_len(p),seq_len(r)]<-face$beta_map
  T[p+1L,r+1L]<-T[p+2L,r+2L]<-1
  C<-face$public_map%*%T[seq_len(p),free,drop=FALSE]
  beta_rep<-as.numeric(face$public_map%*%face$beta_map%*%full[seq_len(r)])
  Vbeta<-C%*%V%*%t(C);beta<-beta_rep;beta[!face$estimable]<-NA_real_
  Vbeta[!face$estimable,]<-NA_real_;Vbeta[,!face$estimable]<-NA_real_
  Uoriginal<-matrix(0,nrow(U),face$original_donors);Uoriginal[,face$retained_donor_indices]<-U
  list(success=TRUE,reason='certified_identifiable_nuisance_limit',face=face,reduced=ans,refinement=refined,
    reduced_route=route,reduced_full=full,free=free,q=q,objective_error=objective_error,A=A,U=U,U_original=Uoriginal,
    covariance=V,beta=beta,beta_representative=beta_rep,Vbeta=Vbeta,
    derivative_map=T[,free,drop=FALSE],full_beta_map=face$beta_map)
}
