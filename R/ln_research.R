# Isolated positive-dispersion JLN research helpers.
# Ported without estimator changes from the frozen 2026-11-03 prototype.
# Research-only fixed positive dispersion JLN. No package namespace changes.
ln_research_context<-function(y,X,id,offset=rep(1,length(y))) {
  stopifnot(is.numeric(y),all(y>=0),all(is.finite(y)),length(y)==nrow(X),length(id)==length(y),
    all(is.finite(X)),all(is.finite(offset)),all(offset>0))
  id<-match(id,unique(id));stopifnot(!is.unsorted(id))
  count<-Matrix(matrix(y,nrow=1),sparse=TRUE)
  fid<-as.integer(c(which(c(1,diff(id))==1),length(id)+1))
  cy<-call_cumsumy(count,fid-1,max(id),1L)
  pc<-center_m(as.matrix(X));stopifnot(sum(pc$sds==0)==1L,all(pc$sds>=0))
  sc<-as.numeric(pc$sds);sc[sc==0]<-1
  ctx<-list(pred=as.matrix(pc$pred),offset=cv_offset(as.double(offset),1,length(id))$offset,
    id=id,fid=fid,cumsumy=cy,posind=list(which(cy[1,]>0)),nb=ncol(X),nind=length(y),k=max(id),
    min=c(1e-4,1e-4),max=c(10,1000))
  list(ctx=ctx,map=ln_joint_map(ctx$nb),posv=call_posindy(t(count),0L,length(y)),
    pub=diag(1/sc),y=y,X=X,offset=offset,id=id,count=count)
}
ln_research_fd<-function(fun,x,mult=1) {
  h<-mult*.Machine$double.eps^(1/3)*(1+abs(x));base<-as.numeric(fun(x))
  vapply(seq_along(x),function(j){plus<-minus<-x;plus[j]<-plus[j]+h[j];minus[j]<-minus[j]-h[j];
    (as.numeric(fun(plus))-as.numeric(fun(minus)))/(2*h[j])},numeric(length(base)))
}
ln_research_spec<-function(dat,fixed) {
  stopifnot(all(names(fixed)%in%c('v','phi')),!anyDuplicated(names(fixed)),
    all(is.finite(unlist(fixed))),all(unlist(fixed)>0))
  p<-dat$ctx$nb;values<-rep(NA_real_,p+2L)
  if(!is.null(fixed$v))values[p+1L]<-log(fixed$v)
  if(!is.null(fixed$phi))values[p+2L]<- -log(fixed$phi)
  free<-which(is.na(values));list(values=values,free=free,fixed=fixed,
    lower=c(rep(-100,p),log(dat$ctx$min))[free],upper=c(rep(100,p),log(dat$ctx$max))[free])
}
ln_research_state<-function(a,dat,spec) {
  full<-spec$values;full[spec$free]<-a
  st<-ln_joint_derivative_state(full,1L,dat$posv,dat$ctx,dat$map)
  idx<-spec$free
  list(full=full,q=st$raw$value,g=st$q_gradient_full[idx],A=st$A_full[idx,idx,drop=FALSE],
    U=st$U_subject_full[idx,,drop=FALSE],full_state=st)
}


ln_research_fixed_fit<-function(dat,fixed,starts,target=2L,derivatives=TRUE) {
  began<-proc.time()[['elapsed']];spec<-ln_research_spec(dat,fixed);p<-dat$ctx$nb
  counts<-c(derivative_state=0L,fd_scalar=0L);cache_x<-NULL;cache_st<-NULL;stage<-"initial";trace<-list()
  ev<-function(a){if(!is.null(cache_x)&&identical(as.numeric(a),cache_x))return(cache_st)
    counts['derivative_state']<<-counts['derivative_state']+1L
    st<-ln_research_state(a,dat,spec)
    flags<-c(q=is.finite(st$q),g=all(is.finite(st$g)),A=all(is.finite(st$A)),U=all(is.finite(st$U)))
    if(!all(flags))trace[[length(trace)+1L]]<<-list(stage=stage,coordinates=a,full=st$full,flags=flags,
      q=st$q,gradient=st$g,A=st$A,eta_range=range(dat$ctx$offset+as.numeric(dat$ctx$pred%*%st$full[seq_len(p)])))
    stopifnot(is.finite(st$q),all(is.finite(st$g)))
    cache_x<<-as.numeric(a);cache_st<<-st;st}
  fail<-function(reason,extra=list())c(list(success=FALSE,reason=reason,fixed=fixed,
    trace=trace,elapsed=proc.time()[['elapsed']]-began,counts=counts),extra)
  if(length(starts)!=2L)return(fail('exactly_two_starts_required'))
  initials<-lapply(starts,function(a)as.numeric(a)[spec$free])
  if(any(vapply(initials,function(a)length(a)!=length(spec$free)||any(!is.finite(a)),logical(1))))return(fail('invalid_start'))
  initials<-lapply(initials,function(a)pmin(pmax(a,spec$lower),spec$upper))
  if(max(abs(initials[[1]]-initials[[2]]))<1e-4)initials[[2]][1]<-initials[[2]][1]+.1
  one_attempt<-function(a,scale) {
    stage<<-paste0("L-BFGS-B scale=",scale);opt<-optim(a,function(x)ev(x)$q/scale,function(x)ev(x)$g/scale,method='L-BFGS-B',lower=spec$lower,upper=spec$upper,
      control=list(maxit=180,factr=1e4,pgtol=1e-10))
    stage<<-"endpoint";a<-opt$par;endpoint<-ev(a)
    stopifnot(all(is.finite(endpoint$A)),all(is.finite(endpoint$U)))
    stage<<-"Newton";steps<-0L
    for(it in seq_len(20)) {
      st<-ev(a);curv<-ln_joint_curvature_diagnostics(st$A)
      if(!isTRUE(curv$positive_definite))break
      d<-solve(st$A,st$g);if(max(abs(d))<1e-11)break
      accepted<-FALSE
      for(rate in 2^-(0:15)) {
        proposal<-a-rate*d
        if(any(proposal<spec$lower|proposal>spec$upper))next
        sn<-tryCatch(ev(proposal),error=function(e)NULL)
        if(is.null(sn)||any(!is.finite(sn$A))||any(!is.finite(sn$U)))next
        if(sn$q<=st$q+1e-12*(1+abs(st$q))){a<-proposal;steps<-steps+1L;accepted<-TRUE;break}
      }
      if(!accepted)break
    }
    list(a=a,state=ev(a),optimization_scale=scale,optimizer=opt[c('convergence','message','counts')],polish_steps=steps)
  }
  one<-function(a) {
    first<-tryCatch(one_attempt(a,1),error=function(e)list(error=conditionMessage(e)))
    if(is.null(first$error))return(first)
    cache_x<<-NULL;cache_st<<-NULL
    second<-one_attempt(a,dat$ctx$k)
    second$initial_attempt_error<-first$error
    second
  }
  trials<-lapply(initials,function(a)tryCatch(one(a),error=function(e)list(error=conditionMessage(e))))
  if(any(vapply(trials,function(x)!is.null(x$error),logical(1))))return(fail('evaluation_or_optimizer_error',list(trials=trials)))
  winner<-which.min(vapply(trials,function(x)x$state$q,numeric(1)));st<-trials[[winner]]$state;a<-trials[[winner]]$a
  stage<-"certification"
  if(any(!is.finite(st$A))||any(!is.finite(st$U)))return(fail("nonfinite_endpoint",list(trials=trials)))
  curv<-ln_joint_curvature_diagnostics(st$A)
  metrics<-c(root=max(abs(st$g)/(1+sqrt(rowSums(st$U^2)))),
    additivity=max(abs(rowSums(st$U)+st$g)/(1+abs(st$g))),
    min_eigen=curv$scaled_min_eigenvalue,condition=curv$scaled_condition_number,
    start_agreement=max(abs(trials[[1]]$a-trials[[2]]$a)/(1+abs(a))))
  margin<-sqrt(.Machine$double.eps)*(1+pmax(abs(spec$lower),abs(spec$upper)))
  free_interior<-all(a>spec$lower+margin & a<spec$upper-margin)
  preliminary<-c(root=metrics['root']<=1e-6,additivity=metrics['additivity']<=1e-10,
    curvature=isTRUE(curv$positive_definite)&&metrics['min_eigen']>1e-8,
    condition=metrics['condition']<1e10,interior=free_interior,starts=metrics['start_agreement']<1e-6)
  if(!all(preliminary))return(fail(paste(names(preliminary)[!preliminary],collapse='|'),
    list(state=st,spec=spec,metrics=metrics,trials=trials)))
  IF<-solve(st$A,st$U);V<-tcrossprod(IF);newton<-solve(st$A,st$g)
  metrics<-c(metrics,newton=max(abs(newton)/(1+abs(a))),newton_in_se=max(abs(newton)/sqrt(diag(V))))
  if(any(!is.finite(V))||any(diag(V)<=0)||metrics['newton']>1e-6||metrics['newton_in_se']>1e-4)
    return(fail('covariance_or_newton',list(state=st,spec=spec,metrics=metrics,trials=trials)))
  if(derivatives) {stage<-"derivative_audit"
    qdon<-function(x){counts['fd_scalar']<<-counts['fd_scalar']+1L;theta<-spec$values;theta[spec$free]<-x
      theta[(p+1):(p+2)]<-exp(theta[(p+1):(p+2)])
      rowSums(ln_joint_objective_components(theta,1L,dat$posv,dat$ctx,dat$map))}
    Ug<- -t(ln_research_fd(qdon,a));Ug2<- -t(ln_research_fd(qdon,a,.5))
    Hfd<-ln_research_fd(function(x)ev(x)$g,a);Hfd2<-ln_research_fd(function(x)ev(x)$g,a,.5)
    scaled<-function(x,y)max(abs(x-y)/(1+abs(y)))
    metrics<-c(metrics,donor_gradient_fd=scaled(Ug,st$U),donor_gradient_half=scaled(Ug,Ug2),
      hessian_fd=scaled(Hfd,st$A),hessian_half=scaled(Hfd,Hfd2))
    if(metrics['donor_gradient_fd']>1e-5||metrics['donor_gradient_half']>1e-4||
       metrics['hessian_fd']>1e-4||metrics['hessian_half']>1e-4)
      return(fail('derivative_audit',list(state=st,spec=spec,metrics=metrics,trials=trials)))
  }
  beta<-as.numeric(dat$pub%*%st$full[seq_len(p)]);Vbeta<-dat$pub%*%V[seq_len(p),seq_len(p),drop=FALSE]%*%t(dat$pub)
  target_if<-as.numeric(dat$pub[target,]%*%IF[seq_len(p),,drop=FALSE]);se<-sqrt(Vbeta[target,target])
  list(success=TRUE,reason='certified_fixed_research_root',fixed=fixed,spec=spec,state=st,trials=trials,
    metrics=metrics,beta=beta,covariance=V,Vbeta=Vbeta,influence=IF,target_if=target_if,
    coefficient=beta[target],se=se,p=2*pnorm(-abs(beta[target]/se)),
    fitted_log_mean=dat$ctx$offset+as.numeric(dat$ctx$pred%*%st$full[seq_len(p)]),
    max_donor_variance_share=max(target_if^2)/sum(target_if^2),
    trace=trace,elapsed=proc.time()[['elapsed']]-began,counts=counts,finite_cluster_multiplier=1)
}

# Explicit research warning policy. No blanket waiver of a failed start.
ln_research_warning_fit<-function(primary,dat,target=2L) {
 fail<-function(reason)list(success=FALSE,reason=reason)
 rec<-primary$joint$candidate_records
 cert<-which(vapply(rec,function(r)isTRUE(r$certified)&&all(r$checks),logical(1)))
 if(!length(cert))return(fail('no_regular_anchor'))
 best<-cert[which.min(vapply(rec[cert],`[[`,numeric(1),'objective'))];anchor<-rec[[best]]
 cores<-Filter(function(r)isTRUE(r$core)&&isTRUE(r$required),rec[cert])
 if(length(unique(lapply(cores,`[[`,'start')))<2L||any(vapply(cores,function(r)isTRUE(r$recovery$adopted),logical(1))))return(fail('insufficient_distinct_original_core_starts'))
 for(r in rec[cert])if(max(abs(r$theta_reduced_log-anchor$theta_reduced_log)/(1+abs(anchor$theta_reduced_log)))>1e-6)
   return(fail('distinct_certified_roots'))
 vetoes<-Filter(function(r)isTRUE(r$resolution$veto),rec)
 allowed<-function(r){bad<-names(r$checks)[!r$checks]
  isTRUE(r$required)&&isTRUE(r$finite_endpoint)&&!isTRUE(r$certified)&&!isTRUE(r$interior)&&
    'interior'%in%bad&&all(bad%in%c('root','newton','curvature','condition','interior'))&&
    identical(r$resolution$reason,'required_candidate_not_certified')&&
    identical(r$objective_instance,anchor$objective_instance)&&identical(r$objective_provenance,anchor$objective_provenance)&&
    (is.null(r$recovery)||!grepl('error|nonfinite|unusable',r$recovery$status))}
 if(!length(vetoes)||!all(vapply(vetoes,allowed,logical(1))))return(fail('veto_outside_warning_scope'))
 known<-ln_recovery_known_points(primary,dat)
 if(!isTRUE(known$success))return(fail('invalid_retained_observations'))
 a<-ln_recovery_scalar_evidence(anchor$theta_reduced,dat)
 a$error<-a$error+abs(a$q-anchor$objective)
 audit<-known$observations
 audit$gap<-audit$q-a$q
 audit$margin<-a$error+audit$error+1e-6*a$scale
 audit$passes<-ifelse(audit$boundary,audit$gap>audit$margin,audit$gap>= -audit$margin)
 if(!all(audit$passes))return(fail('better_or_uncertain_competitor'))
 p<-dat$ctx$nb;spec<-list(values=rep(NA_real_,p+2),free=seq_len(p+2))
 st<-ln_research_state(anchor$theta_reduced_log,dat,spec);curv<-ln_joint_curvature_diagnostics(st$A)
 if(!all(is.finite(c(st$q,st$g,st$A,st$U)))||!isTRUE(curv$positive_definite)||
    curv$scaled_min_eigenvalue<=1e-8||curv$scaled_condition_number>=1e10)return(fail('anchor_curvature'))
 IF<-solve(st$A,st$U);V<-tcrossprod(IF);step<-solve(st$A,st$g)
 if(any(!is.finite(V))||any(diag(V)<=0))return(fail('anchor_covariance'))
 if(max(abs(st$g)/(1+sqrt(rowSums(st$U^2))))>1e-6||max(abs(rowSums(st$U)+st$g)/(1+abs(st$g)))>1e-10||
    max(abs(step)/(1+abs(st$full)))>1e-6||max(abs(step)/sqrt(diag(V)))>1e-4)return(fail('anchor_root_or_newton'))
 # Recheck derivatives at the exact reported full joint root.
 qdon<-function(t){t[(p+1):(p+2)]<-exp(t[(p+1):(p+2)]);rowSums(ln_joint_objective_components(t,1L,dat$posv,dat$ctx,dat$map))}
 Ug<--t(ln_research_fd(qdon,st$full));Ug2<--t(ln_research_fd(qdon,st$full,.5))
 H<-ln_research_fd(function(t)ln_research_state(t,dat,spec)$g,st$full);H2<-ln_research_fd(function(t)ln_research_state(t,dat,spec)$g,st$full,.5)
 if(max(abs(Ug-st$U)/(1+abs(st$U)))>1e-5||max(abs(H-st$A)/(1+abs(st$A)))>1e-4||
    max(abs(Ug-Ug2)/(1+abs(Ug2)))>1e-4||max(abs(H-H2)/(1+abs(H2)))>1e-4)return(fail('anchor_derivatives'))
 beta<-as.numeric(dat$pub%*%st$full[seq_len(p)]);Vb<-dat$pub%*%V[seq_len(p),seq_len(p),drop=FALSE]%*%t(dat$pub)
 se<-sqrt(Vb[target,target]);if(!is.finite(se)||se<=0)return(fail('target_variance'))
 list(success=TRUE,reason='corroborated_interior_with_dominated_unresolved_boundary',coefficient=beta[target],se=se,
  p=2*pnorm(-abs(beta[target]/se)),beta=beta,Vbeta=Vb,state=st,covariance=V,observations=audit,core_starts=length(cores),
  warning='An alternative boundary-start trajectory remains unresolved. Observed endpoints are worse; unvisited competing branches have not been excluded.')
}


ln_research_boundary_spec<-function(primary,alternative=FALSE) {
  if(isTRUE(primary$success))return(list(route='primary',reason='native_available'))
  rec<-primary$joint$candidate_records
  if(!length(rec))return(list(route='withheld',reason='no_candidates'))
  ok<-vapply(rec,function(r)isTRUE(r$finite_endpoint)&&is.finite(r$objective),logical(1))
  if(!any(ok))return(list(route='withheld',reason='no_finite_candidates'))
  idx<-which(ok)[which.min(vapply(rec[ok],`[[`,numeric(1),'objective'))];r<-rec[[idx]]
  if(is.null(r$checks)||anyNA(r$checks)||!all(names(r$checks)[!r$checks]%in%c('root','newton','interior')))
    return(list(route='withheld',reason='failure_outside_scope'))
  d<-utils::tail(r$theta_reduced,2);lo<-r$bounds$lower;hi<-r$bounds$upper
  tol<-sqrt(.Machine$double.eps)*(1+pmax(abs(lo),abs(hi)))
  al<-d<=lo+tol;au<-d>=hi-tol
  if(au[1]||al[2]||!any(c(al[1],au[2])))return(list(route='withheld',reason='no_eligible_active_boundary'))
  cert<-Filter(function(x)isTRUE(x$certified),rec)
  if(length(cert)){
    qc<-min(vapply(cert,`[[`,numeric(1),'objective'))
    if(qc-r$objective<=1e-6*max(1,abs(qc),abs(r$objective)))return(list(route='withheld',reason='no_clear_boundary_advantage'))
  }
  fixed<-list();if(al[1])fixed$v<-if(alternative)1e-6 else 1e-4
  if(au[2])fixed$phi<-if(alternative)1e-5 else 1e-3
  list(route=paste0('fixed_',paste(names(fixed),collapse='_')),fixed=fixed,candidate=idx,reason='eligible_boundary')
}


ln_research_starts <- function(primary, dat) {
  rec <- primary$joint$candidate_records
  finite <- which(vapply(rec,function(r)isTRUE(r$finite_endpoint)&&is.finite(r$objective),logical(1)))
  if(length(finite)) {
    best <- finite[which.min(vapply(rec[finite],`[[`,numeric(1),'objective'))]
    selected <- primary$joint$selected_candidate_index
    if(length(selected)!=1L || is.na(selected) || !selected%in%finite) selected<-best
    return(list(rec[[best]]$theta_reduced_log,rec[[selected]]$theta_reduced_log))
  }
  st <- c(log(mean(dat$y/dat$offset)+.01),rep(0,dat$ctx$nb-1L),log(.1),0)
  list(st,st+c(.1,rep(0,length(st)-1L)))
}
ln_research_boundary_check <- function(fit, dat) {
  if(!isTRUE(fit$success))return(list(success=FALSE,reason='restricted_certificate_failed'))
  p<-dat$ctx$nb;st<-fit$state;full<-st$full
  idx<-which(!is.na(fit$spec$values));g<-st$full_state$q_gradient_full
  U<-st$full_state$U_subject_full
  # Lower log-v: nonnegative minimizing gradient. Upper log-gamma: nonpositive.
  direction<-ifelse(idx==p+1L,1,-1)
  outward<-direction*g[idx]/(1+sqrt(rowSums(U[idx,,drop=FALSE]^2)))
  delta<-vapply(seq_along(idx),function(k){t<-full;t[idx[k]]<-t[idx[k]]+direction[k]*log(2)
    natural<-t;natural[(p+1):(p+2)]<-exp(t[(p+1):(p+2)])
    sum(ln_joint_objective_components(natural,1L,dat$posv,dat$ctx,dat$map))-st$q},numeric(1))
  ok<-all(is.finite(outward))&&all(is.finite(delta))&&all(outward>= -1e-6)&&all(delta>= -1e-10*(1+abs(st$q)))
  list(success=ok,reason=if(ok)'local_outward_boundary_supported'else'boundary_direction_not_supported',
    indices=idx,scaled_outward_gradient=outward,inward_objective_change=delta)
}
