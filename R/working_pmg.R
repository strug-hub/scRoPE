# Certified free-coordinate PGMM working problems. Fixed v is positive, not
# a zero-variance model; all covariance uses the actual free equations.
working_scaled_error <- function(x, y) {
  if (!identical(dim(x), dim(y)) || any(!is.finite(c(x, y)))) return(Inf)
  max(abs(x-y)/(1+abs(y)))
}

working_pmg_problem <- function(y, X, id, offset, fixed = numeric(0)) {
  prep <- hl_prepare_data(y, id, X, offset, cpc=0, mincp=1, verbose=FALSE)
  ctx <- prep$ctx; p <- ctx$nb
  posv <- call_posindy(prep$count, 0L, ctx$nind)
  free <- if (length(fixed)) seq_len(p) else seq_len(p+1L)
  compose <- function(x) {
    z <- if (length(fixed)) c(x, log(fixed[[1L]])) else x
    c(z[seq_len(p)], exp(z[p+1L]))
  }
  evaluate <- function(x) {
    par <- compose(x)
    e <- pmg_ll_der_hes(par, ctx$pred, ctx$offset, posv$Y, ctx$fid,
      ctx$cumsumy[1,], ctx$posind[[1]], posv$posindy, p, ctx$nind, ctx$k)
    T <- diag(c(rep(1,p), par[p+1L]))
    A <- T %*% e$observed_info %*% T
    # Recovery is rooted in the donor estimating equations. Sum within donors
    # first, avoiding the native aggregate's cancellation across all cells.
    natural_score<-rowSums(e$per_subject_gradients)
    U<-(T %*% e$per_subject_gradients)[free,,drop=FALSE]
    A[p+1L,p+1L] <- A[p+1L,p+1L] + par[p+1L]*natural_score[p+1L]
    terms<-scalar_terms(x)
    # A scale-aware consistency tolerance, not a proved libm error bound.
    scalar_tolerance<-128*.Machine$double.eps*(1+terms$absolute_terms)
    list(q=sum(terms$qi), native_q=e$objective, scalar_tolerance=scalar_tolerance,
         g=rowSums(U), native_g=as.numeric(T %*% e$gradient)[free],
         A=A[free,free,drop=FALSE], U=U,
         x=x, par=par, evaluation=e)
  }
  # Independent scalar formula, not differentiation of the score implementation.
  totals <- ctx$cumsumy[1,]
  yy <- as.numeric(prep$count)
  scalar_terms <- function(x) {
    par <- compose(x); v <- par[p+1L]; alpha <- 1/expm1(v)
    log_rate <- log(alpha)-v/2
    eta <- as.numeric(ctx$offset + ctx$pred %*% par[seq_len(p)])
    mu <- as.numeric(rowsum(exp(eta),ctx$id,reorder=FALSE))
    linear <- as.numeric(rowsum(yy*eta,ctx$id,reorder=FALSE))
    increment <- ln_joint_rising_normalized_cpp(alpha, totals)
    qi<- -increment - totals*v/2 + (totals+alpha)*log1p(mu*exp(-log_rate)) - linear
    # Scale by the absolute terms of both equivalent expressions before their
    # cancellations. This diagnostic never changes the objective or root gates.
    absolute_terms<-sum(abs(lgamma(totals+alpha)))+ctx$k*abs(lgamma(alpha))+
      sum(abs(yy*eta))+ctx$k*abs(alpha*log_rate)+
      sum(abs((totals+alpha)*log(mu+exp(log_rate))))+sum(abs(qi))
    list(qi=qi,absolute_terms=absolute_terms)
  }
  donor<-function(x)scalar_terms(x)$qi
  list(ctx=ctx, prep=prep, p=p, free=free, fixed=fixed, compose=compose,
       evaluate=evaluate, donor=donor,
       lower=c(rep(-100,p),log(1e-4))[free],
       upper=c(rep(100,p),log(10))[free])
}

working_pmg_certificate <- function(st, problem, derivatives=TRUE) {
  fail <- function(reason, metrics=list()) list(point=FALSE, covariance=FALSE,
                                               reason=reason, metrics=metrics)
  if (any(!is.finite(c(st$q,st$g,st$A,st$U)))) return(fail("nonfinite_free_state"))
  curvature <- ln_joint_curvature_diagnostics(st$A)
  margin <- sqrt(.Machine$double.eps)*(1+pmax(abs(problem$lower),abs(problem$upper)))
  metrics <- list(root=max(abs(st$g)/(1+sqrt(rowSums(st$U^2)))),
    additivity=max(abs(rowSums(st$U)-st$g)/(1+rowSums(abs(st$U)))),
    native_gradient=st$native_g,donor_gradient=st$g,
    native_gradient_absolute=max(abs(st$native_g-st$g)),
    native_gradient_agreement=max(abs(st$native_g-st$g)/(1+rowSums(abs(st$U)))),
    curvature=curvature,
    interior=all(st$x>problem$lower+margin & st$x<problem$upper-margin),
    scalar_identity=abs(st$q-st$native_q)/st$scalar_tolerance,
    stable_objective=st$q,native_objective=st$native_q,scalar_tolerance=st$scalar_tolerance)
  if (!isTRUE(metrics$interior) || !isTRUE(curvature$positive_definite) ||
      curvature$scaled_min_eigenvalue<=1e-8 || curvature$scaled_condition_number>=1e10 ||
      metrics$root>1e-6 || metrics$additivity>1e-10 || metrics$native_gradient_agreement>1e-10 || metrics$scalar_identity>1)
    return(fail("free_root_or_curvature",metrics))
  delta <- solve(st$A,st$g)
  metrics$newton <- max(abs(delta)/(1+abs(st$x)))
  if (metrics$newton>1e-6) return(fail("free_newton_displacement",metrics))
  if (derivatives) {
    d<-length(st$x);m<-ncol(st$U)
    Ufd <- t(matrix(ln_research_fd(problem$donor,st$x),m,d))
    Uhalf <- t(matrix(ln_research_fd(problem$donor,st$x,.5),m,d))
    Hfd <- matrix(ln_research_fd(function(x)problem$evaluate(x)$g,st$x),d,d)
    Hhalf <- matrix(ln_research_fd(function(x)problem$evaluate(x)$g,st$x,.5),d,d)
    metrics$donor_fd <- working_scaled_error(Ufd,st$U)
    metrics$donor_fd_stability <- working_scaled_error(Ufd,Uhalf)
    metrics$hessian_fd <- working_scaled_error(Hfd,st$A)
    metrics$hessian_fd_stability <- working_scaled_error(Hfd,Hhalf)
    if (metrics$donor_fd>1e-5 || metrics$donor_fd_stability>1e-4 ||
        metrics$hessian_fd>1e-4 || metrics$hessian_fd_stability>1e-4)
      return(fail("derivative_audit",metrics))
  }
  influence <- solve(st$A,st$U); V <- tcrossprod(influence)
  metrics$newton_in_se <- max(ifelse(diag(V)>0,abs(delta)/sqrt(diag(V)),
                                    ifelse(delta==0,0,Inf)))
  cov_ok <- all(is.finite(V)) && metrics$newton_in_se<=1e-4
  list(point=TRUE,covariance=cov_ok,reason=if(cov_ok)"certified" else "covariance_not_certified",
       metrics=metrics,V=V,influence=influence)
}

working_pmg_trial <- function(start, problem) {
  x <- pmin(pmax(start,problem$lower),problem$upper)
  initial <- x
  trace <- list()
  answer <- tryCatch({
    opt <- nlminb(x,function(a)problem$evaluate(a)$q,
      function(a)problem$evaluate(a)$g, lower=problem$lower,upper=problem$upper,
      control=list(iter.max=500L,eval.max=500L,rel.tol=1e-12,x.tol=1e-10))
    x <- opt$par; optimizer_endpoint <- problem$evaluate(x)
    for (it in seq_len(20L)) {
      st <- problem$evaluate(x)
      if (!isTRUE(ln_joint_curvature_diagnostics(st$A)$positive_definite)) break
      step <- solve(st$A,st$g)
      trace[[length(trace)+1L]] <- list(x=x,q=st$q,step=step)
      if (max(abs(step)/(1+abs(x)))<1e-11) break
      accepted <- FALSE
      for (rate in 2^-(0:15)) {
        proposal <- x-rate*step
        if(any(proposal<=problem$lower|proposal>=problem$upper)) next
        newer <- problem$evaluate(proposal)
        if (is.finite(newer$q) && newer$q<=st$q+1e-12*(1+abs(st$q))) {
          x<-proposal;accepted<-TRUE;break
        }
      }
      if (!accepted) break
    }
    st <- problem$evaluate(x)
    list(success=TRUE,state=st,optimizer_endpoint=optimizer_endpoint,
         optimizer=opt[c("convergence","message","iterations","evaluations")])
  },error=function(e)list(success=FALSE,reason="working_evaluation_failed",error=conditionMessage(e)))
  answer$start<-initial;answer$trace<-trace;answer
}

working_pmg_refit <- function(y,X,id,offset,native=NULL,candidates=NULL,allow_fixed=TRUE) {
  full <- working_pmg_problem(y,X,id,offset)
  p<-full$p
  if(is.null(native)) {
    posv<-call_posindy(full$prep$count,0L,full$ctx$nind)
    native<-fit_gene_pmg_unconstrained(1L,posv,full$ctx,control=pmg_resolve_control(),keep_candidates=TRUE)
    candidates<-native$working_candidates;native$working_candidates<-NULL
  }
  fail<-function(reason,...)list(success=FALSE,point_success=FALSE,reason=reason,
                                fit=native,native_candidates=candidates,...)
  par<-c(native$theta$beta,native$theta$subVar)
  if(length(par)!=p+1L || any(!is.finite(par)) || par[p+1L]<=0)
    return(fail("invalid_native_endpoint"))
  fixed<-numeric(0)
  lower_hit<-isTRUE(native$diagnostics$at_lower[p+1L])
  if(lower_hit && allow_fixed) fixed<-c(sigma2=1e-4)
  if(isTRUE(native$diagnostics$at_upper[p+1L]) ||
     any(native$diagnostics$at_lower[seq_len(p)] | native$diagnostics$at_upper[seq_len(p)]))
    return(fail("unsupported_boundary"))
  problem<-if(length(fixed))working_pmg_problem(y,X,id,offset,fixed) else full
  x<-c(par[seq_len(p)],log(par[p+1L]))[problem$free]
  other<-x;other[1L]<-other[1L]+.1
  if(!length(fixed))other[p+1L]<-log(max(2e-4,min(5,par[p+1L]*.7)))
  trials<-lapply(list(x,other),working_pmg_trial,problem=problem)
  if(any(!vapply(trials,function(t)isTRUE(t$success),logical(1))))
    return(fail("working_start_failed",trials=trials,fixed=fixed))
  certificates<-lapply(trials,function(t)working_pmg_certificate(t$state,problem))
  if(any(!vapply(certificates,function(c)isTRUE(c$point),logical(1))))
    return(fail("working_certificate_failed",trials=trials,certificates=certificates,fixed=fixed))
  winner<-which.min(vapply(trials,function(t)t$state$q,numeric(1)))
  st<-trials[[winner]]$state;cert<-certificates[[winner]]
  agreement<-max(vapply(trials,function(t)max(abs(t$state$x-st$x)/(1+abs(st$x))),numeric(1)))
  if(agreement>1e-6)return(fail("working_starts_disagree",trials=trials,certificates=certificates,fixed=fixed))
  # Reevaluate every retained native endpoint in its original problem, then in
  # this fixed mask if applicable. Failed attempts remain recorded and veto.
  if(!length(candidates))return(fail("native_candidate_evidence_missing"))
  competitor<-lapply(candidates,function(c) {
    if(!is.null(c$error)||is.null(c$par))return(list(success=FALSE,reason="native_attempt_failed"))
    natural<-as.numeric(c$par)
    if(length(natural)!=p+1L||any(!is.finite(natural))||natural[p+1L]<=0)
      return(list(success=FALSE,reason="invalid_native_candidate"))
    z<-c(natural[seq_len(p)],log(natural[p+1L]))
    original<-full$evaluate(z)$q
    q<-problem$evaluate(z[problem$free])$q
    list(success=is.finite(q)&&q>=st$q-1e-8*(1+abs(q)),original_q=original,
         working_q=q,delta=q-st$q,par=natural)
  })
  if(any(!vapply(competitor,function(a)isTRUE(a$success),logical(1))))
    return(fail("retained_competitor_not_resolved",trials=trials,competitors=competitor,fixed=fixed))
  # Every retained working optimizer/polish endpoint is also a competitor.
  # These values all use this identical free-coordinate criterion and mask.
  working_competitors<-list()
  for(i in seq_along(trials)) {
    points<-c(list(trials[[i]]$optimizer_endpoint$x),
              lapply(trials[[i]]$trace,`[[`,"x"),list(trials[[i]]$state$x))
    for(j in seq_along(points)) {
      q<-problem$evaluate(points[[j]])$q
      ok<-is.finite(q)&&q>=st$q-1e-10*(1+abs(st$q))
      working_competitors[[length(working_competitors)+1L]]<-list(
        trial=i,endpoint=j,x=points[[j]],q=q,delta=q-st$q,success=ok)
      if(!ok)return(fail("working_competitor_not_resolved",trials=trials,
        working_competitors=working_competitors,fixed=fixed))
    }
  }
  boundary<-NULL
  if(length(fixed)) {
    z<-c(st$par[seq_len(p)],log(st$par[p+1L]));at<-full$evaluate(z)
    inside<-z;inside[p+1L]<-inside[p+1L]+log(2)
    qin<-full$evaluate(inside)$q
    fd<-ln_research_fd(function(w)full$evaluate(w)$q,z)[p+1L]
    outward<-at$g[p+1L]/(1+sqrt(sum(at$U[p+1L,]^2)))
    boundary<-list(outward=outward,inward_difference=qin-at$q,derivative_fd=fd,
      derivative_error=abs(fd-at$g[p+1L])/(1+abs(at$g[p+1L])))
    if(outward< -1e-6 || boundary$inward_difference< -1e-10*(1+abs(at$q)) ||
       boundary$derivative_error>1e-5)
      return(fail("working_mask_direction_not_supported",boundary=boundary,fixed=fixed,trials=trials))
  }
  list(success=TRUE,point_success=TRUE,inference_available=cert$covariance,
       reason=cert$reason,fit=native,par=st$par,beta=st$par[seq_len(p)],
       Vbeta=cert$V[seq_len(p),seq_len(p),drop=FALSE],A=st$A,U=st$U,V=cert$V,
       fixed=fixed,free=problem$free,objective=st$q,ctx=problem$ctx,
       trials=trials,certificates=certificates,competitors=competitor,
       native_candidates=candidates,working_competitors=working_competitors,
       start_agreement=agreement,boundary=boundary)
}
