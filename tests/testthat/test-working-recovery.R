working_fixture <- function(combined=FALSE) {
  # A preceding test may use a different generator. Freeze the intended data
  # locally, without changing the caller's RNG kind, seed or seed absence.
  kind<-RNGkind()
  old<-if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) .Random.seed else NULL
  on.exit({
    do.call(RNGkind,as.list(kind))
    if(is.null(old))rm('.Random.seed',envir=.GlobalEnv) else
      assign('.Random.seed',old,.GlobalEnv)
  })
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  set.seed(if(combined)4003 else 90211)
  m<-if(combined)16L else 60L;k<-if(combined)6L else 8L
  id<-rep(seq_len(m),each=k);x<-rnorm(length(id));z<-rep(rnorm(m),each=k)
  y<-if(combined)rbinom(length(id),4,plogis(.1+.3*x+.15*z)) else
    rpois(length(id),exp(.3+.2*x+.1*z+rep(rnorm(m,0,.5),each=k)))
  X<-cbind(Intercept=1,x=x,region=0,z=z)
  list(y=c(y,0,0),id=c(id,m+1:2),X=rbind(X,c(1,0,1,0),c(1,0,1,0)))
}
working_call <- function(d,model,...) {
  notices<-list()
  out<-withCallingHandlers(scrope_working(d$y,d$id,d$X,model=model,
    cpc=0,mincp=1,verbose=FALSE,...),scrope_recovery_warning=function(w) {
      notices[[length(notices)+1L]]<<-w;invokeRestart("muffleWarning")
    })
  expect_lte(length(notices),1L)
  if(length(notices))expect_identical(notices[[1]]$code,"scrope_recovery_summary")
  attention<-any(!out$routes$native_regular|out$routes$wald_coefficients<out$routes$total_coefficients)
  expect_length(notices,as.integer(attention))
  out
}
working_expect_native_core <- function(actual,expected) {
  exact<-c("convergence","algorithm","random_effect")
  expect_identical(actual[exact],expected[exact])
  expect_recovery_refit_table(actual$summary,expected$summary)
  expect_recovery_refit_table(actual$overdispersion,expected$overdispersion,
    numerical=intersect(c("Subject","Cell"),names(expected$overdispersion)))
}

test_that('working fixtures are invariant to caller RNG kind and restore state', {
  local({
    kind<-RNGkind()
    old<-if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) .Random.seed else NULL
    on.exit({
      do.call(RNGkind,as.list(kind))
      if(is.null(old))rm('.Random.seed',envir=.GlobalEnv) else
        assign('.Random.seed',old,.GlobalEnv)
    })
    for(combined in c(FALSE,TRUE)) {
      reference<-working_fixture(combined)
      for(generator in c('Mersenne-Twister',"L'Ecuyer-CMRG")) {
        RNGkind(generator,'Inversion','Rejection');set.seed(17)
        before_kind<-RNGkind();before_seed<-.Random.seed
        expect_identical(working_fixture(combined),reference)
        expect_identical(RNGkind(),before_kind)
        expect_identical(.Random.seed,before_seed)
        rm('.Random.seed',envir=.GlobalEnv)
        expect_identical(working_fixture(combined),reference)
        expect_identical(RNGkind(),before_kind)
        expect_false(exists('.Random.seed',.GlobalEnv,inherits=FALSE))
      }
    }
  })
})

test_that("opt-in limiting outputs preserve native objects and disclose estimability", {
  d<-working_fixture()
  for(model in c("HL","PGMM")) {
    a<-working_call(d,model)
    b<-working_call(d,model,keep_diagnostics=TRUE)
    native_warnings<-list()
    native<-withCallingHandlers(scrope(d$y,d$id,d$X,model=model,cpc=0,mincp=1,
      verbose=FALSE,keep_diagnostics=FALSE,recovery="none"),warning=function(w) {
        native_warnings[[length(native_warnings)+1L]]<<-w
        invokeRestart("muffleWarning")
      })
    attention<-native$availability$status %in% c("partial","point_only","unavailable")
    expect_length(native_warnings,as.integer(any(attention)))
    for(w in native_warnings)expect_s3_class(w,"scrope_recovery_warning")
    expect_null(a$primary)
    working_expect_native_core(b$primary,native)
    # These are independent recovery fits; discrete output and missingness
    # stay exact while only the named fitted values receive a roundoff budget.
    expect_recovery_refit_table(a$summary,b$summary,
      numerical=c("beta","robust_se","p_value","model_se","model_p_value"))
    expect_true(all(a$summary$route=="identifiable_nuisance_limit"))
    expect_false(a$summary$point_available[a$summary$coefficient=="region"])
    expect_false(any(a$summary$selection_adjusted))
    expect_null(a$diagnostics)
    if(model=="HL") {
      expect_equal(sum(a$summary$point_available),3L)
      expect_false(any(a$summary$available))
      expect_false(any(a$summary$model_available))
      expect_true(all(is.na(a$summary$model_se)))
      expect_true(all(is.na(a$summary$robust_se)))
      expect_true(all(is.na(a$summary$p_value)))
      expect_identical(unique(a$summary$failure[a$summary$estimable]),"hl_limit_sampling_target_unresolved")
      expect_identical(unique(a$summary$failure[!a$summary$estimable]),"nonestimable_coefficient")
    } else expect_equal(sum(a$summary$available),2L)
    changed<-d;changed$X<-d$X[,c("z","Intercept","x","region")]
    changed$X<-sweep(changed$X,2,c(-3,-2,5,2),`*`)
    transformed<-working_call(changed,model)
    s<-transformed$summary[match(a$summary$coefficient,transformed$summary$coefficient),]
    factors<-c(Intercept=if(model=="HL")-2 else 1,x=5,region=2,z=-3)[a$summary$coefficient]
    expect_equal(a$summary$beta,unname(s$beta*factors),tolerance=1e-7)
    expect_equal(a$summary$robust_se,unname(s$robust_se*abs(factors)),tolerance=1e-7)
    expect_identical(a$summary$estimable,s$estimable)
  }
})

test_that("combined PGMM recovery actually fixes and refits the free system", {
  d<-working_fixture(TRUE);a<-working_call(d,"PGMM",keep_diagnostics=TRUE)
  expect_true(all(a$summary$route=="identifiable_limit_fixed_dispersion"))
  expect_equal(a$routes$fixed_donor_parameter,1e-4)
  f<-a$diagnostics[[1]]
  expect_equal(f$fixed,c(sigma2=1e-4))
  expect_equal(nrow(f$A),3L)
  expect_equal(f$V,tcrossprod(solve(f$A,f$U)),tolerance=1e-12)
  expect_equal(f$model_covariance,solve(f$A),tolerance=1e-10)
  expect_equal(f$model_V_report,f$mapping%*%solve(f$A)%*%t(f$mapping),tolerance=1e-10)
  expect_identical(a$summary$model_available,a$summary$available)
  expect_true(all(vapply(f$certificates,function(c)c$point&&c$covariance,logical(1))))
  expect_true(all(vapply(f$working_competitors,`[[`,logical(1),"success")))
  expect_true(all(vapply(f$full_competitors,`[[`,logical(1),"success")))
  expect_true(all(vapply(f$certificates,function(c)c$metrics$root<=1e-6,logical(1))))
})

test_that("the stable PGMM scalar equals the finite-product marginal formula", {
  d<-working_fixture(TRUE);keep<-d$id<=16L;X<-d$X[keep,c("Intercept","x","z")]
  problem<-scRoPE:::working_pmg_problem(d$y[keep],X,d$id[keep],rep(1,sum(keep)))
  ctx<-problem$ctx;yy<-as.numeric(problem$prep$count)
  for(v in c(1e-4,.02,.8,3)) {
    z<-c(.2,-.1,.3,log(v));a<-1/expm1(v);lambda<-a*exp(-v/2)
    eta<-as.numeric(ctx$offset+ctx$pred%*%z[1:3]);S<-ctx$cumsumy[1,]
    mu<-as.numeric(rowsum(exp(eta),ctx$id,reorder=FALSE))
    normalized<-vapply(S,function(s)if(s==0)0 else sum(log1p((0:(s-1))/a)),numeric(1))
    expected<--normalized-S*v/2+(S+a)*log1p(mu/lambda)-
      as.numeric(rowsum(yy*eta,ctx$id,reorder=FALSE))
    expect_equal(problem$donor(z),expected,tolerance=1e-12)
    st<-problem$evaluate(z)
    expect_equal(st$q,sum(expected),tolerance=1e-12)
    expect_identical(st$g,rowSums(st$U))
    if(v>.001) {
      expect_lte(max(abs(st$native_g-st$g)/(1+rowSums(abs(st$U)))),1e-10)
      # Away from a root, this independently checks the nonzero score term
      # in the log-v Hessian transformation as well as the donor-summed score.
      Hfd<-scRoPE:::ln_research_fd(function(w)problem$evaluate(w)$g,z)
      expect_lte(scRoPE:::working_scaled_error(Hfd,st$A),1e-4)
    } else {
      # A boundary is not eligible as a free-v root. The actual fixed route
      # excludes that noisy held score from its free equations and covariance.
      certificate<-scRoPE:::working_pmg_certificate(st,problem,FALSE)
      expect_false(certificate$point)
      expect_false(certificate$metrics$interior)
      fixed<-scRoPE:::working_pmg_problem(d$y[keep],X,d$id[keep],rep(1,sum(keep)),c(sigma2=1e-4))
      sf<-fixed$evaluate(z[1:3])
      expect_lte(max(abs(sf$native_g-sf$g)/(1+rowSums(abs(sf$U)))),1e-10)
    }
    expect_lte(abs(st$q-st$native_q),st$scalar_tolerance)
  }
})

test_that("unsupported face geometry and no outcomes are not reported", {
  d<-working_fixture();d$id[length(d$id)]<-1L
  f<-scRoPE:::working_face_detect(d$y,d$X,d$id)
  expect_false(f$success)
  expect_identical(f$reason,"partial_donor_face_outside_scope")
  d<-working_fixture();d$X<-cbind(d$X,extra=0)
  d$X[nrow(d$X),"region"]<-0;d$X[nrow(d$X),"extra"]<-1
  f<-scRoPE:::working_face_detect(d$y,d$X,d$id)
  expect_false(f$success)
  expect_error(scrope_working(rep(0,20),rep(1:5,each=4),model="PGMM",verbose=FALSE))
})

test_that("regular native results stay exact and compact audits discard donor arrays", {
  d<-working_fixture();keep<-d$id<=60L
  d<-list(y=d$y[keep],id=d$id[keep],X=d$X[keep,c("Intercept","x","z")])
  for(model in c("HL","PGMM")) {
    a<-working_call(d,model)
    b<-working_call(d,model,keep_diagnostics=TRUE)
    native<-scrope(d$y,d$id,d$X,model=model,cpc=0,mincp=1,keep_diagnostics=FALSE,verbose=FALSE,recovery="none")
    expect_null(a$primary)
    working_expect_native_core(b$primary,native)
    expect_null(a$audit[[1]]$native_certificate$influence)
    expect_null(a$audit[[1]]$native_certificate$V)
    expect_false(a$audit[[1]]$retained_details)
    expect_null(a$audit[[1]]$history_audit$observations)
    expect_null(a$audit[[1]]$history_evidence$events)
    expect_null(a$diagnostics)
    if(a$routes$native_regular) {
      # Exact preservation is checked against the native fit from this call.
      expect_identical(b$summary$beta,as.numeric(b$primary$summary[1,grep('^logFC_',names(b$primary$summary))]))
      expect_identical(b$summary$robust_se,as.numeric(b$primary$summary[1,grep('^se_robust_',names(b$primary$summary))]))
    }
  }
})

test_that("HL combined fixed and limiting route audits the augmented free system", {
  d<-working_fixture(TRUE);a<-working_call(d,"HL",keep_diagnostics=TRUE)
  expect_true(all(a$summary$route=="identifiable_limit_fixed_dispersion"))
  expect_equal(a$routes$fixed_donor_parameter,1e-4)
  expect_equal(a$routes$fixed_nb_precision,1000)
  expect_false(a$summary$point_available[a$summary$coefficient=="region"])
  expect_equal(sum(a$summary$available),3L)
  f<-a$diagnostics[[1]]
  expect_length(f$fit$active,0L)
  expect_true(f$model_inference_available)
  expected_model<-f$mapping%*%f$fit$K%*%solve(-f$fit$augmented$B)%*%t(f$fit$K)%*%t(f$mapping)
  expect_equal(f$model_V_report,expected_model,tolerance=1e-10)
  expect_identical(a$summary$model_available,a$summary$available)
  expect_equal(f$certificate$V,tcrossprod(solve(f$fit$augmented$A,f$fit$augmented$psi_subject)),tolerance=1e-12)
  expect_lte(f$certificate$metrics$summed_jacobian_fd,1e-4)
  expect_true(all(vapply(f$terminal_audit,`[[`,logical(1),"success")))
})

test_that("a failed retained PGMM attempt cannot silently certify a recovery", {
  d<-working_fixture(TRUE);keep<-d$id<=16L
  X<-d$X[keep,c("Intercept","x","z")];y<-d$y[keep];id<-d$id[keep]
  problem<-scRoPE:::working_pmg_problem(y,X,id,rep(1,length(y)))
  posv<-scRoPE:::call_posindy(problem$prep$count,0L,problem$ctx$nind)
  native<-scRoPE:::fit_gene_pmg_unconstrained(1L,posv,problem$ctx,
    control=scRoPE:::pmg_resolve_control(),keep_candidates=TRUE)
  candidates<-native$working_candidates
  candidates[[length(candidates)+1L]]<-list(start=rep(0,4),error="retained backend failed")
  a<-scRoPE:::working_pmg_refit(y,X,id,rep(1,length(y)),native,candidates)
  expect_false(a$point_success)
  expect_identical(a$reason,"retained_competitor_not_resolved")
  expect_identical(tail(a$native_candidates,1)[[1]]$error,"retained backend failed")
})

test_that("covariance agreement is relative even for tiny variances", {
  a<-diag(c(4e-16,2));b<-a;b[1,1]<-b[1,1]*1.01506
  expect_equal(scRoPE:::working_covariance_discrepancy(b,a),.01506,tolerance=1e-12)
  expect_equal(scRoPE:::working_covariance_discrepancy(diag(c(0,2)),diag(c(0,2))),0)
  expect_identical(scRoPE:::working_covariance_discrepancy(diag(c(1e-40,2)),diag(c(0,2))),Inf)
})


test_that("intercept-only PGMM recovery retains one-dimensional derivative shapes", {
  id<-rep(1:30,each=5)
  target_totals<-rep(8:12,6)
  y<-unlist(lapply(target_totals,function(s)rep(s%/%5,5)+(seq_len(5)<=s%%5)))
  expect_warning(f<-scrope_working(y,id,model="PGMM",cpc=0,mincp=1,verbose=FALSE),class="scrope_recovery_warning")
  expect_equal(nrow(f$summary),1L)
  expect_identical(f$routes$route,"fixed_positive_dispersion")
  expect_equal(f$routes$fixed_donor_parameter,1e-4)
  expect_true(f$summary$point_available)
  expect_true(f$summary$available)
  totals<-as.numeric(rowsum(y,id));average<-mean(totals)
  expect_equal(f$summary$beta,log(mean(y))-1e-4/2,tolerance=1e-8)
  expect_equal(f$summary$robust_se,sqrt(sum((totals-average)^2))/(30*average),tolerance=1e-8)
  expect_true(all(vapply(f$audit[[1]]$certificates,function(c)c$metrics$hessian_fd<=1e-4,logical(1))))
})

test_that("recovered coefficient refusals do not inherit a fit certificate label", {
  d<-working_fixture();normal<-working_call(d,"PGMM")
  recover<-scRoPE:::working_recover
  scenarios<-list(
    zero_variance=list(variance=0,reason="mapped_variance_not_positive_finite"),
    negative_variance=list(variance=-1,reason="mapped_variance_not_positive_finite"),
    infinite_variance=list(variance=Inf,reason="mapped_variance_not_positive_finite"),
    missing_variance=list(variance=NA_real_,reason="mapped_variance_not_positive_finite"),
    missing_point=list(beta=NA_real_,reason="nonfinite_coefficient_estimate"),
    overflowing_wald=list(beta=1e200,variance=1,reason="nonfinite_wald_statistic"))
  run_scenario<-function(change) {
    testthat::local_mocked_bindings(working_recover=function(...) {
      ans<-recover(...)
      if(!is.null(change$variance))ans$V_report[2,2]<-change$variance
      if(!is.null(change$beta))ans$beta_report[2]<-change$beta
      ans
    },.package="scRoPE")
    working_call(d,"PGMM")
  }
  for(change in scenarios) {
    altered<-run_scenario(change);x<-altered$summary[altered$summary$coefficient=="x",]
    expect_false(x$available)
    expect_identical(x$failure,change$reason)
    expect_true(is.na(x$p_value))
    if(identical(change$reason,"nonfinite_wald_statistic"))expect_equal(x$robust_se,1) else
      expect_true(is.na(x$robust_se))
    expect_identical(altered$primary,normal$primary)
    expect_recovery_refit_table(altered$summary[altered$summary$coefficient!="x",],
      normal$summary[normal$summary$coefficient!="x",],
      numerical=c("beta","robust_se","p_value","model_se","model_p_value"))
    expect_identical(x$point_available,is.null(change$beta)||is.finite(change$beta))
    if(identical(change$reason,"nonfinite_wald_statistic"))expect_true(x$covariance_available)
  }
  for(fit_reason in c("certified","covariance_not_certified")) {
    altered<-local({
      testthat::local_mocked_bindings(working_recover=function(...) {
        ans<-recover(...);ans$inference_available<-FALSE;ans$reason<-fit_reason;ans
      },.package="scRoPE")
      working_call(d,"PGMM")
    })
    expected<-if(fit_reason=="certified")"working_covariance_unavailable" else fit_reason
    expect_false(any(altered$summary$available))
    expect_identical(unique(altered$summary$failure[altered$summary$estimable]),expected)
    expect_recovery_refit_numeric(altered$summary$beta,normal$summary$beta,absolute=1e-12)
    expect_identical(altered$primary,normal$primary)
  }
})
