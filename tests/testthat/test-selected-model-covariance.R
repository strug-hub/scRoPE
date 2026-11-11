# Independent observed-information algebra at the very same selected point.
# HL formula is master (2026-11-06), profile p'' = h'' - h_eta h_eta'/h_etaeta;
# pinned NEBULA 70593f9 opt_pml computes vb2 = vb - vwb' diag(1/vw) vwb.
.model_hl_information <- function(f) {
  d<-f$data;v<-exp(f$rho_full[1]);phi<-exp(f$rho_full[2])
  lambda<-exp(-v/2)/expm1(v)
  mu<-exp(d$log_offset+as.numeric(d$X%*%f$beta)+f$eta[d$subject])
  weight<-phi*(d$y+phi)*mu/(mu+phi)^2
  Ibb<-crossprod(d$X,d$X*weight)
  Ibw<-rowsum(d$X*weight,d$subject,reorder=FALSE)
  Iww<-as.numeric(rowsum(weight,d$subject,reorder=FALSE))+lambda*exp(f$eta)
  Ibb-crossprod(Ibw,Ibw/Iww)
}

test_that('HL model information profiles modes and plugs in the selected dispersion', {
  x<-hl_fo_fixture('R1')
  d<-scRoPE:::hl_fo_make_data(x$y,x$X,x$log_pi,x$subject_index,x$intcol)
  for(mask in list(list(active=c('tau','xi'),fixed=numeric()),
                   list(active='xi',fixed=c(sigma2=.4)),
                   list(active=character(),fixed=c(sigma2=.4,phi=8)))) {
    f<-scRoPE:::hl_fo_fit(d,active=mask$active,fixed_natural=mask$fixed)
    expect_true(f$point_success);expect_true(f$model_covariance_success)
    information<-.model_hl_information(f)
    expect_equal(-f$augmented$B,crossprod(f$K,information%*%f$K),tolerance=1e-10)
    expected<-f$K%*%solve(crossprod(f$K,information%*%f$K))%*%t(f$K)
    expect_equal(f$model_covariance$V_beta,expected,tolerance=1e-10)
    expect_equal(f$covariance$sandwich$V,tcrossprod(solve(f$augmented$A,f$augmented$psi_subject)),tolerance=1e-10)
    if(length(mask$active))expect_gt(max(abs(solve(f$augmented$A)[1:d$p,1:d$p]-expected)),1e-6)
    if(!length(mask$active)) {
      # Independent scalar NB + gamma h likelihood; re-solve every mode for
      # each perturbed coefficient, without the package profile derivatives.
      v<-exp(f$rho_full[1]);phi<-exp(f$rho_full[2]);alpha<-1/expm1(v);lambda<-alpha*exp(-v/2)
      profiled<-function(beta)sum(vapply(seq_len(d$m),function(i) {
        rows<-d$subject_rows[[i]];y<-d$y[rows]
        base<-d$log_offset[rows]+as.numeric(d$X[rows,,drop=FALSE]%*%beta)
        score<-function(eta)sum(y-(y+phi)*plogis(base+eta-log(phi)))+alpha-lambda*exp(eta)
        eta<-uniroot(score,c(f$eta[i]-15,f$eta[i]+15),tol=1e-12)$root
        sum(dnbinom(y,size=phi,mu=exp(base+eta),log=TRUE))+dgamma(exp(eta),shape=alpha,rate=lambda,log=TRUE)+eta
      },numeric(1)))
      expect_equal(unname(-numDeriv::hessian(profiled,f$beta)),unname(information),tolerance=2e-5)
    }
  }
  # A beta constraint changes the free coefficient information before inversion.
  C<-matrix(c(0,1,0),1L);f<-scRoPE:::hl_fo_fit(d,C=C,rhs=.35)
  expect_true(f$point_success);expect_true(f$model_covariance_success)
  I<-.model_hl_information(f);K<-f$K
  expect_equal(f$model_covariance$V_beta,K%*%solve(crossprod(K,I%*%K))%*%t(K),tolerance=1e-10)
  expect_equal(C%*%f$model_covariance$V_beta,matrix(0,1L,d$p),tolerance=1e-12)
})

test_that('HL model covariance is independent of unavailable donor sandwich', {
  x<-hl_fo_fixture('R1');d<-scRoPE:::hl_fo_make_data(x$y,x$X,x$log_pi,x$subject_index,x$intcol)
  f<-scRoPE:::hl_fo_fit(d)
  expected<-f$model_covariance$V_beta
  f$augmented$J[,]<-0
  g<-scRoPE:::hl_fo_finalize_state(f)
  expect_true(g$point_success);expect_false(g$covariance_success)
  expect_true(g$model_covariance_success);expect_equal(g$model_covariance$V_beta,expected,tolerance=0)
})

test_that('HL held cell precision leaves donor dispersion free in the selected fit', {
  x<-hl_fo_fixture('R1')
  d<-scRoPE:::hl_fo_make_data(x$y,x$X,x$log_pi,x$subject_index,x$intcol)
  # Validation mask, not an additional automatic boundary setting.
  f<-scRoPE:::hl_fo_fit(d,active='tau',fixed_natural=c(phi=8))
  expect_true(f$point_success);expect_true(f$model_covariance_success)
  expect_true(f$covariance_success);expect_identical(f$active,1L)
  expect_equal(exp(f$rho_full[2]),8,tolerance=1e-12)
  information<-.model_hl_information(f);K<-f$K
  expect_equal(-f$augmented$B,crossprod(K,information%*%K),tolerance=1e-10)
  expect_equal(f$model_covariance$V_beta,
    K%*%solve(crossprod(K,information%*%K))%*%t(K),tolerance=1e-10)
  expect_equal(f$covariance$sandwich$V,
    tcrossprod(solve(f$augmented$A,f$augmented$psi_subject)),tolerance=1e-10)
  expect_equal(nrow(f$augmented$A),2L*d$p+1L)
})

test_that('HL public model-only partial output retains the same native point', {
  x<-hl_fo_fixture('R1');fitter<-scRoPE:::hl_fo_fit
  testthat::local_mocked_bindings(hl_fo_fit=function(...) {
    f<-fitter(...);f$covariance_success<-FALSE;f$covariance_status<-'controlled_sandwich_failure';f
  },working_recover=function(...)list(success=FALSE,reason='controlled_recovery_failure'),.package='scRoPE')
  f<-suppressWarnings(scrope(x$y,x$subject_index,x$X,offset=x$offset,model='HL',
    cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE))
  expect_identical(f$availability$status,'partial')
  expect_identical(f$availability$selected_route,'native_partial')
  expect_true(all(f$coefficient_availability$model_available))
  expect_false(any(f$coefficient_availability$available))
  expect_identical(f$summary,f$recovery$primary$summary)
})

test_that('JLN recovered model covariance inverts the selected free curvature', {
  fx<-readRDS(test_path('fixtures','jln-research.rds'))
  for(nm in c('ordinary','zero','warning','corner')) {
    x<-fx[[nm]]$input
    a<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE))
    expect_true(all(a$coefficient_availability$model_available))
    if(nm=='ordinary') {
      st<-a$diagnostics$raw[[1]]$ln_joint$state;H<-st$A_full
      map<-a$provenance$coefficient_map_internal_to_public
    } else {
      f<-a$recovery$diagnostics[[1]]$result;H<-f$state$A
      map<-a$provenance$coefficient_map_internal_to_public
      expect_equal(f$model_covariance,solve(H),tolerance=1e-10)
      expect_equal(f$covariance,tcrossprod(solve(H,f$state$U)),tolerance=1e-10)
      expect_identical(nrow(H),if(is.null(f$spec$free))ncol(x$X)+2L else length(f$spec$free))
    }
    p<-ncol(x$X);V<-map%*%solve(H)[seq_len(p),seq_len(p),drop=FALSE]%*%t(map)
    expect_equal(a$coefficient_availability$model_se,sqrt(diag(V)),tolerance=1e-10)
    expect_equal(a$coefficient_availability$model_statistic,a$coefficient_availability$beta^2/diag(V),tolerance=1e-10)
    expect_true(all(a$test_status$success[a$test_status$test=='wald_model']))
  }
})

test_that('fixed information must be restricted before inversion', {
  # A nonzero beta/nuisance cross-block distinguishes fixed from fitted nuisance.
  H<-matrix(c(5,1,2,1,4,1,2,1,3),3L)
  fixed<-scRoPE:::scrope_model_inverse(H[1:2,1:2])$covariance
  fitted<-scRoPE:::scrope_model_inverse(H)$covariance[1:2,1:2]
  expect_equal(fixed,solve(H[1:2,1:2]),tolerance=1e-14)
  expect_gt(max(abs(fixed-fitted)),.01)
  expect_false(scRoPE:::scrope_model_inverse(diag(c(1,-1)))$success)
  expect_false(scRoPE:::scrope_model_inverse(matrix(NA_real_,2,2))$success)
})

test_that('JLN held cell overdispersion leaves donor dispersion free in the selected fit', {
  x<-readRDS(test_path('fixtures','jln-research.rds'))$ordinary$input
  native<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,
    verbose=FALSE,recovery='none',keep_diagnostics=TRUE))
  full<-native$diagnostics$raw[[1]]$ln_joint$state$full_log;p<-ncol(x$X)
  dat<-scRoPE:::ln_research_context(x$y,x$X,x$id,x$offset)
  # Public JLN phi is overdispersion, the reciprocal of NB precision.
  fixed<-list(phi=exp(-full[p+2L]))
  f<-scRoPE:::ln_research_fixed_fit(dat,fixed,list(full,full+c(rep(1e-4,p),0,0)))
  expect_true(f$success);expect_true(f$model_inference_available)
  expect_identical(f$spec$free,seq_len(p+1L))
  expect_equal(exp(-f$state$full[p+2L]),fixed$phi,tolerance=1e-12)
  model<-solve(f$state$A);robust<-tcrossprod(solve(f$state$A,f$state$U))
  expect_equal(f$model_covariance,model,tolerance=1e-10)
  expect_equal(f$model_Vbeta,dat$pub%*%model[seq_len(p),seq_len(p),drop=FALSE]%*%t(dat$pub),tolerance=1e-10)
  expect_equal(f$covariance,robust,tolerance=1e-10)
  expect_equal(f$Vbeta,dat$pub%*%robust[seq_len(p),seq_len(p),drop=FALSE]%*%t(dat$pub),tolerance=1e-10)
})

test_that('mapped model covariance and overflowing Wald have separate output status', {
  primary<-list(summary=data.frame(gene_id=1L,gene='controlled',optimizer_success=TRUE,logFC_x=0,
    se_x=1,p_naive_x=1,model_wald_stat_x=0,model_wald_success_x=TRUE,
    model_wald_failure_x=NA_character_,se_robust_x=1,p_robust_x=1,
    wald_stat_x=0,wald_success_x=TRUE,wald_failure_x=NA_character_),
    overdispersion=data.frame(Subject=.2),convergence=1L,
    diagnostics=list(godambe=list(NULL),raw=list(NULL)),provenance=list(estimator='PGMM'))
  a<-list(success=TRUE,point_success=TRUE,inference_available=TRUE,
    model_inference_available=TRUE,model_inference_reason=NA_character_,
    route='certified_free_refit',reason='certified',beta_report=1e200,
    V_report=matrix(1),model_V_report=matrix(2),estimable=TRUE,
    sampling_target='regular_working_equation_target',elapsed_seconds=0,
    par=c(1e200,.2),dimensions=c(donors=30))
  for(full in c(FALSE,TRUE)) {
    record<-if(full)a else scRoPE:::working_minimal_record(a)
    long<-scRoPE:::working_selected_output(primary,list('1'=record),'PGMM',full)$summary
    expect_true(long$point_available);expect_true(long$covariance_available)
    expect_equal(long$robust_se,1);expect_false(long$available)
    expect_true(is.na(long$p_value));expect_identical(long$failure,'nonfinite_wald_statistic')
    f<-scRoPE:::scrope_apply_working(primary,list('1'=record),'PGMM',full)
    f<-scRoPE:::scrope_finish_recovery(f,'PGMM','auto',keep_diagnostics=full,notify=FALSE)
    ca<-f$coefficient_availability
    expect_true(ca$covariance_available);expect_equal(ca$robust_se,1)
    expect_true(ca$interval_available);expect_false(ca$available)
    expect_false(ca$statistic_available);expect_false(ca$reference_available)
    expect_identical(ca$failure_reason,'nonfinite_wald_statistic')
    expect_equal(f$summary$se_robust_x,1)
    expect_true(is.na(f$summary$wald_stat_x));expect_true(is.na(f$summary$p_robust_x))
    expect_false(f$summary$model_wald_success_x)
    expect_identical(f$summary$model_wald_failure_x,'nonfinite_model_wald_statistic')
    expect_equal(f$summary$se_x,sqrt(2));expect_true(is.na(f$summary$p_naive_x))
    expect_true(ca$model_covariance_available);expect_true(ca$model_interval_available)
    expect_false(ca$model_available)
    if(full) {
      expect_equal(f$diagnostics$godambe[[1]]$model_cov,matrix(2))
      expect_equal(f$diagnostics$godambe[[1]]$robust_cov,matrix(1))
      # The retained native result gains status metadata during finalization;
      # its original estimates, inference and dispersion must remain intact.
      expect_identical(f$recovery$primary$summary[names(primary$summary)],primary$summary)
      expect_identical(f$recovery$primary$overdispersion,primary$overdispersion)
    } else {
      expect_null(f$recovery$primary);expect_null(f$recovery$diagnostics)
      expect_null(record$V_report);expect_null(record$model_V_report)
    }
  }
  a$beta_report<-1;a$model_V_report[,]<-0
  f<-scRoPE:::scrope_apply_working(primary,list('1'=a),'PGMM',TRUE)
  expect_true(f$summary$wald_success_x)
  expect_false(f$summary$model_wald_success_x)
  expect_identical(f$summary$model_wald_failure_x,'mapped_model_variance_not_positive_finite')
  expect_null(f$diagnostics$godambe[[1]]$model_cov)
})


test_that('older LN reporting schemas do not create phantom model tests', {
  sm<-data.frame(gene_id=1L,gene='old',logFC_x=.2,se_x=.1,p_naive_x=.04,
    se_robust_x=.12,p_robust_x=.1,wald_stat_x=(.2/.12)^2,wald_success_x=TRUE,
    wald_failure_x=NA_character_)
  f<-list(summary=sm,convergence=1L,provenance=list(model_covariance='legacy_post_pml_curvature',
    per_gene=data.frame(inference_available=TRUE)),test_status=data.frame(gene_id=1L,gene='old',
    contrast='x',test='wald',success=TRUE,failure_reason=NA_character_,method='old_robust'))
  a<-scRoPE:::scrope_finish_recovery(f,'LN','none',notify=FALSE)
  expect_equal(nrow(a$test_status),1L);expect_false(anyNA(a$test_status$gene_id))
  expect_identical(a$coefficient_availability$model_status,'not_reported')
  expect_identical(a$provenance$model_covariance_convention,'legacy_post_pml_curvature')
  f$test_status<-NULL
  a<-scRoPE:::scrope_finish_recovery(f,'LN','none',notify=FALSE)
  expect_equal(nrow(a$test_status),1L);expect_false(anyNA(a$test_status$gene_id))
})
