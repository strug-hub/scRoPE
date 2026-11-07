test_that("default recovery preserves successful native estimates and supports strict opt-out", {
  fx<-readRDS(test_path('fixtures','jln-research.rds'))
  for(nm in c('ordinary','zero','warning','corner')){
    x<-fx[[nm]]$input
    a<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE))
    n<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,ln_recovery='none'))
    expect_recovery_refit_table(a$recovery$primary$summary,n$summary)
    expect_recovery_refit_table(a$recovery$primary$overdispersion,n$overdispersion,
      numerical=c("Subject","Cell"))
    expect_true(all(a$recovery$coefficients$available))
    expect_false(any(a$recovery$routes$selection_adjusted))
    if(nm=='ordinary')expect_recovery_refit_table(a$summary,n$summary)
    else {
      expect_true(a$recovery$routes$added)
      expect_true(all(is.na(a$summary[grep('^se_',names(a$summary),value=TRUE)[!grepl('robust',grep('^se_',names(a$summary),value=TRUE))]])))
      f<-a$recovery$diagnostics[[1]]$result
      expect_equal(f$covariance,tcrossprod(solve(f$state$A,f$state$U)),tolerance=1e-12)
    }
  }
})

test_that("whole-zero-donor face reports only identifiable coefficients", {
  x<-readRDS(test_path('fixtures','jln-face.rds'))
  a<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE))
  expect_identical(a$recovery$routes$route,'identifiable_nuisance_limit')
  f<-a$recovery$diagnostics[[1]]$result
  expect_identical(f$face$estimable,c(FALSE,TRUE,FALSE,TRUE))
  expect_identical(a$recovery$coefficients$available,c(FALSE,TRUE,FALSE,TRUE))
  expect_true(all(is.na(a$summary[c('logFC_(Intercept)','logFC_region','p_robust_region')])) )
  expect_false(a$provenance$per_gene$certified_root)
  expect_identical(a$availability$availability,'partial_inference')
  expect_true(all(f$U_original[,!seq_len(f$face$original_donors)%in%f$face$retained_donor_indices]==0))
  expect_equal(tcrossprod(f$U_original),tcrossprod(f$U),tolerance=1e-12)
  expect_equal(f$covariance,tcrossprod(solve(f$A,f$U)),tolerance=1e-12)
  # The reduced/oracle comparison checks coefficient and SE mapping across centering.
  keep<-x$id<=52
  o<-suppressWarnings(scrope(x$y[keep],x$id[keep],x$X[keep,c(1,2,4)],cpc=0,mincp=1,verbose=FALSE))
  expect_equal(a$summary$logFC_x,o$summary$logFC_x,tolerance=1e-10)
  expect_equal(a$summary$se_robust_x,o$summary$se_robust_x,tolerance=1e-10)
  compact<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=FALSE))
  expect_recovery_refit_table(compact$summary,a$summary)
  expect_null(compact$diagnostics);expect_null(compact$recovery$diagnostics)
  expect_null(compact$recovery$primary$diagnostics)
})

test_that("face detection rejects mixed donors, approximate relations and nonunique nullspaces", {
  x<-readRDS(test_path('fixtures','jln-face.rds'));detect<-scRoPE:::ln_recovery_face_detect
  expect_true(detect(x$y,x$X,x$id)$success)
  mixed<-x$id;mixed[length(mixed)-1L]<-52L
  expect_identical(detect(x$y,x$X,mixed)$reason,'partial_donor_face_outside_scope')
  near<-x$X;near[1,3]<-1e-8
  expect_false(detect(x$y,near,x$id)$success)
  multi<-cbind(x$X,second_face=c(rep(0,length(x$y)-2),-1,-2))
  expect_identical(detect(x$y,multi,x$id)$reason,'no_unique_one_dimensional_face')
  positive<-x$y;positive[length(positive)]<-1
  expect_false(detect(positive,x$X,x$id)$success)
})

test_that("limiting derivatives and covariance match independent scalar differences", {
  x<-readRDS(test_path('fixtures','jln-face.rds'))
  f<-scRoPE:::ln_recovery_face_fit(x$y,x$X,x$id,x$offset)
  expect_true(f$success)
  keep<-f$face$retained;J<-f$face$columns
  dat<-scRoPE:::ln_research_context(x$y[keep],x$X[keep,J],x$id[keep],x$offset[keep])
  r<-length(J);free<-f$free;full<-f$reduced_full
  qdon<-function(a){t<-full;t[free]<-a;t[(r+1):(r+2)]<-exp(t[(r+1):(r+2)])
    rowSums(scRoPE:::ln_joint_objective_components(t,1L,dat$posv,dat$ctx,dat$map))}
  a<-full[free]
  # numDeriv uses a separate Richardson implementation, not the fitting helper.
  Ufd<- -t(numDeriv::jacobian(qdon,a));Hfd<-numDeriv::hessian(function(a)sum(qdon(a)),a)
  expect_equal(Ufd,f$U,tolerance=1e-5)
  expect_equal(Hfd,f$A,tolerance=1e-4)
  expect_equal(rowSums(f$U),rep(0,nrow(f$U)),tolerance=1e-5)
  full_dat<-scRoPE:::ln_research_context(x$y,x$X,x$id,x$offset);p<-ncol(x$X)
  t<-c(as.numeric(f$full_beta_map%*%full[seq_len(r)]),tail(full,2))
  step<-50/min(-f$face$direction_predictor[f$face$removed]);t[seq_len(p)]<-t[seq_len(p)]+step*f$face$direction
  st<-scRoPE:::ln_joint_derivative_state(t,1L,full_dat$posv,full_dat$ctx,full_dat$map)
  expect_equal(st$raw$value,f$q,tolerance=1e-9)
  expect_equal(crossprod(f$derivative_map,st$A_full%*%f$derivative_map),f$A,tolerance=1e-8)
  expect_equal(crossprod(f$derivative_map,st$U_subject_full),f$U_original,tolerance=1e-8)
})

test_that("zero donor objective identity is positive and converges without changing criterion", {
  x<-readRDS(test_path('fixtures','jln-face.rds'));dat<-scRoPE:::ln_research_context(x$y,x$X,x$id,x$offset)
  p<-ncol(x$X);theta<-c(rep(0,p),.2,1);alpha<-1/expm1(.2);lambda<-alpha*exp(-.1)
  for(b in c(-2,-10,-30)){
    theta[1]<-b;mu<-exp(dat$ctx$offset+dat$ctx$pred%*%theta[seq_len(p)])
    q<-rowSums(scRoPE:::ln_joint_objective_components(theta,1L,dat$posv,dat$ctx,dat$map))
    for(i in c(53,54)){
      m<-mu[dat$id==i];S<-sum(m);rr<-S/(S+lambda)
      exact<-alpha*(-log1p(-rr)-rr)+sum(log1p(alpha*rr*m/S))
      expect_gte(exact,0);expect_equal(unname(q[i]),exact,tolerance=1e-12)
    }
  }
})


test_that("recovery never presents native optional tests or random effects as recovered", {
  x<-readRDS(test_path('fixtures','jln-research.rds'))$zero$input
  a<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,
    additional_tests=c('score','lrt'),lrt_details=TRUE,output_re=TRUE))
  expect_true(a$recovery$routes$added)
  expect_true(all(is.na(a$random_effect)))
  ts<-a$test_status[a$test_status$test!='wald',]
  expect_true(all(!ts$success));expect_true(all(ts$failure_reason=='not_supported_for_recovery'))
  expect_true(all(a$lrt_details$status=='not_supported_for_recovery'))
  expect_true(all(!a$lrt_details$success))
  expect_true(all(is.na(a$lrt_details$p_value)))
})

test_that("local refinement resolves an all-free root and rejects unrelated failures", {
  x<-readRDS(test_path('fixtures','jln-research.rds'))$ordinary$input
  n<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,ln_recovery='none',keep_diagnostics=TRUE))
  primary<-list(joint=n$diagnostics$raw[[1]]$ln_joint)
  dat<-scRoPE:::ln_research_context(x$y,x$X,x$id,x$offset)
  r<-primary$joint$candidate_records[[1]]
  r$theta_reduced_log[1]<-r$theta_reduced_log[1]+1e-5
  r$certified<-FALSE;r$checks['newton']<-FALSE
  primary$joint$candidate_records[[1]]<-r
  f<-scRoPE:::ln_recovery_refine(primary,dat)
  expect_true(f$success)
  expect_length(f$fixed,0L)
  expect_identical(f$spec$free,seq_len(ncol(x$X)+2L))
  expect_equal(f$coefficient,n$summary$logFC_x,tolerance=1e-7)
  expect_equal(f$se,n$summary$se_robust_x,tolerance=1e-7)
  expect_lt(f$metrics['newton_in_se'],1e-4)
  primary$joint$candidate_records[[1]]$interior<-FALSE
  expect_false(scRoPE:::ln_recovery_refine(primary,dat)$success)
})

test_that("a known better full competitor prevents a misleading limiting report", {
  x<-readRDS(test_path('fixtures','jln-face-competitor.rds'))
  f<-suppressWarnings(scrope(x$y,x$id,x$X,offset=x$offset,cpc=0,mincp=1,verbose=FALSE,keep_diagnostics=TRUE))
  expect_identical(f$recovery$routes$route,'withheld')
  expect_identical(f$recovery$routes$reason,'reduced_root_worse_than_known_full_candidate')
  face<-f$recovery$diagnostics[[1]]$face_audit
  expect_gt(face$q,face$known_objective_min)
  expect_false(f$summary$wald_success_x)
  expect_true(is.na(f$summary$p_robust_x))
})
