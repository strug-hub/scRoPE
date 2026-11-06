.hsr <- function(name) getFromNamespace(name, "scRoPE")
.hsr_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    kind <- RNGkind()
    old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
    on.exit({do.call(RNGkind,as.list(kind))
      if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else
        assign(".Random.seed", old, .GlobalEnv)})
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(8761); m <- 40L; ni <- 6L; id <- rep(seq_len(m), each = ni)
    X <- cbind(Intercept = 1, x = rep(rnorm(m), each = ni), w = rnorm(m*ni))
    off <- .1*sin(seq_along(id)); a <- 1/expm1(.3)
    u <- rgamma(m, a, rate = a/exp(.15))
    y <- rnbinom(m*ni, mu = exp(1+off+.2*X[,3])*u[id], size = 4)
    cached <<- list(data = .hsr("hl_fo_make_data")(y, X, off, id),
                   R = matrix(c(0, 1, 0), 1))
    cached
  }
})
.hsr_record <- local({
  cache <- new.env(parent = emptyenv())
  function(active = 1:2) {
    key <- paste0("a", paste(active, collapse = ""))
    if (!exists(key, cache, inherits = FALSE)) {
      f <- .hsr_fixture()
      cache[[key]] <- .hsr("hl_fo_score_record")(f$data, f$R, active = active,
        fixed_natural = c(sigma2 = .3, phi = 4)[setdiff(1:2, active)], keep_diagnostics = FALSE)
    }
    cache[[key]]
  }
})

test_that("installed retry preserves native successes exactly for every active mask", {
  for (active in list(1:2, 1L, 2L, integer(0))) {
    rec <- .hsr_record(active); native <- rec$payload$native; ctx <- rec$payload$context
    expect_true(native$success)
    ans <- .hsr("hl_fo_score_retry")(rec)
    expect_identical(ans$route, "preserved_native")
    expect_identical(ans$value, native)
    expect_identical(ctx$active, active)
    expect_identical(ctx$rho_template[setdiff(1:2,active)],
                     log(unname(rec$payload$request$fixed_natural)))
    # Exercise the actual stable engine even though the public policy bypasses
    # it for these native successes. This verifies fixed-coordinate transport.
    retry <- .hsr(".hl_score_run_retry")(ctx, list(list(state = native$null_fit)),
                                        keep_diagnostics = FALSE)
    expect_true(retry$success, info = retry$status)
    expect_identical(retry$value$null_fit$active, active)
    expect_equal(retry$value$null_fit$rho_full[setdiff(1:2,active)],
                 ctx$rho_template[setdiff(1:2,active)], tolerance = 0)
    expect_equal(retry$value$p, native$p, tolerance = 2e-7)
    expect_true(all(vapply(retry$original_refinements, function(x) x$success, logical(1))))
  }
})

test_that("stable arithmetic is isolated and agrees with high precision references", {
  native <- .hsr("hl_fo_gamma_normalizer_derivatives")
  opts <- options(); before <- body(native); env <- environment(native)
  b <- .hsr(".hl_score_backend")()
  # Independent 256-bit values at the switch and at a cancellation-prone shape.
  expect_equal(unlist(b$hl_fo_gamma_normalizer_derivatives(32,log(32))),
    c(value=.81132533627593884, first=.015706372264739307,
      second=-.00049336652030209011, third=.000030994260214090742), tolerance=1e-12)
  for (a in c(.01,1,31.9)) expect_identical(b$hl_fo_gamma_normalizer_derivatives(a,log(a)), native(a,log(a)))
  q <- b$hl_fo_gamma_normalizer_derivatives(1e8,log(1e8))
  expect_equal(q$first/5e-9, 1 + 1/(6e8), tolerance=2e-15)
  expect_equal(q$second/(-5e-17), 1 + 1/(3e8), tolerance=2e-15)
  expect_identical(body(.hsr("hl_fo_gamma_normalizer_derivatives")), before)
  expect_identical(environment(.hsr("hl_fo_gamma_normalizer_derivatives")), env)
  expect_identical(options(), opts)
})

test_that("record consistency binds data, null, configuration, result and implementation", {
  rec <- .hsr_record(); retry <- .hsr("hl_fo_score_retry")
  mutations <- list(
    function(z) {z$payload$request$data$y[1] <- z$payload$request$data$y[1]+1;z},
    function(z) {z$payload$request$data$X[1,2] <- 1;z},
    function(z) {z$payload$request$data$log_offset[1] <- .2;z},
    function(z) {z$payload$request$data$subject[1] <- 2L;z},
    function(z) {z$payload$request$R[1,2] <- 2;z},
    function(z) {z$payload$request$rhs <- .1;z},
    function(z) {z$payload$context$active <- 1L;z},
    function(z) {z$payload$context$rho_template[1] <- .2;z},
    function(z) {z$payload$context$upper[1] <- 3;z},
    function(z) {z$payload$request$include_laplace_constant <- FALSE;z},
    function(z) {z$payload$native$p <- 0;z})
  for (change in mutations) expect_identical(retry(change(rec))$status, "score_recovery_identity_mismatch")
  path <- tempfile(); on.exit(unlink(path)); saveRDS(rec,path)
  expect_identical(retry(readRDS(path))$value,rec$payload$native)
  rec$payload$implementation <- "old implementation"
  rec$identity <- .hsr(".hl_score_fingerprint")(rec$payload)
  expect_identical(retry(rec)$status,"score_recovery_implementation_mismatch")
  for (bad in list(NULL, 3, structure(list(payload=1),class="hl_score_recovery_record")))
    expect_identical(retry(bad)$status,"score_recovery_identity_mismatch")
})

test_that("mismatched supplied fits and unsupported fixed settings are rejected", {
  rec <- .hsr_record(); args <- rec$payload$request
  fit <- rec$payload$native$constrained_fit
  expect_true(fit$point_success)
  args$constrained_fit <- fit
  expect_silent(.hsr(".hl_score_validate_request")(args))
  bad <- args; bad$constrained_fit$data$y[1] <- 99
  expect_error(.hsr(".hl_score_validate_request")(bad),"does not match.*data")
  bad <- args; bad$constrained_fit$constraint$beta_star[2] <- .1
  expect_error(.hsr(".hl_score_validate_request")(bad),"does not match.*constraint")
  bad <- args; bad$constrained_fit$active <- 1L
  expect_error(.hsr(".hl_score_validate_request")(bad),"active dispersion")
  bad <- args; bad$constrained_fit$guard_bounds$upper[1] <- Inf
  expect_error(.hsr(".hl_score_validate_request")(bad),"actual guards")
  bad <- .hsr_record(1L)$payload$request
  bad$constrained_fit <- .hsr_record(1L)$payload$native$constrained_fit
  bad$constrained_fit$rho_full[2] <- log(8)
  expect_error(.hsr(".hl_score_validate_request")(bad),"fixed dispersions")
  bad <- args; bad$active <- 1L
  expect_error(.hsr(".hl_score_validate_request")(bad),"explicitly fixed")
  bad <- args; bad$data$zero_donors <- 1L
  expect_error(.hsr(".hl_score_validate_request")(bad),"nuisance-limit")
})

test_that("actual native guard expansion and objective constants reach the retry", {
  f <- .hsr_fixture()
  rec <- .hsr("hl_fo_score_record")(f$data,f$R,active=1L,fixed_natural=c(phi=4),
    guard_min=c(.005,.1),guard_max=c(.01,20),include_laplace_constant=FALSE,
    keep_diagnostics=FALSE)
  cf <- rec$payload$native$constrained_fit; ctx <- rec$payload$context
  expect_true(cf$point_success)
  expect_true(length(cf$guard_expansions)>0)
  expect_identical(ctx$lower,cf$guard_bounds$lower)
  expect_identical(ctx$upper,cf$guard_bounds$upper)
  expect_gt(ctx$upper,log(.01))
  expect_false(ctx$include_laplace_constant)
  retry <- .hsr(".hl_score_run_retry")(ctx,list(list(state=rec$payload$native$null_fit)),keep_diagnostics=FALSE)
  expect_true(retry$success,info=retry$status)
  expect_gt(retry$value$null_fit$rho,log(.01))
  expect_equal(retry$value$p,rec$payload$native$p,tolerance=2e-7)
})

test_that("cold failures, wrong fixed states and actual guards cannot be bypassed", {
  calls <- 0L
  backend <- list(hl_fo_score_certify_candidate=function(candidate,...) {
    calls <<- calls+1L; list(success=TRUE,state=candidate$state)
  })
  ctx <- list(active=1L,lower=10,upper=12,rho_template=c(0,log(4)))
  candidate <- list(state=list(rho=11))
  fn <- .hsr("hl_fo_score_refine_candidate")
  fail <- list(evaluate=function(...) list(success=FALSE,status="deliberate_failure"))
  ans <- fn(candidate,fail,ctx,backend)
  expect_false(ans$success); expect_identical(calls,0L)
  expect_identical(ans$reason,"score_recovery_cold_evaluation_failed")
  ev <- list(evaluate=function(rho,cold) {
    expect_true(cold); expect_gt(rho,ctx$lower); expect_lt(rho,ctx$upper)
    list(success=TRUE,rho=rho,rho_full=c(rho,log(4)),active=1L,g0=rho-11.5,G0=matrix(1))
  })
  ans <- fn(candidate,ev,ctx,backend)
  expect_true(ans$success); expect_equal(ans$state$rho,11.5)
  expect_identical(calls,1L)
  failed_trials <- list(evaluate=function(rho,cold) {
    if (rho != 11) return(list(success=FALSE,status="failed_trial"))
    list(success=TRUE,rho=rho,rho_full=c(rho,log(4)),active=1L,g0=-.5,G0=matrix(1))
  })
  expect_identical(fn(candidate,failed_trials,ctx,backend)$reason,"score_recovery_cold_evaluation_failed")
  expect_identical(calls,1L)
  singular <- list(evaluate=function(rho,cold)
    list(success=TRUE,rho=rho,rho_full=c(rho,log(4)),active=1L,g0=-.5,G0=matrix(0)))
  expect_identical(fn(candidate,singular,ctx,backend)$reason,"score_recovery_singular_newton_step")
  expect_identical(calls,1L)
  candidate$state$rho <- 12
  expect_identical(fn(candidate,ev,ctx,backend)$reason,"score_recovery_root_outside_actual_guard")
  fixedctx <- list(active=integer(0),lower=numeric(0),upper=numeric(0),rho_template=c(-1,2))
  wrong <- list(evaluate=function(...) list(success=TRUE,active=integer(0),rho=numeric(0),rho_full=c(-1,3)))
  expect_identical(fn(list(state=list(rho=numeric(0))),wrong,fixedctx,backend)$reason,
                   "score_recovery_state_context_mismatch")
  real <- .hsr_record()$payload; calls <- 0L
  # The initial state succeeds, but certification's independent cold replay fails.
  ev <- list(evaluate=function(...) {calls <<- calls+1L
    if(calls==1L) real$native$null_fit else list(success=FALSE,status="cold_failure")})
  ans <- fn(list(state=real$native$null_fit),ev,real$context,.hsr(".hl_score_backend")())
  expect_false(ans$success)
  expect_identical(ans$reason,"score_cold_replay_failure")
  throws <- .hsr(".hl_score_budget_evaluator")(list(evaluate=function(...)stop("injected"),trace=function()NULL),1L)
  expect_identical(throws$evaluate(0)$status,"score_recovery_evaluation_error")
  expect_identical(throws$evaluate(0)$status,"score_recovery_evaluation_cap")
  expect_length(throws$errors(),1L)
  expect_true(throws$exhausted())
  capped <- .hsr(".hl_score_run_retry")(real$context,list(list(state=real$native$null_fit)),max_evaluations=1L)
  expect_false(capped$success); expect_identical(capped$status,"score_recovery_evaluation_cap")
})

# Requested numerical-recovery controls are newly reported even when the
# native result is preserved. Compare every other field without discarding
# numerical results, statuses or diagnostic objects.
.hsr_public_without_requested_score_policy <- function(x) {
  x$provenance$recovery$controls$hl_score_recovery <- NULL
  x
}

test_that("public opt-in preserves accepted outputs and rejects non-HL use", {
  f <- .hsr_fixture(); d <- f$data
  args <- list(count=matrix(d$y,1),id=d$subject,pred=d$X,offset=exp(d$log_offset),
               cpc=0,mincp=1,verbose=FALSE,additional_tests="score",keep_diagnostics=FALSE,recovery="none")
  native <- do.call(scrope_hl,args)
  optin <- do.call(scrope_hl,c(args,list(hl_score_recovery="stable")))
  expect_identical(optin$provenance$recovery$controls$hl_score_recovery, "stable")
  expect_identical(native$provenance$recovery$controls$hl_score_recovery, "none")
  expect_identical(.hsr_public_without_requested_score_policy(optin),
                   .hsr_public_without_requested_score_policy(native))
  wrapper <- do.call(scrope,c(args,list(model="HL",hl_score_recovery="stable")))
  expect_identical(wrapper,optin)
  for(model in c("LN","PGMM")) expect_error(do.call(scrope,c(args,list(model=model,
    hl_score_recovery="stable"))),"only to model='HL'")
  rec <- .hsr_record()
  native <- do.call(.hsr("hl_fo_score_test"),rec$payload$request)
  optin <- do.call(.hsr("hl_fo_score_test"),c(rec$payload$request,list(recovery="stable")))
  expect_identical(optin,native)
})

test_that("full diagnostic native preservation and contrast bases remain valid", {
  f <- .hsr_fixture(); fixed <- c(sigma2=.3,phi=4)
  rec <- .hsr("hl_fo_score_record")(f$data,f$R,active=integer(0),fixed_natural=fixed)
  expect_true(rec$payload$native$success)
  native <- do.call(.hsr("hl_fo_score_test"),rec$payload$request)
  expect_identical(.hsr("hl_fo_score_retry")(rec)$value,native)
  for(R in list(rbind(f$R,2*f$R),rbind(f$R,c(0,0,1)))) {
    rec <- .hsr("hl_fo_score_record")(f$data,R,active=integer(0),fixed_natural=fixed,
                                      keep_diagnostics=FALSE)
    expect_true(rec$payload$native$success)
    retry <- .hsr(".hl_score_run_retry")(rec$payload$context,
      list(list(state=rec$payload$native$null_fit)),keep_diagnostics=FALSE)
    expect_true(retry$success)
    expect_equal(retry$value$p,rec$payload$native$p,tolerance=1e-12)
    expect_equal(retry$value$df,qr(R)$rank)
  }
})

test_that("permuted and rescaled predictors map to the same tested null", {
  f <- .hsr_fixture(); d <- f$data
  args <- list(count=matrix(d$y,1),id=d$subject,pred=d$X,offset=exp(d$log_offset),
    cpc=0,mincp=1,verbose=FALSE,additional_tests="score",keep_diagnostics=FALSE,recovery="none",
    hl_score_recovery="stable")
  base <- do.call(scrope_hl,args)
  args$pred <- cbind(w=3*d$X[,3],Intercept=2.5,x=d$X[,2]/4)
  mapped <- do.call(scrope_hl,args)
  expect_equal(mapped$summary$p_hl_score_x,base$summary$p_hl_score_x,tolerance=2e-7)
  expect_equal(mapped$summary$p_hl_score_w,base$summary$p_hl_score_w,tolerance=2e-7)
  expect_equal(mapped$summary$logFC_x/4,base$summary$logFC_x,tolerance=2e-7)
  expect_equal(mapped$summary$logFC_w*3,base$summary$logFC_w,tolerance=2e-7)
  args$hl_score_recovery <- "none"
  native <- do.call(scrope_hl,args)
  expect_identical(native$provenance$recovery$controls$hl_score_recovery, "none")
  expect_identical(mapped$provenance$recovery$controls$hl_score_recovery, "stable")
  expect_identical(.hsr_public_without_requested_score_policy(native),
                   .hsr_public_without_requested_score_policy(mapped))
})

test_that("every evaluated failure is recorded and a rejected cold trial is terminal", {
  n <- 0L
  raw <- list(evaluate=function(rho,cold=FALSE){n <<- n+1L
    if(n==1L)return(list(success=FALSE,status="mode_failure",reason="returned"))
    if(n==2L)stop("thrown")
    list(success=TRUE)
  },trace=function()NULL)
  ev <- .hsr(".hl_score_budget_evaluator")(raw,10L)
  expect_false(ev$evaluate(1,cold=FALSE)$success)
  expect_false(ev$cold_failed())
  expect_false(ev$evaluate(2,cold=TRUE)$success)
  expect_true(ev$cold_failed())
  expect_true(ev$evaluate(3,cold=TRUE)$success)
  expect_length(ev$failures(),2L)
  expect_identical(vapply(ev$failures(),function(x)x$kind,character(1)),c("returned_failure","exception"))
  expect_identical(vapply(ev$failures(),function(x)x$cold,logical(1)),c(FALSE,TRUE))
  expect_identical(vapply(ev$failures(),function(x)x$evaluation,integer(1)),1:2)
  expect_length(ev$errors(),1L)

  # This trial would succeed at a shorter step, but the strict policy must stop.
  visits <- numeric(0);certified <- FALSE
  ev <- list(evaluate=function(rho,cold=FALSE){visits <<- c(visits,rho)
    if(rho==1)return(list(success=FALSE,status="mode_failure"))
    list(success=TRUE,rho=rho,active=1L,rho_full=c(rho,log(4)),g0=rho-1,G0=matrix(1))
  })
  b <- list(hl_fo_score_certify_candidate=function(...){certified <<- TRUE;list(success=TRUE)})
  ctx <- list(active=1L,lower=-2,upper=2,rho_template=c(0,log(4)))
  ans <- .hsr("hl_fo_score_refine_candidate")(list(state=list(rho=0)),ev,ctx,b)
  expect_false(ans$success)
  expect_identical(ans$reason,"score_recovery_cold_evaluation_failed")
  expect_identical(visits,c(0,1))
  expect_false(certified)
})

# Inspect the exact inputs to the implementation seal without maintaining a
# second copy of its schema. These diagnostics are also available in disposable
# child processes, where the loaded implementation may differ from the parent.
.hsr_implementation_probe <- function() {
  implementation <- getFromNamespace(".hl_score_implementation", "scRoPE")
  e <- new.env(parent = environment(implementation))
  e$.hl_score_fingerprint <- identity
  environment(implementation) <- e
  parts <- implementation()
  fingerprint <- getFromNamespace(".hl_score_fingerprint", "scRoPE")
  list(implementation = fingerprint(parts),
       components = vapply(parts, fingerprint, character(1)),
       functions = vapply(parts$functions, fingerprint, character(1)),
       numerical_libraries = parts$numerical_libraries,
       execution = parts$execution)
}

.hsr_verify_replay <- function(rec, parent) {
  retry <- getFromNamespace("hl_fo_score_retry", "scRoPE")(rec)
  if (identical(retry$route, "preserved_native")) return(invisible(retry))
  fingerprint <- getFromNamespace(".hl_score_fingerprint", "scRoPE")
  child <- tryCatch(.hsr_implementation_probe(), error = conditionMessage)
  cat("HL_REPLAY_DIAGNOSTICS\n")
  print(list(status = retry$status, route = retry$route,
    recorded_payload = rec$identity, replayed_payload = fingerprint(rec$payload),
    recorded_implementation = rec$payload$implementation,
    parent_implementation = parent$implementation,
    child_implementation = if (is.list(child)) child$implementation else child))
  if (is.list(child)) {
    print(cbind(parent = parent$components, child = child$components))
    changed <- union(names(parent$functions), names(child$functions))
    changed <- changed[!vapply(changed, function(n)
      identical(parent$functions[n], child$functions[n]), logical(1))]
    print(list(changed_functions = cbind(parent = parent$functions[changed],
                                         child = child$functions[changed]),
      parent_numerical_libraries = parent$numerical_libraries,
      child_numerical_libraries = child$numerical_libraries,
      recorded_execution = rec$payload$execution,
      parent_execution = parent$execution, child_execution = child$execution))
  }
  stop("Saved HL record refused before native-library mutation: ", retry$status)
}

test_that("implementation ignores source metadata but binds executable expressions", {
  implementation <- .hsr(".hl_score_implementation")
  ns <- new.env(parent = environment(implementation))
  e <- new.env(parent = environment(implementation))
  environment(implementation) <- e
  lines <- c("function(x = 1) {", "  x + 1", "}")
  src <- srcfilecopy("identity-fixture.R", lines)
  fn <- eval(parse(text = lines, srcfile = src, keep.source = TRUE))
  environment(fn) <- ns
  ns$hl_fo_score_test <- e$hl_fo_score_test <- fn
  # Reproduce installed keep-source laziness without depending on how this
  # particular test installation was built. Forcing lines changes serialization
  # of source metadata, but neither the executable body nor its defaults.
  delayedAssign("lines", c("function(x = 1) {", "  x + 1", "}"),
                eval.env = baseenv(), assign.env = src)
  source_before <- serialize(body(fn), NULL, version = 2)
  before <- implementation()
  expect_identical(serialize(body(fn), NULL, version = 2), source_before)
  expect_identical(getSrcLines(src, 1L, 3L), lines)
  expect_false(identical(serialize(body(fn), NULL, version = 2), source_before))
  expect_identical(implementation(), before)
  src$filename <- "another source display path.R"
  expect_identical(implementation(), before)

  changed <- fn
  body(changed) <- quote({ x + 2 })
  ns$hl_fo_score_test <- changed
  expect_false(identical(implementation(), before))
  changed <- fn
  formals(changed)$x <- 2
  ns$hl_fo_score_test <- changed
  expect_false(identical(implementation(), before))
  changed <- fn
  attr(body(changed), "semantic_marker") <- "must remain bound"
  ns$hl_fo_score_test <- changed
  expect_false(identical(implementation(), before))
  ns$hl_fo_score_test <- fn
  expect_identical(implementation(), before)
})

test_that("identity covers the executed native backend and native file contents", {
  exec <- .hsr(".hl_score_execution_identity")()
  expect_identical(exec$inner,"compiled")
  dll <- getLoadedDLLs()[["scRoPE"]]
  expect_identical(exec$native_libraries$scRoPE$md5,unname(tools::md5sum(dll[["path"]])))
  expect_identical(exec$native_libraries$nloptr$md5,
    unname(tools::md5sum(getLoadedDLLs()[["nloptr"]][["path"]])))
  p <- tempfile();on.exit(unlink(p),add=TRUE)
  writeBin(as.raw(1:10),p)
  id1 <- .hsr(".hl_score_native_library")(list(path=p,dynamicLookup=FALSE))
  writeBin(as.raw(2:11),p)
  id2 <- .hsr(".hl_score_native_library")(list(path=p,dynamicLookup=FALSE))
  expect_identical(id1$path,id2$path);expect_false(identical(id1$md5,id2$md5))
  implementation <- .hsr(".hl_score_implementation")
  e <- new.env(parent=environment(implementation));environment(implementation)<-e
  e$.hl_score_execution_identity <- function() exec
  before <- implementation()
  exec$native_libraries$scRoPE$md5 <- "different native build under identical R code"
  expect_false(identical(before,implementation()))

  # Source inspection is a real trigger in retained-source installations. The
  # full retention tests also exercise trace()/untrace() before these replays.
  before <- .hsr(".hl_score_implementation")()
  src <- attr(body(.hsr("hl_fo_fit")), "srcfile")
  if (!is.null(src)) invisible(getSrcLines(src, 1L, 1L))
  expect_identical(.hsr(".hl_score_implementation")(), before)

  # Unload only in a disposable R process: it must invalidate a saved record,
  # while a fresh all-fixed R-backend record is valid and preserves its native result.
  rec <- .hsr_record(integer(0));rds<-tempfile(fileext=".rds");script<-tempfile(fileext=".R")
  on.exit(unlink(c(rds,script)),add=TRUE)
  saveRDS(list(record=rec,parent=.hsr_implementation_probe()),rds)
  lib <- dirname(getNamespaceInfo(asNamespace("scRoPE"),"path"))
  quote_r <- function(x)paste(deparse(x),collapse="")
  child_setup <- c(sprintf('library(scRoPE,lib.loc=%s)',quote_r(lib)),
    paste('.hsr_implementation_probe <-',paste(deparse(.hsr_implementation_probe),collapse="\n")),
    paste('.hsr_verify_replay <-',paste(deparse(.hsr_verify_replay),collapse="\n")),
    sprintf('bundle<-readRDS(%s);rec<-bundle$record',quote_r(rds)),
    '.hsr_verify_replay(rec,bundle$parent)')
  writeLines(c(child_setup,
    'before<-scRoPE:::.hl_score_execution_identity()',
    'stopifnot(before$inner=="compiled",scRoPE:::hl_fo_score_retry(rec)$route=="preserved_native")',
    'dll_path<-getLoadedDLLs()[["scRoPE"]][["path"]];dyn.unload(dll_path)',
    'after<-scRoPE:::.hl_score_execution_identity()',
    'stopifnot(after$inner=="R",is.null(after$native_libraries$scRoPE))',
    'stopifnot(scRoPE:::hl_fo_score_retry(rec)$status=="score_recovery_implementation_mismatch")',
    'fresh<-do.call(scRoPE:::hl_fo_score_record,rec$payload$request)',
    'stopifnot(fresh$payload$execution$inner=="R",fresh$payload$native$success)',
    'stopifnot(identical(scRoPE:::hl_fo_score_retry(fresh)$value,fresh$payload$native))',
    'dyn.load(dll_path)',
    'stopifnot(scRoPE:::hl_fo_score_retry(rec)$status=="score_recovery_implementation_unavailable")',
    'stopifnot(scRoPE:::hl_fo_score_retry(fresh)$status=="score_recovery_implementation_unavailable")',
    'cat("BACKEND_SWITCH_VERIFIED\\n")'),script)
  out <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(script)),stdout=TRUE,stderr=TRUE)
  expect_null(attr(out,"status"),info=paste(out,collapse="\n"))
  expect_true(any(grepl("BACKEND_SWITCH_VERIFIED",out)),info=paste(out,collapse="\n"))
  writeLines(c(child_setup,
    'stopifnot(scRoPE:::hl_fo_score_retry(rec)$route=="preserved_native")',
    'path<-getLoadedDLLs()[["nloptr"]][["path"]];dyn.unload(path);dyn.load(path)',
    'stopifnot(scRoPE:::hl_fo_score_retry(rec)$status=="score_recovery_implementation_unavailable")',
    'cat("NLOPT_STALE_BINDING_VERIFIED\\n")'),script)
  out <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(script)),stdout=TRUE,stderr=TRUE)
  expect_null(attr(out,"status"),info=paste(out,collapse="\n"))
  expect_true(any(grepl("NLOPT_STALE_BINDING_VERIFIED",out)),info=paste(out,collapse="\n"))
})

.hsr_portable <- local({
  cache <- new.env(parent=emptyenv())
  function(which=c("recovery","multiplicity")) {
    which<-match.arg(which)
    if(exists(which,cache,inherits=FALSE))return(cache[[which]])
    input<-dget(testthat::test_path("fixtures",paste0("hl-score-",which,".R")))
    args<-c(input,list(model="HL",additional_tests="score",cpc=0,mincp=1,verbose=FALSE,recovery="none"))
    native<-do.call(scrope,c(args,list(keep_diagnostics=TRUE)))
    stable<-do.call(scrope,c(args,list(keep_diagnostics=TRUE,hl_score_recovery="stable")))
    compact<-do.call(scrope,c(args,list(keep_diagnostics=FALSE,hl_score_recovery="stable")))
    score<-native$diagnostics$hl_adj[[1]]$x$score
    rec<-.hsr("hl_fo_score_record")(score$root_candidates[[1]]$state$data,
      score$contrast$R,rhs=score$contrast$rhs,
      unrestricted_fit=native$diagnostics$hl[[1]],constrained_fit=score$constrained_fit,
      keep_diagnostics=FALSE)
    retry<-.hsr("hl_fo_score_retry")(rec)
    cache[[which]]<-list(input=input,args=args,native=native,stable=stable,compact=compact,
                        record=rec,retry=retry)
    cache[[which]]
  }
})

test_that("portable public ambiguity recovers while genuine multiplicity stays withheld", {
  for(which in c("recovery","multiplicity")) {
    f<-.hsr_portable(which); n<-f$native$diagnostics$hl_adj[[1]]$x$score
    s<-f$stable$diagnostics$hl_adj[[1]]$x$score
    expect_true(f$native$summary$fit_success)
    expect_identical(n$status,"score_null_multiple_roots")
    expect_identical(f$record$payload$native$status,"score_null_multiple_roots")
    expect_gte(length(n$root_candidates),2L)
    expect_identical(s$numerical_recovery$primary,n)
    expect_true(all(vapply(s$numerical_recovery$attempt$original_refinements,function(x)x$success,logical(1))))
    expect_length(s$numerical_recovery$attempt$original_refinements,length(n$root_candidates))
    expect_false(any(vapply(s$numerical_recovery$attempt$evaluation_failures,function(x)x$cold,logical(1))))
    cols<-setdiff(names(f$native$summary),grep("hl_score|test_status|test_failure_detail",names(f$native$summary),value=TRUE))
    expect_identical(f$stable$summary[,cols],f$native$summary[,cols])
    expect_identical(f$stable$overdispersion,f$native$overdispersion)
    expect_identical(f$stable$convergence,f$native$convergence)
    expect_identical(f$stable$diagnostics$hl,f$native$diagnostics$hl)
    expect_identical(f$compact$summary,f$stable$summary)
    expect_identical(f$compact$hl_score_recovery,f$stable$hl_score_recovery)
    expect_null(f$compact$diagnostics)
    route<-subset(f$compact$hl_score_recovery,contrast=="x")
    expect_equal(nrow(route),1L);expect_equal(route$cold_failures,0L)
    expect_identical(route$backend,"compiled")
    if(which=="recovery") {
      expect_true(s$success);expect_true(f$retry$success)
      expect_identical(route$route,"stable_gamma_cold_score_refinement")
      expect_equal(s$p,.7862556,tolerance=2e-6)
      expect_equal(f$retry$value$p,s$p,tolerance=1e-8)
      expect_length(s$root_candidates,1L)
    }else{
      expect_false(s$success);expect_false(f$retry$success)
      expect_identical(route$route,"withheld")
      expect_identical(route$retry_status,"score_null_multiple_roots")
      expect_identical(f$retry$value,f$record$payload$native)
      expect_identical(s[setdiff(names(s),"numerical_recovery")],n)
      expect_length(s$root_candidates,2L)
    }
    args<-f$args;args$model<-NULL
    via_hl<-do.call(scrope_hl,c(args,list(keep_diagnostics=FALSE,hl_score_recovery="stable")))
    expect_identical(via_hl,f$compact)
  }
})

test_that("a returned cold failure vetoes the entire real retry and warm failures stay auditable", {
  rec<-.hsr_portable("recovery")$record
  run<-.hsr(".hl_score_run_retry");env<-new.env(parent=environment(run));environment(run)<-env
  factory<-.hsr(".hl_score_backend")
  inject<-function(cold_failure){
    env$.hl_score_backend<-function(){
      b<-factory();make<-b$hl_fo_score_make_evaluator
      b$hl_fo_score_make_evaluator<-function(...){
        raw<-make(...);evaluate<-raw$evaluate;n<-0L;injected<-FALSE
        raw$evaluate<-function(rho,cold=FALSE){
          n<<-n+1L
          if(!injected && ((cold_failure && cold && n==2L) || (!cold_failure && !cold))){
            injected<<-TRUE
            return(list(success=FALSE,status="mode_failure",reason="injected_returned_failure"))
          }
          evaluate(rho,cold=cold)
        }
        raw
      }
      b
    }
    run(rec$payload$context,rec$payload$native$root_candidates,keep_diagnostics=FALSE)
  }
  cold<-inject(TRUE)
  expect_false(cold$success);expect_identical(cold$status,"score_recovery_cold_evaluation_failed")
  expect_length(cold$evaluation_failures,1L);expect_length(cold$errors,0L)
  expect_identical(cold$evaluation_failures[[1]]$status,"mode_failure")
  expect_identical(cold$evaluation_failures[[1]]$kind,"returned_failure")
  expect_true(cold$evaluation_failures[[1]]$cold)
  warm<-inject(FALSE)
  expect_true(warm$success,info=warm$status)
  expect_length(warm$evaluation_failures,1L);expect_length(warm$errors,0L)
  expect_false(warm$evaluation_failures[[1]]$cold)
  expect_equal(warm$value$p,.7862556,tolerance=2e-6)
})

test_that("implementation identity uses locale-independent ordering", {
  old<-Sys.getlocale("LC_COLLATE");on.exit(Sys.setlocale("LC_COLLATE",old),add=TRUE)
  before<-.hsr(".hl_score_implementation")()
  for(locale in unique(c("C", "C.UTF-8", "en_US.UTF-8"))) {
    selected<-suppressWarnings(Sys.setlocale("LC_COLLATE",locale))
    if(nzchar(selected)) expect_identical(.hsr(".hl_score_implementation")(),before)
  }
})
