# Optional CI observers; restore every binding after the test process.
# Original decisions/returned values are preserved. No retry is performed.
.hl_history_diagnostic_install <- function(directory) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  directory <- normalizePath(directory,mustWork=TRUE)
  ns <- asNamespace('scRoPE')
  original <- lapply(c('working_hl_history_audit','working_hl_profile','working_hl_refine'),
    function(n)get(n,ns,inherits=FALSE))
  names(original) <- c('working_hl_history_audit','working_hl_profile','working_hl_refine')
  counter <- 0L
  context <- NULL
  metadata <- list(R=R.version,session=sessionInfo(),numerical_libraries=extSoftVersion(),
    machine=.Machine,matprod=getOption('matprod'),
    environment=Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
      'VECLIB_MAXIMUM_THREADS','R_KEEP_PKG_SOURCE','R_ENABLE_JIT')),
    backend=get('.hl_score_execution_identity',ns)())
  profile_observer <- function(...) {
    withCallingHandlers(original$working_hl_profile(...),error=function(e) {
      if(is.null(context))return(invisible(NULL))
      # A calling handler sees the original profile's local failed result before
      # the audit catches its error. Diagnostics themselves cannot replace it.
      tryCatch({
        frames <- sys.frames()
        indices <- which(vapply(seq_along(frames),function(i)
          identical(sys.function(i),original$working_hl_profile),logical(1)))
        frame <- if(length(indices))frames[[tail(indices,1L)]] else NULL
        values <- if(is.null(frame))NULL else mget(
          intersect(c('rho','active','source','ans'),ls(frame,all.names=TRUE)),
          envir=frame,inherits=FALSE)
        context$profile_errors[[length(context$profile_errors)+1L]] <-
          list(error=conditionMessage(e),original_locals=values)
      },error=function(ignored)NULL)
      invisible(NULL)
    })
  }
  refine_observer <- function(...) {
    ans <- original$working_hl_refine(...)
    if(!is.null(context)&&!isTRUE(ans$success))
      tryCatch({
        context$refinement_failures[[length(context$refinement_failures)+1L]] <- ans
      },error=function(ignored)NULL)
    if(!is.null(context)&&identical(ans$reason,'coalesced_stalled_history_root'))
      tryCatch({
        context$refinement_resolutions[[length(context$refinement_resolutions)+1L]] <-
          ans[c('success','reason','trace','stationarity','coalescence')]
      },error=function(ignored)NULL)
    ans
  }
  audit_observer <- function(history,fit,...) {
    previous <- context
    context <<- new.env(parent=emptyenv())
    context$profile_errors <- context$refinement_failures <- context$refinement_resolutions <- list()
    on.exit(context <<- previous,add=TRUE)
    ans <- original$working_hl_history_audit(history,fit,...)
    # Failed serialization must not introduce a warning/error into the tests.
    withCallingHandlers(tryCatch({
      data <- fit$data
      focal_fixture <- identical(data$n,480L)&&identical(data$m,60L)&&
        identical(as.integer(data$zero_donors),2L)
      if(!isTRUE(ans$success)||focal_fixture||length(context$refinement_resolutions)) {
        counter <<- counter+1L
        path <- file.path(directory,sprintf('hl-history-%s-%03d.rds',Sys.getpid(),counter))
        record <- list(schema=1L,environment=metadata,history=history,selected_fit=fit,audit=ans,
          original_profile_errors=context$profile_errors,
          original_refinement_failures=context$refinement_failures,
          original_refinement_resolutions=context$refinement_resolutions,
          keep_details=ans$details_retained,calls=lapply(tail(sys.calls(),12L),deparse))
        saveRDS(record,path)
        cat('HL_HISTORY_DIAGNOSTIC ',basename(path),' success=',isTRUE(ans$success),
          ' reason=',ans$reason,' processed=',ans$processed,' retained=',ans$retained,
          ' resolved_stalls=',length(context$refinement_resolutions),'\n',sep='')
        if(!isTRUE(ans$success)) {
          w <- ans$failure_witness$observation
          cat('HL_HISTORY_FAILURE_WITNESS\n');dput(ans$failure_witness)
          cat('HL_HISTORY_FAILURE_CONFIGURATION\n')
          key <- if(length(w$configuration)==1L)w$configuration else
            ans$failure_witness$configuration
          dput(if(length(key)==1L)ans$configurations[[key]] else NULL)
          cat('HL_HISTORY_ORIGINAL_PROFILE_ERRORS\n');dput(context$profile_errors)
          cat('HL_HISTORY_ORIGINAL_REFINEMENT_FAILURES\n');dput(context$refinement_failures)
        }
      }
    },error=function(e)cat('HL_HISTORY_DIAGNOSTIC_WRITE_ERROR ',conditionMessage(e),'\n',sep='')),
      warning=function(w)invokeRestart('muffleWarning'))
    ans
  }
  restore <- function()for(n in names(original))assignInNamespace(n,original[[n]],ns='scRoPE')
  installed <- FALSE
  on.exit(if(!installed)restore(),add=TRUE)
  assignInNamespace('working_hl_profile',profile_observer,ns='scRoPE')
  assignInNamespace('working_hl_refine',refine_observer,ns='scRoPE')
  assignInNamespace('working_hl_history_audit',audit_observer,ns='scRoPE')
  installed <- TRUE
  restore
}
