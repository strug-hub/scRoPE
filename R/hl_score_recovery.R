# Opt-in numerical recovery of an ambiguous finite full-chi Score null.
# Records are consistency-sealed, not cryptographic authentication objects.
.hl_score_fingerprint <- function(x) {
  path <- tempfile("scrope-hl-score-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

.hl_score_native_library <- function(dll) {
  path <- normalizePath(dll[["path"]], mustWork = TRUE)
  digest <- unname(tools::md5sum(path))
  if (length(digest) != 1L || is.na(digest)) stop("Cannot fingerprint loaded native library: ", path)
  list(path = path, md5 = digest, dynamic_lookup = dll[["dynamicLookup"]])
}

.hl_score_execution_identity <- function() {
  symbol <- "_scRoPE_hl_fo_inner_eval_cpp"
  compiled <- isTRUE(tryCatch(is.loaded(symbol, PACKAGE = "scRoPE"), error = function(e) FALSE))
  loaded <- getLoadedDLLs()
  native <- list()
  # nloptr supplies the recession check even when inner evaluation uses R.
  for (name in c(if (compiled) "scRoPE", "nloptr")) {
    dll <- loaded[[name]]
    if (is.null(dll)) stop("Required native library is not loaded: ", name)
    native[[name]] <- .hl_score_native_library(dll)
  }
  nlopt_binding <- get("NLoptR_Optimize", asNamespace("nloptr"), inherits = FALSE)
  nlopt_bound <- tryCatch(getNativeSymbolInfo("NLoptR_Optimize", PACKAGE = nlopt_binding$dll),
                          error = function(e) NULL)
  if (is.null(nlopt_bound) ||
      !identical(nlopt_bound$dll[["path"]], loaded[["nloptr"]][["path"]]))
    stop("The registered nloptr binding is stale or invalid.")
  if (compiled) {
    binding <- get(symbol, environment(hl_fo_score_test), inherits = FALSE)
    resolved <- getNativeSymbolInfo(symbol, PACKAGE = "scRoPE")
    # A name lookup can succeed after dyn.unload()/dyn.load() even though the
    # namespace's original registered binding has a dead DLLInfo pointer.
    bound <- tryCatch(getNativeSymbolInfo(symbol, PACKAGE = binding$dll),
                      error = function(e) NULL)
    if (is.null(bound)) stop("The registered HL binding is stale or invalid.")
    if (!identical(binding$dll[["path"]], resolved$dll[["path"]]) ||
        !identical(binding$name, resolved$name) ||
        !identical(binding$numParameters, resolved$numParameters))
      stop("The compiled HL binding differs from the loaded dispatch symbol.")
  }
  list(inner = if (compiled) "compiled" else "R", symbol = if (compiled) symbol else NULL,
       native_libraries = native)
}

.hl_score_implementation <- function() {
  ns <- environment(hl_fo_score_test)
  names <- sort(grep("^\\.?(hl_fo_|hl_score_)", ls(ns, all.names = TRUE, sorted = FALSE),
                     value = TRUE), method = "radix")
  functions <- lapply(names, function(n) {
    f <- get(n, ns, inherits = FALSE)
    if (is.function(f)) list(formals = formals(f), body = body(f)) else NULL
  })
  .hl_score_fingerprint(list(schema = 2L, R = R.version$version.string,
                            execution = .hl_score_execution_identity(),
                            platform = R.version$platform,
                            packages = lapply(c("scRoPE", "nloptr", "tools"), getNamespaceVersion),
                            numerical_libraries = extSoftVersion(),
                            functions = setNames(functions, names)))
}

.hl_score_validate_request <- function(args) {
  d <- args$data
  rebuilt <- hl_fo_make_data(d$y, d$X, d$log_offset, d$subject, d$intcol)
  if (!identical(d[names(rebuilt)], rebuilt))
    stop("Score recovery requires an intact canonical data object.")
  if ((!is.null(d$zero_donors) && d$zero_donors != 0) ||
      (!is.null(d$inference_m) && d$inference_m != d$m))
    stop("Score recovery does not support nuisance-limit donor augmentation.")
  active <- .hl_fo_score_active(args$active)
  fixed <- args$fixed_natural
  idx <- match(names(fixed), c("sigma2", "phi"))
  if (length(fixed) && (is.null(names(fixed)) || anyNA(idx) ||
      anyDuplicated(idx) || any(!is.finite(fixed)) || any(fixed <= 0)))
    stop("Recovery fixed_natural must have unique positive sigma2/phi entries.")
  if (!setequal(idx, setdiff(1:2, active)))
    stop("Recovery requires every inactive dispersion to be explicitly fixed.")
  for (g in c("guard_min", "guard_max")) {
    if (length(args[[g]]) != 2L || any(!is.finite(args[[g]])) || any(args[[g]] <= 0))
      stop("Recovery guards must contain two finite positive values.")
  }
  if (any(args$guard_min >= args$guard_max)) stop("Recovery guard_min must be below guard_max.")
  if (!is.logical(args$include_laplace_constant) ||
      length(args$include_laplace_constant) != 1L || is.na(args$include_laplace_constant))
    stop("include_laplace_constant must be TRUE or FALSE.")
  for (name in c("unrestricted_fit", "constrained_fit")) {
    fit <- args[[name]]
    if (is.null(fit)) next
    # A supplied fit is used only as a start, but must still belong to this data.
    if (is.null(fit$data) || !identical(fit$data[names(rebuilt)], rebuilt))
      stop("Supplied ", name, " does not match the score data.")
    if (name == "constrained_fit" && isTRUE(fit$point_success)) {
      affine <- hl_fo_affine_constraint(args$R, args$rhs, p = d$p)
      if (!isTRUE(affine$success) || !isTRUE(fit$constraint$success) ||
          !isTRUE(all.equal(fit$constraint$beta_star, affine$beta_star, tolerance = 1e-12)) ||
          !isTRUE(all.equal(tcrossprod(fit$constraint$K), tcrossprod(affine$K), tolerance = 1e-12)))
        stop("Supplied constrained_fit does not match the score constraint.")
      bounds <- fit$guard_bounds
      if (length(active) && (length(bounds$lower) != length(active) ||
          length(bounds$upper) != length(active) ||
          any(!is.finite(c(bounds$lower, bounds$upper))) || any(bounds$lower >= bounds$upper)))
        stop("Supplied constrained_fit has invalid actual guards.")
      if (!identical(as.integer(fit$active), active))
        stop("Supplied constrained_fit has a different active dispersion mask.")
      if (length(idx) && !isTRUE(all.equal(as.numeric(fit$rho_full[idx]), log(as.numeric(fixed)), tolerance = 1e-12)))
        stop("Supplied constrained_fit has different fixed dispersions.")
    }
  }
  invisible(TRUE)
}

hl_fo_score_record <- function(data, R, rhs = NULL, unrestricted_fit = NULL,
                               guard_min = c(1e-4, 1e-4), guard_max = c(10, 1000),
                               active = c("tau", "xi"), fixed_natural = numeric(0),
                               constrained_fit = NULL, include_laplace_constant = TRUE,
                               keep_diagnostics = TRUE) {
  args <- as.list(environment())
  .hl_score_validate_request(args)
  implementation <- .hl_score_implementation()
  execution <- .hl_score_execution_identity()
  context <- NULL
  capture <- function(x) context <<- x
  native <- do.call(hl_fo_score_test, c(args, list(recovery = "none", .context_callback = capture)))
  if (!identical(implementation, .hl_score_implementation()))
    stop("The HL implementation changed during native fitting.")
  payload <- list(request = args, context = context, native = native,
                  implementation = implementation, execution = execution)
  structure(list(payload = payload, identity = .hl_score_fingerprint(payload)),
            class = "hl_score_recovery_record")
}

.hl_score_gamma_stable <- function(alpha, log_alpha) {
  if (alpha < 32) return(hl_fo_gamma_normalizer_derivatives(alpha, log_alpha))
  B <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
  n <- 2 * seq_along(B)
  inv <- 1 / alpha
  list(value = .5 * (log_alpha - log(2*pi)) - sum(B/(n*(n-1))*inv^(n-1)),
       first = inv/2 + sum(B/n*inv^n),
       second = -inv^2/2 - sum(B*inv^(n+1)),
       third = inv^3 + sum((n+1)*B*inv^(n+2)))
}

.hl_score_backend <- function() {
  # Lexically isolated copies of the installed first-order functions. No source
  # checkout, global option, shared binding or namespace mutation is involved.
  ns <- environment(hl_fo_score_test)
  backend <- new.env(parent = ns)
  names <- grep("^\\.?hl_fo_", ls(ns, all.names = TRUE), value = TRUE)
  for (n in names) {
    f <- get(n, ns, inherits = FALSE)
    if (is.function(f)) {
      environment(f) <- backend
      assign(n, f, backend)
    }
  }
  # Keep this closure's namespace environment, so its small-shape branch calls
  # the original normalizer rather than recursively calling itself.
  backend$hl_fo_gamma_normalizer_derivatives <- .hl_score_gamma_stable
  backend
}

.hl_score_budget_evaluator <- function(evaluator, limit) {
  count <- 0L
  exhausted <- FALSE
  errors <- list()
  failures <- list()
  evaluate <- function(rho, cold = FALSE) {
    if (count >= limit) {
      exhausted <<- TRUE
      return(list(success = FALSE, status = "score_recovery_evaluation_cap"))
    }
    count <<- count + 1L
    kind <- "returned_failure"
    answer <- tryCatch(evaluator$evaluate(rho, cold = cold), error = function(e) {
      kind <<- "exception"
      errors[[length(errors) + 1L]] <<- list(evaluation = count, rho = rho,
        cold = isTRUE(cold), message = conditionMessage(e))
      list(success = FALSE, status = "score_recovery_evaluation_error", reason = conditionMessage(e))
    })
    if (!is.list(answer)) answer <- list(success = FALSE,
      status = "score_recovery_invalid_evaluation", reason = "evaluator_returned_non_list")
    if (!isTRUE(answer$success)) failures[[length(failures) + 1L]] <<- list(
      evaluation = count, rho = as.numeric(rho), cold = isTRUE(cold), kind = kind,
      status = answer$status, reason = answer$reason)
    answer
  }
  list(evaluate = evaluate, count = function() count,
       exhausted = function() exhausted, errors = function() errors,
       failures = function() failures,
       cold_failed = function() any(vapply(failures, function(x) x$cold, logical(1))),
       trace = evaluator$trace)
}

hl_fo_score_refine_candidate <- function(candidate, evaluator, context, backend,
                                          max_iterations = 12L, keep_diagnostics = TRUE) {
  fail <- function(reason, state = NULL, trace = list())
    list(success = FALSE, status = "score_null_no_convergence", reason = reason,
         state = state, refinement_trace = trace)
  q <- length(context$active)
  rho <- candidate$state$rho
  valid_rho <- function(r) is.numeric(r) && length(r) == q && all(is.finite(r)) &&
    all(r > context$lower) && all(r < context$upper)
  if (!valid_rho(rho)) return(fail("score_recovery_root_outside_actual_guard"))
  current <- evaluator$evaluate(rho, cold = TRUE)
  trace <- list()
  for (iteration in 0:max_iterations) {
    if (!isTRUE(current$success)) return(fail("score_recovery_cold_evaluation_failed", current, trace))
    inactive <- setdiff(1:2, context$active)
    if (!valid_rho(current$rho) || length(current$rho_full) != 2L ||
        any(!is.finite(current$rho_full)) ||
        !identical(as.numeric(current$rho_full[context$active]), as.numeric(current$rho)) ||
        !identical(as.integer(current$active), context$active) ||
        !identical(as.numeric(current$rho_full[inactive]), as.numeric(context$rho_template[inactive])))
      return(fail("score_recovery_state_context_mismatch", current, trace))
    if (!q) break
    step <- tryCatch(solve(current$G0, current$g0), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) return(fail("score_recovery_singular_newton_step", current, trace))
    if (isTRUE(keep_diagnostics)) trace[[length(trace) + 1L]] <- list(iteration = iteration, rho = current$rho,
                                     score = current$g0, step = step)
    if (max(abs(step)) <= 1e-9) break
    if (iteration == max_iterations) return(fail("score_recovery_refinement_cap", current, trace))
    accepted <- FALSE
    for (rate in 2^-(0:10)) {
      proposed <- current$rho - rate * step
      if (!valid_rho(proposed)) next
      trial <- evaluator$evaluate(proposed, cold = TRUE)
      if (!isTRUE(trial$success))
        return(fail("score_recovery_cold_evaluation_failed", trial, trace))
      next_step <- tryCatch(solve(trial$G0, trial$g0), error = function(e) NULL)
      if (!is.null(next_step) && all(is.finite(next_step)) && max(abs(next_step)) < max(abs(step))) {
        current <- trial
        accepted <- TRUE
        break
      }
    }
    if (!accepted) return(fail("score_recovery_no_cold_improvement", current, trace))
  }
  candidate$state <- current
  certified <- backend$hl_fo_score_certify_candidate(
    candidate, evaluator, context$R, context$rhs, context$C_R, context$K,
    keep_diagnostics = keep_diagnostics)
  certified$refinement_trace <- trace
  certified
}

.hl_score_run_retry <- function(context, originals, max_evaluations = 3000L,
                                keep_diagnostics = TRUE) {
  backend <- .hl_score_backend()
  raw <- backend$hl_fo_score_make_evaluator(
    data = context$data, rho_template = context$rho_template, active = context$active,
    beta_star = context$beta_star, K = context$K, seed_beta = context$seed_beta,
    include_laplace_constant = context$include_laplace_constant,
    constraint_R = context$R, constraint_rhs = context$rhs, keep_diagnostics = keep_diagnostics)
  ev <- .hl_score_budget_evaluator(raw, max_evaluations)
  retained <- lapply(originals, hl_fo_score_refine_candidate,
                     evaluator = ev, context = context, backend = backend, keep_diagnostics = keep_diagnostics)
  fail <- function(reason, runs = list(), attempts = list()) list(
    success = FALSE, status = reason, original_refinements = retained,
    runs = runs, certification_attempts = attempts, evaluations = ev$count(),
    errors = ev$errors(), evaluation_failures = ev$failures(), trace = ev$trace())
  # An old certified root may not disappear merely because the retry failed to
  # refine it. Ambiguity is resolved only when every original root is accounted for.
  if (ev$exhausted()) return(fail("score_recovery_evaluation_cap"))
  if (ev$cold_failed()) return(fail("score_recovery_cold_evaluation_failed"))
  if (length(ev$errors())) return(fail("score_recovery_evaluation_error"))
  if (any(!vapply(retained, function(x) isTRUE(x$success), logical(1))))
    return(fail("score_recovery_original_root_unresolved"))
  runs <- lapply(context$starts, backend$hl_fo_score_optimize_start,
                 lower = context$lower, upper = context$upper, evaluator = ev,
                 keep_diagnostics = keep_diagnostics)
  algebraic <- Filter(function(x) isTRUE(x$success), runs)
  attempts <- lapply(algebraic, hl_fo_score_refine_candidate,
                     evaluator = ev, context = context, backend = backend, keep_diagnostics = keep_diagnostics)
  if (ev$exhausted()) return(fail("score_recovery_evaluation_cap", runs, attempts))
  if (ev$cold_failed()) return(fail("score_recovery_cold_evaluation_failed", runs, attempts))
  if (length(ev$errors())) return(fail("score_recovery_evaluation_error", runs, attempts))
  if (any(!vapply(attempts, function(x) isTRUE(x$success), logical(1))))
    return(fail("score_recovery_candidate_unresolved", runs, attempts))
  selection <- backend$hl_fo_score_select_candidates(c(retained, attempts))
  if (!isTRUE(selection$success)) return(c(fail(selection$status, runs, attempts), list(selection = selection)))
  state <- selection$selected$state
  state$object_type <- "full_chi_score_null"
  value <- backend$hl_fo_score_postroot(state$augmented$psi_subject, state$augmented$A,
                                       context$R, state$beta, context$rhs, length(context$active))
  if (identical(value$status, "score_adjusted_sensitivity_singular")) {
    value$algebra_status <- value$status
    value$status <- "score_nuisance_singular"
  }
  value$statistic <- value$stat
  value$p_value <- value$p
  value$null_fit <- state
  value$null_fit_status <- "ok"
  value$root_certification <- state$certification
  value$root_candidates <- if (isTRUE(keep_diagnostics)) selection$candidates else
    lapply(selection$candidates, .hl_fo_score_compact_run)
  list(success = isTRUE(value$success), status = value$status, value = value,
       original_refinements = retained, runs = runs, certification_attempts = attempts,
       selection = selection, evaluations = ev$count(), errors = ev$errors(), evaluation_failures = ev$failures(), trace = ev$trace())
}

hl_fo_score_retry <- function(record, max_evaluations = 3000L) {
  failure <- function(status) list(success = FALSE, status = status, route = "withheld",
    value = .hl_fo_score_root_failure(status), primary = if (is.list(record) && is.list(record$payload)) record$payload$native else NULL)
  if (!is.list(record) || !inherits(record, "hl_score_recovery_record") || !is.list(record$payload) ||
      !identical(record$identity, .hl_score_fingerprint(record$payload)))
    return(failure("score_recovery_identity_mismatch"))
  current <- tryCatch(.hl_score_implementation(), error = function(e) NULL)
  if (is.null(current)) return(failure("score_recovery_implementation_unavailable"))
  if (!identical(record$payload$implementation, current))
    return(failure("score_recovery_implementation_mismatch"))
  if (length(max_evaluations) != 1L || !is.finite(max_evaluations) ||
      max_evaluations < 1 || max_evaluations > .Machine$integer.max ||
      max_evaluations != as.integer(max_evaluations))
    stop("max_evaluations must be one positive integer.")
  native <- record$payload$native
  if (isTRUE(native$success)) return(list(success = TRUE, status = native$status,
    route = "preserved_native", value = native, primary = native))
  if (!identical(native$status, "score_null_multiple_roots"))
    return(list(success = FALSE, status = native$status, route = "unsupported_native_failure",
                value = native, primary = native, execution = record$payload$execution))
  ctx <- record$payload$context
  if (is.null(ctx) || !length(native$root_candidates)) return(failure("score_recovery_context_missing"))
  answer <- tryCatch(.hl_score_run_retry(ctx, native$root_candidates, max_evaluations,
                                           record$payload$request$keep_diagnostics),
    error = function(e) list(success = FALSE, status = "score_recovery_internal_error", reason = conditionMessage(e)))
  if (isTRUE(answer$success)) {
    answer$value$constrained_fit <- native$constrained_fit
    answer$value$constrained_fit_status <- native$constrained_fit_status
    answer$value$null_runs <- answer$runs
    answer$value$starts <- ctx$starts
    answer$value$root_trace <- answer$trace
  }
  list(success = isTRUE(answer$success), status = answer$status,
       route = if (isTRUE(answer$success)) "stable_gamma_cold_score_refinement" else "withheld",
       value = if (isTRUE(answer$success)) answer$value else native,
       primary = native, context = ctx, identity = record$identity,
       execution = record$payload$execution, attempt = answer)
}
