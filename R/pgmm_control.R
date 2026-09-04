# Public controls apply only to the main unconstrained PGMM fit. Internal
# constrained/profile fits retain their existing control and retry policies.
pmg_resolve_control <- function(control = NULL) {
  defaults <- list(iter.max = 500L, eval.max = 500L)
  if (is.null(control)) return(defaults)
  if (!is.list(control) || is.object(control) ||
      (length(control) && (is.null(names(control)) ||
       anyNA(names(control)) || any(!nzchar(names(control))) ||
       anyDuplicated(names(control)) || any(!names(control) %in% names(defaults))))) {
    stop("pgmm_control must be a named list containing only iter.max and eval.max.")
  }
  for (name in names(control)) {
    value <- control[[name]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0 || value != trunc(value) ||
        value > .Machine$integer.max) {
      stop("pgmm_control$", name, " must be a positive finite integer.")
    }
    defaults[[name]] <- as.integer(value)
  }
  defaults
}
