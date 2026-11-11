\name{scrope_jln_research}
\alias{scrope_jln_research}
\title{Opt-in JLN boundary and search-warning research workflow}
\usage{
scrope_jln_research(count, id, pred, target, offset = NULL,
  boundary = TRUE, dominated_warning = TRUE, cpc = 0.005,
  mincp = 5L, ncore = 1L, keep_diagnostics = TRUE, verbose = TRUE)
}
\arguments{
\item{count}{Gene-by-cell count matrix, sparse matrix, or one response vector.}
\item{id}{Complete donor identifiers, with cells contiguous within each donor.}
\item{pred}{Finite numeric design matrix with uniquely named columns and an intercept of ones in the first column.}
\item{target}{One prespecified named non-intercept coefficient. This interface reports its zero-null Wald test, not general linear contrasts.}
\item{offset}{Positive cell exposures, not log exposures; defaults to ones.}
\item{boundary}{Enable the fixed-positive-dispersion fallback.}
\item{dominated_warning}{Enable corroborated interior-root reporting with an unresolved-search warning.}
\item{cpc, mincp}{Native expression-filter thresholds, as in \code{scrope}. Only analyzed genes appear in the summary.}
\item{ncore}{Positive integer number of workers for native JLN. Research recovery is serial.}
\item{keep_diagnostics}{Retain per-gene research and native diagnostics. FALSE removes diagnostic objects after computing the same summary.}
\item{verbose}{Print native fit progress.}
}
\description{
An explicit research entry point using ordinary JLN first, a narrowly qualified
interior-root warning policy, and a fixed-positive-dispersion fallback.
This compatibility wrapper preserves the earlier opt-in rule and explicitly
uses strict native \code{scrope()} fits. The main interface now supports
automatic recovery, including identifiable limiting faces.
}
\details{
The native fit uses joint-LN, L-BFGS-B, the default dispersion box, Wald inference,
and no operational fallback. Accepted native results are retained. The warning
route requires at least two distinct original core starts at compatible certified
interior roots and clear objective dominance over retained unresolved boundary
observations. It retains the full joint sandwich and warns that an unvisited
competing branch is not excluded.

The boundary route requires an eligible best finite boundary candidate, no
unrelated certificate failure, and clear superiority to any certified interior
competitor. Only implicated coordinates are fixed: donor variance at 1e-4 and/or
cell overdispersion at 1e-3 (NB size 1000). All remaining parameters are jointly
refitted using the original positive-domain JLN criterion. The free Hessian and
donor scores are selected before covariance inversion; no extra donor-count
factor, ridge or generalized inverse is used. Endpoint derivative, curvature,
root, Newton, start-agreement and local boundary-direction checks remain required.
Unrelated rank and nuisance-limit failures remain unavailable.

Fixed-route inference concerns a regular fixed-mask working target under the
donor M-estimation assumptions; the scientific target need not equal it.
Warning-route inference assumes selection of the intended regular branch.
Automatic route selection is not selectively calibrated. Empirical undercoverage
has been observed in finite-donor studies. No general coverage or FDR guarantee
is claimed. P-values and intervals are research outputs, not a validation
certificate. This interface does not provide adjusted Score, ratio tests or
research random effects.

Native availability warning text is retained in \code{primary_warnings} and
replaced by one combined research-availability warning. Other native warnings
are still emitted. Disabling both research flags retains only native inference.
Gene IDs, rather than potentially duplicated display names, identify rows.
}
\value{
A list with \code{summary}, \code{primary}, \code{diagnostics},
\code{primary_warnings}, and \code{policy}. The summary contains one row per
analyzed gene, with gene ID, target, route, availability, coefficient, robust SE,
p-value, normal 95 percent interval, attempted fixed values, inferential scope,
reason and warning. Unavailable research coefficients and inference fields are
NA; native diagnostic endpoints remain in \code{primary}. Routes are
\code{ordinary_jln}, \code{ordinary_jln_with_search_warning}, \code{fixed_v},
\code{fixed_phi}, \code{fixed_v_phi}, or \code{withheld}.
The native \code{scrope} result is retained in \code{primary}; diagnostic
retention follows \code{keep_diagnostics}. \code{se_model}, \code{p_model} and
\code{model_available} use the selected free inverse curvature, following the
JLN working-model convention. Optional adjusted tests remain unsupported.
}
\seealso{\code{\link{scrope}}, \code{\link{group_cell}}}
\examples{
\dontrun{
fit <- scrope_jln_research(count, id, pred, target = "CFTR_centered",
                          offset = exposure)
fit$summary[, c("gene", "route", "available", "coefficient", "se_robust")]
fit$summary$warning
}
}
