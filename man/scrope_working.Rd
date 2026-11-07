\name{scrope_working}
\alias{scrope_working}
\title{Explicit working-target HL and PGMM recovery}
\usage{
scrope_working(count, id, pred = NULL, offset = NULL,
  model = c("HL", "PGMM"), cpc = 0.005, mincp = 5,
  keep_diagnostics = FALSE, verbose = TRUE)
}
\arguments{
\item{count}{Gene-by-cell count matrix, dense or sparse.}
\item{id}{Subject identifiers aligned with columns of count.}
\item{pred}{Cell-level design matrix; defaults to an intercept.}
\item{offset}{Optional positive exposure offsets.}
\item{model}{Either \code{"HL"} or \code{"PGMM"}; no switching between estimators.}
\item{cpc}{Minimum counts-per-cell threshold for gene filtering.}
\item{mincp}{Minimum number of nonzero cells required to keep a gene.}
\item{keep_diagnostics}{Retain full native and recovery objects. The default
retains compact acceptance/failure summaries and a bounded failure witness,
without complete event ledgers or donor arrays. Computation and acceptance
are unchanged; native fitting is sequential in this interface.}
\item{verbose}{Print progress messages.}
}
\description{
Runs the requested native estimator first, preserving its result in
\code{primary}. An explicit call opts into separately labelled estimates
from certified positive-dispersion working fits or identifiable limiting
equations. Ordinary \code{scrope()} and \code{scrope_hl()} now assess eligible
working routes automatically; \code{recovery="none"} selects native-only behavior.
}
\details{
Recovery uses the default dispersion guards and fixed donor parameter
\code{sigma2=1e-4} and, for eligible HL cell boundaries, NB precision
\code{phi=1000}. These positive values define working models. A nuisance
limit must be one-dimensional and remove only whole all-zero donors.
Nonestimable coefficients remain missing. Free-donor-dispersion HL limits
return points only, without SEs or p-values, including when NB precision is fixed.

Routes are \code{native}, \code{certified_free_refit},
\code{fixed_positive_dispersion}, \code{identifiable_nuisance_limit},
\code{identifiable_limit_fixed_dispersion}, or \code{withheld}.
Inspect \code{routes$fixed_donor_parameter} and \code{routes$fixed_nb_precision};
\code{NA} means not fixed, not zero. Reported HL coefficients use original-design
coordinates; PGMM preserves its centered unit-intercept convention with input-scale
slopes. Their intercept estimability can differ. See
\code{provenance$coefficient_coordinates}.

Recovered inference is coefficient-wise donor-robust Wald only. No model-based
covariance, optional tests, arbitrary contrasts or random effects are supplied
for recovery. Outcomes select the route. Working-target equality to the
scientific target, selective coverage, finite-donor calibration and FDR control
are not established by these computational checks. No finite-cluster multiplier
or pseudoinverse is used. Successful regular native fits are retained exactly.
The raw sandwich requires an identified regular free root, suitable independent-
donor moment/LLN/CLT conditions and no dominating donor. Aggregate expected-score
zero alone does not justify raw donor outer products under arbitrary heterogeneous
fixed designs. Numerical certification does not establish these sampling conditions.

HL recovery captures and audits outer profile observations across guard rounds,
trust/polish searches, reduced-face fits and working refinement independently
of optional diagnostic retention. It validates original scalar evidence and
compares cold profiles under the selected configuration; raw objectives from
different targets are not compared. Tied root-like profiles require complete-root
and covariance agreement. Failed no-profile exploratory attempts remain visible;
missing finite evidence, inconsistent provenance, required replay failures,
better competitors or budget exhaustion refuse recovery. The limits are 10,000
captured events, 64 MiB of configuration/event payload accounting and 5,000
unique selected-configuration replays. This is not a process-memory or wall-time
limit. The installed guide describes scope and differences from JLN.
Ordinary accepted native fits retain their existing preservation policy.
}
\value{
A list with long coefficient \code{summary}, per-gene \code{routes}, untouched
native \code{primary}, compact per-gene \code{audit}, optional full
\code{diagnostics}, and \code{provenance}. Point, covariance and Wald
availability are separate. \code{native_available} is the original reported Wald
flag; \code{native_regular} records the preservation checks. Use \code{available}
together with \code{scope} and \code{sampling_target} for the selected output.
A reported coefficient variance does not establish joint covariance rank or
calibration.

Use \code{summary} for reported coefficients and uncertainty. Unmasked
\code{beta_report}/\code{V_report} representatives are retained only in full
recovery diagnostics; nonestimable entries are not reportable estimates.
For recovered PGMM fits, selected full-system parameters and equations are at
the top level of \code{diagnostics[[gene]]}; its \code{fit} is the native input
endpoint. HL's recovered \code{fit} contains its selected augmented state.
Native public results remain in \code{primary}; full fits are optional.
On withheld routes, fixed values may describe an attempted mask rather than
an available selected result. \code{failure} describes coefficient-level
withholding separately from the fit's certificate.

Compact HL audits retain identities, checked/generated/pending/dropped counts,
limits, comparison summaries and one bounded failure witness. With
\code{keep_diagnostics=TRUE}, \code{history_audit$observations} contains detailed
processed observations and \code{history_audit$history_evidence} retains the
augmented refusal ledger; top-level \code{history_evidence} is the original
collector snapshot. Full event/replay ledgers are absent in compact mode.

}
\seealso{
\code{\link{scrope}}, \code{\link{scrope_hl}}.
The installed guide is available at
\code{system.file("doc", "hl-pgmm-working-recovery.md", package="scRoPE")}.
}
