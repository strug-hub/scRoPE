> Retention update: detailed audit observations, complete histories and full fits
> are returned only with `keep_diagnostics=TRUE`. Compact mode retains resource
> and acceptance summaries plus a bounded failure witness. Mandatory evidence is
> still checked. Historical descriptions of retained event ledgers below refer
> to detailed mode. The default public `scrope()`/`scrope_hl()` API now uses
> automatic eligible recovery and compact retention.

# HL and PGMM working fits

`scrope_working()` is the explicit long-format working-fit entry point. It first
runs the requested estimator and preserves that native result in `primary`.
Supported recovery can supply an estimable coefficient from a refitted working
model or a limiting face. It never switches between HL and PGMM. Ordinary
`scrope()` now assesses the same eligible working routes automatically; use
`recovery="none"` for native-only behavior. See [the common guide](automatic-recovery.md).

## Calling the interface

```r
working <- scrope_working(count, id, pred, offset = exposure,
                          model = "PGMM", # or "HL"
                          cpc = 0.005, mincp = 5,
                          keep_diagnostics = FALSE, verbose = TRUE)
working$summary[, c("gene", "coefficient", "beta", "robust_se", "p_value",
                    "estimable", "point_available", "covariance_available",
                    "available", "route", "scope", "sampling_target", "failure")]
working$routes
working$audit[[1]]
native <- working$primary
```

Counts, donors and design follow `scrope()` conventions: genes by cells, one
donor ID and design row per cell, and positive exposure offsets (not log offsets).
Gene filtering still uses `cpc` and `mincp`. Specify the model explicitly.
There is no `ncore`, custom contrast or optional-test argument in this initial
wrapper. Calls are sequential; compact diagnostics are the default. There is
no external wall-time limit built into this API, and HL native searches can be
slow. Full histories with `keep_diagnostics=TRUE` can consume substantial memory.

## Reading the returned output

| Component/field | Meaning |
|---|---|
| `primary` | Native result from this package build, before working recovery. |
| `summary` | One row per retained gene and reported coefficient. |
| `point_available` | A finite selected point exists for an estimable coefficient. |
| `covariance_available` | The coefficient covariance passes the applicable output gate. |
| `available` | The selected coefficient Wald result is reported; read its target and route too. |
| `native_available` | Original native Wald flag, not the preservation or recovery decision. |
| `native_regular` | Native preservation decision under the recorded numerical policy. |
| `estimable` | The reported functional annihilates any verified escaping direction. |
| `failure` | Coefficient-level withholding reason; it is separate from fit-level certification. |
| `routes` | Per-gene route, fixed parameters, native preservation, scope and timing. |
| `audit` | Compact numerical acceptance/refusal ledgers and diagnostic objects. |
| `diagnostics` | Full internal fitting/estimating-equation objects only when requested. |
| `provenance` | Model, implementation identity, coefficient-coordinate convention and limitations. |

`available=TRUE` is a computational result. It does not establish scientific-target
equality, valid selection-adjusted coverage, finite-donor calibration or FDR control.
The selected scope is `native_estimator`, `data_selected_working_target`,
`point_only_hl_limit_sampling_unresolved`, or `unavailable`.
`sampling_target="regular_working_equation_target"` denotes the declared working
system; it does not turn a fixed-mask or selected-face fit into the native target.

## Route labels and fixed values

| Route | Interpretation |
|---|---|
| `native` | Preserve an accepted regular native endpoint and its coefficient output. |
| `certified_free_refit` | Numerical refit with all dispersions free, passing recovery checks. |
| `fixed_positive_dispersion` | Actual refit with an eligible positive dispersion held fixed. |
| `identifiable_nuisance_limit` | Verified one-dimensional face removes whole all-zero donors; only estimable coefficients can be reported. |
| `identifiable_limit_fixed_dispersion` | A verified face and fixed-positive-dispersion refit are both used. |
| `withheld` | The supported numerical/geometry/estimability contract was not established. |

Fixed values are explicit in `routes$fixed_donor_parameter` and
`routes$fixed_nb_precision`; `NA` means the parameter was not fixed, not that
it equals zero. Details also appear in `audit[[gene]]$fixed`.
For a `withheld` route, nonmissing fixed values can describe the attempted mask;
they do not imply an available selected result.

- HL or PGMM donor parameter: `sigma2=1e-4` when its lower-face conditions pass.
- HL NB precision: `phi=1000` when its upper-face conditions pass, possibly
  together with fixed donor dispersion. PGMM has no cell-level dispersion.

These are positive working-model constants, not exact zero-variance/Poisson laws.
All remaining coordinates are refitted and their own estimating equations supply
the covariance. No covariance is obtained by deleting rows from an old inverse.
The Gamma donor parameter is the model's `sigma2` parameterization; it should not
be silently relabelled as the exact variance of log-Gamma effects. HL's `phi` is
NB precision. JLN's reported cell-overdispersion `phi` is an inverse-size quantity,
so equal parameter names do not imply equal numeric values across methods.

## Coefficient coordinates and estimability

HL reports original-design coefficients. PGMM preserves its centered unit-intercept
convention, with slopes in input predictor units. Read
`provenance$coefficient_coordinates`; intercepts need not be directly comparable.
For example, if group=0 donors remain and group=1 donors are wholly zero, an
original-design intercept can remain identifiable while a centered intercept
containing a multiple of the escaping group coefficient does not. This can give
different reported estimability patterns across methods without a coding error.

The face detector canonicalizes a unique nonzero constant intercept, including
permuted or rescaled intercept columns, then maps back to the model's reporting
coordinates. Multiple recession directions, partial-donor removal, unverified
near-identities and unsupported designs remain refusals. A finite slope does
not make the escaping intercept or nuisance coefficient identifiable.

Use `summary` for reported points, SEs and p-values. Compact `beta_report` and
`V_report` are **unmasked internal diagnostic representatives**: their entries
for nonestimable coefficients must not be used as estimates or zero-variance
inference. Public masking and the `estimable` field are authoritative.

## Why the three HL limits remain point-only

PGMM's vanished-donor limiting contribution is zero. For HL, with donor parameter
v and a=1/expm1(v), each vanished donor retains the Laplace contribution

    r(v) = (a-1/2)log(a) - log Gamma(a) - a + log(2*pi)/2.

When v is estimated, its nonzero dispersion score and sensitivity must stay in
the HL equations. A raw donor outer product is not automatically an exact
fixed-design covariance when individual donor-score means are nonzero, and an
outcome-selected face introduces a separate target/selection question. Therefore
**every recovered HL face with free donor dispersion is point-only**, even if
the NB precision happens to be fixed. Its recovered SEs and p-values are missing.
When donor dispersion is fixed, r(v) is constant in the remaining free system;
this permits the separately qualified fixed-working-model route if all checks pass.

In the matched HLCA panel, REQ038/NHERF1, REQ039/NHERF2 and REQ040/SLC26A9 are
the three point-only cases. Subsequent research derived an explicitly hypothetical
fixed-face covariance-upper-bound reference and studied limited symmetric null
simulations. That reference, its intervals and any finite-donor adjustment are
**not part of this package**. The research does not establish general
selection-adjusted coverage or scientific-target equality for those HLCA cases.

## Numerical checks and their scope

The feature follows the same principles as JLN: explicit eligibility, actual
free-coordinate refitting, root/curvature checks, a sandwich from the same
equations, coefficient estimability and working-target disclosure. The estimators
and audit implementations are distinct.

PGMM recovery uses stable algebra for the exact marginal criterion, donor-summed
scores, independent finite-difference checks, agreeing deterministic starts and
same-mask comparisons of retained native and working endpoints. Its covariance
is `tcrossprod(solve(A,U))`, with no inverse score meat, ridge, pseudoinverse or
finite-donor multiplier. Singular meat can still support a positive scalar
variance; arbitrary joint recovery tests are not supplied.

HL keeps the full augmented system, including auxiliary chi and all free
dispersions; its sensitivity need not be symmetric. Independent perturbations
cold-solve donor modes and check donor/summed derivatives and root displacement.
Held-coordinate boundary support is checked on the reprofiled criterion.
Certified candidate roots are cold-refined and must agree as complete roots;
retained native/boundary-round/selected terminal endpoints are cold-profiled
under the selected mask. Failed required cold evaluations remain refusals.

**HL search-history audit:** the opt-in wrapper captures every outer profile
call across native guard rounds, trust proposals, gradient/cold polishing,
reduced-face searches and working refinement. Capture is independent of
`keep_diagnostics`. A record identifies its full log-dispersion coordinates,
profile beta and F (when available), success/status, source and warm/cold origin.
Shared configurations bind data, vanished donors, affine restrictions, held
coordinates, the named HL criterion and the loaded backend identity.

For finite evidence the audit cold-solves donor modes at the original beta/rho
and independently reconstructs the original scalar. It then cold-profiles the
recorded dispersion coordinates under the selected target and fixed mask.
Original F and cold F are both compared when the scalar targets and full
coordinates agree; an active-mask difference alone does not make those scalars
incomparable. Changed faces or held values instead generate explicitly labelled
projected dispersion probes. Raw F from different targets is never compared.
Contradictory duplicate records cannot be hidden by replay caching.

A tied, root-like cold probe must coalesce as a complete beta/rho/eta/A/J root;
its covariance must also agree to support inference. Already certified complete
matches need no redundant polishing. Other root-like ties undergo bounded
refinement, with every generated trial appended to the same audit queue.
Original certified-candidate refinement is unchanged. Selected-objective comparisons
use the existing `1e-10*(1+abs(F))` scale; scalar reconstruction/cold agreement and
complete-root/covariance comparisons use their declared `1e-8` scales. These
are numerical tolerances, not proved floating-point error bounds.

A failed exploratory attempt without a finite HL profile remains visible but
is not automatically a competing root. Finite evidence with missing coordinates,
incompatible provenance, scalar disagreement, a required replay failure or a
better profile refuses recovery. Existing stricter terminal-endpoint checks
remain in force. Ordinary accepted native results bypass recovery vetoes, with
that policy recorded in `history_capture`.

Capture is bounded at 10,000 events and a 64 MiB accounting budget for unique
configuration/event payloads; this is not a bound on process memory. At most
5,000 distinct selected-configuration replays are attempted. Audit-generated
refinement calls additionally obey the event/byte limits and the existing bounded
refinement procedure. Exhaustion is a visible refusal. No silent truncation or
rounded-coordinate merging is used. Refusals retain captured pending evidence,
processed/pending counts, and observed/retained/dropped counts. Inspect
`audit[[gene]]$history_audit` and `history_capture`. On an audit refusal,
`history_audit$history_evidence` includes the complete augmented ledger with
generated pending probes; top-level `history_evidence` is the original collector
snapshot before audit-generated trials.

This closes the former terminal-only coverage gap for retained HL profile
observations. HL and JLN still use different criteria and audit algorithms:
HL's nested criterion requires the inner p-score root, while JLN audits its
enumerated joint-objective candidate/recovery/basin records. Inner mode/beta
solver iterations without a valid HL profile, unevaluated guard proposals and
derivative perturbations are not competing HL profile observations. Held-coordinate
boundary checks remain separate. Neither audit proves global optimality,
scientific-target equality, selection-adjusted inference or general calibration.

Native preservation also has a recorded policy. PGMM keeps qualified native
fits under its existing 1e-3-SE Newton displacement limit; recovery uses the
stricter 1e-4-SE and scaled-root checks. A stricter recovery diagnostic alone
does not replace an accepted native result.
HL preserves an accepted native point when there is no unresolved face or
recession. It retains the native inference flags; it does not impose the new
recovery finite-difference audit on an already accepted native result.

Selected public coordinates are in the masked `summary`. Unmasked
`beta_report`/`V_report` and full system objects are available only with full
diagnostics. For PGMM, top-level `par`, `beta`, `A`, `U`, `V`, `mapping` and `fixed`
describe the selected recovered system; its nested `fit` is the native input
endpoint. HL's recovered `fit` is its selected augmented state. Preserved native
routes have no selected recovery state; full native fits are in `primary`
diagnostics when retained. Compact audit output contains summaries, not these
full fit objects or complete event/replay ledgers.
Do not interchange these objects or coordinate systems.

## Inference assumptions and unsupported output

The expected free-equation root defines a working target on a fixed/stably
selected face and fixed mask. A conventional raw sandwich needs an identified
regular root, independent donors, appropriate iid/superpopulation or separately
justified individual-score-mean conditions, sufficient moments, a donor LLN/CLT
and no dominating donor. Aggregate expected-score zero alone does not establish
the raw sandwich for arbitrary heterogeneous fixed designs. Unequal cell counts
must satisfy the corresponding scaling and moment conditions.

Outcomes select routes and faces here. Numerical certification does not prove
these sampling assumptions, stable face selection, target equality or finite-
donor coverage. The subsequent limited HL null study found 50/59 coverage at
12 donors, with similarly low coverage in regular controls; larger controls were
compatible with 95% but imprecise. It is research evidence, not an API calibration
guarantee or a validated threshold based only on donor count.

Recovered output supplies coefficient-wise donor-robust Wald inference only.
Recovered model-based SEs, Score/LRT, arbitrary joint contrasts and random effects
are not provided. Genuine nonexistence, nonidentifiability, unresolved competitors
or failed numerical checks remain explicit outcomes.

## Numerical prerequisites and review evidence

The combined patch includes two previously reviewed numerical components:

1. PGMM direct-sandwich calculation avoids an unnecessary inverse of donor-score
   meat. It preserves the historical object when all requested blocks agree in
   relative variance/SD-product units, including tiny variances. Joint rank gates
   prevent singular contrast covariance from being re-enabled by a ridge or an
   optional-test fallback. This can affect numerically defective ordinary PGMM
   output without changing its estimator or public defaults.
2. `hl_score_recovery="stable"` is a separate optional Score retry through
   `scrope()`/`scrope_hl()`. Its default remains `"none"`. It preserves successful
   native Score results, checks data/backend identity and withholds genuine
   multiplicity or failed cold evaluations. It adds no Wald recovery and is not
   an option of `scrope_working()`.

In the completed matched49-case HLCA assessment, the focal coefficient had
49/49 HL points and46/49 HL Wald outputs (43 native,3 fixed-working,3 point-only),
and49/49 PGMM Wald outputs (43 native,3 limits,1 fixed-working,2 free numerical
refits). Legacy PGMM already emitted49/49 raw Wald flags; the candidate qualifies
and refits nonregular cases rather than claiming49 newly available p-values.
Nonestimable coefficients remain missing. These are computational availability
counts, not calibration evidence.

The sealed working-fit archive passed928 focused assertions. An earlier feature
build passed6746 assertions with two parallel skips; it preceded the last PGMM
changes and is not a full-suite test of the final archive. The numerical
prerequisite separately passed6673 assertions under Accelerate and reference
BLAS on the same macOS system, with the same two parallel skips. The mandatory
two-worker gate could not execute because only one allocated core was detected.
No second-OS validation is claimed. The known pre-existing multi-contrast LRT
zero-quadratic issue remains separate from this patch. The review bundle records
exact source hashes, validation stages, current focused results and remaining
limitations; earlier evidence is not relabelled as a new full-suite pass.
