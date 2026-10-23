# Methods, inference availability and migration

This guide accompanies the scRoPE 0.3.0 local release candidate. Start with
[the package quick start](../../README.md) for the 24-gene omitted-slope example.
The additional route examples below use the separate, smaller `sample_data`
teaching object (8 genes, 30 donors); they do not use `toy_data`.

```r
library(scRoPE)
data(sample_data)
design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
```

## Supported tests and output fields

`model` chooses the fitting engine; `additional_tests` requests supported
inferential tests. The default engine is joint-LN. Coefficient-wise Wald
results require no additional-test option.

| Engine | Model-based output | Donor-robust output | Optional requests |
|:--|:--|:--|:--|
| Joint-LN | `se_*`, `p_naive_*` | `se_robust_*`, `p_robust_*` | Experimental scalar-zero Score: `score_adj_*`, `p_score_adj_*`; Pace ratio: `lrt_adj_*`, `p_lrt_adj_*` |
| Corrected first-order HL | `se_*` and `p_naive_*` intentionally `NA` | `se_robust_*`, `p_robust_*` | `additional_tests="score"`: `hl_score_adj_*`, `p_hl_score_*` |
| Exact PGMM | `se_*`, `p_naive_*` | `se_robust_*`, `p_robust_*` | Adjusted Score: `score_adj_*`, `p_score_adj_*`; scaled adjusted LRT: `lrt_adj_*`, `p_lrt_adj_*` |

For example, using the separate `sample_data` and `design` defined above:

```r
fit_hl_score <- scrope(sample_data$count, sample_data$sid,
  pred = design, offset = sample_data$offset, model = "HL", ncore = 1,
  additional_tests = "score")
fit_pgmm_tests <- scrope(sample_data$count, sample_data$sid,
  pred = design, offset = sample_data$offset, model = "PGMM", ncore = 1,
  additional_tests = c("score", "lrt"), lrt_details = TRUE)
```

Inspect `fit$provenance$test_capabilities` and per-test success/failure fields.
A capability is not a guarantee that every fit or null constraint succeeds.
An LRT labelled `wald` or `wald_ridge` is a fallback Wald result, not an accepted
adjusted LRT. That fallback behavior belongs to PGMM; JLN optional tests do
not substitute a Wald result. Raw objective/likelihood quantities in diagnostics
are not a separate public model-based Score/LRT test API. Corrected HL LRT is
unsupported. JLN's optional tests have the experimental scope described below.

Public JLN Wald inference is per coefficient. `lrt_contrasts` must not be
presented as a supported JLN custom-Wald interface. It can request eligible
experimental Score/Pace tests or broader diagnostic-only calculations.
Corrected HL accepts linear
contrasts for robust Wald/Score. PGMM accepts adjusted-test contrasts,
including nonzero nulls and robust Wald diagnostic/fallback calculations.
Coordinate conventions and test-specific details are described below.

### JLN experimental Score and Pace tests

JLN offers two explicitly requested **experimental, asymptotic** tests:
the donor-robust efficient Score and the Pace-adjusted working-objective ratio.
The default call still returns coefficient-wise model-based and donor-robust
Wald inference only. This candidate exposes existing diagnostic statistics;
the request does not change the fitted coefficients or Wald calculation.

- `additional_tests = "score"` requests robust Score.
- `additional_tests = "lrt"`, or its alias `lrt = TRUE`, requests ratio-related
  inference. For JLN, this is the **Pace-adjusted working-objective ratio**,
  not an ordinary exact-likelihood LRT. Use `c("score", "lrt")` to request both.
- `lrt_contrasts` specifies hypotheses in reported coefficient coordinates.
  With LRT requested, `lrt_details = TRUE` returns a tidy diagnostic table.
  If neither optional test was requested, either `lrt_details = TRUE` or a
  non-`NULL` `lrt_contrasts` also activates the LRT request. A Score-only
  request plus `lrt_details = TRUE` does not additionally request LRT or
  produce that table.
- There is no `additional_tests = "ratio"` or `"pace"` option.

Initial public scope is **one non-intercept coefficient tested against zero**.
This includes cell- and donor-level predictors. A custom contrast must have
exactly one row and exactly one finite nonzero entry selecting a non-intercept
coefficient, with `b` exactly zero. Any finite nonzero rescaling of that selector
expresses the same null; fitting and numerical checks still apply. Otherwise
valid nonzero nulls, mixed-coefficient contrasts, intercept-containing contrasts
and multiple-row hypotheses remain diagnostic-only. Their public statistics
and p-values are `NA`, with reason
`joint_ln_optional_scope_unsupported`; this scope limit is not a claim that
their internal formulas are invalid.

Both public optional tests require a certified unrestricted primary JLN fit
and a certified constrained fit, preserving the acceptance policy used in the
retained study. A native Score diagnostic can exist without an accepted
primary fit, but that larger set is not exposed here. Score additionally
requires valid nuisance-adjusted donor-score variability. Pace also requires
matching objective provenance at both endpoints, effective curvature and
valid ratio components. Passing the primary Wald checks does not guarantee
optional-test availability. No ridge, clipping, Wald substitution or legacy
fallback is used to rescue a failed JLN optional test.

```r
fit_joint_tests <- scrope(sample_data$count, sample_data$sid,
  pred = design, offset = sample_data$offset, ncore = 1,
  additional_tests = c("score", "lrt"), lrt_details = TRUE)
fit_joint_tests$summary[, c(
  "gene", "p_score_adj_X1", "score_status_X1", "score_failure_X1",
  "p_lrt_adj_X1", "lrt_status_X1", "lrt_failure_X1"
)]
fit_joint_tests$provenance$test_capabilities
fit_joint_tests$test_status
```

For an eligible, accepted request, `score_method_*` is
`joint_ln_robust_score` and `lrt_method_*` is `joint_ln_pace_adjusted_ratio`.
The corresponding `score_status_*` / `lrt_status_*` value is
`experimental_asymptotic`. Success records numerical availability, not a
finite-sample calibration guarantee. Failed or out-of-scope requests have
unavailable status, typed `NA` statistics/tails and explicit per-test reasons.
These public statuses remain available with `keep_diagnostics = FALSE`.
`test_capabilities` labels the optional scope as
`scalar_zero_null_single_nonintercept` and availability as
`experimental_asymptotic_when_requested_and_eligible`.

With `keep_diagnostics = TRUE`, inspect a constrained entry under
`fit$diagnostics$raw[[gene_index]]$lrt[[contrast_name]]$constrained$ln_joint`:
`score_diagnostic` contains the robust efficient Score calculation;
`pace_ratio_diagnostic` contains the adjusted working-objective ratio.
Each has its own `diagnostic_available` and `failure_reason` fields, while
`public_available` remains `FALSE` and `public_status` is `diagnostic_only`.
Those fields describe the internal helper; the public wrapper independently
applies the request, scope and acceptance gates before generating a tail.
Keeping a diagnostic object does not broaden the public scope.

The raw ratio `wP` is twice the difference between unrestricted and constrained
JLN working log objectives. The Pace adjustment multiplies it by
`wU / qP`, using the robust efficient Score quadratic `wU` and its curvature-based
counterpart `qP`. Public Score uses `wU`; public Pace uses `wP * wU / qP`.
Their p-values are upper tails of chi-square with one degree of freedom.
A finite raw `wP` is not a public model-based LRT statistic. `p_naive_*`
continues to denote model-based **Wald** inference, not an unadjusted ratio test.

#### What the retained calibration study supports

The retained joint-LN 0.2.5.9004 study contains 15 settings with 4,000 planned
datasets each. Each dataset supplies a scalar zero-null cell predictor `x`
and donor predictor `z` for both tests: 120,000 requests per test, not 240,000
independent datasets. Settings include Gamma reference donors, moment-matched
lognormal donors and strong independent donor/cell hidden normal variation
(variance 2 at each level). Donor counts are 30, 50, 100 and 200 with about
400 cells per donor, plus 100 donors with about 40 cells. The two generating
slopes are zero. This study is separate from the bundled teaching data and
the later omitted-slope Wald illustration.

At nominal 5%, the following rates are **conditional on available results**
for 200 donors and about 400 cells per donor:

| Setting | Cell `x`: Score / Pace | Donor `z`: Score / Pace |
|:--|--:|--:|
| Gamma reference | 4.78% / 4.78% | 6.34% / 6.01% |
| Lognormal donors | 5.12% / 5.12% | 8.62% / 7.92% |
| Strong hidden variation | 5.21% / 5.21% | 16.39% / 14.60% |

Across all 15 settings, each cell-`x` pointwise Wilson 95% interval includes
5%, while every donor-`z` interval lies above 5% for both tests. These are
pointwise summaries, not simultaneous calibration certification. Smaller
cell-`x` tails are conservative in several settings. For strong hidden
variation at 200 donors, Score rejects 639/3,899 donor nulls (95% interval
15.26–17.58%) and Pace rejects 569/3,898 (13.52–15.74%). The donor-level
inflation is substantial and is not explained by Monte Carlo uncertainty or
by counting every missing result as a nonrejection. Do not treat an eligible
donor-level p-value as demonstrated nominal finite-donor inference.

Score has 118,514 available results and Pace 118,432 of 120,000 requests each.
For each test, 1,436 planned requests had held input identities and 50 failed
root checks; Pace had 82 additional effective-curvature failures. Held input
records are not optimizer failures. Retain all planned requests and report
test-specific failures separately; conditional calibration does not establish
unconditional performance when outcomes are unavailable. Do not replace
missing p-values with 1 or silently discard them from a planned analysis.

The chi-square reference requires regularity, independent donors, valid
nuisance-adjusted variability and the appropriate working target; Pace also
requires its local objective-ratio approximation. Numerical checks do not
establish these population assumptions. Null-tail evidence does not establish
non-null power, bias correction, arbitrary-contrast calibration or family-level
FDR control. Earlier pilot and deterministic checks exercised nonzero and
multivariate hypotheses, but the finalized study above covers scalar zero
nulls only. Intercept-containing contrasts also require a clearly defined
population target because JLN reports a sample-centered intercept.

## Specialist and compatibility interfaces

- Internal PGMM point/grid **relative-profile likelihood** functions support
  scalar fixed-effect targets, with ordinary and robust profiles. They are
  implemented development tools in `R/pgmm_profile.R`, not exported package
  functions or an extra `scrope()` model. They are distinct from the optional
  public adjusted LRT.
- `scrope_hl()` wraps the same corrected first-order HL engine. This route is
  sequential. JLN and PGMM have multicore execution routes.
- Explicit `nebula(method="HL")` retains a **historical model-based** interface,
  separate from scRoPE's corrected HL inference. Its `model` choices are
  `"NBGMM"` (default), `"NBLMM"` and `"PMM"`; outputs use historical model-based
  Wald uncertainty, not scRoPE donor-robust uncertainty. `reml=1` is supported
  for NBLMM; other choices reset it to zero with a warning. Do not confuse
  these compatibility models with maintained `scrope(model=...)` engines.
- Default or explicit `nebula(method="LN")`, `legacy_post_pml` and legacy
  operational fallback remain retired. Explicit requests error; none is
  silently converted to JLN. Per-gene model switching is also retired.
- `scToNeb()` converts Seurat/SingleCellExperiment inputs, `group_cell()`
  organizes donor blocks, and `nbresidual()` reconstructs optional Pearson
  residuals from compatible retained fits. These are utilities, not extra
  inferential engines. Conditional residuals need retained donor predictions.

## Interpreting the teaching illustration

The 24-gene example compares two uncertainty calculations around the same
JLN coefficient. Its diamonds are known generating population log-mean
contrasts. Mean preservation alone does not establish equality with the
approximate JLN population target under omitted slopes, especially for
non-null genes. Observed interval inclusion is not a coverage estimate,
and sandwich variance does not correct coefficient bias or a target mismatch.
Robust inference depends on independent donors, adequate donor information,
regularity and a working estimating equation that identifies the intended
target. See the [complete example design](../examples/toy24/README.md) for
parameters, RNG chronology, all outcomes and runtime qualifications.
Historical null calibration is separate from this illustrative dataset.

## Joint-LN estimation and availability

The default and the explicit request below use the same estimator:

```r
fit_joint <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  offset = sample_data$offset,
  ncore = 1,
  ln_estimator = "joint_ln"
)

fit_joint$availability
```

The retained joint-LN numerical policy can
resolve an extra nearby numerical endpoint as belonging to the same solution
basin. When that happens, model-based and donor-robust Wald inference is
reported from the unchanged certified selected root. The original strict
decision is retained in provenance. Boundary, distinct-basin, hard-stop, and
unresolved cases continue to withhold inference. Numerical certification is
distinct from finite-sample calibration; this remains a local release candidate.

`fit_joint$availability` has exactly one row per analyzed gene/request, rather
than one row per coefficient or contrast. `availability = "inference"` means
the coefficient, covariance, robust standard error, Wald statistic, and p-value
are all evaluated at the same certified joint-LN fit. When primary inference is
unavailable, the affected standard errors, Wald statistics, and p-values remain
typed `NA`; sc-RoPE does not store `p = 1` or silently substitute results from
another estimator. Technical route and failure details remain in `provenance`
and, when requested, `diagnostics`.

The joint-LN estimator estimates coefficients and dispersions together. It is
not expected to reproduce the legacy estimator exactly, and its verified
solution is local rather than a claim of global uniqueness.

Technical route labels are retained in `provenance` for reproducibility and
debugging, but users do not need them to select results. Use `inference_status`,
`resolution_status`, and `inference_source` for the primary joint-LN fit.
`strict_inference_status`, `enabled_inference_status`, and `newly_enabled`
record how the calibration-candidate rule differs from the original strict
decision. The earlier `route`, `inference_available`, and `failure_reason`
fields remain available for compatibility and expert review.

With `output_re = TRUE`, joint-LN donor modes are auxiliary fixed-beta
conditional PML predictions evaluated at the plug-in joint-LN parameters. They
cannot change coefficient inference.

Model-based and donor-robust Wald are the default JLN tests. Experimental
Score and Pace-adjusted ratio p-values require an explicit optional request
and the scope and numerical checks described above.
Switching from the legacy estimator to `joint_ln` may change coefficients,
standard errors, p-values, convergence statuses, and downstream results.

### Migration from legacy LN

`ln_estimator = "legacy_post_pml"` and
`ln_operational_fallback = "legacy_post_pml"` are retired from `scrope()` and
produce actionable errors. They are not silently translated into joint-LN.
`ln_operational_fallback = "none"` remains accepted. A rejected joint-LN fit
continues to report unavailable inference; the candidate does not obtain a
coefficient or p-value from another estimator. Per-gene model switching is
also unavailable.

Historical reproduction requires the exact original version and settings,
including its source changes. The retained simulations used multiple versions
(.9001, .9004 and .9005); a version label or the older public Git checkout alone
is not sufficient provenance. Switching an earlier default call to this
candidate intentionally changes its estimator and can change its results.

For current results, use `availability == "inference"` for certified primary
inference and retain `unavailable` outcomes in the planned denominator. Do not
replace missing p-values with 1. A finite uncertified endpoint remains a
diagnostic quantity. Historical `estimate_only` sidecars are not generated by
the maintained workflow.

The separate exported `nebula()` compatibility function also rejects its
historical default and explicit `method = "LN"` requests. Use `scrope(...,
model = "LN", ln_estimator = "joint_ln")` for the maintained LN workflow,
or the exact historical package for reproduction. Explicit
`nebula(method = "HL")` remains available as a historical compatibility API;
it is distinct from corrected first-order HL via `scrope(model = "HL")`.
No call is silently redirected from one estimator to another.

## PGMM optional adjusted Score and LRT inference

The example below explicitly selects PGMM. Its adjusted Score/LRT outputs,
eligibility and fallback behavior differ from the experimental JLN tests
described above.

Request the adjusted Score and scaled adjusted likelihood-ratio tests through
`additional_tests`:

```r
fit_tests <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  offset = sample_data$offset,
  model = "PGMM",
  ncore = 1,
  additional_tests = c("score", "lrt"),
  lrt_details = TRUE
)

head(fit_tests$summary[, c(
  "gene",
  "p_robust_X1",
  "p_score_adj_X1",
  "p_lrt_adj_X1",
  "lrt_method_X1"
)])
```

On this fixed `sample_data` example, 23 of the 24 PGMM adjusted LRTs are available.
For G08 / X2, the constrained fit passes optimization but its profiled curvature
is nonpositive. The call warns and records `lrt_method_X2 = "wald"`; that finite
fallback p-value is a Wald result, **not an adjusted LRT**. All 24 adjusted
Scores are available. Check `lrt_details` and method/availability fields before
interpreting or selecting optional tests; do not drop gated contrasts from a
planned denominator. These are example-specific outcomes, not failure-rate or
calibration estimates.

The robust Wald test is the primary default because it is directly available
after the unconstrained fit. Score and LRT inference require constrained fits
for every gene and contrast, so they are more computationally intensive. The
PGMM LRT additionally requires accepted unrestricted and constrained fits;
a positive likelihood difference alone is insufficient. Score inference
requires an accepted constrained null fit, independently of the unrestricted
Wald fit. Unavailable tests retain explicit failure reasons and missing tails.

- `score_adj_*` and `p_score_adj_*` contain adjusted Score statistics and
  p-values.
- `lrt_adj_*` and `p_lrt_adj_*` contain scaled adjusted likelihood-ratio
  statistics and p-values.
- `lrt_method_*` records whether the reported result used the adjusted LRT,
  ridge regularisation, or a documented fallback.
- `lrt_details = TRUE` returns a tidy `lrt_details` data frame containing the
  raw likelihood-ratio quantities and diagnostic information.

## Working Likelihoods within sc-RoPE

LN, HL, and PGMM are alternative working-likelihood routes within the same
sc-RoPE robust-inference framework. They are not separate inferential
frameworks: each route provides fixed-effect estimates, curvature, and
donor-level score contributions to the Godambe adjustment.

| Working-likelihood route | Call | Variance structure | Typical use |
|---|---|---|---|
| NBGMM, LN approximation | `scrope()` | Donor heterogeneity and cell-level extra-Poisson variation | Default scalable association analysis |
| NBGMM, first-order nested HL | `scrope(model = "HL")` | Donor heterogeneity and cell-level extra-Poisson variation | Explicit HL fitting or comparison with the LN approximation |
| PGMM, exact marginal likelihood | `scrope(model = "PGMM")` | Donor heterogeneity without a separate cell-level overdispersion parameter | Simpler exact working likelihood for fixed-effect inference |

The NBGMM working likelihood is the more complete variance model. PGMM is a
targeted simplification: it retains the donor-level gamma random effect but
omits the cell-level gamma overdispersion component.

For example, the corrected first-order HL route is requested through the same
front door:

```r
fit_hl <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  offset = sample_data$offset,
  model = "HL",
  ncore = 1
)
```

The existing `scrope_hl()` function remains available as a compatibility
wrapper. The current first-order HL implementation runs sequentially; asking
for `ncore > 1` produces a warning rather than silently claiming parallel
execution.

Test availability remains route-specific. PGMM supports optional adjusted
Score and LRT tests; certified joint-LN exposes model-based and donor-robust
Wald, with opt-in experimental scalar-zero Score and Pace tests; first-order HL
exposes robust Wald and Score while its objective-ratio/LRT remains unsupported.
The `provenance$test_capabilities` table records these distinctions in every
fit.

## PGMM Working Likelihood

The PGMM route integrates out its donor-level gamma random effect exactly. It
can be used with both subject-level and cell-level fixed-effect predictors.

### Subject-level predictor example

In `sample_data`, `X2` is constant within donor:

```r
design_subject <- model.matrix(~ X2, data = sample_data$pred)

fit_pmg_subject <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design_subject,
  offset = sample_data$offset,
  model = "PGMM",
  ncore = 1
)

head(fit_pmg_subject$summary)
```

### Cell-level predictor example

`X1` and `cc` vary within donors, so the same design defined at the top of this guide can
also be fitted under the PGMM working likelihood:

```r
design_cell <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

fit_pmg_cell <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design_cell,
  offset = sample_data$offset,
  model = "PGMM",
  ncore = 1
)

head(fit_pmg_cell$summary[, c(
  "gene",
  "logFC_X1",
  "se_robust_X1",
  "p_robust_X1"
)])
```

This second call prints a short message identifying the cell-level predictors
and recommending the robust results when cell-level overdispersion may be
present. These predictors are supported; the message reflects that PGMM omits
a separate cell-level overdispersion component.

When donors are the independent clusters, the mean model preserves the target,
and the number of donors is sufficiently large, the donor-level Godambe
adjustment provides large-sample robust inference for both subject-level and
cell-level fixed effects. In particular, it reduces sensitivity of the
standard errors to omitted cell-level overdispersion. This adjustment does not
make PGMM the true data-generating model, correct mean-model misspecification,
or rescue an analysis that ignores the donor as the biological unit of
replication.

PGMM can therefore be useful when fixed-effect inference is the main target and
the detailed decomposition of cell-level variability is secondary. Prefer an
LN or HL NBGMM working likelihood when cell-level overdispersion is itself of
scientific interest, when residual overdispersion is strong, or when the number
of donors is modest and the data are sparse, highly unbalanced, or strongly
heterogeneous.

The PGMM main unconstrained fit uses `nlminb` with `iter.max = 500L` and
`eval.max = 500L`. Override these with, for example,
`pgmm_control = list(iter.max = 150L, eval.max = 200L)` to reproduce the earlier
limits. The control applies to the main fit only; constrained/profile fits keep
their existing controls. It does not add retries or change starts, tolerances,
objective scaling, or acceptance rules. Existing recovery after an optimizer
error remains available, so this argument is not a total work cap. The 500/500
choice is a development default and does not guarantee convergence.

PGMM reports `optimizer_success` separately from its composite convergence code.
When optimization is unsuccessful or its native status is missing, model and
robust Wald statistics, p-values, and standard errors are unavailable. A retained
`logFC_*` or dispersion in that case describes the failed endpoint for diagnosis;
it is not an accepted estimate. Full endpoint and covariance values remain in
`diagnostics$raw` when diagnostics are requested. Model and sandwich covariance
availability are checked independently: failure of one does not automatically
suppress an otherwise valid test based on the other. `model_wald_*` fields
report the model-Wald statistic, success, and failure reason alongside the
existing robust `wald_*` fields.

Effective main limits are retained in `provenance$optimizer`, including when
`keep_diagnostics = FALSE`. Per-fit diagnostics record the actual backend,
status, controls, native iterations/evaluations and attempted optimizer calls.
Information conditioning and a scaled Newton correction describe the endpoint;
they do not move it or introduce a new acceptance threshold.

Adjusted Score and likelihood-ratio tests can be added to either PGMM example
using the same `additional_tests` and `lrt_details` arguments shown above.

## Contrast coordinates and residuals

For maintained LN and PGMM calls, custom `lrt_contrasts` use the same coefficient
units as the reported `logFC_*` values, with columns ordered as in `pred`.
The fitting design may be standardized internally; both the contrast matrix
and its null hypothesis are mapped consistently. A nonzero slope null is in
the original predictor units; for JLN it remains diagnostic-only in this
initial experimental release scope. The reported intercept for LN/PGMM is
sample-centered; an intercept-containing contrast does not automatically have
a fixed population target or the same interpretation as the intercept at
unobserved zero-valued predictors. Corrected first-order HL retains its
original-design coefficient coordinates.

`nbresidual()` recognizes corrected first-order HL output and reconstructs its
means using the original design coordinates. Supply the original `pred`,
positive `offset`, and matching cell/donor ordering. JLN and PGMM retain their
centered-intercept convention. Objects without estimator provenance retain
the historical reconstruction convention; the function cannot infer an
unrecorded historical fitting route. Conditional residuals require the
corresponding retained donor predictions.

## Understanding the Main Outputs

For each modeled predictor, `summary` includes:

- `logFC_*`: estimated fixed-effect coefficient;
- `se_*` and `p_naive_*`: model-based standard error and Wald p-value for
  JLN/PGMM; corrected HL intentionally returns `NA` in these fields;
- `se_robust_*` and `p_robust_*`: robust standard error and Wald p-value,
  returned by default;
- `wald_stat_*`, `wald_success_*`, and `wald_failure_*`: robust Wald statistic
  and diagnostics;
- optional Score results use engine-specific names: JLN/PGMM `score_adj_*` /
  `p_score_adj_*`, corrected HL `hl_score_adj_*` / `p_hl_score_*`; and
- optional `lrt_adj_*`, `p_lrt_adj_*`, and `lrt_method_*` contain JLN's
  experimental Pace-adjusted working-objective-ratio results or PGMM's scaled
  adjusted likelihood-ratio results, according to the selected engine.

For JLN optional tests, inspect `score_status_*`, `lrt_status_*`, the per-test
success/failure fields and `test_status`. Out-of-scope or failed requests remain
unavailable; corrected HL LRT remains unsupported.

The returned object also contains:

- `overdispersion`: fitted working-model overdispersion parameters;
- `convergence` and `algorithm`: per-gene fitting status;
- `availability`: one joint-LN status row per analyzed gene/request;
- `test_status`: JLN per-test availability and failure reasons, retained even
  when detailed diagnostic objects are omitted;
- `provenance`: estimator definition and technical information about each fit;
- `diagnostics`: detailed per-gene fitting and robust-inference information;
- `lrt_details`: optional tidy likelihood-ratio diagnostics;
- `random_effect`: retained donor predictions for JLN/HL when
  `output_re = TRUE`; PGMM warns and returns `NULL` for this request.

## Advanced PGMM Profile-Likelihood Diagnostics

The exact PGMM route includes internal tools for ordinary and robust relative
profile likelihoods for scalar fixed-effect targets. These tools are intended
primarily for methodological development, exact-likelihood simulations, and
likelihood-paradigm analyses.

The profile-likelihood infrastructure includes:

- equality-constrained PGMM fitting for scalar contrasts;
- pointwise and grid-based ordinary and robust relative profiles;
- warm-started constrained fits with optimizer retries;
- raw likelihood-ratio and boundary safeguards; and
- information-scaled convergence checks for high-information profile fits.

For routine PGMM association analyses, the summary-level robust Wald, adjusted
Score, and scaled adjusted likelihood-ratio results are the main user-facing
outputs.

