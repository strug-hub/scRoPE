# scRoPE

**0.2.5.9006 — local release candidate, not a published release.** This
checkout starts from the preserved 0.2.5.9005 source used for the newer
simulation campaigns. Historical results retain their original versions;
this candidate does not replace or relabel those runs.

sc-RoPE (single-cell Robust Parameter Estimation) is a framework for robust
fixed-effect inference in multi-donor single-cell RNA-seq studies. It supports
population-level differential expression and regression-based co-expression
analyses while retaining cell-level observations and accounting for donor
heterogeneity.

The framework combines a hierarchical working likelihood with donor-level
score information. The working likelihood supplies fixed-effect estimates and
observed curvature; sc-RoPE uses the variability of the donor-level score
contributions within a Godambe framework to obtain robust uncertainty
quantification. This protects inference against broad variance-model
misspecification when the mean model still identifies the fixed-effect target.

Robust Wald inference is returned by default. Adjusted Score and scaled
adjusted likelihood-ratio tests are optional and depend on the selected
estimator; they remain unavailable for the default joint-LN estimator pending
calibration.

## Local candidate and dependencies

Use this exact source for candidate validation. The GitHub package has not been
updated by this work. R and package requirements are listed in `DESCRIPTION`;
no dependency installation is performed automatically.

With those requirements already available, build and check in a separate local
output directory:

```sh
R CMD build /path/to/this/source
R CMD check --no-manual scRoPE_0.2.5.9006.tar.gz
```

The examples below assume that this local candidate has been loaded. Check
`packageVersion("scRoPE")` before comparing outputs with another version.

## Quick start

The default `scrope()` call uses the scalable LN approximation to the
negative-binomial gamma mixed-model (NBGMM) working likelihood. The design can
contain both subject-level and cell-level predictors. Joint-LN (`joint_ln`)
is now the default and the maintained LN estimator.

The `model` argument explicitly selects one complete fitting route: `"LN"`,
`"HL"`, or `"PGMM"`. sc-RoPE does not silently replace the requested route with
another route.

```r
library(scRoPE)

data(sample_data)

design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

fit <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  ncore = 1
)
```

A compact view of the default robust Wald results for `X1` is:

```r
head(fit$summary[, c(
  "gene",
  "logFC_X1",
  "se_robust_X1",
  "p_robust_X1"
)])
```

Here, `logFC_X1` is the estimated coefficient, `se_robust_X1` is its robust
standard error, and `p_robust_X1` is the robust Wald p-value. The complete
`fit$summary` also contains model-based results and results for the other
predictors.

## Joint-LN estimation and availability

The default and the explicit request below use the same estimator:

```r
fit_joint <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
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
unresolved cases continue to withhold inference. This behavior is for
simulation calibration and is not yet a calibrated public release.

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

Adjusted Score, Pace scaling, and scaled adjusted likelihood-ratio tests remain
unavailable for `joint_ln`; model-based and donor-robust Wald are the supported
tests for this option.
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

## Optional Adjusted Score and LRT Inference

The example below explicitly selects PGMM. Adjusted Score and LRT tests are
unavailable for default or explicit joint-LN pending calibration.

Request the adjusted Score and scaled adjusted likelihood-ratio tests through
`additional_tests`:

```r
fit_tests <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
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
Wald, with Score, Pace, and LRT inference calibration-pending; first-order HL
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
  model = "PGMM",
  ncore = 1
)

head(fit_pmg_subject$summary)
```

### Cell-level predictor example

`X1` and `cc` vary within donors, so the same design used in the quick start can
also be fitted under the PGMM working likelihood:

```r
design_cell <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

fit_pmg_cell <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design_cell,
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
the original predictor units. The reported intercept for LN/PGMM is
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
- `se_*` and `p_naive_*`: model-based standard error and Wald p-value;
- `se_robust_*` and `p_robust_*`: robust standard error and Wald p-value,
  returned by default;
- `wald_stat_*`, `wald_success_*`, and `wald_failure_*`: robust Wald statistic
  and diagnostics;
- optional `score_adj_*` and `p_score_adj_*`: adjusted Score results, when
  supported by the selected estimator; and
- optional `lrt_adj_*`, `p_lrt_adj_*`, and `lrt_method_*`:
  scaled adjusted likelihood-ratio results and method diagnostics.

The returned object also contains:

- `overdispersion`: fitted working-model overdispersion parameters;
- `convergence` and `algorithm`: per-gene fitting status;
- `availability`: one joint-LN status row per analyzed gene/request;
- `provenance`: estimator definition and technical information about each fit;
- `diagnostics`: detailed per-gene fitting and robust-inference information;
- `lrt_details`: optional tidy likelihood-ratio diagnostics;
- `random_effect`: donor-level random effects when `output_re = TRUE`.

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

For routine association analyses, the summary-level robust Wald, adjusted
Score, and scaled adjusted likelihood-ratio results are the main user-facing
outputs.

## Attribution

scRoPE is derived from the `nebula` package (GPL-2).

Upstream source:

https://github.com/lhe17/nebula

High-level changes in this repository include:

- Added Godambe-based robust Wald inference for LN and HL NBGMM pipelines.
  Adjusted Score and scaled adjusted likelihood-ratio tests remain
  calibration-pending for `joint_ln`; legacy LN is retired from `scrope()`.
- Implemented contrast-level testing and expanded robust diagnostics.
- Improved HL estimation and stability, including analytic Hessians, cached
  scores, and negative-LR safeguards.
- Added an exact PGMM working-likelihood route with robust Wald, adjusted Score,
  and adjusted likelihood-ratio inference.
- Added exact PGMM profile-likelihood diagnostics, including constrained
  profiling, warm-started profile grids, raw-LR safeguards, and
  information-scaled convergence checks.
