# scRoPE

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
estimator; they are currently unavailable for the opt-in joint-LN estimator.

## Installation

```r
# latest development version
install.packages("remotes")
remotes::install_github("strug-hub/scRoPE")
```

## Quick start

The default `scrope()` call uses the scalable LN approximation to the
negative-binomial gamma mixed-model (NBGMM) working likelihood. The design can
contain both subject-level and cell-level predictors. The released
`legacy_post_pml` estimator remains the default.

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

```text
  gene     logFC_X1 se_robust_X1 p_robust_X1
1    A -0.016755225   0.03156698   0.5955689
2    B -0.002679074   0.03738061   0.9428644
3    C  0.017954707   0.03697886   0.6272935
4    D -0.005698984   0.02849276   0.8414687
5    E  0.011557090   0.03633498   0.7504314
6    F  0.013483039   0.03827191   0.7246163
```

Here, `logFC_X1` is the estimated coefficient, `se_robust_X1` is its robust
standard error, and `p_robust_X1` is the robust Wald p-value. The complete
`fit$summary` also contains model-based results and results for the other
predictors.

## Opt-in joint-LN estimation

The public default remains unchanged. Users who want the estimator-aligned
joint-LN fit can request it explicitly:

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

In the 0.2.5.9000 calibration-candidate build, a narrowly bounded review can
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

### Optional legacy fallback

For operational sensitivity analyses, sc-RoPE can try the released legacy
post-PML method when a joint-LN fit cannot provide inference because a
dispersion estimate reached its allowed boundary:

```r
fit_joint_fallback <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  ncore = 1,
  ln_estimator = "joint_ln",
  ln_operational_fallback = "legacy_post_pml"
)

fallback <- fit_joint_fallback$operational_fallback

fit_joint_fallback$availability

fallback$per_gene[, c(
  "gene",
  "fallback_status",
  "fallback_route_class",
  "fallback_boundary"
)]

subset(fallback$coefficients, coefficient_available %in% TRUE)
fallback$overdispersion
```

The bundled quick-start data do not need fallback, so the filtered fallback
table in this example is empty.

The default is `ln_operational_fallback = "none"`, which performs no fallback
fit and returns no separate fallback results. Joint-LN emits one consolidated
availability warning when one or more genes lack certified primary inference;
it does not emit one availability warning per contrast. The
`"legacy_post_pml"` option is available only for `model = "LN"` with
`ln_estimator = "joint_ln"` and `allow_per_gene_switch = FALSE`.

Each fallback attempt is a new, independent legacy post-PML fit. Its
coefficient and dispersions come from that same fit, while the primary joint-LN
result remains unchanged. The sidecar deliberately contains no covariance,
standard error, test statistic, p-value, q-value, or inferential-success field;
`inference_supported` is always `FALSE`.

For routine use, interpret the three values in `fit$availability` as follows:

- **`inference`:** use the certified primary `summary` result.
- **`estimate_only`:** an explicitly requested fallback has
  `fallback_status = "coefficient_available"`; its separate coefficient can be
  reported operationally, without uncertainty or hypothesis testing.
- **`unavailable`:** neither result is available; do not assign a p-value.

This precedence is fixed: certified primary inference, then a successful
explicitly requested coefficient-only fallback, then unavailable. A finite
diagnostic endpoint from an uncertified joint-LN fit is never promoted to
`estimate_only`. The `inference_status`, `resolution_status`,
`inference_source`, `reason_group`, `status_message`, and `technical_code`
columns provide a short explanation while preserving the complete expert
records in `provenance`, `diagnostics`, and `operational_fallback`. The table
also records `estimator`, `certified_root`, `strict_inference_status`,
`enabled_inference_status`, `newly_enabled`, `availability_transition`,
`resolution_policy_version`, `gene_id`, `gene`, `inference_supported`,
`estimate_source`, and `attention_flag`.

The remaining `fallback_status` values distinguish a failed fallback attempt
from a gene that did not need or qualify for fallback. They are diagnostic and
do not provide inference.

An available fallback coefficient is still a point estimate from a different,
potentially rescued or boundary, estimator route. Inspect
`fallback_route_class`, `fallback_boundary`, and the distance-to-boundary
diagnostics before using it descriptively. It cannot be combined with primary
p-values in an FDR procedure because no fallback p-value is produced.

## Optional Adjusted Score and LRT Inference

The example below uses the default estimator. Adjusted Score and LRT tests are
unavailable for `joint_ln` pending calibration.

Request the adjusted Score and scaled adjusted likelihood-ratio tests through
`additional_tests`:

```r
fit_tests <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
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
scaled adjusted LRT can also be more sensitive to approximation error in the LN
or HL objective.

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

Test availability remains route-specific. Legacy LN and PGMM retain their
existing optional tests; certified joint-LN exposes model-based and donor-robust
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

Adjusted Score and likelihood-ratio tests can be added to either PGMM example
using the same `additional_tests` and `lrt_details` arguments shown above.

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
- `operational_fallback`: separate inference-free legacy coefficients,
  dispersions, and diagnostics, returned only when
  `ln_operational_fallback = "legacy_post_pml"`; and
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
  Adjusted Score and scaled adjusted likelihood-ratio tests remain available
  on retained legacy routes; they are calibration-pending for `joint_ln`.
- Implemented contrast-level testing and expanded robust diagnostics.
- Improved HL estimation and stability, including analytic Hessians, cached
  scores, and negative-LR safeguards.
- Added an exact PGMM working-likelihood route with robust Wald, adjusted Score,
  and adjusted likelihood-ratio inference.
- Added exact PGMM profile-likelihood diagnostics, including constrained
  profiling, warm-started profile grids, raw-LR safeguards, and
  information-scaled convergence checks.
