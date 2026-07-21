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
adjusted likelihood-ratio tests are route-dependent optional tests; the opt-in
joint-LN route currently returns them as unavailable pending calibration.

## Installation

```r
# latest development version
install.packages("remotes")
remotes::install_github("strug-hub/scRoPE")
```

## Quick Start: Current Default LN Route

The default `scrope()` call uses the scalable LN approximation to the
negative-binomial gamma mixed-model (NBGMM) working likelihood. The design can
contain both subject-level and cell-level predictors. The current default is
`ln_estimator = "legacy_post_pml"` for legacy output reproduction.

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

## Opt-in Joint-LN Estimator

The estimator-aligned LN route is available explicitly while the public
default remains unchanged:

```r
fit_joint <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  ncore = 1,
  ln_estimator = "joint_ln"
)

fit_joint$provenance$per_gene[, c(
  "gene",
  "route",
  "start_provenance",
  "inference_available",
  "failure_reason",
  "mode_status"
)]
```

For genes classified `R0`, `joint_ln` reports a stationary local maximizer of
the NEBULA-LN approximate marginal criterion that passed the declared root,
curvature, derivative, and start-agreement gates. This is not a claim of global
uniqueness. Model covariance is evaluated from the full joint-LN curvature,
and robust covariance is evaluated from the full donor-level sandwich at the
same root.

sc-RoPE 0.2.3 operationally reported `repml$beta` after the raw-LN fit. The
NEBULA paper describes its LN approximate-marginal-likelihood M-estimator and a
separate h-likelihood step for standard-error computation, but does not clearly
state this software-level coefficient replacement. Accordingly, `joint_ln` is
defined by the NEBULA-LN approximate marginal criterion; it is not claimed to
reproduce legacy NEBULA coefficients.

The structured `provenance` object distinguishes:

- `R0`: regular joint-LN root with model and robust Wald inference available;
- `R1`: legacy post-PML output reproduction, without a claim of
  estimator-aligned joint-LN inference; and
- `R2`: unsupported analytic joint-LN inference. Covariance, standard errors,
  Wald statistics, and p-values are unavailable rather than filled by a
  fallback statistic.

With `output_re = TRUE`, joint-LN donor modes are auxiliary fixed-beta
conditional PML predictions evaluated at the plug-in joint-LN parameters. They
cannot change coefficients, dispersions, covariance, tests, convergence, or
route classification.

Adjusted Score and scaled adjusted likelihood-ratio tests remain unavailable
for `joint_ln`. Requests through `additional_tests` return unavailable
statistics with status `joint_ln_calibration_pending`; same-objective
constrained fits are retained for diagnostics, not as calibration evidence.
Switching from the legacy route to `joint_ln` may change coefficients, standard
errors, p-values, convergence statuses, FDR calls, and downstream results.

## Optional Adjusted Score and LRT Inference

The example below uses the current `legacy_post_pml` default. It does not
describe `joint_ln`, for which adjusted Score and LRT results remain
unavailable with status `joint_ln_calibration_pending`.

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
| NBGMM, HL/APHL approximation | `scrope_hl()` | Donor heterogeneity and cell-level extra-Poisson variation | HL fitting or comparison with the LN approximation |
| PGMM, exact marginal likelihood | `scrope(model = "PGMM")` | Donor heterogeneity without a separate cell-level overdispersion parameter | Simpler exact working likelihood for fixed-effect inference |

The NBGMM working likelihood is the more complete variance model. PGMM is a
targeted simplification: it retains the donor-level gamma random effect but
omits the cell-level gamma overdispersion component.

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
- optional, route-dependent `score_adj_*` and `p_score_adj_*`: adjusted Score
  results; and
- optional, route-dependent `lrt_adj_*`, `p_lrt_adj_*`, and `lrt_method_*`:
  scaled adjusted likelihood-ratio results and method diagnostics.

The returned object also contains:

- `overdispersion`: fitted working-model overdispersion parameters;
- `convergence` and `algorithm`: per-gene fitting status;
- `provenance`: estimator definition, parameter maps, and per-gene route,
  start, inference-availability, failure, and optional-mode status;
- `diagnostics`: detailed per-gene fitting and robust-inference information;
- `lrt_details`: optional tidy likelihood-ratio diagnostics; and
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
