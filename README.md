# scRoPE

scRoPE provides robust inference for multi-donor single-cell RNA-seq
association analyses using mixed models. It extends the `nebula` framework
with Godambe-based robust Wald inference, adjusted Score tests, and adjusted
likelihood-ratio tests for negative-binomial and Poisson-gamma mixed-model
working likelihoods.

The package supports both scalable approximate NBGMM fitting paths and an
exact PGMM likelihood path. The exact PGMM path is also used for
profile-likelihood diagnostics and likelihood-paradigm analyses, including
relative ordinary and robust profile likelihoods for scalar fixed-effect
targets.

## Installation

```r
# latest development version
install.packages("remotes")
remotes::install_github("strug-hub/scRoPE")
```

## Fitting Paths

scRoPE currently exposes three main fitting paths:

### `scrope()`

Default LN-based NBGMM pipeline. This path uses the large-number approximation
from `nebula` and adds Godambe-based robust inference.

### `scrope_hl()`

HL/APHL-based NBGMM pipeline. This path uses the hierarchical-likelihood
fitting strategy and adds the same robust inference framework.

### `scrope(model = "PGMM")`

Exact Poisson-gamma mixed model path. This path fits the exact PGMM marginal
likelihood and supports robust Wald, adjusted Score, and adjusted
likelihood-ratio inference. It also supports exact PGMM profile-likelihood
diagnostics for scalar fixed-effect targets.

## Quick Start

```r
library(scRoPE)

data(sample_data)

design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

fit <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  ncore = 1,
  additional_tests = c("score", "lrt"),
  lrt_details = TRUE
)

head(fit$summary)
```

## Exact PGMM Path

`scrope(model = "PGMM")` fits an exact Poisson-gamma mixed model. This path is
inherited from the legacy `nebula(model = "PMM")` model class, but has been
extended in scRoPE with robust inference and profile-likelihood diagnostics.

```r
library(scRoPE)

data(sample_data)

design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)

fit_pmg <- scrope(
  sample_data$count,
  sample_data$sid,
  pred = design,
  model = "PGMM",
  ncore = 1,
  additional_tests = c("score", "lrt"),
  lrt_details = TRUE
)

head(fit_pmg$summary)
```

The PGMM working likelihood includes a subject-level gamma random effect and
uses the exact marginal likelihood obtained after integrating out the
subject-level random effect. It can be used for both subject-level and
cell-level fixed-effect targets.

A key purpose of the robust inference framework is to protect fixed-effect
inference when the working variance structure is imperfect. Therefore, the
PGMM path should not be interpreted as requiring the complete variance model to
be correct. In particular, even though the PGMM working likelihood does not
include a separate cell-level overdispersion parameter, the Godambe-based
robust adjustment is designed to account for variance misspecification when
estimating uncertainty for fixed-effect targets.

The PGMM path is especially useful when:

- an exact marginal likelihood is desired;
- the analysis focuses on likelihood-based or profile-likelihood diagnostics;
- robust inference is needed for subject-level or cell-level fixed-effect
  targets under a simpler working mixed model;
- simulation studies require an exact likelihood benchmark.

## Profile Likelihood Diagnostics

The exact PGMM path includes internal helpers for ordinary and robust relative
profile likelihoods for scalar fixed-effect targets. These tools are useful for
likelihood-paradigm analyses, including simulations of the probability of
misleading evidence and Royall-style bump-function diagnostics.

The profile-likelihood infrastructure includes:

- equality-constrained PGMM fitting for scalar contrasts;
- ordinary relative profile likelihoods;
- robust adjusted relative profile likelihoods;
- pointwise and grid-based profiling;
- warm-started constrained profile grids;
- raw likelihood-ratio safeguards;
- optimizer diagnostics and retry logic;
- information-scaled convergence checks for high-information constrained
  profile fits.

These diagnostics are primarily intended for methodological development and
simulation studies. For routine differential expression or association testing,
the summary-level Wald, Score, and adjusted likelihood-ratio outputs are the
main user-facing results.

## Choosing a Path

- Use `scrope()` as the default scalable NBGMM path for multi-donor scRNA-seq
  association testing when both subject-level and cell-level overdispersion are
  explicitly modeled.
- Use `scrope_hl()` when the HL/APHL fitting pipeline is desired explicitly,
  for example when the LN approximation is insufficient or when comparing
  fitting strategies.
- Use `scrope(model = "PGMM")` when an exact Poisson-gamma marginal likelihood
  is desired, including exact-likelihood simulations, likelihood-paradigm
  analyses, and robust inference under a simpler working mixed model.

The LN/HL NBGMM paths and the exact PGMM path serve different purposes. The
NBGMM paths model both subject-level and cell-level overdispersion explicitly.
The PGMM path uses a simpler exact likelihood and relies on robust inference to
protect fixed-effect uncertainty under variance misspecification.

## Attribution

scRoPE is derived from the `nebula` package (GPL-2).

Upstream source:

https://github.com/lhe17/nebula

High-level changes in this repository include:

- Added Godambe-based adjusted Wald, Score, and scaled adjusted
  likelihood-ratio tests for LN and HL NBGMM pipelines.
- Implemented contrast-level testing and expanded robust diagnostics.
- Improved HL estimation and stability, including analytic Hessians, cached
  scores, and negative-LR safeguards.
- Added an exact PGMM fitting path with robust Wald, adjusted Score, and
  adjusted likelihood-ratio inference.
- Added exact PGMM profile-likelihood diagnostics, including constrained
  profiling, warm-started profile grids, raw-LR safeguards, and
  information-scaled convergence checks.
