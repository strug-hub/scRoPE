# scRoPE

scRoPE provides robust inference for multi-donor single-cell RNA-seq
association analyses with mixed models. It extends the `nebula` framework with
Godambe-based robust Wald inference, adjusted Score tests, and gated adjusted
likelihood-ratio tests for both approximate NBGMM fits and the exact PGMM path.

## Installation

```r
# latest development version
install.packages("devtools")
devtools::install_github("strug-hub/scRoPE")
```

## Fitting Paths

scRoPE currently exposes three package-facing fitting paths:

- `scrope()`
  Default LN-based NBGMM pipeline, with optional adjusted Score and adjusted
  LRT inference.
- `scrope_hl()`
  HL/APHL pipeline for the same NBGMM mean/variance structure.
- `scrope(model = "PGMM")`
  Exact Poisson-gamma mixed model path inherited from legacy
  `nebula(model = "PMM")`.

## Quick Start

```r
library(scRoPE)
data(sample_data)
design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
fit <- scrope(sample_data$count, sample_data$sid, pred = design, ncore = 1,
              additional_tests = c("score", "lrt"), lrt_details = TRUE)
head(fit$summary)
```

## Exact PGMM Path

`scrope(model = "PGMM")` fits the exact Poisson-gamma mixed model inherited
from legacy `nebula(model = "PMM")`. This path supports the exact PGMM
likelihood together with robust Wald, adjusted Score, and gated adjusted LRT
inference.

```r
library(scRoPE)
data(sample_data)
design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
fit_pmg <- scrope(sample_data$count, sample_data$sid, pred = design,
                  model = "PGMM", ncore = 1)
head(fit_pmg$summary)
```

Important caveat:

- `PGMM` allows cell-level fixed predictors, matching legacy `PMM`.
- `PGMM` only models subject-level overdispersion and does not include a
  cell-level random-effect / overdispersion term.
- Because of that, `PGMM` should not be the default choice for testing a
  cell-level predictor when cell-level overdispersion matters. Use the LN/HL
  NBGMM path for that setting.

## Choosing A Path

- Use `scrope()` as the default package entry point for multi-donor scRNA-seq
  association testing when both subject-level and cell-level overdispersion
  matter.
- Use `scrope_hl()` when you want the HL/APHL fitting pipeline explicitly.
- Use `scrope(model = "PGMM")` when you want the exact Poisson-gamma marginal
  likelihood and only a subject-level overdispersion component is appropriate.

## Attribution

scRoPE is derived from the nebula package (GPL-2). Upstream source:
https://github.com/lhe17/nebula

High-level changes in this repo include:
- Added Godambe-based adjusted Wald, Score, and scaled adjusted
  likelihood-ratio tests for both LN and HL pipelines.
- Implemented contrast-level testing and expanded robust diagnostics for the
  adjusted tests.
- Improved HL estimation and stability (analytic Hessians, cached scores, and
  negative LRT safeguards) for reliable inference.
