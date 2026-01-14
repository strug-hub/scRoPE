# scRoPE

Robust inference for multi-donor single-cell RNA-seq association analyses using
negative binomial mixed models.

## Installation

```r
# latest development version
install.packages("devtools")
devtools::install_github("strug-hub/scRoPE")
```

## Example

```r
library(scRoPE)
data(sample_data)
design <- model.matrix(~ X1 + X2 + cc, data = sample_data$pred)
fit <- scrope(sample_data$count, sample_data$sid, pred = design, ncore = 1,
              additional_tests = c("score", "lrt"), lrt_details = TRUE)
head(fit$summary)
```

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
