# 24-gene omitted-slope example

Load `data(toy_data)` in scRoPE 0.3.0 for a 24 by 8,000 integer count matrix
and its complete metadata. This is a teaching dataset, not paper/HLCA output.
The earlier `sample_data` / `sample_seurat` objects remain separate.

## Run the installed example

In a fresh R session with the intended local package installed:

```r
library(scRoPE)
source(system.file("examples", "toy24", "quick_start.R", package = "scRoPE"))
results             # all genes, both intervals, truth, availability and reasons
fit$availability
fit_time
```

The script uses one worker, JLN and no estimator fallback, with `cpc = 0` and
`mincp = 1` so all planned genes are analyzed. It writes no files, installs
nothing, and exposes the complete `fit`, `results` and `fit_time` in the caller.
An entire-call error remains an R error; individual unavailable genes are
retained by the result extractor, with missing uncertainty and explicit reasons.

To regenerate the data without fitting or writing files:

```r
e <- new.env()
sys.source(system.file("examples", "toy24", "generate.R", package = "scRoPE"), e)
sys.source(system.file("examples", "verify_regeneration.R", package = "scRoPE"), e)
regenerated <- e$generate_toy24()
data(toy_data)
e$verify_toy24_regeneration(regenerated, toy_data)
```

The verifier requires exact counts, design metadata, identifiers, RNG settings
and stored streams, including their types and attributes. Only the continuous
donor intercept multipliers and slopes allow last-bit numerical differences:
each absolute difference must be at most
`64 * .Machine$double.eps * max(1, abs(actual), abs(expected))`.
Continuous values can differ at that scale across R versions and numerical
environments even when counts and streams agree. The bundled object retains
the reviewed reference unchanged; same-environment replay remains exact.

The generator restores `RNGkind()` and `.Random.seed`. Base R cannot restore
a caller's hidden Box–Muller normal cache; the declared Inversion generator
is unaffected. No random seed needs to be chosen by the user.

## Generating model and interpretation

200 independent donors; 40 cells per donor, balanced at `x = -0.5` and `0.5`;
unit multiplicative offsets. Each gene has independent latent donor draws:

```text
U ~ Gamma(shape = rate = 1 / 0.35)
b ~ Normal(mean = 0, variance = tau_squared)
Y | U, b ~ Negative Binomial(size = 4, mean = mu)
mu = a * U * exp(beta * x + b * x - tau_squared * x^2 / 2)
```

Reference genes have `tau_squared = 0`; omitted-slope genes have
`tau_squared = 1`. In each scenario, six genes have `beta = 0`, three
`beta = log(1.5)` and three `beta = -log(1.5)`. Each six-gene category has
three genes at baseline `a = 1` and three at `a = 2`. Sign and baseline are
not a fully balanced factorial; genes do not share paired latent draws.

Because `E(U) = 1` and `E(exp(b*x - tau_squared*x^2/2)) = 1`,
`E(Y | x) = a * exp(beta*x)`. Thus `beta` is exactly the generating population
log-mean contrast between the two conditions. Baseline `a` is the mean at
`x = 0`, not the pooled non-null mean. A particular donor's conditional ratio
is `exp(beta + b)`; the finite sampled donors need not reproduce the population
contrast exactly.

The fitted model retains a donor random intercept and a common fixed x
coefficient, but omits donor-specific slopes. Both uncertainty calculations
use the same fitted coefficient. Equality of the approximate JLN population
target and the generating non-null contrast under omission is unproved.
Comparing fitted intervals with that reference does not certify coverage or
bias correction. Twenty-four genes do not establish calibration, power or FDR.
Historical .9004 null simulations used 400 cells per donor; those separate
results do not calibrate this 40-cell 0.3.0 example and are not bundled here.

## Fixed design and provenance

The master seed is the fixed date-based `20261021`, using L'Ecuyer-CMRG,
Inversion and Rejection. Gene k uses stream k, advanced with
`parallel::nextRNGStream()` from the preceding initial state. Draw order is
200 Gamma intercept multipliers, 200 Normal slopes (including zero SD), then
8,000 negative-binomial counts. Gene metadata and all latent draws/streams
are retained in `toy_data`.

Eight genes were first reviewed separately. After reviewing those outcomes,
16 appended genes and unused streams 9–24 were fixed before their generation
or fitting on 21 October 2026. Original rows/metadata/latent draws remain exact.
There was no favorable-seed search, significance-based selection or
post-outcome parameter tuning. The complete reviewed object has RDS SHA-256
`a4eef333bf34f6dd8a3df6e91385b9991365db4c360fc379664d734f63b12be8`;
that is the input RDS identifier, not the differently serialized package RDA.

`reference_results.csv` retains all 24 reviewed outputs from scRoPE 0.3.0
source commit `55899156dedb14841c79ae0d36f56be8ecb17c0a` (R 4.5.2, macOS,
Accelerate BLAS). Integration changes data/docs/tests only, not the fitter.
Exact floating-point results can vary with platform and numerical libraries.
A fit comparison is a software reproducibility check, not new calibration.

All 24 reviewed fits supplied inference without warnings. Model/robust
interval inclusion was 6/6 vs 6/6 in each reference category, 3/6 vs 6/6 for
omitted null genes, and 2/6 vs 6/6 for omitted non-null genes. Three omitted
null genes reject zero with model-based tests, none with robust tests.
All 12 non-null genes reject zero with both; zero testing is distinct from
including their nonzero generating contrast. All outcomes remain visible.

## Save a table or redraw the figure

After running `quick_start.R`, save `results` to a new output path:

```r
write.csv(results, "toy-results.csv", row.names = FALSE, na = "NA")
system.file("examples", "toy24", "render.R", package = "scRoPE")
```

Use the printed script path in a terminal, with existing ggplot2, svglite and
ragg packages (not needed for generating or fitting):

```sh
Rscript /path/to/render.R toy-results.csv NEW_FIGURE_DIRECTORY
```

The renderer reads saved results only, shows all 24 genes on common axes,
retains unavailable genes, and refuses to overwrite a figure directory.
To redraw the reviewed figure without fitting, use the CSV returned by
`system.file("examples", "toy24", "reference_results.csv", package="scRoPE")`.
For source maintainers, `Rscript data-raw/rebuild_toy24.R` regenerates only
`data/toy_data.rda`; it never refits or changes the reference result table.

## Recorded runtime

Benchmark configuration: M3 Max MacBook, R 4.5.2 / scRoPE 0.3.0.
0.047 seconds to generate; 10.964 seconds for the first fit; 11.957 seconds
for the complete one-fit R process including startup and exports; 10.409
seconds for an identical-input warm fit. Final PNG/SVG rendering took 0.988
seconds, excluding R startup. These are individual local measurements, not
cold-cache benchmarks or estimates for typical laptops.

The measured baseline used `ncore = 1`, one foreach worker, Accelerate BLAS,
and a BLAS thread query of one. Numerical thread environment caps were one;
OpenMP thread count was unavailable. The simple R quick start uses one fit
worker but does not modify the user's native-thread settings. Hardware
core/RAM probing was unavailable, so no independent hardware specification
or universal speed claim is made. A second first-fit timing briefly overlapped
plotting and was not used as the isolated baseline.
