# Reproducing the synthetic teaching examples

These examples were designed for this unpublished local 0.3.0 preparation.
They contain no NEBULA/pancreas/HLCA observations or paper simulation results.
Parameters were committed before counts or fitted outputs were examined.
Seed 20261020 is fixed; there was no seed search or gene selection.

## Regenerate from source

From the package source root, with the intended scRoPE build and existing
Matrix, Seurat, ggplot2, svglite and ragg installations available:

    Rscript data-raw/rebuild_examples.R
    Rscript data-raw/fit_readme_examples.R
    Rscript data-raw/render_readme.R

The first command recreates both teaching data files and the generating
coefficient table. The second fits all eight genes in each scenario and writes
the complete illustrative result CSV. An optional path argument saves a full
fit/session receipt outside the source tree. The third renders SVG/PNG figures
from that CSV; it does not fit models. These scripts install nothing. The
source-only data-raw directory is excluded from the built package.

Installed users can generate the same inputs without fitting:

    e <- new.env()
    sys.source(system.file("examples", "synthetic_parameters.R", package="scRoPE"), e)
    sys.source(system.file("examples", "generate_synthetic.R", package="scRoPE"), e)
    examples <- e$generate_synthetic_examples()
    sample_data <- examples$reference
    sample_seurat <- e$make_synthetic_seurat(sample_data)

The generated list also holds the paired misspecified scenario, coefficients,
conditional means and donor effects. The generator restores RNGkind and
.Random.seed; a caller's hidden Box-Muller cache cannot be restored by base R.
The declared Inversion generator is unaffected by that limitation.

## Recorded environment and interpretation

Initial generation: R 4.5.2, Matrix 1.7.4, Seurat 5.3.0, SeuratObject 5.2.0.
Plots: ggplot2 4.0.0, svglite 2.1.3, ragg 1.4.0.
Initial fits: scRoPE 0.3.0 local source commit
62b05bc8c5ccd93e499d93dd5391a5bf3e06243d. The synthetic-example changes do not
alter the JLN fitter. Exact results can depend on numerical libraries and R
versions; regenerate with the intended version rather than trusting a version
label alone. Tests check data regeneration, conversion, and public inference
algebra. A Seurat serialization can vary across dependency versions even when
its counts and metadata agree.

All 16 fits supplied inference in the fixed run. The complete CSV retains the
availability column and all genes. The figure scripts show unavailable fits
if encountered in a future explicitly documented regeneration; do not change
seeds or silently delete them. Robust/model SE ratios are descriptive, not a
ranking of accuracy. Both intervals share one estimate, and the sandwich does
not correct bias. No repeated simulation was run, so neither calibration nor
power improvement is established. Historical examples remain in historical
sources; changing maintained fixtures does not alter historical results or
resolve the inherited code-license conflict.
