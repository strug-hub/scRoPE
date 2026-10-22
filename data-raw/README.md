# Teaching-example source scripts

The main 24-gene quick start is documented in
[../inst/examples/toy24/README.md](../inst/examples/toy24/README.md).
`rebuild_toy24.R` regenerates its data without fitting. The scripts below
continue to maintain the separate original synthetic/Seurat examples.

Read SYNTHETIC_DESIGN.md before regenerating. Its parameters and seed were
frozen in a separate commit before outcomes were inspected. The pure generator
and parameters are installed under inst/examples so users can inspect and
regenerate them. The wrappers here are source-only development utilities.

Run rebuild_examples.R, fit_readme_examples.R, then render_readme.R from the
source root. Full instructions and dependency versions are in
inst/examples/README.md. No simulation campaigns, HLCA workflows, dependency
installations or external jobs are invoked.
