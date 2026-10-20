# Installed-package software release gates

These small checks exercise package execution and public contracts. They are
not calibration studies and do not support type-I error, coverage, power, FDR,
or method-comparison claims. Paper simulation and HLCA workflows are separate.

`run_workable_version_smoke.R` fits two bundled genes through default JLN,
explicit JLN, and PGMM, and one through first-order HL. It checks route-specific
provenance schemas (JLN version 3; HL and PGMM version 2), agreement between the
default and explicit JLN calls, test availability, complete saved-result
round-tripping, and elapsed times. Five error checks verify retirement of the
legacy estimator, legacy fallback, default `nebula()`, explicit `nebula()` LN,
and `nebula()` PMM with LN. An unrelated input error does not satisfy them.

`run_installed_parallel_gate.R` requires at least two cores under normal
`parallelly::availableCores()` resource detection. It starts two distinct PSOCK
processes and checks their installed namespace and compiled DLL identity. It
then runs exactly these installed package regressions:

- `installed joint-LN PSOCK fits use the running package`
- `joint LN single-core and multicore results agree`

Both tests must execute and pass. Missing installed tests, insufficient cores,
a skipped test, a warning, or an error makes this gate fail. Do not override
resource detection to make it pass; use a host with the required allocation.
The tests use deterministic one-gene fixtures and request `ncore = 2L` for the
parallel fits. The saved receipt records resources, worker PIDs, DLL checksums,
test outcomes, and the R session.

Both scripts load the explicit installed package path and verify that its DLL
is loaded from that installation. Source trees and source-loading fallbacks
cannot count as compiled-package validation. Use a fresh R process for each
command after installing the exact source archive with its tests:

```sh
mkdir -p /tmp/scrope-rc-library
R CMD INSTALL --install-tests --library=/tmp/scrope-rc-library scRoPE_0.3.0.tar.gz
Rscript /tmp/scrope-rc-library/scRoPE/tests/validation/run_workable_version_smoke.R \
  /tmp/scrope-workable-smoke.rds /tmp/scrope-rc-library/scRoPE
Rscript /tmp/scrope-rc-library/scRoPE/tests/validation/run_installed_parallel_gate.R \
  /tmp/scrope-parallel-gate.rds /tmp/scrope-rc-library/scRoPE
```

The Linux GitHub Actions workflow builds one source archive, runs `R CMD check`
with tests enabled, and installs that same archive into a fresh library for
these gates. Declared dependencies, including Suggests, are required; check
errors and warnings fail the workflow. No package publication is performed.
Passing this workflow covers Linux with the current R release; it does not
establish macOS, Windows, older R, or scientific calibration results.
