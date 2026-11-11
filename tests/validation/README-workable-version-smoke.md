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
- `PGMM two-worker observation preserves optional output and gene RNG streams`

All three tests must execute and pass. Missing installed tests, insufficient cores,
a skipped test, a warning, or an error makes this gate fail. Do not override
resource detection to make it pass; use a host with the required allocation.
The LN tests use deterministic one-gene fixtures and request `ncore = 2L` for the
parallel fits. The PGMM test uses two retained genes, requires two distinct
gene-executing worker PIDs, and checks exact serial/parallel optional-test output,
provenance, observer records (apart from PIDs), and gene RNG draws. The gate sets
`SCROPE_RUN_AUTO_PGMM_PARALLEL=true` and `NOT_CRAN=true` for this required run;
the PGMM test remains opt-in during ordinary package checks.
The saved receipt records resources, probe-worker PIDs, DLL checksums,
test outcomes, R/package versions and platform, without installation paths or
host/user names. All three regressions are attempted before reporting an aggregate
failure, so a failure in one does not hide the others' outcomes. The PGMM result
row records whether its assertions passed; its gene-worker PIDs are asserted
inside that test, rather than saved in the gate receipt.

Each LN regression saves its original serial/parallel synthetic fit pair before
assertions, alongside the gate receipt. Pace diagnostic witness indices must
point to the first certified, interior candidate matching their own endpoint's
exact coordinates, provenance and existing objective-agreement threshold.
Those positions are local to a fit's candidate list: only validated positions
are normalized in comparison copies. Original indices and candidate records
remain in the receipts; numerical tolerances, runtime keys, statuses and all
other integer comparisons are unchanged.

Both scripts load the explicit installed package path and verify that its DLL
is loaded from that installation. Source trees and source-loading fallbacks
cannot count as compiled-package validation. Use a fresh R process for each
command after installing the exact source archive with its tests:

```sh
mkdir -p /tmp/scrope-rc-library
R CMD INSTALL --install-tests --library=/tmp/scrope-rc-library scRoPE_0.4.1.tar.gz
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
