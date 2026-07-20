# First-order HL Wald/Score convergence validation

This directory contains the standalone common-random-number Monte Carlo gate
for the corrected first-order nested HL pipeline. It is intentionally separate
from ordinary `testthat` tests because the acceptance run fits 250 replicates at
each of four subject counts.

The locked regular generator uses five cells per subject,
`beta = (-1.2, 0.35, -0.25)`, `sigma2 = 0.4`, `phi = 8`, and tests
`R = (0, 1, 0)` against `rhs = 0.35`. Each replicate is generated once at
`m = 2000`; smaller cases use the first `m` subjects and their cells, giving
common random numbers across `m = 250, 500, 1000, 2000`.

Both dispersions are fixed at their generating truth (`q = 0`). This is an
explicitly supported pipeline edge case. It isolates the required first-order
equivalence of robust Wald, based on the unrestricted affine engine and its
augmented sandwich, and Score, based on the separate constrained full-chi null
root and adjusted subject meat. The choice is recorded in every artifact's
design metadata; it is not a dummy outer fit.

Run the two-replicate benchmark/smoke from the package root:

```sh
Rscript tests/validation/run_hl_first_order_wald_score_mc.R \
  --mode=smoke --replicates=2 --workers=2
```

The runner is resumable. It saves one complete generated `m = 2000` array and
one result object per replicate, then rebuilds combined CSV/RDS summaries from
those objects. Replicate RNG streams are derived deterministically from the
locked L'Ecuyer-CMRG master stream and do not depend on worker count or
scheduling. Re-running the command reuses artifacts only when their generated
file checksum and locked design match.

The runner first loads the exact worktree through `pkgload`. If build tools are
temporarily unavailable, it records that reason and uses the solver's exact R
reference path from the corrected first-order source files. This is slower but
mathematically identical; the manifest records which loader was used.

After focused mathematical and solver gates are green and the full run is
explicitly approved, use:

```sh
Rscript tests/validation/run_hl_first_order_wald_score_mc.R \
  --mode=full --replicates=250 --workers=8 \
  --confirm-full=RUN_LOCKED_250
```

Full mode is guarded: it requires exactly 250 replicates and the confirmation
token. Smoke mode cannot run more than two replicates.

Artifacts are written beneath the ignored directory
`temp/hl_first_order_validation/`:

- `design_metadata.rds`: locked generator, inference, `q = 0`, and threshold
  metadata;
- `generated/replicate_NNN_m2000.rds`: all maximum-size observed and latent
  arrays used for subject subsetting;
- `per_replicate/replicate_NNN_results.rds`: fit/covariance/Wald/Score statuses,
  statistics, p-values, ranks, discrepancies, and timings for all four sizes;
- `fit_results.{rds,csv}`: combined per-fit records;
- `wald_score_summary.{rds,csv}`: success counts/rates, `D_m`, and the fitted
  log-log slope; and
- `run_manifest_{smoke,full}.rds`: runtime, resume counts, package version, and
  session information.

For successful pairs the discrepancy is
`abs(W - Q_EE) / (1 + W + Q_EE)`, and `D_m` is its median. The full gate is
evaluated only if all 250 Wald/Score pairs are successful at every size. It
requires `D_2000 <= 0.05` and the least-squares slope of `log(D_m)` on `log(m)`
to be at most `-0.25`. Smoke output is always labeled
`smoke_only_not_acceptance`.
