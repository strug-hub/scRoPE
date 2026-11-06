# Automatic recovery, availability and retained evidence

This candidate uses `recovery="auto"` and `keep_diagnostics=FALSE` by default
for public `scrope()` calls. `scrope_hl()` has the same defaults. These are
intentional changes from the prior HL/PGMM API and prior retention default.

```r
fit <- scrope(count, id, pred, offset=exposure, model="PGMM")
native <- scrope(count, id, pred, offset=exposure, model="PGMM", recovery="none")
detailed <- scrope(count, id, pred, offset=exposure, model="PGMM",
                   keep_diagnostics=TRUE)
```

The requested estimator is retained. Recovery does not switch between LN, HL
and PGMM. Accepted native results are preserved under each method's reviewed
qualification rule. Otherwise an eligible numerical refit, positive fixed
dispersion working fit, or identifiable whole-zero-donor limit is assessed.
Passing the numerical gates does not establish calibration after selection.

## Compatibility and controls

An explicit LN `ln_recovery="none"` remains a native-only request when the new
`recovery` argument is omitted. Both arguments may be supplied if they agree;
conflicting explicit LN values are an error. `ln_recovery` remains ignored by
HL/PGMM. The separate `hl_score_recovery="stable"` option acts only on requested
native HL Score tests. It remains usable with `recovery="none"`.

Working recovery requires the reviewed default dispersion box. A custom box
retains the native result and reports `recovery_requires_default_box` (the LN
reason also names its L-BFGS-B requirement). Native HL/PGMM already ignore
`opt="trust"` with a warning; their actual solver is unchanged. A custom
`pgmm_control` applies unchanged to the main native optimizer. It is not a total
work limit: recovery uses its separately recorded, reviewed refit/audit policy.

The HL/PGMM fixed donor parameter is `sigma2=1e-4`. For eligible HL cell boundaries,
the fixed NB precision is `phi=1000` (inverse precision 0.001). These are positive
working parameters, not a zero-variance reference distribution.

## Reading results

`summary` retains the familiar coefficient columns. The common compact tables are:

| Field | Meaning |
| --- | --- |
| `availability` | One row per analyzed gene, selected route/source, separate point/covariance/Wald counts, recovery eligibility, target, reason, and native convergence code. |
| `coefficient_availability` | One row per coefficient, estimates, estimability, separate point/covariance/statistic/reference/interval availability, and exact failure reason. |
| `test_status` | One row per requested/ordinary reported test; `availability_status` is the common status while method-specific `status` remains intact. |
| `recovery$primary` | Native-only public result before whole-fit recovery selection, including its original estimates and requested tests. Absent when recovery was not assessed. |
| `recovery$routes` | Method-specific route and fixed-parameter disclosure, when recovery was assessed. |
| `fit_provenance` | PGMM native optimizer backend/code, controls and compact convergence/covariance evidence; available without full raw diagnostics. |
| `provenance$recovery` | Effective recovery mode, controls and retention policy. |

Common gene/coefficient statuses are `available`, `available_with_limitations`,
`partial` (gene-level), `point_only`, and `unavailable`. Test rows additionally
identify `unsupported`; unrequested optional tests are not fabricated as failures.
Native JLN experimental-asymptotic labels remain in the original test status.
An optional-test failure does not erase available coefficient Wald results.

The legacy `convergence` vector is deliberately not redefined. Native JLN uses
1/-70 for available/withheld inference, and accepted JLN recovery sets 1 even if
only some coefficients are estimable. HL uses 1 for an accepted native point
fit, even when its covariance fails, and -25 otherwise. PGMM retains its existing
optimizer/covariance codes. HL/PGMM recovery leaves that native code unchanged.
Use the common availability tables for the selected result.

Selection is for the whole gene fit. PGMM qualification includes all parameter
variances, information and Newton displacement. A native scalar previously
reported as available can therefore be replaced or withheld. Its original
result remains under `primary`; coefficients from different targets are not
combined in the selected summary.

Recovered output supports coefficient-wise donor-robust Wald only. Model-based
SEs/p-values, arbitrary joint-contrast tests, optional Score/LRT and random
effects are unavailable for selected recovered rows. Their native counterparts
remain under `primary`. For the free-donor-dispersion HL limiting face,
estimable points can be reported but every recovered SE and p-value is missing;
the reason is `hl_limit_sampling_target_unresolved`. Numerical covariance stored
in detailed research evidence does not override that reporting restriction.

## Retention and memory

With `keep_diagnostics=FALSE`, each LN/PGMM gene task computes its required
native tests and recovery before returning a compact payload to the coordinator.
HL remains sequential. Completed donor arrays, full candidate/root objects and
optional null fits are not accumulated across genes. Compact HL audit output
keeps identities, dimensions, limits, checked/generated/pending/dropped counts,
classification/error summaries and one bounded failure witness, not complete
event or replay lists. Failure paths obey the same retention contract.

With `TRUE`, full native/recovered fits and detailed audit histories remain
available and can be large. Retention does not change the mandatory acceptance
checks, their limits, numerical order or inference policy. Required working
history and caches are still held while a gene is assessed. This does not promise
a fixed peak RSS or a genome-wide memory bound; parallel workers also hold inputs
and active fits. Object size, uncompressed serialization, compressed files and
process high-water RSS measure different quantities.

The warning class `scrope_recovery_warning` has code
`scrope_recovery_summary` and named counts. It summarizes selected recovery and
coefficient availability; detailed reasons stay in the tables. Native control
and optional-test warnings retain their own meanings.

## Reproducible simulations

Record the exact package/archive and adapter identity, input/generation seeds
or hashes, explicit model/recovery policy, optimizer controls, requested tests,
and endpoint-specific unavailable denominators. These dataset-level records
belong in the simulation manifest rather than repeated cell arrays inside each
fit. Existing pinned simulation adapters are not automatically migrated to this
candidate. In particular, future PGMM adapters can use compact `fit_provenance`
instead of requiring raw diagnostics. Limiting JLN recovery now isolates nested
fit RNG state; reproducing historical post-call RNG side effects requires the
historical version even when fitted numerical outputs agree.

Working-target equality to a scientific target, selective coverage, general
finite-donor calibration and FDR control are not established by numerical
availability. See the method-specific JLN and HL/PGMM guides for their distinct
criteria and limitations.
