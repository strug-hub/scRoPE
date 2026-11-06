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

The common contract applies to **`scrope()` for all three models and
`scrope_hl()`**, including `recovery="none"` and `keep_diagnostics=FALSE`.
The separate `scrope_working()` and `scrope_jln_research()` interfaces retain
their own long-format schemas and reporting behavior.

| Public field | Allowed status | Meaning |
| --- | --- | --- |
| `availability$status` | `available` | All reported coefficient Wald results are available from the native fit. |
| | `available_with_limitations` | All reported coefficient Wald results are available from a selected recovery route. |
| | `partial` | Some, but not all, coefficient Wald results are available. |
| | `point_only` | At least one accepted coefficient estimate is available, but no coefficient Wald result is available. |
| | `unavailable` | Neither accepted coefficient estimates nor coefficient Wald results are available. |
| `coefficient_availability$status` | `available`, `available_with_limitations`, `point_only`, `unavailable` | The corresponding classification for a single coefficient; there is no coefficient-level `partial` value. |
| `test_status$availability_status` | `available`, `available_with_limitations`, `unsupported`, `unavailable` | Per-test success, successful recovered inference with limitations, an unsupported requested test, or another failure. |

Gene `selected_source` and coefficient `source` use `native`, `recovery` or
`unavailable`. The unavailable source denotes a withheld recovery route; a
native-only failed result can have source `native` and status `unavailable`.
Route names, target descriptions and detailed failure reasons remain
method-specific. `selection_adjusted` is FALSE; common status names do not make
the statistical targets or calibration guarantees identical.

Gene rows include `selected_point_coefficients`, `selected_covariance_coefficients`,
`selected_wald_coefficients` and `total_coefficients`. Coefficient rows separate
`point_available`, `covariance_available`, `statistic_available`,
`reference_available`, `interval_available`, `available` and `estimable`.
`available` requires the coefficient Wald success flag and a finite p-value;
`statistic_available` additionally checks a finite statistic. Interval
availability requires an accepted point, covariance and finite positive reported
SE; it is a reporting flag, not a coverage guarantee. A finite raw estimate in
`summary` alone does not establish any of these qualifications.

Test rows use `success` and `failure_reason`. The original `test_status$status`
is retained separately: for example, a successful native LN Score/LRT can have
`status="experimental_asymptotic"` while its common `availability_status` is
`"available"`. Unrequested optional tests are not fabricated as failures. An
optional-test failure does not erase available coefficient Wald results or
change the gene's coefficient-Wald classification. `summary$test_status` is
`"ok"` or a comma-separated list of failed test families, not a convergence code.

LN also retains the compatibility column `availability$availability`, with
values such as `"inference"` and `"partial_inference"`, plus its existing
technical/status-message columns. Use `availability$status` for cross-model
classification. Common reporting has `recovery$report_schema_version=1L` and
`provenance$result_schema_version=4L`.

## Numeric convergence and native evidence

The legacy `fit$convergence` vector is deliberately not redefined:

| Model | Native numeric codes | Effect of selecting recovery |
| --- | --- | --- |
| LN / joint-LN | `1`: native inference available; `-70`: native inference withheld. | Accepted recovery sets `1`, even if only some coefficients are estimable. |
| HL | `1`: accepted native point fit, including a fit whose covariance fails; `-25`: native point fit unavailable. | The original native code is unchanged, even when recovered Wald inference is available. |
| PGMM | `1`: the native composite code remains positive; `-50`: optimizer not accepted; `-25`: model covariance wholly unavailable; `-27`: robust covariance unavailable when the code would otherwise be positive. | The original native code is unchanged. `-25` can override `-50`; a positive code alone does not qualify every coefficient. |

PGMM's raw optimizer codes have another meaning. In `fit$fit_provenance`,
`raw_convergence` is interpreted together with `optimizer`: `nlminb` accepts
`0`, whereas `bobyqa` accepts `1:4`. NLopt codes `5/6` are budget stops, not
accepted convergence. `optimizer_success` records this backend-specific check;
it is not a selected-recovery success flag. These native optimizer fields remain
available without full diagnostics.

`fit$availability$native_convergence` records the original code for each gene.
When recovery is assessed, `fit$recovery$primary` preserves the native estimates,
tests and convergence; `fit$summary` and the common availability tables describe
the selected output. With `recovery="none"`, there is no separate `primary`
result, and the shared tables describe native output directly. Do not use
`fit$convergence == 1` as a cross-model inference filter.

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

## Filtering results

Select the coefficient of interest from the common table so that an unavailable
intercept or nuisance coefficient does not discard an available focal contrast:

```r
coefficient_results <- fit$coefficient_availability
wald <- subset(coefficient_results,
               available & statistic_available & reference_available)
wald[, c("gene_id", "gene", "coefficient", "beta", "robust_se", "p_value",
         "source", "route", "sampling_target")]

native_wald <- subset(wald, source == "native")
point_only <- subset(coefficient_results, point_available & !available)
complete_genes <- subset(fit$availability,
                         selected_wald_coefficients == total_coefficients)

# Empty if optional tests were not requested; keep the original status label.
optional <- subset(fit$test_status, test %in% c("score", "lrt") & success)
optional[, c("gene_id", "gene", "contrast", "test", "status", "availability_status")]
```

These are numerical-availability filters. Keep source, route, target and
method-specific status with the values used downstream; filtering to native
rows does not itself establish unconditional or selection-adjusted coverage.
Recovered working-target rows must not be silently pooled with native rows as
though the targets were identical. Retain missing-result denominators in summaries.

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

## Recovery warnings and other messages

The common aggregate is a **warning condition**, not a `message()` or a numeric
convergence code. Its classes are `scrope_recovery_warning`, `warning` and
`condition`; its `code` is `"scrope_recovery_summary"`. These are condition
fields, not additional columns in the returned fit. The message starts with
`"scRoPE LN recovery:"`, `"scRoPE HL recovery:"` or `"scRoPE PGMM recovery:"`.

At most one such aggregate is emitted per public call, across analyzed genes,
when any selected source is recovery or any gene is `partial`, `point_only` or
`unavailable`. It is silent when all genes have native coefficient Wald results.
It still runs with **`verbose=FALSE`**, `keep_diagnostics=FALSE`, or
`recovery="none"`. In native-only mode its recovered count is zero: the warning
also reports availability, and does not imply that recovery was attempted.
Optional-test failure alone does not trigger this common warning when native
coefficient Wald is complete. Normalizing `recovery$primary` is silent.

The condition's `counts` vector has names `genes`, `recovered`, `all_wald`,
`partial_wald`, `point_only`, and `unavailable`. `recovered` counts selected
recovery routes, including point-only recovery; it is not a count of successful
Wald tests and overlaps the availability categories. `all_wald` and
`partial_wald` count genes by coefficient availability. Genes removed by the
input filters are outside these counts. A recovered result also appends the
working-target/selection/calibration caution to the message.

Capture and selectively muffle only this classified summary when a pipeline
records it elsewhere; other warnings still propagate:

```r
recovery_notices <- list()
fit <- withCallingHandlers(
  scrope(count, id, pred, offset=exposure, model="HL", verbose=FALSE),
  scrope_recovery_warning=function(w) {
    recovery_notices[[length(recovery_notices) + 1L]] <<-
      list(code=w$code, counts=w$counts, message=conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
```

This captures the condition without parsing its human-readable text and does
not change the returned statuses. Native ignored-control/rank warnings and
optional LRT fallback/ridge warnings are separate and can occur even with
`verbose=FALSE`, including during the native stage before recovery selection.
Verbose-only output includes progress, the PGMM cell-predictor notice and HL
optional-test messages. The separate research wrappers keep their own messages;
the shared condition is not a replacement for every method-specific diagnostic.

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
