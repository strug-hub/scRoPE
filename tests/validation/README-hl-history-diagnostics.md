# Temporary HL history diagnostics

Set `SCROPE_HL_HISTORY_DIAGNOSTICS=true` when running `tests/testthat.R` to observe the original history-audit decisions. The Linux workflow enables this flag. Ordinary package calls do not install the observers. There are no additional fits, retries, changed acceptance rules or changed identity checks.

The existing `linux-r-release-gates` CI artifact retains:

```
scrope-check/scRoPE.Rcheck/tests/testthat/_diagnostics/hl-history-<pid>-<counter>.rds
```

A record is written for every audit returning failure and for successful audits whose selected data have 480 cells, 60 nonempty donors and two zero donors (the failing nuisance-limit fixture and its transformed counterpart). Tests deliberately exercising refusal also produce records; inspect `calls`, the configuration and source labels to identify the focal test. An unexpected exception escaping the audit itself has no returned result to record.

Each schema-1 record contains:

- `history`: the full history supplied to the audit, including event coordinates, source labels and full configurations/data.
- `selected_fit`: the unchanged selected fit used by required replay.
- `audit`: the original result, including failure witness, selected replay coordinates, counts, configuration summaries and any retained details. Events generated inside the audit are not added to the saved input history; their available evidence is in the audit result, failure witness and refinement records.
- `original_profile_errors`: errors raised by the original cold profile, including its local `rho`, `active`, `source` and failed `ans`, where those locals exist. An error thrown before `ans` is assigned cannot provide that result.
- `original_refinement_failures`: unchanged failed refinement results, including their reasons and traces.
- `original_refinement_resolutions`: returned traces, stationarity checks and full-state comparisons when a stalled historical root satisfies the coalescence fallback. Such audits are saved even outside the focal fixture; the log's `resolved_stalls` count identifies actual use of the fallback.
- `environment`, `keep_details`, `calls`: R/session/numerical-library information, thread-related environment settings, execution/DLL identity, retention mode and calling context.

`HL_HISTORY_DIAGNOSTIC` lines identify saved files. Failed witnesses and raw failure details are printed to `testthat.Rout.fail` (or `testthat.Rout` when the tests pass), retained by the same artifact upload. R CMD check may show only an excerpt in the Actions console. A diagnostic write error is printed as `HL_HISTORY_DIAGNOSTIC_WRITE_ERROR` and does not replace the audit result or add a warning. Missing records must not be interpreted as successful replay.

The observers wrap only `working_*` functions in the test process, call each original function once, and restore the previous bindings when the runner returns or errors. They do not replace the implementation fingerprint or native-library checks. Test mocks may temporarily replace an observer; mocked failures can therefore lack raw profile/refinement details.

Run the installed diagnostic behavior gate with:

```sh
Rscript --vanilla /path/to/library/scRoPE/tests/validation/run_hl_history_diagnostic_gate.R \
  /path/to/diagnostic-gate.rds /path/to/library/scRoPE
```

The gate verifies unchanged real audit results, exact synthetic refusal results and captured inputs, both raw failure paths, write-error isolation, and binding/identity restoration. The synthetic failures are not evidence of the unresolved Linux cause.

These observers add allocations, call frames and file I/O; they may affect timing or reproducibility of an intermittent numerical failure. A passing diagnostic run is evidence to compare, not a solver repair. Artifacts contain synthetic test data and runtime/library paths. Remove the temporary workflow flag after the Linux witness has been investigated.
