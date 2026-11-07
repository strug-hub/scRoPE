`hl-inner-event106.rds` contains numerical inputs from PR #3 Linux run
42346009411, artifact 13564027939, record `hl-history-8403-016.rds`:
event 106, `nuisance_face:round0:trust_proposal:96`. The underlying nuisance-limit
fixture is generated with seed 90211 in `test-working-recovery.R` and
`test-auto-public-recovery.R`. Exact reduced data and the selected starting beta
are retained because tiny changes can alter the original solver trajectory.

The additional `false_convergence` beta and donor modes come from a local
controlled experiment using those same inputs and a cancellation-first line
search with the old donor accuracy. Its reported score is below `1e-10`, while
independently polished donor roots give a score above `2e-10`. The regression
uses independently bracketed scalar roots and directly summed analytic profile
scores and curvatures to check that discrepancy. Its oracle does not call a
production mode solver, evaluator or primitive. Both the compiled and literal R
implementations must resolve the false-convergence state under the existing
acceptance limits.

The file contains no functions, execution identities, audit decisions or history
seals. It is numerical input for a new local solve, never a cross-platform replay
record. Regressions require accurate results, not a particular old failure or
exact iteration count on every platform.
