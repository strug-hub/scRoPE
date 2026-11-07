## Automatic recovery and shared availability (local review candidate)

- Refine converged HL donor modes when their combined residuals can affect the inner profile score. Bound this extra work, preserve the score and step limits, and require accurately reprofiled Armijo evidence when objective subtraction is dominated by roundoff. Required history replay and covariance safeguards remain in force.
- Make `recovery="auto"` the ordinary `scrope()` policy for LN, HL and PGMM; `recovery="none"` selects native-only behavior. `scrope_hl()` has the same policy. This supersedes the earlier opt-in-only HL/PGMM stage below; the candidate is not a published release.
- Default detailed fit/history retention to `keep_diagnostics=FALSE`. Common gene, coefficient and test availability tables remain available in compact and native-only modes.
- Unify `availability$status`, `coefficient_availability$status` and `test_status$availability_status`, with separate point, covariance, statistic, reference and interval flags. Preserve method-specific test labels and the original native numeric convergence and optimizer-code meanings; use the shared tables to select inference.
- Add the classified warning `scrope_recovery_warning`, code `scrope_recovery_summary`, with aggregate analyzed-gene counts. It reports selected recovery or incomplete coefficient Wald availability even with `verbose=FALSE`; other native warnings/messages remain separate.
- Preserve native evidence in `recovery$primary` when recovery is assessed and select whole fits. Recovered working-target inference remains coefficient-wise robust Wald only; three free-dispersion HL limits remain point-only. Numerical availability does not establish selection-adjusted coverage or general calibration.
- Document filtering recipes, raw-code meanings, native-versus-selected results and selective warning handling in the common recovery guide and public help.

## Isolated joint LRT zero-quadratic repair

- Fix the shared multi-contrast adjusted-LRT error when its profile quadratic is zero or tiny. Normalize nonzero score directions before evaluating the homogeneous adjustment ratio.
- Extend the statistic to zero only at a verified coincident, accepted null fit with positive full-rank unregularized profile matrices. Report an undefined multiplier (`scale = NA`) and `gate$small_quadratic$mode = "coincident_zero_extension"`; no direction-independent multiplier is claimed.
- Withhold undefined new joint cases, including deficient rank or zero direction with positive LR, without a numeric fallback in the LRT column. Preserve scalar and ordinary positive-quadratic behavior. This optional-test repair changes neither fitting nor Wald covariance and does not establish finite-sample calibration.

## Earlier isolated HL search-history revision

- Extend opt-in HL recovery to an identified, bounded ledger of native/working outer profile and refinement probes, independent of diagnostic verbosity. Validate original scalar evidence before selected-configuration comparison; retain failed, pending and overflow evidence explicitly.
- Require complete-root/covariance agreement for tied root-like probes and audit every generated refinement trial. Preserve ordinary native results and the free-donor-dispersion HL point-only policy. That revision changed neither PGMM fitting nor public defaults; the later automatic policy is described above.
- Document capture/replay limits, runtime trade-offs and the method-specific differences from JLN; the audit establishes neither global optimality nor calibration.

## Earlier opt-in HL and PGMM working-recovery candidate

- Initially add `scrope_working(..., model="HL" or "PGMM")` for explicit working-target fits. At that stage ordinary `scrope()` defaults were unchanged; the later automatic policy is described above. The explicit long-format interface remains available and preserves accepted regular native endpoints in both the coefficient output and `primary`.
- Refit eligible fixed-positive-dispersion systems and one-dimensional whole-zero-donor limits, with free-equation sandwiches, independent derivative checks, retained competitor audits and coefficient estimability. Free-donor-dispersion HL limits remain point-only.
- Report fixed parameters, native preservation policy, route, computational certification, failures and selection limitations. Compact audits are the default; recovered optional tests and arbitrary contrasts are not supported. Numerical availability does not establish calibration or scientific-target equality.
- Distinguish coefficient-level nonfinite estimates, invalid mapped variances and overflowing Wald statistics from a fit-level certificate in recovered-output failure labels; numerical gates and accepted results are unchanged.
- Add an installed HL/PGMM working-fit guide covering route labels, fixed masks, coefficient coordinates, diagnostic provenance and the three point-only HL limits.

## Isolated numerical integration candidate

- Compute PGMM robust covariance directly from solved donor influences when the historical calculation fails or disagrees in variance-scaled units. Preserve agreeing native results.
- Add opt-in `hl_score_recovery="stable"` for certified numerical refinement of ambiguous HL Score roots; the default remains `"none"`. Preserve point/Wald results and native successful Score tests; withhold genuine multiplicity and failed cold evaluations.
- These numerical prerequisites alone do not establish small-donor calibration or enable the separately described working-fit recovery. No permutation method is included.

## scRoPE 0.3.1.9000 (draft patch candidate)

- Make JLN-first recovery available through `scrope(..., ln_recovery = "auto")`, the new LN default. `ln_recovery = "none"` preserves strict native behavior. Accepted ordinary fits are retained unchanged.
- Support bounded all-free numerical refinement, the narrowly audited dominated-start warning, and fixed-positive-dispersion working fits with the exact free-parameter donor sandwich.
- Add a general identifiable limiting-face route for a unique separating design direction affecting only whole all-zero donors. Report only estimable coefficients; nonestimable intercept/nuisance coefficients remain missing. The reduced JLN criterion is the exact full-objective limit.
- Preserve native results in `recovery$primary`; expose route/working-target/selection warnings and coefficient-level availability. Recovered rows support donor-robust coefficient Wald tests only. Their model-based SE, optional Score/LRT and random effects are unavailable.
- Audit retained native recovery endpoints, including unadopted attempts, before promoting a local fit. Keep active failure metadata consistent with the reported route and label both standardized and public covariance coordinates.
- Keep the earlier `scrope_jln_research()` wrapper for reproducibility of the old opt-in rule. It explicitly uses strict native fits and does not include the new limiting-face route.
- Data-selected route/face coverage, scientific-target equality and FDR calibration are not guaranteed. Prior finite-donor undercoverage remains documented. This development candidate is under review and is not a released version.

## scRoPE 0.3.0.9003 (isolated research candidate)

- Add `scrope_jln_research()` as an explicit opt-in combined JLN research workflow with labelled native, dominated-start warning, fixed-dispersion and withheld routes. Existing entry points and defaults are unchanged.
- Retain the original positive-domain criterion and free-parameter donor sandwich, with a bounded evaluation-guard repair. Research selection/working-target limitations and observed finite-donor undercoverage are explicit; no general selective or FDR calibration is claimed.

# scRoPE 0.3.0 — merged source

The 0.3.0 source was merged into `main` in PR #1. The development-history notes
below describe its preparation.

- Expose JLN donor-robust Score and Pace-adjusted working-objective-ratio
  p-values as experimental, asymptotic outputs through the existing opt-in
  `additional_tests` argument. Initial public scope is a single non-intercept
  coefficient tested against zero; Wald remains the default. Preserve the
  studied accepted-primary requirement and all native per-test numerical
  gates. Broader contrasts retain diagnostic-only results with explicit
  reasons; unavailable tests are not replaced with Wald or legacy inference.
  The linked methods guide distinguishes favorable cell-level calibration
  from observed donor-level inflation in the retained validation study.
- Refresh the README with supported fitting/test choices and acknowledgments;
  keep technical interpretation and migration detail in the linked guide.
- Apply MIT to the maintained 0.3.0 source under the maintainer's instruction,
  retaining verbatim upstream NEBULA MIT notices and separate dependency terms.
  Historical sources and their original license records are unchanged.


- Add the reviewed 24-gene omitted-donor-slope `toy_data` quick start, with
  exact regeneration, complete reference results and an all-gene interval plot.
  Keep existing `sample_data` and `sample_seurat` interfaces unchanged.
  Shorten the README and retain advanced methods/migration in a linked guide.
  This illustration does not establish calibration or non-null target equality.

- Integrate the preserved .9005 baseline and the reviewed .9006 candidate as
  present-day local commits. The version is a proposed release target, not an
  approved or published release; no historical campaign is relabeled.
- Update the installed-package smoke for JLN schema 3 and HL/PGMM schema 2,
  default/explicit JLN parity, and actionable legacy-request rejection.
- Add a mandatory installed-package two-worker gate and a proposed Linux CI
  workflow. A skipped parallel test does not satisfy that gate.
- Document source hashes, recover the NEBULA 1.5.3 development ancestry, and
  restore example-data attribution. Initially record conflicting license
  metadata and missing generation details without changing terms; the later
  maintained-source MIT change is described above.

- Replace inherited teaching fixtures with paired, reproducibly generated
  synthetic multi-donor examples and a matching Seurat object. Supply generator,
  fixed seed/parameters, conversion tests and illustrative README figures.
  Historical objects and scientific results retain their original sources.
  These examples do not establish calibration, power or code-license clearance.

The changes below first formed the unpublished .9006 review candidate and are
carried into this release preparation. The remaining version entries describe
source development history; they do not imply GitHub releases or tags.

# scRoPE 0.2.5.9006 — local release candidate

- Imported the exact preserved 0.2.5.9005 scientific source into an isolated
  development worktree. This candidate is unpublished and does not replace
  the original source versions or outputs of historical simulations.
- Made joint-LN the maintained `scrope()` LN default. Explicit legacy estimator
  and operational-fallback requests are rejected with migration guidance;
  rejected JLN fits are not silently replaced by another estimator.
- Retired default and explicit `nebula(method = "LN")` calls with actionable
  guidance. Explicit `nebula(method = "HL")` remains a separate historical
  compatibility API, with no automatic method substitution.
- Require accepted PGMM endpoints for optional inference: LRT needs accepted
  unrestricted and constrained fits; Score needs an accepted constrained null
  fit. Positive LR values do not override optimizer rejection.
- Map PGMM custom contrasts from reported coefficient units to the internal
  fitting coordinates, including nonzero null values. Retain the existing
  sample-centered intercept interpretation.
- Reconstruct corrected first-order HL residual means in original-design
  coordinates. Existing JLN/PGMM centered-coordinate behavior is retained.
- Add focused regressions and bounded comparisons against preserved .9005.
  The .9006 snapshot kept Joint-LN Score/Pace/LRT calibration status unchanged;
  the later scoped 0.3.0 exposure is described above. No new inferential
  calibration or final release approval was implied by this snapshot.

# scRoPE 0.2.5.9005

- Corrected-HL now passes `keep_diagnostics` through its nested and Score
  solvers. With `FALSE`, optional iteration/evaluation histories are not
  built, duplicate candidate records are compact, and full per-gene fits
  and test objects are not accumulated across genes. Required warm starts,
  objectives, derivatives, certification, covariance and test calculations
  are unchanged. The default `TRUE` retains the existing detailed output.
  Public summaries, status/reason fields and requested random effects remain
  available under either setting. Retained-object reductions do not establish
  a representative runtime or peak-memory improvement.

# scRoPE 0.2.5.9004

- Required unrestricted joint-LN core starts can receive one bounded L-BFGS
  recovery from their own failed interior endpoint, followed by ordinary
  polishing and a new complete certificate. Recovery is accepted only within
  the existing objective allowance. Original failures and recovery decisions
  remain in candidate diagnostics; required-start resolution and all inference
  gates are unchanged. Boundary, constrained, optional and PQL starts bypass
  recovery. This development candidate awaits full-fit validation.

# scRoPE 0.2.5.9002 (development)

- Joint-LN evaluates the same working objective with stable scalar arithmetic
  during optimization, polishing, certification, and donorwise audits. This
  reduces cancellation near stationary points while preserving derivatives,
  scientific certification thresholds, bounds, and inference restrictions.

- Normalized gamma ratios avoid storage proportional to the largest count.
  This development implementation requires independent review and subsequent
  cohort validation before scientific adoption.

# scRoPE 0.2.5.9001 (development)

- Fixed the PGMM BOBYQA fallback argument to use `x0`. Optimizer acceptance
  now respects each backend's return codes: `nlminb` code 0 and BOBYQA codes
  1 through 4. Budget stops and missing or malformed statuses do not count
  as convergence. Existing fallback triggers and objective ordering remain
  unchanged.

- Public PGMM main fits now use explicit `nlminb` limits of 500 iterations
  and 500 objective evaluations. `pgmm_control` permits bounded overrides,
  including the earlier 150/200 limits. Tolerances, starts, constrained-fit
  controls, and retry policy are unchanged. These limits are a development
  choice, not a guarantee of convergence or statistical calibration.

- PGMM model and robust Wald statistics, p-values, and standard errors are
  unavailable when the native optimizer has not converged. New independent
  optimizer and model-Wald status fields distinguish optimizer failure from
  covariance failure. Failed endpoints and covariance remain in raw diagnostics.

- PGMM diagnostics retain selected optimizer controls, native iterations and
  evaluations, attempted optimizer calls, covariance availability, information
  conditioning, and a scaled Newton correction. The correction is diagnostic
  only and never updates or certifies an endpoint.

# scRoPE 0.2.5.9000 (development)

- Joint-LN can now use one bounded review to recognize an additional optimizer
  endpoint as the same numerical solution. Model-based and donor-robust Wald
  inference then continues from the unchanged certified fit. Boundary,
  distinct, hard-stop, or unresolved solutions still withhold inference.

- Added plain-language per-gene fields for the estimator, inference status,
  numerical resolution, and inference source. The original strict decision and
  technical route codes remain in provenance for reproducibility.

This development build is intended for calibration simulations and is not a
public release.

# scRoPE 0.2.5 (2026-08-24)

- Improved diagnostics when joint-LN optimization finds more than one plausible
  solution. scRoPE now distinguishes near-duplicate solutions from boundary or
  genuinely competing solutions, making withheld results easier to investigate.
  In 0.2.5, this assessment is diagnostic-only and does not change which results
  are eligible for inference.

- Added an `availability` table for joint-LN fits, with one row per analyzed
  gene. It reports whether inference is available, whether only an explicitly
  requested operational estimate is available, or whether the result is
  unavailable. Concise status messages explain unavailable results, while full
  technical diagnostics remain available for advanced users.

- Improved the performance of corrected first-order HL fits by reusing
  recession checks when the design is unchanged. This reduces repeated work
  without changing the acceptance criteria or statistical results.

# scRoPE 0.2.4 (2026-07-23)

- Added one explicit `scrope(model = ...)` front door for the LN, corrected
  first-order HL, and PGMM routes. `scrope_hl()` remains available as a
  compatibility wrapper, and structured provenance advertises route-specific
  Wald, Score, and LRT availability.
- Fixed installed-package PSOCK execution for joint-LN by exposing fitter
  dependencies directly and starting workers with the exact loaded candidate
  library. Serial and two-worker fits now agree across public output,
  provenance, and internal Score/Pace diagnostics.
- Made `nbresidual()` select coefficients by exact `logFC_*` design names, so
  residuals remain valid when a route interleaves coefficient and inference
  columns. Legacy unnamed leading coefficient blocks remain supported.
- Route-specific no-op arguments now fail clearly or warn independently of
  `verbose`; legacy optional tests that are skipped record the fit-specific
  reason in `provenance$test_capabilities`.
- Added the explicit opt-in
  `ln_operational_fallback = "legacy_post_pml"`. When joint-LN inference is
  unavailable because a dispersion estimate reached its allowed boundary,
  sc-RoPE can independently fit the released legacy estimator. The separate
  sidecar exposes only operational coefficients, dispersions, boundary status,
  and diagnostics; it contains no covariance, standard errors, statistics,
  p-values, q-values, or FDR-eligible output. The primary joint-LN result is
  unchanged, and the default remains `"none"`.
- Added the opt-in `ln_estimator = "joint_ln"`. For genes that pass its
  numerical checks, it evaluates model and donor-sandwich covariance at the
  same joint fit. Otherwise, standard errors and Wald tests are unavailable.
  The public default remains `legacy_post_pml`.
- Added structured estimator and per-gene provenance, including numerical
  details, start provenance, inference availability, normalized failure reason,
  and optional-mode status.
- Added auxiliary fixed-beta conditional donor modes for joint-LN fits. These
  modes do not alter coefficient inference.
- Kept adjusted Score and scaled adjusted likelihood-ratio tests unavailable
  for `joint_ln`, with status `joint_ln_calibration_pending`, pending separate
  calibration work.
- Corrected LN donor-score allocation for zero-based positive-observation
  indices and stabilized donor/aggregate score summation. This intentionally
  changes affected legacy robust covariance and dependent tests while retaining
  the legacy coefficient-estimator architecture.
- Documented that sc-RoPE 0.2.3 operationally reported `repml$beta` after the
  raw-LN fit, whereas `joint_ln` targets the NEBULA-LN approximate marginal
  criterion. Coefficients, standard errors, p-values, convergence statuses,
  and downstream results may therefore change.

# scRoPE 0.2.3 (2026-07-14)

- Restored adjusted LRT inference for the LN NBGMM route by limiting the exact
  PGMM scalar-profile calculation to fits that expose its required derivative
  state.
- Corrected constrained LN convergence and covariance calculations to use the
  free-parameter curvature, avoiding false nonconvergence from the expected
  rank deficiency in the full constrained parameter space.
- Reorganized the README around default robust Wald inference, with optional
  adjusted Score/LRT tests and separate subject- and cell-level PGMM examples.
- Replaced the PGMM cell-level predictor warning with a brief informational
  message directing users to the robust results when cell-level overdispersion
  may be present.

# scRoPE 0.2.2 (2026-05-04)

- Strengthened exact PGMM profile-likelihood optimization with constrained and
  unconstrained retry starts, optional average-objective scaling, warm-started
  profile grids, and richer optimizer diagnostics.
- Added reduced-Hessian information-scaled convergence diagnostics for exact
  PGMM constrained profile fits, including scaled Newton-step and
  gradient-quadratic LR-bound checks.
- Enabled the information-scaled convergence check by default for exact PGMM
  profile point/grid helpers while keeping hard validity checks for finite
  fits, satisfied constraints, nonnegative raw LR, positive `Hpsi`, and positive
  `Jpsi`.
- Improved exact PGMM profile output diagnostics, including raw LR clipping
  indicators, acceptance flags, optimizer-warning indicators, retry metadata,
  gradient norms, boundary flags, and reduced-information diagnostics.
- Added regression tests for tiny negative raw-LR clipping, warning-code
  acceptance, information-scaled convergence acceptance, and rejection of
  profiles with nonpositive robust profile information.

# scRoPE 0.2.1 (2026-04-16)

- Added an internal exact PGMM profile-likelihood layer for scalar fixed-effect
  targets, including unconstrained and constrained helpers plus pointwise and
  grid evaluation utilities.
- Added ordinary relative profile likelihood and robust adjusted relative
  profile likelihood calculations for exact PGMM profile points and profile
  grids.
- Refactored shared scalar-profile reparameterization and bookkeeping in the
  PGMM LRT helpers so the new profile layer and the existing robust PGMM LRT
  path use aligned constrained-fit machinery.
- Added regression tests covering profile-point outputs, profile-grid
  structure, standardized profile-curve maxima, custom scalar contrasts, and
  agreement between the pointwise robust adjusted profile LR quantity and the
  existing robust PGMM LRT at matched null points.

# scRoPE 0.2.0 (2026-04-09)

- Added a top-level exact `scrope(model = "PGMM")` path with robust Wald,
  adjusted Score, and gated adjusted LRT support.
- Added PGMM regression coverage against legacy `nebula(model = "PMM")` and
  residual tests for the one-column PGMM overdispersion output.
- Clarified the PGMM caveat: cell-level fixed predictors are allowed, but PGMM
  only models subject-level overdispersion, so PGMM should be used with caution
  for testing cell-level predictors.
- Fixed HL Laplace tau-tau block to include the `b_i f_{i,tau,tau}` term.
- Fixed constrained HL Hessian transform from `log(phi)` to `phi` by adding the
  nonlinear diagonal correction in the `phi-phi` block.
- Fixed constrained LN score/Godambe evaluation point to use the rescued
  constrained beta vector, and aligned LN state diagnostics (`theta_ln`) with
  the pre-rescue LN objective.
- Added regression tests covering HL constrained refresh mappings, HL tau-tau
  Laplace curvature, and LN constrained/state consistency.
- Clarified `use_betas` behavior in LN/HL documentation.

# scRoPE 0.1.0

- First public release of scRoPE, derived from NEBULA with LN/HL estimation.
- Added robust Godambe-based adjusted Wald, Score, and scaled adjusted LRTs for
  linear contrasts, with optional detailed diagnostics.
- Strengthened HL pipeline stability (analytic Hessians, cached scores, and
  safeguards for negative or unstable LRT statistics).
