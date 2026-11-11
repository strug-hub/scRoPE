# scRoPE 0.4.1 (2026-11-11)

- Return a complete, qualified native fit directly after its existing optimizer,
  geometry, curvature and inference checks. Remove the extra PGMM post-acceptance
  0.001-SE recovery trigger and redundant native recovery certificate.
- Attempt eligible recovery for genuinely incomplete inference. A complete,
  qualified recovered result supplies estimates and covariance from the same fit;
  its working target is explicit. If recovery fails or remains incomplete, retain
  usable native partial results. Full superseded native fits are retained only
  when diagnostics are requested; compact original status and reasons remain.
- Supply model-based and donor-robust coefficient Wald inference from the same
  selected fit on supported PGMM, JLN and corrected HL routes. PGMM/JLN use inverse
  free information/curvature; HL uses inverse profile h-likelihood coefficient
  information at fitted dispersions. The HL augmented robust sandwich is unchanged.
- Keep covariance, standard-error, interval and Wald-test availability distinct.
  A valid standard error remains available when a squared Wald statistic overflows;
  the statistic and p-value are withheld with an explicit reason.
- Unify native qualification, recovery attempted/outcome and availability messages
  across methods while preserving method-specific reasons and convergence codes.
  Existing public tables gain compact fields; result/report schemas advance to 6/3.
- Compensate positive-cell PGMM objective summation to reduce cancellation and
  summation error. This arithmetic change is separately reviewed from recovery policy.

Supported recovered inference now includes both model-based and donor-robust
coefficient Wald results. Free-donor-variance HL limiting faces remain point-only;
optional recovered Score/LRT, arbitrary joint tests and random effects are not
added. Working targets may differ from scientific targets. Availability does not
establish general finite-donor calibration, selection-adjusted coverage or FDR
control. Historical analyses retain their actual package versions and policies.
See the [selected-model inference guide](inst/doc/selected-model-inference.md)
and [automatic recovery guide](inst/doc/automatic-recovery.md).

The development entries below preserve the unpublished candidate lineage. The
0.4.1 summary above describes the final policy and supersedes interim restrictions.

# scRoPE 0.4.0.9009 (isolated status-reporting candidate)

- Add compact native qualification, inference flags/counts and original reasons
  to the existing common tables. Report actual recovery attempted/outcome,
  caught errors, selected source and a consistent availability message.
- Keep valid selected JLN recovery SEs when the Wald statistic overflows; withhold
  the statistic/p-value and report the exact coefficient refusal reason.
- Preserve existing failure normalization prefixes, method-specific codes,
  convergence codes and optional-test metadata. Extend aggregate warning counts
  without changing native acceptance or recovery selection.
- Correct coefficient-partial and warning-count documentation. Full superseded
  native fits remain diagnostics-only; no inference capability is broadened.

# scRoPE 0.4.0.9008 (isolated reporting repair)

- Retain a valid robust SE in working and common output when only the Wald
  statistic overflows. Covariance, SE and interval availability remain distinct
  from Wald statistic and p-value availability; fitting and selection are unchanged.
- Permanently test held-cell/free-donor HL and JLN model and robust covariance
  reconstruction at the same selected fit using existing deterministic fixtures.

# scRoPE 0.4.0.9007 (isolated model-inference candidate)

- Report model and robust coefficient Wald inference from the same selected fit
  for supported PGMM, JLN and corrected HL native/recovery routes. PGMM/JLN use
  inverse free information/curvature. HL uses inverse profile h-likelihood
  coefficient information at fitted dispersions; its augmented sandwich is unchanged.
- Expose both families separately in compact availability and warning summaries.
  Genuine missing model inference may now trigger recovery; complete recovered
  replacement requires both families. Retain native partial inference when recovery
  is incomplete, and retain full superseded objects only with diagnostics.
- Preserve unsupported free-donor-variance HL limiting faces as point-only.
  No optional recovered Score/LRT, selection correction or calibration is added.
- Mathematical reference: standalone 2026-11-06 master, source e087a14; NEBULA
  comparison pinned to 70593f95342e623f721ec5919adc53868e25a2b0.

# scRoPE 0.4.0.9006 (local review candidate)

- Report native PGMM model-only inference failure as partial with the same
  classified warning under recovery="none" and recovery="auto". The repair
  reads existing coefficient flags and does not refit, certify, change estimates,
  or add any model-based inference to recovered fits.

# scRoPE 0.4.0.9005 (isolated candidate)

- Qualify complete native HL/PGMM once. Remove PGMM post-acceptance 0.001-SE
  recovery trigger and redundant native recovery certificate. Reuse native HL
  geometry and existing coefficient inference flags.
- Attempt recovery for incomplete robust inference; select complete certified
  recovered robust output as one fit, otherwise retain qualified native partial
  results. Missing model-only inference does not trigger robust-only recovery.
- Store superseded native results only with `keep_diagnostics=TRUE`, including
  explicit working output and JLN retention; preserve compact native status.
- Use one aggregate warning contract across ordinary and explicit working output.
- Separately: compensate positive-cell PGMM objective summation and identify
  actual recovery-certificate failure gates without changing their tolerances.

# scRoPE 0.4.0 (2026-11-07)

- Make automatic recovery the default for LN, HL and PGMM through
  `scrope()` and `scrope_hl()`. Qualified native fits are preserved; eligible
  numerical refits, fixed-positive-dispersion working fits and identifiable
  whole-zero-donor limits disclose their routes, fixed parameters and targets.
  Set `recovery="none"` for native-only fitting. Existing explicit LN
  `ln_recovery` policies and the explicit working/research interfaces remain
  supported.
- Default to compact results with `keep_diagnostics=FALSE`. Shared gene,
  coefficient and test availability tables distinguish point, covariance,
  statistic, reference and interval availability. Original native evidence is
  preserved in `recovery$primary` when recovery is assessed; selection applies
  to a whole fit. A classified `scrope_recovery_warning` reports aggregate
  recovery and availability counts even with `verbose=FALSE`.
- Improve PGMM robust covariance using solved donor influences when the
  historical sandwich fails or disagrees in variance-scaled units. Native
  qualification can replace or withhold previously reported scalar results;
  agreeing qualified native results are preserved.
- Improve HL donor-mode polishing and profile refinement with bounded work,
  retained-history audits, accurate Armijo evidence near roundoff, and cold
  replay checks. Covariance and unresolved-root safeguards remain active.
  Opt-in `hl_score_recovery="stable"` refines eligible ambiguous Score roots;
  its default remains `"none"` and genuine multiplicity stays withheld.
- Handle zero or tiny profile quadratics in the shared joint adjusted LRT.
  An exact zero extension requires coincident accepted null states and
  full-rank positive profile matrices; its undefined multiplier remains `NA`.
  Unsupported joint cases remain withheld.
- Isolate caller RNG state during nested JLN limiting recovery. Require
  installed-package smoke checks and LN/PGMM two-worker regressions in CI.

Recovered inference supports coefficient-wise donor-robust Wald only; recovered
model-based SEs, optional Score/LRT, arbitrary joint tests and random effects
remain unavailable. Free-donor-dispersion HL limiting fits remain point-only,
without recovered SEs or p-values. Working targets may differ from scientific
targets, and selection-adjusted coverage, general finite-donor calibration and
FDR control are not established. Historical simulations and HLCA results retain
the versions and policies that produced them; this release does not relabel or
revalidate those analyses. Recovery audits can be expensive, and compact output
does not establish a genome-wide peak-memory bound.

## Development history

The following entries document earlier source snapshots and their preparation;
they do not imply published releases or tags. The 0.4.0 summary above supersedes
the earlier opt-in and retention defaults.

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
