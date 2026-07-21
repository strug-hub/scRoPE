# scRoPE (development version)

- Added the opt-in `ln_estimator = "joint_ln"` route. For `R0` genes it
  reports a stationary local maximizer of the NEBULA-LN approximate marginal
  criterion and evaluates model and donor-sandwich covariance at that same
  root. The public default remains `legacy_post_pml`.
- Added structured estimator and per-gene provenance, including numerical
  route, start provenance, inference availability, normalized failure reason,
  and optional-mode status. `R2` genes receive unavailable covariance and Wald
  results rather than a fallback statistic.
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
  FDR calls, and downstream results may therefore change.

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
