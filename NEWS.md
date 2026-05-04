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
