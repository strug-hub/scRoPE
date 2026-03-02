# scRoPE 0.1.1 (2026-03-02)

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
