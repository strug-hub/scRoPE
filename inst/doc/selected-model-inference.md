# Model and robust coefficient inference from the selected fit

Candidate 0.4.0.9007 extends model coefficient Wald output to supported recovery
routes and corrected native HL. It preserves the requested estimator, selected
point, free-parameter mask, coordinate map and estimability mask. `se_*`,
`p_naive_*` and `model_wald_*` identify model inference. `se_robust_*`,
`p_robust_*` and `wald_*` identify donor-robust inference. Both test the same
reported coefficient against zero, using the usual one-degree-of-freedom
chi-square reference. Neither family receives a selection or finite-cluster
correction.

## Primary references and version scope

The primary mathematical reference is the 61-page standalone
`scRoPE_Master_Mathematical_Reference_20261106_e087a14.pdf`, source revision
`e087a14fb90ed124381b867b99261ac07e597dc7`. It supersedes the October 22 reference.
The November 6 reference describes the joint LN working criterion, corrected
nested HL estimator, PGMM and working-recovery routes. Its descriptions of
then-current public availability are historical; this candidate adds the model
output authorized here without changing those estimators.

For HL, the reference's profile h-likelihood derivatives (Section 6.1, equations
82–87 on PDF page 24) supply the coefficient information. The conventional
covariance also matches the information convention in NEBULA revision
`70593f95342e623f721ec5919adc53868e25a2b0`: `src/optimization.cpp` computes
`vb2 = vb - vwb' diag(1/vw) vwb` in `opt_pml`, returns it as `var`, and
`R/nebula.R` inverts that matrix for coefficient SEs. This comparison concerns
the algebra at matching coordinates and fitted nuisance values. It does not
assert that the two packages' optimizers or dispersion estimators return
numerically identical fits. In particular, the canonical conditional HL intercept
and a marginal or mean-centered intercept must not be interchanged across
parameterizations. The standalone reference separately derives the
augmented HL sensitivity used for robust inference.

## PGMM: free marginal information

For the negative marginal log likelihood q(a), a contains all free coefficients
and free log donor variance. Let H = d²q/da da' at the selected root, and U contain
one free-score column per donor. The model covariance is H^(-1); the robust
covariance is H^(-1) U U' H^(-T). A public coefficient derivative map D gives
D H^(-1) D' and D H^(-1) U U' H^(-T) D'. There is no extra donor-count multiplier.

When variance is fixed, restrict H and U to the actual free coordinates *before*
inversion. When variance is estimated, invert the full free system and then
extract/map the coefficient block. Inverting a coefficient block with estimated
nuisance variance would instead give a different conditional covariance.
The existing native PGMM calculation and optimizer qualification are unchanged.

## JLN: established working-model inverse curvature

For the joint LN negative working criterion, the same H^(-1) and donor sandwich
construction applies to the selected free system, including any free log donor
variance and log cell precision. This is the established native JLN model
convention, extended to refined, search-warning, fixed-mask and supported
identifiable-face fits. The criterion is an approximation; the label makes that
convention explicit and does not claim exact lognormal marginal likelihood.

Cell overdispersion phi in JLN public output is reciprocal NB precision gamma.
A fixed public phi=0.001 means gamma=1000 in the criterion. Fixed values do not
receive score columns or covariance rows in the free inverse. Public coefficients
use the existing centered-design scaling map. A limiting-face map is applied
only after inversion of the regular reduced free system. Nonestimable
coefficients and corresponding reported covariance rows and columns stay missing.

## Corrected HL: profile h information, not augmented sensitivity

Write p_i(beta,d) = h_i(beta, eta_hat_i(beta,d), d), where eta_hat solves the
individual mode equation and d denotes the selected dispersion values. The
profile coefficient Hessian is

    p_beta_beta = h_beta_beta - h_beta_eta h_eta_beta / h_eta_eta.

For beta = beta_star + K gamma, the existing inner Jacobian is
B = sum_i p_{i,gamma gamma}. The conventional model coefficient covariance is

    V_gamma_model = (-B)^(-1),
    V_beta_model  = K V_gamma_model K'.

The original-design map M_X then gives M_X V_beta_model M_X'. Dispersions are
plugged in at their selected values whether they were estimated or held fixed.
This is a conventional profile h-likelihood coefficient-information covariance;
it is not a covariance for all fitted dispersions and does not assert an
information identity for the full nested estimator. It profiles the individual
modes, rather than treating them as fixed.

Robust HL inference continues to use the existing generally nonsymmetric
augmented estimating-equation sensitivity A and donor contributions Psi:
A^(-1) Psi Psi' A^(-T), followed by the existing scientific/public map.
Neither A^(-1) nor its symmetric part is used as HL model covariance. The native
inner/outer curvature, root and robust certification rules are unchanged.

## Selection, failure and storage

Complete accepted native fits return directly with no recovery certificate,
extra root audit or refit. A genuinely incomplete native family can trigger
existing eligible recovery. Replacing a point-qualified native partial fit
requires a certified recovered result with both families complete for every
estimable coefficient. Incomplete/failed recovery retains the native fit and
each supported family, without splicing points or covariance between fits.
Failure of an added model calculation does not erase valid robust inference or
invalidate an otherwise qualified point. Existing strict native JLN qualification
is unchanged; this extension does not revive native candidates rejected by that
policy. This candidate adds no condition-number threshold.

`coefficient_availability` reports `model_available`, `model_covariance_available`,
`model_se`, `model_statistic`, `model_p_value`, `model_interval_available`,
`model_failure_reason`, `model_status`, and the covariance convention.
The existing `available`, `covariance_available`, `interval_available`,
`statistic` and `p_value` remain robust-family fields. `robust_status` makes that
family explicit. A valid SE follows covariance availability even when its Wald
statistic overflows; the statistic and p-value remain unavailable. This applies
to both families, including compact and diagnostic working output.
`availability` counts both families separately. Overall
`point_only` means neither family supplies Wald inference; one available family
and an unavailable supported family produces `partial`. Optional-test failures
are recorded separately and do not trigger whole-fit recovery.

`scrope_working()` adds model SE/p-value/availability columns to its long table.
The compatibility `scrope_jln_research()` wrapper adds `se_model`, `p_model` and
`model_available` for its prespecified target; existing robust intervals retain
their meaning. Full covariance and superseded native objects require diagnostics.
Compact results retain coefficient-scale values and short status/provenance,
without donor arrays, a second full native fit, or full model covariance matrices.

## Limits that remain

Fixed-positive and data-selected limiting-face fits have explicitly changed
working targets. Neither covariance construction proves equality with the
scientific generating target, coverage after route/face selection, finite-donor
calibration or FDR control. The search-warning route retains its unresolved
search caution. Supported fixed-v HL limiting faces use their regular reduced
profile information; free-v HL limiting faces continue to report points only
because the sampling interpretation of the vanished-donor contribution remains
unresolved. Both families are withheld there, even if a numerical matrix exists
in full diagnostics. Unsupported geometry, competing roots and failed required
checks remain unsupported. Recovered optional Score/LRT, arbitrary joint
contrasts and random effects are not added.
