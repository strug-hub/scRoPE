These literal inputs contain synthetic counts and designs only. They require no
research files, RNG implementation, external download or cached fitted object.

* Recovery: selected draw 35, seed 910035, from the bounded adapter5 regression
  search (30 donors, 8 cells each, NB size 10, baseline mean 10, Gamma donor
  dispersion parameter sigma2=0.01, predictor equally spaced from -1 to 1). The native full-chi
  solver has numerically split roots near log(sigma2)=-7.6817; stable cold
  refinement resolves them. The search stopped at its first successful public
  multiple-root case. This selection is regression discovery, not calibration.
* Multiplicity: previously frozen simulation draw 93, seed 202615093
  (12 donors, 6 cells each, Gamma-NB size 4, dispersion parameter sigma2=0.3, binary donor type).
  Two separated full-chi roots persist after stable refinement.

The tests assert the public native/recovery transition, lower-level record retry,
compact reporting, full diagnostics and preservation of all point/Wald fields.
The cases are arithmetic regressions; no inferential calibration claim follows.
