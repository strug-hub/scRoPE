# First-order HL difficult-case fixture contract v2

The authoritative `R1/C1/C2/K0/Q0/Q1/Z1/D1/W1/X0/N0/B1/S1` arrays are
preserved.  In particular, the `Z1` seed and all four `D1` seeds are not
rewritten after observing estimator behavior.

The implementation contract records two distinct facts:

- `spec_expected` retains the required result from the original fixture table;
- `expected` records the result supported by the locked generated array under
  the corrected first-order estimator.

`Z1` and the two low-variance `D1` arrays pass the public defaults
`cpc=0.005` and `mincp=5`, but seek a dispersion face.  They are therefore
boundary regressions, not low-expression-filter examples and not successful
interior fits.  The high-variance `D1` arrays remain active interior gates.

Supplemental fixtures provide estimator-testable interior cases without
concealing those boundary results:

- `Z2` uses the first certified seed in sequential order for the declared
  zero-heavy design.  Seed 104773 is retained as the first boundary attempt;
  seed 104774 is the selected interior array.
- `D2/low_v_high_phi` uses the first sequential candidate at `sigma2=0.1`.
- `D2/low_v_low_phi` first tries five sequential `n_i=5` arrays at the same
  `sigma2=0.1` and `phi=0.3`.  Only after those remain boundary seeking does
  it increase repeated cells per subject to `n_i=10`; it does not clip a
  fitted variance, add a ridge, or choose a seed by proximity to the truth.

The complete declared search path, including every rejected array and its
public-filter and fit status, is stored in
`hl_first_order_fixture_contract_v2_selection.csv`.  Unit tests regenerate
every row and lock its count and subject summaries.

The cheap fixture and public-filter checks run with the package tests.  The
complete multistart fits and cold replays are an explicit validation job:

```sh
Rscript tests/validation/save_hl_first_order_seeded_fixtures.R .
Rscript tests/validation/run_hl_first_order_fixture_contract_v2.R .
```

For a focused replay, the optional second argument is an output directory and
the optional third argument is a comma-separated case list, for example:

```sh
Rscript tests/validation/run_hl_first_order_fixture_contract_v2.R . \
  temp/hl_first_order_validation/fixture_contract_v2_targeted \
  D1_high_v_high_phi
```

Artifacts are written below the ignored directory
`temp/hl_first_order_validation/`.  The active gate requires every declared
boundary status and, for each interior replacement, point and covariance
success, stationary negative curvature, full sensitivity/meat ranks, distance
from all numerical guard faces, and a maximum cold-replay discrepancy no
larger than `1e-8`.
