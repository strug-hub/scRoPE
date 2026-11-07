# Numerical test policy

Decide what an assertion promises before choosing its comparator. An independent
optimization can follow a slightly different arithmetic path on a different R,
compiler or BLAS build, even with identical data and controls. Bitwise agreement
between those fits is not the package's numerical accuracy guarantee.

## Exact contracts

Use `expect_identical()` for a value copied or preserved from the **same fit**,
serialized identity records, loaded backend identity, status and refusal reasons,
route selection, masks, names, dimensions, types, fixed parameters, and RNG
restoration. Keep explicitly required mathematical zeroes exact. A status change
or an available value becoming missing is a behavior change, never roundoff.

For example, a successful native Score record must replay its recorded result
exactly. Comparing two separately fitted compact and detailed reports is a
different assertion, even when both calls use the same data.

## Independent refits

Use `expect_recovery_refit_table()` with an explicit list of the continuous
columns being compared. It keeps table structure and all unlisted columns exact,
and checks numeric types, attributes, `NA`, `NaN` and signed infinities exactly.
It compares **every finite element** using

```
abs(actual - expected) <= absolute + 1e-10 * max(abs(actual), abs(expected))
```

The implementation scales before subtraction to avoid overflow. There is no
averaging across a vector: a large coefficient cannot conceal an inaccurate
small value. The relative budget is 10,000 times smaller than the existing
`1e-6` root/Newton certificates. The reported Linux differences in the limiting
PGMM fixture were at most approximately `3e-12` in absolute terms and `2e-12`
relative to nonzero fitted values, comfortably inside this budget.

The absolute budget is `1e-12` only for signed coefficients (`beta`,
`coefficient`, `logFC_*`) and interval endpoints (`ci_lower`, `ci_upper`) that
can cross zero under roundoff. This floor is a test precision requirement in
the reported coefficient units, not a universal scientific equivalence margin.
Use a separately justified scale if testing radically rescaled predictors.
SEs, variance diagonals, statistics, dispersions and p-values have
**no absolute floor**. In particular, `1e-100` and `1.01e-100` must fail, and a
positive value changing to zero must fail. Keep separate invariants for variance
positivity, estimability, finite inference and numerical certificates.

For independent covariance or information matrices, use
`expect_recovery_refit_matrix()`. Its off-diagonal budget uses the product of
the coordinate scales `sqrt(abs(diag_i)) * sqrt(abs(diag_j))` (taking the larger
diagonal magnitude from the two matrices). Off-diagonals can cancel to nearly
zero even when their coordinates have a substantial scale. Diagonal errors
remain relative to each diagonal itself: a 1% error in a `4e-16` variance still
fails beside a variance of 2. Check known mathematical zeroes separately with
exact assertions; this comparator does not assert that an empirical covariance
must vanish. Names, dimensions and nonfinite states remain exact.

For example:

```r
expect_recovery_refit_table(compact$coefficient_availability,
  detailed$coefficient_availability,
  numerical = c("beta", "robust_se", "statistic", "p_value"))
expect_recovery_refit_table(compact_working$summary, detailed_working$summary,
  numerical = c("beta", "robust_se", "p_value"))
```

The JLN helper's default field selection retains its existing public-summary
naming contract; new schemas should list their columns explicitly. Do not pass
all numeric columns automatically: numeric identifiers, counters and fixed
settings are exact contracts too. `test-refit-comparison.R` checks that the
helper accepts roundoff but rejects status changes, missingness changes, schema
changes, and material relative errors even for tiny entries.

## Mathematical accuracy and future failures

Nested solvers need an accuracy check below their parent solver's acceptance
scale. Individually admissible donor residuals can combine into an inadmissible
HL profile-score error. `test-hl-inner-accuracy.R` retains the numerical inputs
of a failed required cold profile, checks both implementations against independent
bracketed donor roots, and includes a state whose unpolished score falsely suggests
convergence. It also checks that positive objective roundoff cannot override
negative or unresolved line-integral evidence, and that polishing, backtracking
and iteration budgets still refuse unsuccessful solves. These are production
accuracy regressions; a refit comparator cannot replace them.

For accuracy against an independent formula, derivative, solution or frozen
reference, use a problem-specific error budget. Refit agreement alone cannot
prove accuracy. Existing examples include `helper-pmg-portability.R`, which
uses independently solved donor-total and conditional-slope equations together
with Newton corrections, and the tiny-variance sandwich tests, which use
variance-scaled discrepancies. Preserve these checks rather than replacing them
with the refit helper.

When a platform failure occurs:

1. Save full-precision actual/expected values, statuses and missingness, plus R,
   dependency, BLAS and platform information. Identify whether the assertion
   compares preserved data, separate fits, or independent mathematics.
2. Confirm that model choice, inference availability and accuracy certificates
   agree. Investigate any change to these before modifying a comparator.
3. For independent fits only, select the continuous fields and justify their
   scale and budget. Add a negative control showing that the intended error or
   refusal would still fail. Never round outputs, loosen production criteria,
   remove identity checks, or blanket-replace exact assertions to make CI pass.
4. Run the changed complete test files, relevant comparison-helper users, and
   the installed-package and actual two-worker gates. Require the Linux PR job
   to pass before considering cross-platform validation complete.

The current repair covers the demonstrated independent limiting-PGMM refits
and separates same-fit native preservation from independent public refits in
the two affected files. Other exact comparisons in retention and streaming tests
still need classification when those areas change: use same-fit fixtures to
test exact preservation and the explicit comparator for independent refits.
Do not infer that every remaining exact numeric assertion is erroneous.
