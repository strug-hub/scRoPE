# Synthetic teaching example: frozen design

Declared on 2026-10-20 before generating outcomes or inspecting fitted results.
The machine-readable specification is inst/examples/synthetic_parameters.R.
Seed 20261020; RNGversion 4.0.0; Mersenne-Twister / Inversion / Rejection.
One fixed realization per scenario, no seed search, no gene selection, no
repeated calibration or power experiment. All eight genes and unavailable fits
must be reported. Changes require an explicit dated explanation.

Pre-outcome design review changed an initial treatment-dependent variance idea
to a common lognormal donor distribution. That avoids confounding the example
with a treatment-dependent variance assumption and its coefficient-target
complications. No counts or fitted results existed at the time of this change.
No mixture/outlier scenario or alternative seed was evaluated.

## Population and generating mean

30 independent donors, 40 cells each: D01-D20 controls and D21-D30 treated.
8 genes, G01-G08, have generating treatment log coefficients
(-0.6, -0.3, 0, 0, 0, 0.3, 0.6, 0.9). X2 is the donor's 0/1 treatment indicator.
X1 is a standard-normal cell predictor. Each donor has 20 cells in each cc
category, case/control; positions are randomly permuted, case is the reference.
These are artificial labels, not real cell types. Offsets are
exp(Uniform(-0.25, 0.25)). All predictors/offsets are independent of latent donor
effects. X1 coefficients are 0.2 and cccontrol coefficients are -0.25; intercepts
and all parameters are listed in the specification.

For gene g, donor i and cell j:

    mu_gij = offset_ij * exp(beta0_g + 0.2*X1_ij + betaX2_g*X2_i
                            - 0.25*I(cc_ij == "control")) * omega_gi
    Y_gij | omega_gi ~ NegativeBinomial(mean = mu_gij, size = 1 / 0.25)

Thus conditional variance is mu + 0.25*mu^2. In both scenarios E(omega)=1,
Var(omega)=0.12, the marginal mean has the stated log-linear form, and the donor
effect distribution is the same in both treatment groups. Genes have
independent donor effects and conditional count draws.

## Two paired scenarios

Reference: omega ~ Gamma(shape=1/0.12, rate=1/0.12).
Misspecified: log(omega) ~ Normal(-log(1.12)/2, log(1.12)).
Only donor-distribution shape changes; the first two moments agree. At variance
0.12 the difference is modest, and a robust benefit is not guaranteed.

The gamma reference is representable by the JLN working family. Internally,
JLN uses donor dispersion v with E(omega)=exp(v/2) and relative variance
expm1(v); mean-one gamma is represented by v=log1p(0.12) and an intercept
shift -v/2. Therefore 0.12 is NOT the package's Subject parameter. Slopes are
unchanged by this intercept reparameterization. The public intercept also
centers covariates and is not compared to beta0 here. JLN still uses a
large-cell approximation at 40 cells per donor; reference does not mean exact
finite-sample likelihood or calibrated Wald intervals.

Both scenarios share covariates, offsets, donor-by-gene uniforms and
cell-by-gene uniforms. The appropriate inverse CDF transforms those uniforms.
Pairing reduces unrelated differences; it creates neither replicates nor a
guaranteed direction of the SE change. Draw order: X1, within-donor cc
permutations, offsets, gene-by-donor uniforms, gene-by-cell uniforms.

## Analysis and display plan

Fit both full datasets with scrope(model="LN", ln_estimator="joint_ln",
ln_operational_fallback="none"), supplied offsets, model.matrix(~X1+X2+cc),
ncore=1, keep_diagnostics=TRUE, cpc=0, mincp=1, otherwise ordinary defaults.
No optional calibrated Score/LRT claim. Target: reported X2 slope, null=0;
the slope is unaffected by intercept centering.

Report every gene's fit/availability, generating betaX2, estimate, model and
robust SE, robust/model SE ratio, z-based 95% Wald intervals and zero-null
Wald p-values. Publicly unavailable inference stays unavailable. Display both
scenarios side by side in fixed gene order with comparable axes and generating
coefficient markers. Both uncertainty calculations use the SAME estimate.

Primary illustrative metric: robust/model SE ratio for every gene. A larger
ratio is not inherently better. Intervals are nominal, pointwise and not
multiplicity-adjusted. Truth inclusion in one draw is not estimated coverage;
p-value differences are not evidence of power or Type I error control.
Generating coefficients are not asserted to be the exact pseudo-targets of a
misspecified approximate objective. Sandwich uncertainty does not repair bias,
mean-model errors or dependence between donors. Thirty donors leaves
finite-sample uncertainty. Narrower, wider or unavailable robust intervals must
be displayed honestly. If no clear benefit appears, say so.

## Package objects and boundaries

sample_data is the reference scenario's sparse count matrix, donor IDs,
predictors and positive offsets. sample_seurat wraps exactly those counts and
cells with replicate/celltype/tech compatibility metadata and X1/X2/cc/offset.
These are synthetic teaching objects, not HLCA, biological discoveries or
revised-paper simulation results. No historical object, snapshot, test receipt
or analysis output changes. The generator and schema/round-trip tests belong
to this maintained version. Code attribution and license declarations remain.
