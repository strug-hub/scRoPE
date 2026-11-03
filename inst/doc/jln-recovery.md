# Identifiable JLN inference on a whole-zero-donor face

## Objective identity

Use the existing JLN parameterization: donor lognormal variance v > 0, cell size gamma > 0 (cell overdispersion phi = 1/gamma), and mu_ij = exposure_ij exp(x_ij beta). For an all-zero donor, write S = sum mu, alpha = 1/expm1(v), lambda = alpha exp(-v/2), r = S/(S+lambda), and w_j = mu_j/S. Its negative JLN working log likelihood is exactly

q_i = alpha { -log(1-r) - r } + gamma sum_j log(1 + alpha r w_j/gamma).

Both terms are nonnegative. As all means tend to zero, r tends to zero and q_i tends to zero. This expression includes the original donor constants; it is not a replacement NB/Poisson objective.

Suppose a design direction d satisfies X_F d = 0 on the retained face and X_R d < 0 elsewhere, all removed counts are zero, and every donor is wholly retained or wholly removed. Along beta(t) = beta_0 + t d, retained contributions are unchanged and removed contributions tend to zero. Every finite full objective is at least its retained objective. Conversely any retained linear predictor is attained on this path when the retained design has rank p−1. Consequently the infima agree in the model closure. Numerical branch/global optimization remains a separate issue; the implementation requires certified retained roots and rejects a reduced objective worse than a known finite full candidate.

On a compact neighborhood with positive finite v and gamma and regular retained parameters, the projected first and second derivatives and donor scores also converge. Numerical finite-path checks verify this using the existing full JLN evaluator. This is not uniform in vanishing or diverging dispersion.

The competitor audit reevaluates retained candidate, optimizer, polishing and native-recovery history endpoints, including unadopted attempts and basin probes, with the same scalar criterion. It checks objective provenance and recorded-scalar agreement, includes numerical evaluation allowances, and withholds recovery when retained evidence cannot be verified. The warning and all-free refinement routes use the same retained-history audit.

## Coordinates and covariance

The detector uses the same centered/scaled design as JLN. A unique one-dimensional nullspace of positive-count rows must satisfy machine-scale predictor identities, strict negative margins on removed cells, and whole-donor removal. Multiple null directions, partial donor removal and ambiguous near-identities are unsupported. A column basis of the retained design defines reduced coefficients eta with S_F B eta = S_reduced eta. Centering differences are included in B. Both dispersions remain estimated unless the separately certified fixed-positive route applies.

For the actual free vector a (reduced coefficients and all free dispersions), u_i(a) = −d q_i/da, A = sum_i d²q_i/da da', B_score = sum_i u_i u_i'. The covariance is A^−1 B_score A^−T, evaluated as tcrossprod(solve(A,U)); there is no extra 1/G and no finite-cluster multiplier. Free derivative blocks are selected before inversion. Removed donors have exactly zero limiting scores; embedding their zero columns changes neither A nor B_score. No singular full inverse, generalized inverse, ridge or pinned nuisance value is used.

A requested coefficient/contrast c is estimable only if c annihilates d in the corresponding coordinate scale. Map the free covariance by the actual linear derivative map, D V D'. Nonestimable coefficients and their covariance rows/columns are NA; a finite arbitrary representative is retained only internally for objective checks. A finite slope does not make an escaping intercept or region coefficient estimable.

## Inferential scope

The regular sandwich statement concerns a fixed limiting face and, when used, fixed positive dispersion mask/epsilon. Require an identified regular free root, positive nonsingular free sensitivity, independent iid donor draws from an appropriate superpopulation (or separately justified negligible individual score means / centering for nonidentical donors), sufficient moments, a donor CLT and no dominating donor. Aggregate expected-score zero alone is insufficient for raw outer products under arbitrary nonidentical donor misspecification. Unequal donor sizes must satisfy the corresponding scaling/CLT conditions.

The expected free-score root defines the working target. It need not equal the full finite JLN target or the scientific generating coefficient. A structural-zero population face with consistent face identification can support fixed-face asymptotics; a face selected because rare positive-mean donors happened to have zero counts is a selection event. Neither the sandwich identity nor deterministic detection proves selective coverage. The prespecified fresh simulation screen examines this practical transition without claiming universal or FDR calibration. Existing observed fixed-route undercoverage remains relevant. Ordinary Wald inference is withheld if the remaining free root is nonregular or a requested contrast is not estimable.

## Numerical refinement

A separate bounded all-free refinement is allowed only with at least two distinct certified required core starts, agreement of all certified roots, and remaining finite interior candidates failing only root/Newton tolerances. Each such endpoint and the anchor are reoptimized with the original criterion, exact derivatives and full free sandwich; both must certify and agree with the anchor. This resolves a numerical tolerance failure and does not waive an unresolved boundary or competing branch.

## Bounded empirical evidence

The fresh 300-dataset nuisance-face screen used 100 datasets per setting. With 54 total donors (52 retained donors × 18 cells and two singleton donors), automatic recovery reported 97/100 structural-zero datasets, rejecting 7/97 reported (7/100 generated), and 100/100 finite rare-mean datasets, rejecting 5/100. With 24 total donors (22 × 6 cells and two rare singletons), it reported 99/100 and rejected 10/99 reported (10/100 generated). The oracle retained-donor procedure also rejected 10/99 in the sparse setting. The automatic robust RMS SE / empirical SD ratios were 0.961, 0.982 and 0.904 respectively. Thus the sparse setting shows undercoverage; this implementation is not a promise of nominal finite-sample coverage.

Three structural-zero fits were withheld because a reduced local root lost to a known full boundary competitor. One sparse fit was also unavailable to the oracle. Failures are retained in generated-dataset denominators. This small screen does not establish general post-selection or FDR validity. Earlier fixed-dispersion studies also showed finite-donor undercoverage (65/1000 null rejections in the 60-donor zero-variance development setting, including 54/900 in its fresh follow-up). No correction was selected by searching these p-values.
