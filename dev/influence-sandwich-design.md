# Full influence-function inference: implementation audit and workflow

## Objective and scope

Implement full empirical person-level covariance for the actual REPAIR M-MOM
pipeline, rather than transplanting the research OLS-outcome sandwich. Retain
ordinary point estimates and model-based SEs. The new repair_sandwich() result
keeps full and frozen-measurement influences/covariances separate, including
signed cross-stage terms and explicit reasons for unavailable inference.

Initial measurement adapter: a single Gaussian lmer with independent target
people, target-contained nuisance effects, prior weights/offsets and aligned
cohorts. nlme and multiple-model derivatives are not implemented and are rejected.
No covariance independence assumption substitutes for missing cross-model code.

## Mathematical construction

For fractional person frequencies f summing to N, T(f,psi) is the exact canonical
finite-sample structural point estimator at measurement parameter psi. Its
frequency influence is computed along f_j=1+epsilon*(I[j=i]-1/N):

    U_i = N * d T(f(epsilon;i), psi_hat) / d epsilon.

The user-selected regularization policy/custom alpha and N remain fixed,
matching N-person cases sampling. Data-dependent multipliers, preliminary
coefficients, structural variance, case weights and weighted moments all update.

For marginal Gaussian ML, S_i is the cluster likelihood score. For REML, use the
joint beta/covariance objective

    sum_i f_i log p(y_i; beta,V_i) - .5 log|sum_i f_i X_i' V_i^-1 X_i|
       + p/2 log(2*pi).

The derivative with respect to beta gives its GLS equation; covariance scores
include the restricted determinant. Case-frequency derivatives also differentiate
that global determinant, not just the additive ML likelihood contributions.
Their centered values are S_i^*. With H=-(1/N)*d score_total/d psi,

    IF_psi,i = H^-1 S_i^*,
    D = d T(1,psi)/d psi,
    IF_total,i = U_i + D IF_psi,i.

Return V_full = sum_i IF_total,i IF_total,i' / N^2. The decomposition is
V_stage2 + V_measurement + V_cross, with V_cross the two signed empirical
cross-products. There is no assumption of independence between stages and no
rule that the total correction must increase every SE. Normal HC0/Wald inference
is primary; no automatic finite-sample correction is labeled as established.

## Implementation choices and checks

1. Preserve the robust-refit PR (#20) and its recorded bootstrap simulations;
   work on codex/influence-sandwich. Freeze numerical code before the comparison.
2. Add private case frequencies to normalized M-MOM inputs and route every
   point-estimation step through the same functions used by fit_mmom(). All
   public calls retain unit frequencies. Fractional frequencies are an internal
   point-functional interface; they do not define a new public weighted SE API.
3. Build person marginal covariances from native lme4 block designs/Cholesky
   parameters. Log-transform positive diagonal theta and residual SD; off-diagonal
   theta remains unconstrained. Target-contained nuisance covariance remains in
   likelihood scoring; target random-effect prior covariance is excluded from
   the likelihood-only score calculation, exactly as in ordinary REPAIR.
4. Evaluate analytic Gaussian scores and case derivatives. Numerically
   differentiate those scores for observed sensitivity, and differentiate the
   complete structural pipeline for person/parameter influences. This is not a
   fully symbolic closed-form implementation or resampling variance.
5. Check centered native ML/REML scores, marginal likelihood agreement, scaled
   positive-definite curvature, derivative step halving/one-sided agreement and
   unchanged active branches. Refuse full inference at singular/weakly identified
   measurement fits or unstable sensitivities. Stable local structural floors
   remain reported; they do not imply uniform regularity near a boundary.
6. Preserve original coefficients/conditional covariance and caller RNG. Keep
   ID-aligned influences, scores, bread, structural sensitivities, criterion,
   versions and signed decomposition. Failed full inference stays inspectable;
   accessors never quietly substitute stage2 or ordinary conditional inference.
7. Verify using independent finite-difference likelihood gradients, ML/REML
   frequency-score perturbations, explicit integer person duplication, direct
   likelihood reoptimization under fractional perturbations, the OLS HC0 special
   case, nested weighted/offset scoring, saved extraction and failure handling.
   Existing point/SE reference fixtures must continue to pass.

## Simulation comparison

Use all six profiles/100 outer replications from the existing robust screen.
Arm A retains the original two noisy-predictor partial regression, comparing
conditional, stage2/full influences and frozen/full cases bootstrap intervals.
The older OLS-outcome sandwich is not applied to this incompatible estimand.

Arm B adds an observed x_gamma = gamma*b1/.9 + sqrt(1-gamma^2)*z using the same
measurement data and independent original structural standard-normal z. Then
E[b1|x_gamma]=.9*gamma*x_gamma. Analyze re(time)~x_gamma, report coefficients/SEs
on the population-SD scale by dividing by .9, and include the unchanged research
stacked OLS HC0-HC3 helper as a comparator. Balanced schedules give an exact
point-estimator/HC0 anchor. Unequal schedules preserve the same population
estimand with distinct weighted M-MOM and unweighted legacy OLS estimators.
The legacy helper is evaluated with ML as its original caller; REPAIR retains
REML. Its focal OLS slope is invariant to common fixed-effect centering. Retain
the legacy t_(N-2) critical value and HC1 convention; separate those differences
from cross-stage variance propagation.

Reuse recorded measurement draws, convergence decisions and refitted residual
SDs to reconstruct Arm B's focal bootstrap only. In this design the focal slope
is invariant to common score centering; intercept bootstrap inference is not
reconstructed this way. Re-estimate every structural M-MOM equation for every
draw, then verify the focal coefficients and eligibility against complete public
bootstrap calls on fixed prefixes. This reuses actual native optimizations,
not idealized converged fits. Legacy covariance, numerical status, actual
eligibility and derivative-step/B sensitivity remain separately auditable.

The executed protocol, scripts and results live in the research workspace under
vig_hallquist_2026/analysis/repair_bootstrap_pilot/influence_comparison/. The study
is a paired exploratory screen on reused data, not independent confirmatory
coverage evidence. Report denominators, matched uncertainty, bias, covariance
contributions and all analytic-unavailability reasons. Do not count missing
intervals as valid non-rejections or infer that an SE correction fixes point bias.

## Executed findings (2026-10-30)

Numerical baseline: 558bdcead1c15cedd1e332aed6f01e86d49f3721, built on robust-refit
head 7770874. Later pointwise formula/weight scope guards do not change this
study's time models. Source/library hashes and all old checkpoints are retained.
See [the review bundle](influence-simulation/README.md) for independently
verifiable endpoints, detailed report, drivers, diagnostics and source manifests.

The 600 measurement datasets produced 2,400 structural analyses and 17,200
method-specific interval records (172 cells, each 100 planned datasets). Full
influence inference was available for 1,351/2,400; stage2 for all 2,400. Of 1,049
full refusals, 1,048 were singular native covariance boundaries. One was an
unstable parameter derivative (outcome arm, low_larger, rep 68, gamma .40,
coordinate 6). There were no driver/public-API/reconstruction errors. Low-profile
full availability was 37–66/100 versus bootstrap 95–98/100; high-control
availability was 100/100 for both. Different eligibility must not be mistaken
for a same-sample comparison.

| Noisy-predictor profile, truth .40 | Conditional coverage | Full IF coverage (n) | Bootstrap coverage (n) |
| --- | ---: | ---: | ---: |
| Low, balanced | 95% | 78.4% (37) | 93.9% (98) |
| Negative correlation | 90% | 83.3% (48) | 93.7% (95) |
| Small, short | 78% | 73.2% (41) | 79.4% (97) |
| Low, larger N | 99% | 81.8% (66) | 95.8% (96) |
| High reliability | 94% | 92.0% (100) | 91.0% (100) |
| Unequal information | 85% | 82.6% (46) | 91.8% (98) |

All conditional intervals above were available (n=100). These are available-case
coverage rates; the full report supplies Wilson MC intervals and paired contrasts.
On common eligible noisy-predictor datasets, IF minus bootstrap coverage at truth
.40 was −.167 (MC bounds −.290,−.043; n=36) for low balanced,
−.113 (−.192,−.033; n=62) for low larger, and −.109 (−.200,−.018; n=46)
for unequal information. These exploratory contrasts are not multiplicity-adjusted.
The small/short point bias was −.207; none of the SE methods removes it.

The balanced OLS-outcome identity passed all 1,000 paired comparisons: max point
gap 5.55e-16 and stage2/legacy HC0 variance gap 6.09e-8. Full and stage2 focal
variance cancel to the same value in these balanced designs. Legacy HC3/t
nonzero-effect coverage was 95–97% with 100 finite intervals/profile, but its
unweighted estimator, HC adjustment and t critical value must not be attributed
to a general cross-stage improvement. The helper is inapplicable to arm A.

Full/stage2 focal variance ratios across available noisy-predictor analyses ranged
roughly .79–3.62. Signed cross terms can offset the positive measurement term.
The median change is modest in most profiles; the new SE is not uniformly larger.
Many low-information M-MOM fits lie on active moment corrections and a forced
residual-variance floor. A stable local derivative does not capture global branch
switching or remove point bias, helping explain why correct local propagation
is insufficient for this screen's finite-sample coverage.

Validation: 24 complete public B99 outcome bootstraps (2,376 attempts) matched
draws, full-vector eligibility and focal estimates (max gap 6.66e-16). Independent
recalculation verified all 172 rates from endpoints. Covariance decomposition max
relative error was 1.45e-16; ID alignment/centering and original noisy-predictor
bootstrap identity passed. On the fixed first-three/profile subset, halving the
step changed no availability decisions and changed matched full covariance by
at most 1.35e-7 relatively. All 156 fixed-subset B499/B1999 intervals were
available; one coverage decision changed to covered, none to uncovered. B499
draws/estimates are exact prefixes. These subset checks are not a full-grid
stability or precise tail-coverage guarantee.

Decision: expose the complete influence method as experimental and diagnostic,
with inspectable failures and unchanged existing inference. Retain person
bootstrap as the practical alternative for difficult settings, and clearly
document that neither method repairs short-panel bias universally. Further
analytic finite-sample correction work needs a new, separately specified study;
do not choose an HC multiplier after inspecting these results.

## Final package verification

The source-package check ran all 2,099 assertions with zero test failures,
warnings or skips; examples and all vignette rebuilds passed. R CMD check had
zero errors/notes and one installation warning from the external R 4.5.2 header
R_ext/Boolean.h (unknown Clang warning group `-Wfixed-enum-extension`). The
warning is retained in the bundle's check/install logs; no warning suppression
or check-policy exception was introduced. The pkgdown build passed validation
of 57 HTML pages, local links/assets and search destinations.

The later pointwise-recipe scope guard was checked through the public API on
the fixed first-three/profile subset (72 arm/effect calls), with identical
availability, IDs and covariances to the frozen numerical study (relative
tolerance 1e-10). Recomputed `scale()` and sample-wide weight recipes have
explicit rejection tests. The package remains opt-in and experimental.
