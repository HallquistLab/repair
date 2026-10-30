# Influence-function versus bootstrap comparison protocol

Recorded before inspecting full-grid influence-function performance. Preserve all
previous simulation outputs. This extends the 2026-10-30 robust-refit screen;
it is a paired exploratory comparison on reused data, not independent replication.

## Estimator and local inference

Implement a separate repair_sandwich() result retaining the REPAIR point estimate.
Person-frequency derivatives propagate through the canonical preliminary/final
M-MOM moment corrections, structural variance, estimated case weights and active
regularization multiplier. No structural weights/corrections are frozen. Native
Gaussian ML or REML measurement influence is obtained from the exact score and
its numerically differentiated Jacobian. REML includes the global fixed-design
information determinant AND its case-frequency derivative. Simple ML scores at
a REML fit are not substituted. Measurement perturbations reconstruct complete
likelihood scores/error covariances, including target-contained nuisance effects.

Total person influence = frozen-measurement structural influence + propagated
measurement-parameter influence. Covariance is its empirical outer product/N^2,
including the signed cross-stage contribution. Return full, stage2, measurement
and cross contributions. Use normal 95% intervals, HC0 scaling. Validate central
sensitivities at half steps, one-sided consistency and unchanged active branches.
Full inference is unavailable for singular/weakly identified measurement models,
non-centered likelihood scores, or unstable/branch-crossing derivatives. Keep
stage2 availability separate. Stable local zero-variance/correction branches are
reported, not treated as a guarantee of uniform validity near a boundary.

Initial adapter: one independent-person Gaussian lmer, including contained
nuisance effects, offsets and prior weights; coinciding measurement/analysis
cohorts. nlme/multiple-model derivatives are rejected explicitly. Existing
conditional SEs, bootstrap failure policy and point estimator are unchanged.
Numerical differentiation is a local influence calculation, not resampling or
a claimed fully symbolic closed-form variance. No finite-sample HC correction is
silently added to this primary method.

## A: Existing noisy-predictor experiment

Reuse all 600 original measurement datasets / 1200 null/effect structural
analyses from ../robust_refit/results/main (100 outer replications/profile;
gamma=0/.40). Keep calibration, population scale, nuisance intercept effect,
paired streams and B499 ordinary/stratified/frozen bootstraps unchanged. Compare
conditional M-MOM Wald, stage2 influence Wald, full influence Wald, and recorded
full/frozen percentile bootstrap intervals, with actual availability denominators
and matched comparisons. Include signed covariance decomposition and boundary/
derivative rejection causes. The older OLS-outcome sandwich has no valid entry
for this noisy-predictor estimand; do not apply it by reversing regression roles.

## B: Matched random-slope-outcome experiment

Use the SAME measurement data and native first-stage fits/draws. Construct an
observed predictor x_gamma = gamma*b1/.9 + sqrt(1-gamma^2)*z, where z is the
shared original structural standard-normal draw independent of latent b and
measurement errors. Joint normality gives E[b1 | x_gamma]=.9*gamma*x_gamma.
Analyze re(time) ~ x_gamma, reporting the coefficient and SE divided by .9.
This preserves marginal G, measurement reliability, all schedules and the
true standardized slope gamma=0/.40. It is a distinct estimand/design from A.

Compare conditional/stage2/full M-MOM inference and person bootstrap with the
existing stacked_sandwich_for_corrected_scores() HC0-HC3 OLS-outcome method,
using the unmodified research helper and derivative backends with ML as its
original caller does. The REPAIR measurement fit retains REML. OLS slopes and
the HC0/2/3 slope covariance are invariant to the common first-stage centering;
first-stage fitting-criterion differences consequently cancel for that focal
OLS coefficient in this design. Record the criteria and boundary diagnostics.
Legacy finite-covariance outputs follow the existing comparator behavior; do
not apply a new IF boundary rule retroactively to that comparator. HC1 retains
its original total-stack-dimension scaling. All intervals use the original
legacy t_(N-2) reference. Full influence uses normal intervals as prespecified.

In balanced schedules, M-MOM's error-free-predictor equation has constant case
weights and agrees with OLS. In the unequal schedule M-MOM estimates precision
weights, while legacy OLS is unweighted: same population target but different
point estimators. Identify that difference in all tables and interpretation.

For each resample, the focal M-MOM slope is invariant to the common measurement
fixed-effect centering. Thus its full bootstrap can reuse the recorded draw map,
measurement convergence decisions and refitted sigma_hat from the previous
actual native fits, while recomputing B's structural estimator at every draw.
Use original fixed effects solely for that common centering. Do NOT claim this
reconstructs intercept inference. Require public complete-bootstrap prefix
checks (first two datasets/profile, both outcomes, B99) to match source draws,
focal coefficients and eligibility before interpreting this reuse. Use the
same default complete-vector failure gate; no omitted/replaced failed draws.
Record frozen-measurement cases intervals separately for diagnosis.

## Validation and summaries

Freeze numerical code before the main grid; use excluded replication 801 or
package teaching data for development/timing. Package checks must preserve all
existing numerical fixtures. Verify exact native marginal likelihood agreement,
independent objective-score differences, ML/REML case-score derivatives,
integer-frequency versus explicit duplication, covariance decomposition,
centered ID-aligned influences, unchanged original results and RNG, OLS HC0
special case, nested/weighted/offset scoring, unavailable-inference behavior,
and direct end-to-end likelihood reoptimization under fractional case shifts.

For A, coefficients/conditional SEs must reproduce the recorded robust run.
For B, validate the cached focal bootstrap against complete public calls, the
legacy comparator against its unchanged source, and balanced OLS identities.
Inspect numerical errors separately from intended analytic unavailability.
Summarize per arm/profile/gamma/method: eligible counts, coverage, null rejection,
non-null detection, width, bias, empirical SD and SE calibration; pair contrasts
on common eligible datasets. Do not count unavailable intervals as valid
non-rejections. Separate differences in covariance estimation, critical values,
point estimators, and eligibility. Evaluate B1999 on the same preselected first
three datasets/profile and derivative-step sensitivity on that fixed subset.
Archive code/library/reference hashes, sessions, raw analytic results and compact
review evidence. Do not infer that analytic SEs remove small/short-panel bias or
that finite-case agreement establishes nominal coverage universally.
