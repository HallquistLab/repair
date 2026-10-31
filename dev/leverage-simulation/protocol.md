# Fresh full-system leverage / delete-one comparison: frozen protocol

2026-10-31. Recorded before generating or inspecting this study's main results.
Do not replace the earlier robust-refit or influence-comparison datasets/results.

## Estimators

Primary: conditional M-MOM normal Wald, uncorrected full influence normal Wald,
full-system HC3-style normal Wald, full person delete-one jackknife normal Wald,
and full person percentile bootstrap. Stage2 influence and frozen-measurement
percentile intervals diagnose conditioning. Secondary: HC3-style and jackknife
t_(N-1) intervals, and normal Wald intervals using bootstrap SDs. No selected
finite-sample multiplier or bias subtraction is tuned to the results.

The leverage correction explicitly augments Gaussian measurement equations with
unweighted observed/error augmented moments, weighted observed/error moments,
and the final corrected coefficient equations. Preliminary coefficients,
residual variance and inverse case weights are functions of those nuisance
moments and individual error covariances; final correction multipliers depend
on the weighted moments. Differentiate the complete system, retaining REML's
exact deleted fixed-design determinant. Solve each deleted system's Jacobian at
the original parameter vector and sum uncentered coefficient-shift outer products
(HC3-style scaling). Original N and policy/alpha remain fixed in this analytic
correction. In OLS this must equal ordinary HC3. The system representation is
specified because nonlinear leverage corrections are representation dependent.
This is not claimed to be a derived CR2/Satterthwaite or uniform boundary result.

Pre-pilot normalization specification: include log(mean case weight) as an
additional nuisance parameter, with per-person equation w_i/mean(w)-1, and
normalize the final coefficient equation by that parameter. This cancels a
common estimated weight in the OLS special case, avoiding spurious leverage
from an irrelevant common variance scale. Weighted moment equations retain
their original unnormalized weights. Structurally zero error-moment entries
are excluded from the state, preserving the known error-free intercept/design.

Whole-person delete-one uses actual native optimization and the entire M-MOM
estimator on N-1 people. Finite-N corrections therefore update. Jackknife variance
is (N-1)/N times centered delete-one cross-products. Original point estimates
remain unchanged; jackknife bias is diagnostic, not subtracted. Every complete
deletion vector is required; no omission or successes-only jackknife inference.
Singular but converged native deletions follow the existing refit policy and are
recorded. Analytic inference retains its original nonsingularity/stability gate;
all deleted leverage systems must be usable. Do not silently return plain IF.

## Fresh design

Six profiles, each 100 new measurement datasets, following the prior noisy-
predictor design: N=50/n=10/reliability=.25/rho=.5; rho=-.5 variant;
N=30/n=3 variant; N=100/n=10 low reliability; N=100/n=10 reliability=.8;
N=50 with half n=3/half n=17 independently assigned, reliability=.25.
Use the unchanged calibration and data generator from robust_refit/engine.R.
Time is centered/sample-SD standardized, random coefficient SDs .9, Gaussian iid
measurement errors, measurement Y=.5*time+b0+b1*time+eps. Structural outcome is
Z=1.5+.4*b0+(gamma/.9)*b1+eta, Var(Z)=1, gamma=0/.40 with paired Gaussian noise.
Both recovered random coefficients are noisy predictors. Report focal coefficient,
SEs and endpoints multiplied by the population SD .9. No old OLS-outcome helper
is applied to this incompatible structural estimator; OLS HC3 is a package test.

NEW L'Ecuyer-CMRG master seed 202610319 (old study used 202610301), independent
outer streams, shared across profiles for paired precision. Bootstrap seed
340000000+replication. Main uses reps 1:100, B499. First three/profile extend
to B1999 and halve analytic derivative steps. Development/timing pilots use rep
801, B9 and are excluded. Save generated data and stream hashes, verify no input
dataset is identical to any earlier main dataset. No dataset regeneration or
replacement after a convergence failure. No optional stopping based on coverage.

## Validation and computation

Freeze package numerical source/library and study driver hashes before main.
Retain tested robust refit policy. Every bootstrap and delete-one measurement
fit is an actual native fit, with convergence/optimizer/singularity diagnostics.
Immutable per-person OLS/design-inverse inputs may be cached exactly for this
iid matching-design DGP; refitted fixed effects/sigma and complete structural
moments update for every draw/deletion. Verify all coefficient vectors and usable
flags against complete public APIs: bootstrap B99 on first TWO datasets/profile,
both outcomes (24 calls); delete-one on first ONE dataset/profile, both outcomes
(12 calls). No scalar focal-only gate and no failed vector repair/omission.

Before main, test canonical augmented-equation reproduction, full IF/augmented
HC0 agreement, OLS HC3 identity, exact ML/REML deletion scores, cross covariance,
scale/conditioning failures, explicit public N-1 refit identities, jackknife
covariance/bias, strict failure gates, schema/IDs, saved extraction and RNG.
Pilot refinements fix implementation errors only; never choose a correction by
main-grid coverage. Archive pilot failures and any protocol amendments separately.

## Reporting

Per profile/truth/method: planned and returned counts; available-case coverage,
null rejection/non-null detection with Wilson MC intervals; delivery per all
planned data; original and eligible point bias/SD; SE calibration and widths.
Compare methods on common eligible datasets with paired decision MC uncertainty.
Missing intervals are NA, never valid non-rejections. With 100 planned datasets,
coverage near .95 has MC SE about .022; no broad superiority claims from small
differences. Report numerical/boundary refusal reasons and leverage conditioning.

Separate smooth-interior structural branches from active corrections/floors,
and singular-original measurement strata. Report actual delete-one floor and
correction changes. Quantify bootstrap-SD/Wald versus percentile differences,
full/measurement/cross components, leverage versus exact deletion shifts and
prespecified B/half-step sensitivity. Conditioning on availability/branches can
select different samples; label denominators and avoid causal interpretations.
Preserve previous experiments, include verifiable compact evidence in the PR,
and treat fresh results as an exploratory validation screen, not universal proof.
