# Full refitted-workflow validation for Bayesian plug-in REPAIR

Status: prospective plan, not completed statistical validation. Software checks
and native sampling smoke tests do not establish structural interval coverage.

## Questions and targets

Evaluate the current package's M-MOM estimator after estimating all shared
measurement parameters anew in every independent cohort. The primary questions
are structural bias, SE calibration and normal 95% confidence-interval coverage,
including prior sensitivity, imbalance, weak measurement information, failures
and regularization. The primary criterion is frequentist repeated-sampling
performance at fixed generating parameters and a fixed declared prior policy.
No posterior is claimed for plug-in REPAIR.

## Data generation and role-specific estimands

Generate independent people with numeric repeated predictor t and

    y_ij = 0.5 + 0.4*t_ij + b0_i + b1_i*t_ij + epsilon_ij,
    epsilon_ij ~ N(0,sigma^2).

Use t schedules on [-0.4,1.4], optionally jittered independently of effects.
Do not center every schedule at zero: off-diagonal intercept/slope measurement
error must be exercised. Fit `y ~ t + (1 + t | id)` with one correlated block,
estimating fixed effects, residual variance and the full random covariance.
Use zero-mean marginal Gaussian random coefficients and three separate targets:

| Target | Structural DGP | Truth and REPAIR formula |
|---|---|---|
| Random predictors, observed outcome | (b0,b1) ~ N(0,G); D = 1+.25*b0+.60*b1+.30*c+zeta; c independent N(0,1) | .60 for `D ~ re("(Intercept)") + re("t") + c` |
| Random outcome, observed predictor | b0 ~ N(0,.8^2); c ~ N(0,1); b1=.60*c+eta; b0 independent of c/eta | .60 for `re("t") ~ c` |
| Random outcome and predictor | b0 ~ N(0,.8^2); b1=.60*b0+eta | .60 for `re("t") ~ re("(Intercept)")` |

Reference settings for target 1: SDs (.8,.7), correlation .4, Var(zeta)=.8^2.
For target 2 use Var(eta)=.7^2, hence Var(b1)=.60^2+.7^2.
For target 3 use Var(eta)=.7^2-(.60^2)*(.8^2), giving marginal slope SD .7.
Keep all innovations independent unless their dependence is explicitly varied.
Repeat null versions by setting .60=0 and recalculating G. Save latent values
for oracle comparisons; never supply them to the candidate estimator.

The measurement fit marginalizes over structural variables. Do not add c or c*t
to step 1 in target 2: that changes the remaining random-slope deviation and the
estimand. Observed structural variables are measured without error.

## Staged design and resource budget

Freeze a machine-readable profile table, seeds, priors and code commit before
confirmation. Time a small pilot before launching a larger study. Avoid an
unbounded Cartesian product. Proposed anchor profiles for each target:

| Profile | People | Records/person | Residual SD | Change from reference |
|---|---:|---|---:|---|
| Reference | 150 | 8 | 1 | Correct likelihood, weak proper priors |
| Few people | 50 | 8 | 1 | Shared-parameter uncertainty |
| Many people | 500 | 8 | 1 | Large-cohort comparison |
| Sparse records | 150 | 4 | 1 | Weak individual information |
| Dense records | 150 | 20 | 1 | High information |
| Low precision | 150 | 8 | 2 | Low slope reliability |
| Unequal precision | 150 | Half 4, half 20 | 1 | Random assignment independent of effects |
| Correlated effects | 150 | 8 | 1 | Correlation .85 for target 1; vary other targets only with consistent structural truth |
| Regularizing priors | 150 | 8 | 1 | Prespecified stronger scale priors |
| Conflicting priors | 50 | 4 | 1 | Deliberately wrong fixed-effect/variance priors |
| Near-degenerate variation | 150 | 8 | 1 | Small slope SD where compatible with the target DGP |
| Rank failure | 150 | 8; 5% constant t | 1 | Identification and explicit fail/omit policy |

Report actual reliability distributions rather than treating residual SD as an
exact reliability setting. Keep the correctly specified core separate from
intentional prior conflict, rank failures, informative observation counts and
optional residual misspecification. AR(1)-generated errors fitted as independent
are a misspecification experiment, not support for Bayesian AR(1) extraction.
Omission changes the cohort and requires a separate target/selection analysis.

Run 20 outer replicates/profile/backend to estimate runtime/storage and sampling
failure, then 100 for engineering. Resolve code issues before confirmation and
use new held-out seed streams for at least 1,000 outer replicates in the frozen
core and selected stress profiles. At 95% coverage, 1,000 replicates yield MCSE
about .0069; 250 yield .0138. Increase to 2,000+ where greater precision is needed.
Freeze a runtime budget and report curtailed profiles and their MC uncertainty.

## Paired estimators and prior policies

Evaluate the same cohorts with:

* Bayesian plug-in REPAIR for each backend, `minimum` regularization primary;
  `none` and `step_up` are distinct sensitivity analyses with failures and
  selected constants retained.
* Frequentist lmer REPAIR with a fixed convergence policy.
* Oracle regression on the actual latent coefficients.
* Naive regression on likelihood-only scores, isolating ME correction.
* Regression on posterior means, describing its uncertainty method explicitly.
  Predictor calibration is a legitimate competitor; distinguish predictor-only
  use from random outcomes and joint latent roles.
* A joint Bayesian measurement/structural model with the same estimand in selected
  anchor profiles. Specify the structural factorization instead of multiplying
  unrelated distributions over the same random coefficients.

Declare priors on original outcome/predictor scales, including fixed means,
SD/residual scales and covariance priors. Disable automatic scaling where
possible; otherwise record adjusted priors per fit and analyze that adaptive
policy. brms LKJ/SD and rstanarm decov priors are not equivalent: use declared
native policies and do not attribute their differences to numerical errors.
Check numerical equivalence separately at identical supplied parameters.

## Complete outer replicate and sampling policy

1. Generate a fresh cohort and save data-generation/sampling seeds, profile ID
   and source hash before fitting. Use separate deterministic RNG streams.
2. Fit the Bayesian measurement model natively on that cohort. Reusing compiled
   Stan code is allowed; reusing nuisance estimates, draws or another cohort's
   posterior is not.
3. Start with four chains, 1,000 warmup and 1,000 retained iterations/chain,
   adapt_delta=.95 and max_treedepth=12. Record warnings, Rhat, bulk/tail ESS,
   MCSE, divergences, depth hits, BFMI when available, time and versions.
4. Create a new `repair_bayes_plugin()` snapshot and run ordinary `repair()` with
   explicit ID matching. Save coefficient tables, covariance, scopes, failed
   scores, joint error blocks and regularization diagnostics from current APIs.
5. Run paired comparators. Repeat 20 prespecified successful datasets with
   independent sampler seeds and longer chains to quantify MCMC variability
   relative to repeated-cohort structural variability.

Predeclare one diagnostic-triggered rescue (for example twice the iterations,
adapt_delta=.99 and max_treedepth=15). Retain both attempt results and diagnostics;
report original and rescued policies separately. Do not retry until a favorable
coefficient/interval appears or silently exclude unresolved sampling issues.

## Propagation experiments

Public refitting/bootstrap APIs reject Bayesian snapshots in this PR. A separate
research-only whole-person bootstrap must resample complete people, assign new
IDs to duplicated people, retain paired records/structural rows, **refit the
Bayesian model**, snapshot and rerun REPAIR. Preserve prior/scaling policies and
independent recorded sampling streams. Do not bootstrap posterior means alone,
hold nuisance parameters fixed, or append posterior variances to Psi_i.

Pilot its cost on selected outer datasets. Use at least 999 attempts for tail
assessment, report all failures/availability, and assess stability with larger B.
Coverage claims need many independent outer cohorts, not one cohort with many
bootstrap samples. Compare percentile and conditional intervals and the joint
model. Per-draw REPAIR quantiles or Rubin-rule pooling require separate derivation.

## Reporting and prospective decision criteria

Report attempts, usable posteriors, interval availability and failures by stage;
mean estimate, bias/MCSE, RMSE, empirical SD, mean SE, SE/SD ratio, 95% width and
coverage, null rejection and regularization frequencies for every profile/method.
Give successful-analysis conditional coverage and the all-attempt covered
fraction (missing intervals count as not covering) alongside availability. Show
paired-complete comparisons with their denominators; preserve all failed rows.

Screening flags: |coverage-.95| > max(.02,2*MCSE), |bias|/empirical-SD > .10,
SE/SD outside [.9,1.1], or interval availability below .99 in the correct-model
core. Flags trigger diagnosis, not automatic support claims. Report Monte Carlo
uncertainty and multiplicity across profiles. Stratify problems by information
and sampling status. Freeze any refined procedure before held-out confirmation.

Archive profile/replicate/attempt tables, exact priors, diagnostics, RNG state,
source manifests and session information; retain selected fits/data for independent
reproduction and keep large caches/executables out of Git. Review the argument
and independent confirmation results before strengthening inference claims.
Simulation-based calibration additionally checks genuine joint Bayesian
implementations; it does not turn the hybrid estimator into a posterior.

## Existing evidence and immediate software checks

Earlier Gaussian prototypes verified direct likelihood reconstruction against
conditional prior removal and ran simulations holding most shared parameters
fixed, using an earlier structural implementation. They are not the study above.
Current tests check the production scorer/estimator at matched shared parameters,
joint covariance, identification, diagnostics and actual native sampling:

```sh
Rscript -e 'devtools::test(filter = "bayes-plugin", stop_on_failure = TRUE)'
REPAIR_RUN_STAN_TESTS=true Rscript -e 'devtools::test(filter = "bayes-plugin-native", stop_on_failure = TRUE)'
# With CmdStan installed, set REPAIR_BRMS_BACKEND=cmdstanr for the second call.
```

Native tests cover shuffled/missing rows, unused grouping levels, independent
per-person OLS/covariance checks, E(sigma^2) versus E(sigma)^2, structural fitting,
prior/sampler provenance and rejection guards. They are software correctness
checks, not a Monte Carlo coverage study.
