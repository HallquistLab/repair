# Person-level bootstrap and refitting: future-version design draft

Status: proposal for discussion, not an implemented or approved public API.
Recorded: 2026-10-29. Baseline: main at `150e56f` (package checks merged).
This document preserves the design discussed after the stable extractors and
website work. It does not commit the first release to bootstrap inference.
See [the release-readiness audit](release-readiness-audit.md) for priorities.

## Purpose and estimand

Implement a nonparametric bootstrap of whole people that reruns the complete
REPAIR analysis: measurement-model fitting, likelihood-only scoring,
sampling-error covariance assembly, and M-MOM structural estimation.

The purpose is to approximate sampling variability of the complete fitted
procedure, including measurement-model estimation. Existing conditional
model-based standard errors remain available for comparison. Successful
implementation does not establish calibrated coverage, especially near weak
identification, variance-component boundaries, or regularization transitions.
Bootstrap inference should initially be labeled experimental pending simulation
validation; it is not a general correction for misspecification, informative
missingness, dependence between people, or weak identification.

People are independent sampling units. Resampling their complete records varies
observed designs, cluster sizes, outcomes, and covariates together. The target is
repeated sampling of people with their observation histories, not inference
conditional on an unchanged set of covariates/designs. Retaining an originally
selected cohort conditions on that selection and its observed missingness pattern.
The procedure does not simulate new within-person trajectories.

## Proposed interface

These names and arguments are illustrative and do not exist in the current package.

```r
# Proposed API -- not executable with the current package.
refitted <- repair_refit(
  fit,
  measurement_data = list(growth = repeated_data),
  data = people,
  id = "Subject"
)

boot <- repair_bootstrap(
  fit,
  measurement_data = list(growth = repeated_data),
  data = people,
  id = "Subject",
  B = 1999,
  seed = 20261029
)

summary(boot)
repair_coefficients(boot)
confint(boot, level = .95, type = "percentile")
repair_diagnostics(boot, table = "replicates")
```

`repair_refit()` would replay an analysis on explicitly supplied replacement data.
`repair_bootstrap()` would use the same validated refitting machinery repeatedly.
Named measurement datasets must correspond to the named fitted models. Structural
and measurement ID column names may differ, but identity must be explicit.

Requiring data is intentional: current fitted objects save models, assembled
inputs, and an original call, but not necessarily a self-contained, reproducible
fitting recipe. The engine should normalize these into an internal replay
specification; it should not blindly re-evaluate calls in the global environment.
A callback extension can support unusual fits, subject to the same replay checks.

## One bootstrap replicate

1. Sample N people with replacement from the declared sampling cohort.
2. Copy each selected person's structural row and all effective repeated
   observations from every measurement dataset.
3. Assign each occurrence a new bootstrap ID, synchronized across all datasets.
4. Refit every referenced measurement model under its recorded specification.
5. Recompute likelihood-only scores and the complete sampling-error covariances.
6. Reassemble structural inputs and rerun REPAIR, recomputing data-dependent
   moment corrections, preliminary estimates, case variances, and weights.
7. Save coefficients and diagnostics, or a structured failure record.

Sampling source IDs A, C, A must produce three distinct bootstrap people, such as
b001, b002, b003. Both copies of A have the same source records but different
grouping identities. Otherwise they could be treated as one enlarged cluster.
Store a separate mapping from replicate/draw position/bootstrap ID to source ID;
never parse identity back out of display labels. Preserve character IDs such as
`"001"` and the original within-person observation order.

Do not resample observations within a person in the initial implementation.
Do not hold fitted scores, measurement covariance estimates, or structural weights
fixed while claiming to bootstrap the full analysis. Keep the chosen
regularization *policy* fixed while recomputing its sample-dependent quantities.

## Cohort contract

The first implementation should have a deliberately narrow sampling contract.

| Situation | Proposed behavior |
|---|---|
| Measurement models and structural analysis use the same people | Supported |
| People have unequal observation counts or different schedules | Supported; retain each full effective observation block |
| Multiple models use the same people | One shared person draw across every model |
| Measurement models include additional people outside the structural cohort | Stop; require an explicitly revised reference fit or a future sampling-plan interface |
| People share higher-level clusters or crossed random effects | Unsupported by the initial person-level bootstrap |
| Original analysis omitted people | Record that inference conditions on the retained cohort and observed missingness pattern |

The current package allows measurement-only people, so the extra-person case is
important. Silently dropping them before resampling changes measurement fitting
and therefore changes the estimator. Independently sampling each measurement
cohort would also destroy the required cross-model alignment. A future sampling
plan could support overlapping cohorts, measurement-only participants, and strata,
but its joint sampling design must be specified first.

If the user intentionally restricts the measurement cohort, create and report a
new reference fit using that restriction before bootstrapping. Do not attach its
bootstrap distribution to the old point estimate. In the initial design, the
identity replay and cohort audit must agree before any random draws begin.

## Replay specification and preflight checks

Before resampling, replay the original data with every person appearing once.
Check agreement with the supplied fit within documented numerical tolerances.
Compare model/effect identities, effective rows and IDs, score/covariance outputs,
structural coefficients, and regularization settings. A mismatch should stop with
an actionable explanation, not proceed with a changed estimator.

The replay specification must preserve:

- Measurement formulas, ML/REML choice, contrasts, offsets, weights, optimizer
  controls, and supported nlme correlation/variance structures.
- Effective observations, within-person ordering, time variables, and ID mappings.
- Structural formula, coefficient identities, regularization policy/custom alpha,
  and cross-model covariance policy.
- Declared behavior for sample-dependent preprocessing and external constants.

For a cases bootstrap, `lme4::refit()` alone is insufficient: that interface mainly
changes the response while retaining the model representation. Resampling people
changes grouping identities, designs, and possibly row counts. Rebuild the model
on the resampled data through a validated adapter or explicit callback.

Do not silently reset optimizer controls, switch ML/REML, reapply a stale subset
expression, read mutable global variables, or allow omission to change effective
rows unexpectedly. Formula dependencies that cannot be made reproducible require
an explicit fitting recipe. Stored fitted values alone are not enough.

Coefficient meanings must remain stable. Freeze the intended factor coding and
validate the resulting coefficient schema. If resampling removes a factor level
and makes a coefficient unidentifiable, mark that replicate invalid rather than
dropping or renaming it. The current structural assembler drops unused levels;
bootstrap integration therefore needs explicit schema validation.

Scores retain the package's random-deviation interpretation. Re-estimating fixed
effects can change their centering, particularly the structural intercept.
Document this estimand rather than comparing arbitrary coefficient positions.
Sample-dependent standardization/centering requires an explicit frozen-versus-
reestimated rule and matching original-data estimator; automatic standardization
is outside the initial scope.

## Multiple measurement models and cross-model covariance

Synchronized resampling preserves empirical associations across a person's
measurement processes. It does not establish independent scoring errors.

- `cross_model_error = "independent"` remains an explicit maintained assumption
  in every replicate. Within-model covariances are still recomputed and retained.
- Covariance blocks estimated from the data must be recomputed after refitting
  through a declared callback if their estimation variability is to be represented.
- Genuinely external covariance blocks may be copied and relabeled with sampled
  people, but resulting inference remains conditional on those fixed blocks.

Initially support the explicit independence policy for multiple models. Defer
estimated supplied-block covariance until its refitting contract is validated.
Do not silently reuse estimated original blocks and label that full uncertainty
propagation. Callback outputs need the current named coefficient and bootstrap-ID
coverage, symmetry/positive-semidefinite checks, and original-data replay checks.

## Failure and warning policy

B counts attempted resamples, not successful fits. Retain all attempts. Never
replace failed draws with new draws until a requested success count is reached.

| Event | Proposed treatment |
|---|---|
| Measurement optimizer fails | Invalid replicate; retain model and reason |
| Singular random-effect covariance | Record; retain if convergence, scoring, and structural estimation remain usable |
| Score calculation fails or structural design loses rank | Invalid replicate |
| M-MOM fails or returns nonfinite coefficients | Invalid replicate |
| Numerical retry | Retry the same resample under a declared fallback and record both attempts |

Do not equate every warning with a failure: distinguish optimizer failure,
singularity, missing variance-parameter uncertainty diagnostics, and informational
warnings under each backend. Do not silently drop people within a replicate to
rescue estimation. Extreme but valid estimates remain in the distribution; no
automatic trimming or winsorization.

Proposed conservative default: return all draws and diagnostics, but withhold
inferential summaries when replicates fail unless the user explicitly requests
summaries conditional on successful replicates. Such summaries must display
attempted/usable counts and failure reasons. A threshold such as 95% successful
replicates is not a validity guarantee. The precise opt-in interface and minimum
numerical requirements remain design decisions, not implemented defaults.

Use a common set of eligible coefficient vectors for a joint covariance matrix;
avoid unnoticed term-by-term deletion. Keep numerical requirements (enough usable
vectors to calculate an estimate) distinct from adequate Monte Carlo precision
and statistical validity. A result can remain useful for diagnostics even when
inference is unavailable.

## Result object and stable extraction

Return a separate bootstrap result object; leave the original REPAIR fit and its
conditional inference unchanged. Store the original estimate, replicate
coefficient matrix, resampling map, replay specification, seed information,
version metadata, and replicate/model diagnostics. Retaining all fitted models
should be optional because memory use can be substantial.

| Interface | Proposed meaning |
|---|---|
| `coef(boot)` | Original-data REPAIR estimate |
| `vcov(boot)` | Empirical covariance of eligible bootstrap coefficient vectors |
| `confint(boot)` | Explicitly selected bootstrap interval |
| `repair_coefficients(boot)` | Original estimates, bootstrap SEs/intervals, inference label, attempt and success counts |
| `repair_diagnostics(boot, table = "replicates")` | One row per attempted resample, with validity and failure stage |

Preserve the existing coefficient-table contract for ordinary REPAIR fits. The
new result class needs a separately documented schema and inference metadata;
resolve whether a dedicated bootstrap table extractor is clearer than extending
the current function. Do not place bootstrap SEs into the existing conditional
Wald statistic/p-value columns as though their interpretation were unchanged.

Replicate diagnostics should identify model, stage, optimizer/convergence status,
singularity, effective row/person counts, warnings, structural status, elapsed
time, and retries. Keep diagnostic records even when no coefficients are usable.

## Initial inference options

The bootstrap SE is the standard deviation of eligible replicate estimates; their
joint sample covariance supplies `vcov(boot)`. Do not add bootstrap variance to
the existing conditional variance, which would generally double-count variation.

Start with percentile intervals and empirical covariance, retaining the original
point estimate. Label intervals explicitly. Retain the original conditional
Wald output as a comparison, not as a bootstrap p-value calculation.

Defer automatic bias correction, bootstrap p-values, BCa, studentized intervals,
and simultaneous inference. BCa requires validated influence/jackknife machinery;
studentization needs a defensible variance estimate for each replicate. The
current conditional SE is not automatically a valid studentizing quantity for
full-pipeline uncertainty.

`B = 1999` is a proposed practical default, not a precision guarantee. Small B is
useful for exercising code; inference needs adequate tail resolution and checks
of Monte Carlo stability. Increasing B cannot repair model or bootstrap invalidity.

## Reproducibility and execution

Start with a serial engine. Assign random-number streams by replicate, independent
of worker scheduling, so a later parallel implementation can reproduce the same
draws. Preserve the caller's RNG state under a documented seed contract. Store
replicate IDs, draws/seed streams, data/specification fingerprints, and relevant
R/package versions. Never retain data in a log or checkpoint unexpectedly.

Optional checkpoints/resume and model retention should be deliberate interfaces.
Resuming must validate the same data, specification, seeds, and completed draws;
changing B must not silently relabel an already completed run. Parallel support
must not change valid draws, failure classifications, or results within the
same numerical environment. Avoid promising bitwise optimizer equivalence across
platforms or different dependency versions.

## Acceptance gates

### Software correctness

- Identity replay reproduces the original analysis.
- Independently written manual bootstrap loops agree for deterministic draws.
- Duplicate source IDs become distinct bootstrap people with synchronized IDs
  across every model; source records and observation order are retained.
- lmer and supported nlme controls, weights, offsets, contrasts, covariance
  structures, and effective rows are preserved.
- Coefficient schemas remain fixed; absent factor levels produce explicit failure.
- Multi-model pairing and covariance policies remain correct under permutations.
- Attempts, failures, retries, and inferential eligibility are counted correctly.
- Singular-but-usable fits are distinguished from nonconvergence/scoring failure.
- Existing manuscript/simulation fixtures remain unchanged and continue to pass.
- Seeded serial/resumed/parallel runs agree when those modes become available.

### Scientific validation

Compare current conditional Wald intervals and bootstrap intervals across
representative manuscript conditions: sample size, reliability, unbalanced
observations, correlated random coefficients, error in predictors/outcomes,
multiple processes, and weak corrected moments. Include well-behaved controls
and zero-measurement-error/known-parameter reference cases where appropriate.

Report point-estimate bias, empirical variability versus estimated SE, interval
width distributions, coverage, Type I error for a prespecified interval-based
test, failure/eligibility rates, and the fraction of attempted analyses producing
a covering interval. Distinguish within-bootstrap replicate failure from failure
of a whole simulated analysis to yield an interval. Report Monte Carlo uncertainty
and prespecify tolerances; do not select successful subsets after seeing results.

The existing algebra audit specifically motivates comparing full-pipeline cluster
bootstrap and stabilized Wald inference in difficult conditions. Passing software
tests is not evidence of nominal coverage, and no such new coverage study has been
performed by this design draft.

## Implementation sequence and unresolved decisions

1. Validated replay/refitting and cohort audit, with identity tests.
2. Serial person bootstrap, explicit failure records, percentile intervals, and
   a separate result/extractor contract; retain an experimental label.
3. Scientific validation, followed by reproducible parallel execution and
   carefully scoped additional sampling/covariance policies.

Before implementation, decide the initial lmer/nlme adapter scope, callback shape,
cohort restrictions, inference opt-in for incomplete replicates, output schema,
and computational defaults. Do not make this scope a prerequisite for a first
release that accurately reports conditional inference.

## Sources and related material

- Current contracts: [package plan](package_plan.Rmd),
  [extractors](extractor-contracts.md), [multiple models](multimodel-validation.md),
  and [reliability](reliability-validation.md).
- Loy and Korobova, [Bootstrapping Clustered Data in R using lmeresampler](https://journal.r-project.org/articles/RJ-2023-015/).
- lme4 [refit documentation](https://lme4.github.io/lme4/reference/refit.html).
- R boot [confidence interval documentation](https://stat.ethz.ch/R-manual/R-devel/library/boot/html/boot.ci.html).
- Companion repository: `vig_hallquist_2026/analysis/fuller_algebra_audit.Rmd`,
  especially its discussion of stabilization, coverage, and full cluster bootstrap.
  The package's existing fixture provenance identifies the frozen numerical
  reference revision; future coverage work needs its own versioned design/results.


## Refitting foundation implemented (2026-10-30)

The initial single-model replay foundation is now implemented as `repair_refit()`.
Its concrete API distinguishes explicit reference `measurement_data`/`data` from
`new_measurement_data`/`new_data`. It verifies original-data agreement first,
checks effective rows and schemas, supports data-only native lmer/lme recipes and
explicit fitting callbacks, and audits warnings and preserved control settings.
Deterministic source-ID draws can duplicate whole people with fresh target and
eligible nested nuisance identities, without making random draws or computing
bootstrap inference. See `refitting-validation.md` for current limits and the
implemented workflow. The illustrative bootstrap API and scientific validation
plan above remain future work.
