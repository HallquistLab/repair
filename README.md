
# repair <a href="https://hallquistlab.github.io/repair/"><img src="man/figures/logo.png" align="right" height="139" alt="repair website" /></a>

[![R-CMD-check](https://github.com/HallquistLab/repair/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/HallquistLab/repair/actions/workflows/R-CMD-check.yaml)

`repair` implements **Random-Effect Prediction Adjustment for Inference
in Regression (REPAIR)**, a two-stage method for estimating
relationships involving person- or cluster-specific processes measured
with repeated observations.

## Why use REPAIR?

Do people who improve more rapidly during treatment have better
subsequent outcomes? Does a person’s response to stress predict their
symptoms? Are rates of change measured in two different tasks related?
These questions concern individual differences in processes, such as
change over time or within-unit associations, and how those differences
relate to other variables.

Multilevel (mixed) models describe these differences using random
intercepts and slopes. A common next step is to extract the predicted
random effects and use them as predictors or outcomes in another
regression. However, these predictions are estimates, with precision
that can vary across people because of their observation counts,
measurement schedules, and predictor distributions. Empirical Bayes
predictions, often called BLUPs, also shrink uncertain individual values
toward the population mean. Shrinkage changes their variances and
covariances, and treating the predictions as perfectly measured can
distort regression coefficients and their standard errors. The
scientific goal is to learn about relationships among the underlying
processes, rather than just relationships among their predicted values.

REPAIR carries information about measurement precision into that second
regression. It first obtains **likelihood-only predictions** and their
unit-specific sampling-error covariances from fitted Gaussian mixed
models. These predictions remove shrinkage but remain noisy estimates.
It then uses a **modified method-of-moments (M-MOM)** regression to
account for their sampling error when random coefficients serve as
predictors, the outcome, or both.

This workflow is useful when you want to connect processes measured
across different tasks, settings, or timescales while keeping the MLMs
separate from the regression of interest. The same principle applies to
other clusters, such as students nested within classrooms. The
motivation follows Vig and Hallquist’s *REPAIR: Random-effect prediction
adjustment for inference in regression* manuscript; run
`citation("repair")` for the method and software references.

## Supported analyses

Version 1.0.0 is the first release candidate and has not yet been
submitted to CRAN. The supported structural workflow is an additive
regression with an intercept, using random coefficients from **fitted
`lme4::lmer()` or linear `nlme::lme()` models**. Select one target
grouping factor per model; lmer models can also have nuisance groups
contained within each target person. Random coefficients may be
predictors, the outcome, or both. Numeric, logical, and factor observed
predictors are supported. See the [workflow
guide](vignettes/repair-workflow.Rmd) for matching, omission,
interpretation, and current limits.

**Experimental Bayesian inputs:** explicitly wrap a narrowly supported sampled
Gaussian brms/rstanarm fit with `repair_bayes_plugin(fit)`. This uses posterior
summaries of shared measurement parameters with the existing estimator.
Structural inference remains conditional and model-based; step-1 uncertainty
is not propagated, full refitted-workflow coverage is unvalidated, and Bayesian
refitting/bootstrap APIs are unsupported. See the
[Bayesian plug-in guide](vignettes/bayesian-plugin.Rmd) for the precise scope.

## Package website

Browse the [package website](https://hallquistlab.github.io/repair/) for
a guided introduction, searchable function reference, executable
articles, and release notes. The site is rebuilt from `main`; it
documents the current package version.

## Install

``` r
# install.packages("remotes")
remotes::install_github("HallquistLab/repair")
```

Installation from source requires R \>= 4.3 and a C++20 compiler with
the R build toolchain. On Windows, use the Rtools version appropriate
for your R installation.

## Fit, inspect, and summarize

Fit the MLM using its usual `lmer()` or `lme()` interface, or use
`repair_mlm()` to specify fitting inside `repair()`. Its weights,
offsets, contrasts, estimation choices, and effective observations
remain under your control. Supply structural variables in a separate
data frame containing one row per analysis ID.

``` r
library(repair)
data(repair_example)  # Entirely simulated teaching data; see ?repair_example
growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
people <- repair_example$people

repair_effects(list(growth = growth))
#>    model  group        term n_levels
#> 1 growth person (Intercept)       60
#> 2 growth person         day       60
fit <- repair(
  re("day", model = "growth") ~ practice,
  models = list(growth = growth),
  data = people,
  id = "person"
)
summary(fit)
#> REPAIR structural regression
#> re("day", model = "growth") ~ practice
#> Units: 60 analyzed / 60 requested; 0 excluded
#> MLM: growth (lme4; group person; 0 IDs outside scope)
#> Cross-model prediction errors: not applicable (one model)
#> Predictions: likelihood-only random deviations
#>   model  group term    role
#>  growth person  day outcome
#>
#> Posterior reliability (structural scope; equal unit weights):
#>   model term    role      measure reliability n_units n_defined
#>  growth  day outcome     marginal       0.943      60        60
#>  growth  day outcome residualized       0.897      60        60
#>
#>             estimate      se  zstat      pval  ci_low ci_high
#> (Intercept)  -0.7895 0.11473 -6.882 5.905e-12 -1.0144 -0.5647
#> practice      0.1755 0.02142  8.193 2.548e-16  0.1335  0.2175
#>
#> Regularization (preliminary / final): minimum / minimum
#> Error-moment subtraction multipliers: 0.95 / 0.95
#> Largest analysis-weight fraction: 0.0179 ; weight-equivalent count: 59.8
#> Moment conditioning and unit weights: repair_diagnostics(fit).
#> Conditional on fitted MLM parameters; model-based normal SEs and intervals.
```

For an integrated call, supply an unfitted first-stage MLM
specification:

``` r
inline_fit <- repair(
  re(day) ~ practice,
  models = repair_mlm(score ~ day + (day | person),
    data = repair_example$repeated),
  data = people, id = "person"
)
repair_coefficients(inline_fit)
#>          term   estimate std.error statistic      p.value   conf.low  conf.high
#> 1 (Intercept) -0.7895289 0.1147252 -6.881913 5.905407e-12 -1.0143862 -0.5646716
#> 2    practice  0.1754854 0.0214190  8.192979 2.548385e-16  0.1335049  0.2174659
#>   conf.level
#> 1       0.95
#> 2       0.95
```

`repair_mlm()` describes the first-stage MLM; the outer formula
describes the structural model. `engine = "lme"` uses an explicit
fixed-effects formula and a separate `random` argument. Native controls
can be supplied in the specification. The helper defers fitting until
consumed, and each consuming call fits afresh; reuse
`inline_fit@mlm_fit` for further analyses without refitting. Fitted
models and specifications can be mixed in named lists. See `?repair_mlm`
and the [workflow guide](vignettes/repair-workflow.Rmd) for both engines
and replay limits.

The marker `re("day", model = "growth")` selects the random slope for
`day`. It represents a deviation from the fitted fixed slope. Use
`re("(Intercept)", model = "growth")` for the random intercept. When
there is only one matching effect, `model` may be omitted: `re("day")`
and `re(day)` are equivalent. A bare name is read literally, without
looking up a variable named `day`. Quoted names remain supported.
Multiple models can require `model`, and repeated coefficient names
across groups in one MLM can require `group`; use `repair_effects()` to
inspect the available names. Qualifier values remain quoted, for example
`re(day, model = "growth", group = "person")`.

``` r
coef(fit)
#> (Intercept)    practice
#>  -0.7895289   0.1754854
confint(fit)
#>                  2.5 %     97.5 %
#> (Intercept) -1.0143862 -0.5646716
#> practice     0.1335049  0.2174659
nobs(fit)  # people, not repeated observations
#> [1] 60
```

## Inspect what enters the estimator

``` r
analysis <- repair_data(
  re("day", model = "growth") ~ practice,
  models = list(growth = growth), data = people, id = "person"
)
analysis
#> REPAIR analysis data: 60 people, 1 predictor columns
#> Excluded: 0 of 60 requested IDs
#>   model  group term                       label    role
#>  growth person  day re("day", model = "growth") outcome
head(repair_predictions(analysis))
#>     id  model  group term                       label    role prediction
#> 1 P001 growth person  day re("day", model = "growth") outcome  0.4451909
#> 2 P002 growth person  day re("day", model = "growth") outcome  0.4105728
#> 3 P003 growth person  day re("day", model = "growth") outcome -0.4013259
#> 4 P004 growth person  day re("day", model = "growth") outcome  0.6289102
#> 5 P005 growth person  day re("day", model = "growth") outcome  0.6052267
#> 6 P006 growth person  day re("day", model = "growth") outcome -1.2120879
#>   std.error prediction.error.variance      scope analyzed status reason
#> 1 0.1401763                0.01964941 structural     TRUE     ok
#> 2 0.1401763                0.01964941 structural     TRUE     ok
#> 3 0.2268540                0.05146273 structural     TRUE     ok
#> 4 0.1618617                0.02619921 structural     TRUE     ok
#> 5 0.1401763                0.01964941 structural     TRUE     ok
#> 6 0.2268540                0.05146273 structural     TRUE     ok
head(repair_prediction_covariance(analysis))
#>     id model1 group1 term1                      label1   role1 model2 group2
#> 1 P001 growth person   day re("day", model = "growth") outcome growth person
#> 2 P002 growth person   day re("day", model = "growth") outcome growth person
#> 3 P003 growth person   day re("day", model = "growth") outcome growth person
#> 4 P004 growth person   day re("day", model = "growth") outcome growth person
#> 5 P005 growth person   day re("day", model = "growth") outcome growth person
#> 6 P006 growth person   day re("day", model = "growth") outcome growth person
#>   term2                      label2   role2 covariance       source
#> 1   day re("day", model = "growth") outcome 0.01964941 within_model
#> 2   day re("day", model = "growth") outcome 0.01964941 within_model
#> 3   day re("day", model = "growth") outcome 0.05146273 within_model
#> 4   day re("day", model = "growth") outcome 0.02619921 within_model
#> 5   day re("day", model = "growth") outcome 0.01964941 within_model
#> 6   day re("day", model = "growth") outcome 0.05146273 within_model
```

IDs determine alignment; data-frame row order never determines a match.
The unit-level data frame defines the requested scope. Missing/nonfinite
values, unmatched IDs, and failed predictions stop by default. Use
`na_action = "omit"` to explicitly exclude them with a warning and an
ID/reason audit trail.

Use public table extractors for reporting, joins, and export:

``` r
repair_coefficients(fit, conf.level = .95)
#>          term   estimate std.error statistic      p.value   conf.low  conf.high
#> 1 (Intercept) -0.7895289 0.1147252 -6.881913 5.905407e-12 -1.0143862 -0.5646716
#> 2    practice  0.1754854 0.0214190  8.192979 2.548385e-16  0.1335049  0.2174659
#>   conf.level
#> 1       0.95
#> 2       0.95
head(repair_predictions(fit))
#>     id  model  group term                       label    role prediction
#> 1 P001 growth person  day re("day", model = "growth") outcome  0.4451909
#> 2 P002 growth person  day re("day", model = "growth") outcome  0.4105728
#> 3 P003 growth person  day re("day", model = "growth") outcome -0.4013259
#> 4 P004 growth person  day re("day", model = "growth") outcome  0.6289102
#> 5 P005 growth person  day re("day", model = "growth") outcome  0.6052267
#> 6 P006 growth person  day re("day", model = "growth") outcome -1.2120879
#>   std.error prediction.error.variance      scope analyzed status reason
#> 1 0.1401763                0.01964941 structural     TRUE     ok
#> 2 0.1401763                0.01964941 structural     TRUE     ok
#> 3 0.2268540                0.05146273 structural     TRUE     ok
#> 4 0.1618617                0.02619921 structural     TRUE     ok
#> 5 0.1401763                0.01964941 structural     TRUE     ok
#> 6 0.2268540                0.05146273 structural     TRUE     ok
repair_diagnostics(fit, table = "excluded")
#> [1] id     reason
#> <0 rows> (or 0-length row.names)
```

`repair_prediction_covariance(fit)` preserves the assembled within- and
cross-model prediction-error covariance. It differs from `vcov(fit)`,
which describes structural coefficient uncertainty. These ordinary data
frames have explicit ID columns and the same schema for one or multiple
MLMs. See the [extractor guide](vignettes/extracting-results.Rmd) for
column contracts, scope selection, covariance reconstruction, and safe
joins. Internal S7 properties remain available for advanced inspection;
custom estimation uses `fit_mmom()`.

## Include occasions nested within people

A MLM such as `y ~ x + (1 + x | person) + (1 | person:day)` can supply
person coefficients to REPAIR. Use `group = "person"` in `re()` when
names such as `(Intercept)` occur in both groups. Prediction extraction
and posterior reliability include all nested nuisance covariance. Every
fitted nuisance level must belong to one target person; shared school,
item, or calendar-day effects remain rejected, even after structural
omissions. See the executable [nested-nuisance
guide](vignettes/nested-nuisance.Rmd) for discovery, fitting, multiple
nuisance groups, and independent numerical references.

## Verify and refit an analysis

`repair_refit()` first replays an analysis using explicit reference
measurement and unit-level data. Separate `new_mlm_data` and `new_data`
arguments supply changes only after the reference checks pass. One
supported lmer or lme MLM can be rebuilt, retaining the structural
equation and controls:

``` r
replayed <- repair_refit(fit, mlm_data = list(growth = repair_example$repeated),
  data = people)
summary(replayed)$refit$checks
#>       stage n_rows n_requested n_analyzed           controls verified
#> 1 reference    586          60         60 recorded_optimizer     TRUE
```

The [refitting guide](vignettes/refitting.Rmd) covers native fitting
recipes, callbacks, factor-schema checks, revised cohorts, and
deterministic whole-person draws with fresh target/nuisance IDs. This
provides deterministic replay for the experimental bootstrap below;
ordinary refit SEs remain conditional model-based SEs.

## Adaptive step-up regularization

`regularization = "step_up"` enables the manuscript’s M-MOMS search. It
selects the preliminary and final alphas separately and records the
search, conditioning criteria, and any fallback to minimum alpha. Refits
and bootstrap samples rerun the selection. See the [step-up
guide](vignettes/step-up.Rmd).

``` r
adaptive <- repair(practice ~ re("day"), growth, people, "person",
  regularization = "step_up")
repair_diagnostics(adaptive, table = "regularization")
#>          step  policy alpha subtraction_multiplier regularized
#> 1 preliminary step_up     3                   0.95        TRUE
#> 2       final step_up     3                   0.95        TRUE
```

## Experimental person-level bootstrap

`repair_bootstrap()` samples whole people with replacement and reruns
`repair_refit()`, including measurement fitting, scoring covariance, and
structural estimation. Repeated source people receive fresh target and
nuisance IDs.

``` r
# Five attempts exercise the API; use many more and assess Monte Carlo stability.
boot <- repair_bootstrap(fit, repair_example$repeated, people, B = 5, seed = 41)
repair_bootstrap_coefficients(boot, conf.level = .90)
#>          term   estimate  std.error   conf.low  conf.high conf.level     method
#> 1 (Intercept) -0.7895289 0.15787603 -0.9447660 -0.6092946        0.9 percentile
#> 2    practice  0.1754854 0.03381545  0.1337869  0.2111553        0.9 percentile
#>   attempted usable failed conditional_on_success                     inference
#> 1         5      5      0                  FALSE Experimental person bootstrap
#> 2         5      5      0                  FALSE Experimental person bootstrap
```

Results are separate from the original fit. `confint(boot)` gives
percentile intervals; `vcov(boot)` gives empirical covariance. Every
attempt and source-ID map is retained. Failed samples are never
replaced, and inference is withheld after any failure unless explicitly
requested conditional on successful replicates. The [bootstrap
guide](vignettes/person-bootstrap.Rmd) explains diagnostics, RNG
restoration, supported single-model cohorts, and saved results. This
feature is experimental: broader independent coverage validation remains
outstanding.

For repeated measurement convergence problems, `repair_lmer_refitter()`
provides an opt-in, deterministic fitting recipe with likelihood-checked
BOBYQA and bounded optimizer retries. Use it to rebuild the original
measurement model and REPAIR analysis first; saved models then carry the
same strategy into reference and bootstrap refits. Inspect the retained
internal attempts with
`repair_bootstrap_diagnostics(boot)$optimization`. The default
failed-sample gate is unchanged. A six-profile simulation screen
improved default interval availability from 1–53% to 95–100%, with mixed
coverage gains and persistent small/short-panel undercoverage. See the
[audit and simulation
results](https://github.com/HallquistLab/repair/blob/2dda7622ae7dceb24d726ff2c3aeb56db390083f/dev/robust-refit-audit.md)
and `vignette("person-bootstrap")` for scope, usage and limitations.

## Connect two measurement processes

Use model-qualified references when an equation draws on multiple fits:

``` r
fit <- repair(
  re("time", model = "symptoms") ~ re("time", model = "stress") + age,
  models = list(symptoms = symptoms_fit, stress = stress_fit),
  data = people, id = "person",
  cross_model_error = "independent"
)
```

This explicitly assumes zero **cross-model prediction-error
covariance**; latent processes can still be correlated. Separate model
fits alone do not justify the assumption. If the prediction errors are
correlated, supply per-person covariance blocks instead. Omitting the
policy is an error for multiple referenced models. Within-model
covariance is always retained. See the executable [multiple-model
guide](vignettes/multiple-models.Rmd) for design guidance,
`repair_cross_covariance()` to propagate explicitly assumed observation
error covariance, and `repair_sensitivity()` to compare assumptions with
estimates, SEs, intervals, and regularization diagnostics. Automatic
covariance estimation is not included. The guide also reproduces a
dual-process analysis against the paper’s simulation code.

## Reliability and diagnostics

Assess the selected random coefficients using the manuscript’s posterior
reliability definitions. The default uses analyzed units with equal
weights:

``` r
repair_reliability(fit)[, c("term", "role", "measure", "reliability", "n_units", "n_defined")]
#>   term    role      measure reliability n_units n_defined
#> 1  day outcome     marginal   0.9433813      60        60
#> 2  day outcome residualized   0.8967932      60        60
diagnostics <- repair_diagnostics(fit)
diagnostics$weight_summary
#>   n_units max_fraction weight_equivalent_n
#> 1      60   0.01788085            59.82168
diagnostics$conditioning
#>          step dimension min_eigenvalue max_eigenvalue condition_number
#> 1 preliminary         2       17.52838       1747.322         99.68531
#> 2       final         2       74.13875       7551.600        101.85768
```

``` r
person_reliability <- repair_reliability(fit, detail = "unit")
plot(person_reliability)
plot(person_reliability, type = "observations")
```

Posterior reliability is distinct from the likelihood-only
prediction-error covariance used by the estimator. Undefined values
retain reasons and explicit denominators; diagnostics do not impose
automatic adequacy cutoffs. See the [reliability and diagnostics
guide](vignettes/reliability-diagnostics.Rmd) for marginal versus
intercept-residualized reliability, scope selection, and interpretation
of numerical and weight diagnostics.

## Support, troubleshooting, and citation

The [support guide](vignettes/support-and-troubleshooting.Rmd)
consolidates the supported model/inference matrix, a complete
`repair_example` analysis, common diagnostic messages, and saved-fit
guidance. The teaching data include both repeated observations and
meaningful synthetic unit-level variables.

Run `citation("repair")` for separate software and method references.
Zachary A. Vig and Michael N. Hallquist are authors and co-maintainers;
Zach is the primary R package maintainer contact. The package version is
1.0.0; see the release notes for changes.

## Interpretation and limits

- Standard errors are model-based and conditional on estimated MLM
  quantities. They use the modified final predictor moment, consistent
  with the manuscript and primary simulation estimator.
  Measurement-model parameter uncertainty is not propagated; intervals
  and tests use a normal reference.
- `regularization = "minimum"` is the default. `"none"` disables
  stabilization but **still corrects for measurement error**. Custom
  settings require `alpha`.
- Observed structural variables are assumed measured without error. The
  complete error covariance is preserved when random coefficients appear
  on both sides.
- Interactions or nonlinear transformations involving random
  coefficients, models without an intercept, and joint multiple-equation
  fitting are not implemented. Unsupported formulas fail explicitly.
- Single-group `nlme` models support independent/AR(1)/continuous-time
  AR(1) correlation and constant/`varIdent`/`varPower` residual
  variance. Keep the data in the fit (`keep.data = TRUE`). See the [nlme
  reproduction guide](vignettes/nlme-workflow.Rmd).
- lmer nuisance groups must be contained within the selected target
  group on all effective fitted rows. Shared nuisance effects remain
  rejected because their cross-person prediction errors cannot enter
  this structural estimator. Lower-level prediction extraction of shared
  structures remains experimental. Multilevel nlme fits are rejected.

## Development and reproducibility

The [design plan](dev/package_plan.Rmd) describes the implemented
contracts and extension sequence. The [audit and
decisions](dev/design-audit.md) trace findings to corrections and tests.
The test suite includes frozen REPAIR simulation and independent
algebra-audit fixtures for both R and C++ backends, plus three
single-model nlme and two dual-process reproduction datasets; see
[fixture provenance](tests/testthat/fixtures/README.md).

Pull requests and pushes to main run package checks on Linux, macOS, and
Windows, with additional Linux checks for previous, development, and
minimum-supported R. Tests, examples, and all vignettes run with
suggested dependencies installed; errors and warnings fail CI. See the
[CI guide](dev/ci-checks.md) for the matrix, check artifacts, local
reproduction, and the `Package checks` branch-protection target.

``` r
# From a source checkout:
testthat::test_local()
```
