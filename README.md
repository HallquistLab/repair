
# repair

[![R-CMD-check](https://github.com/zacharyvig/repair/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/zacharyvig/repair/actions/workflows/R-CMD-check.yaml)

`repair` implements REPAIR: structural regression involving random
coefficients from a Gaussian mixed model, with a modified
method-of-moments correction for sampling error in the estimated
coefficients. It uses **likelihood-only scores** and their
sampling-error covariance, rather than treating shrunken BLUPs as
error-free observed variables.

This is a development package. The supported structural workflow is an
additive regression with an intercept, using random coefficients from
**fitted `lme4::lmer()` or linear `nlme::lme()` models, each with one
grouping factor**. Random coefficients may be predictors, the outcome,
or both. Numeric, logical, and factor observed predictors are supported.
See the [workflow guide](vignettes/repair-workflow.Rmd) for matching,
omission, interpretation, and current limits.

## Install

``` r
# install.packages("remotes")
remotes::install_github("zacharyvig/repair")
```

Installation from source requires R \>= 4.3 and a C++20 compiler with
the R build toolchain. On Windows, use the Rtools version appropriate
for your R installation.

## Fit, inspect, and summarize

Fit the measurement model using its usual `lmer()` or `lme()` interface.
Its weights, offsets, contrasts, estimation choices, and effective
observations remain under your control. Supply structural variables in a
separate data frame containing one row per analysis ID.

``` r
library(repair)
growth <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)

# Illustrative person-level predictor, not an actual sleepstudy measurement.
people <- data.frame(
  Subject = levels(lme4::sleepstudy$Subject),
  observed = seq_len(18)
)

repair_effects(list(growth = growth))
#>    model   group        term n_levels
#> 1 growth Subject (Intercept)       18
#> 2 growth Subject        Days       18
fit <- repair(
  re("Days", model = "growth") ~ observed,
  models = list(growth = growth),
  data = people,
  id = "Subject"
)
summary(fit)
#> REPAIR structural regression
#> re("Days", model = "growth") ~ observed
#> People: 18 analyzed / 18 requested; 0 excluded
#> Measurement model: growth (lme4; group Subject; 0 IDs outside cohort)
#> Cross-model scoring errors: not applicable (one model)
#> Scores: likelihood-only random deviations
#>   model   group term    role
#>  growth Subject Days outcome
#>
#> Posterior reliability (analyzed cohort; equal person weights):
#>   model term    role      measure reliability n_people n_defined
#>  growth Days outcome     marginal       0.849       18        18
#>  growth Days outcome residualized       0.828       18        18
#>
#>             estimate    se   zstat   pval  ci_low ci_high
#> (Intercept)  -2.9721 3.215 -0.9245 0.3552 -9.2729  3.3287
#> observed      0.3129 0.297  1.0534 0.2922 -0.2692  0.8949
#>
#> Regularization (preliminary / final): minimum / minimum
#> Error-moment subtraction multipliers: 0.8333 / 0.8333
#> Largest analysis-weight fraction: 0.0556 ; weight-equivalent count: 18
#> Moment conditioning and individual weights: repair_diagnostics(fit).
#> Conditional on fitted measurement-model parameters; model-based normal SEs and intervals.
```

The marker `re("Days", model = "growth")` selects the random slope for
`Days`. It represents a deviation from the fitted fixed slope. Use
`re("(Intercept)", model = "growth")` for the random intercept. When
there is only one matching effect, `model` may be omitted.

``` r
coef(fit)
#> (Intercept)    observed
#>  -2.9721078   0.3128535
confint(fit)
#>                 2.5 %    97.5 %
#> (Intercept) -9.272872 3.3286569
#> observed    -0.269238 0.8949449
nobs(fit)  # people, not repeated observations
#> [1] 18
```

## Inspect what enters the estimator

``` r
analysis <- repair_data(
  re("Days", model = "growth") ~ observed,
  models = list(growth = growth), data = people, id = "Subject"
)
analysis
#> REPAIR analysis data: 18 people, 1 predictor columns
#> Excluded: 0 of 18 requested IDs
#>   model   group term                        label    role
#>  growth Subject Days re("Days", model = "growth") outcome
head(repair_scores(analysis))
#>    id  model   group term                        label    role      score
#> 1 308 growth Subject Days re("Days", model = "growth") outcome 11.2974165
#> 2 309 growth Subject Days re("Days", model = "growth") outcome -8.2055005
#> 3 310 growth Subject Days re("Days", model = "growth") outcome -4.3523872
#> 4 330 growth Subject Days re("Days", model = "growth") outcome -7.4592132
#> 5 331 growth Subject Days re("Days", model = "growth") outcome -5.2012672
#> 6 332 growth Subject Days re("Days", model = "growth") outcome -0.9005181
#>   std.error error.variance   cohort analyzed status reason
#> 1  2.817564       7.938667 analyzed     TRUE     ok
#> 2  2.817564       7.938667 analyzed     TRUE     ok
#> 3  2.817564       7.938667 analyzed     TRUE     ok
#> 4  2.817564       7.938667 analyzed     TRUE     ok
#> 5  2.817564       7.938667 analyzed     TRUE     ok
#> 6  2.817564       7.938667 analyzed     TRUE     ok
head(repair_score_covariance(analysis))
#>    id model1  group1 term1                       label1   role1 model2  group2
#> 1 308 growth Subject  Days re("Days", model = "growth") outcome growth Subject
#> 2 309 growth Subject  Days re("Days", model = "growth") outcome growth Subject
#> 3 310 growth Subject  Days re("Days", model = "growth") outcome growth Subject
#> 4 330 growth Subject  Days re("Days", model = "growth") outcome growth Subject
#> 5 331 growth Subject  Days re("Days", model = "growth") outcome growth Subject
#> 6 332 growth Subject  Days re("Days", model = "growth") outcome growth Subject
#>   term2                       label2   role2 covariance       source
#> 1  Days re("Days", model = "growth") outcome   7.938667 within_model
#> 2  Days re("Days", model = "growth") outcome   7.938667 within_model
#> 3  Days re("Days", model = "growth") outcome   7.938667 within_model
#> 4  Days re("Days", model = "growth") outcome   7.938667 within_model
#> 5  Days re("Days", model = "growth") outcome   7.938667 within_model
#> 6  Days re("Days", model = "growth") outcome   7.938667 within_model
```

IDs determine alignment; data-frame row order never determines a match.
The person-level data frame defines the requested cohort.
Missing/nonfinite values, unmatched IDs, and failed scores stop by
default. Use `na_action = "omit"` to explicitly exclude them with a
warning and an ID/reason audit trail.

Use public table extractors for reporting, joins, and export:

``` r
repair_coefficients(fit, conf.level = .95)
#>          term   estimate std.error  statistic   p.value  conf.low conf.high
#> 1 (Intercept) -2.9721078 3.2147349 -0.9245266 0.3552122 -9.272872 3.3286569
#> 2    observed  0.3128535 0.2969909  1.0534109 0.2921527 -0.269238 0.8949449
#>   conf.level
#> 1       0.95
#> 2       0.95
head(repair_scores(fit))
#>    id  model   group term                        label    role      score
#> 1 308 growth Subject Days re("Days", model = "growth") outcome 11.2974165
#> 2 309 growth Subject Days re("Days", model = "growth") outcome -8.2055005
#> 3 310 growth Subject Days re("Days", model = "growth") outcome -4.3523872
#> 4 330 growth Subject Days re("Days", model = "growth") outcome -7.4592132
#> 5 331 growth Subject Days re("Days", model = "growth") outcome -5.2012672
#> 6 332 growth Subject Days re("Days", model = "growth") outcome -0.9005181
#>   std.error error.variance   cohort analyzed status reason
#> 1  2.817564       7.938667 analyzed     TRUE     ok
#> 2  2.817564       7.938667 analyzed     TRUE     ok
#> 3  2.817564       7.938667 analyzed     TRUE     ok
#> 4  2.817564       7.938667 analyzed     TRUE     ok
#> 5  2.817564       7.938667 analyzed     TRUE     ok
#> 6  2.817564       7.938667 analyzed     TRUE     ok
repair_diagnostics(fit, table = "excluded")
#> [1] id     reason
#> <0 rows> (or 0-length row.names)
```

`repair_score_covariance(fit)` preserves the assembled within- and
cross-model scoring-error covariance. It differs from `vcov(fit)`, which
describes structural coefficient uncertainty. These ordinary data frames
have explicit ID columns and the same schema for one or multiple
measurement models. See the [extractor
guide](vignettes/extracting-results.Rmd) for column contracts, cohort
selection, covariance reconstruction, and safe joins. Internal S7
properties remain available for advanced inspection; custom estimation
uses `fit_mmom()`.

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

This explicitly assumes zero **cross-model scoring-error covariance**;
latent processes can still be correlated. Separate model fits alone do
not justify the assumption. If the scoring errors are correlated, supply
per-person covariance blocks instead. Omitting the policy is an error
for multiple referenced models. Within-model covariance is always
retained. See the executable [multiple-model
guide](vignettes/multiple-models.Rmd) for both policies, ID matching,
and a dual-process reproduction against the paper’s simulation code.

## Reliability and diagnostics

Assess the selected random coefficients using the manuscript’s posterior
reliability definitions. The default uses analyzed people with equal
weights:

``` r
repair_reliability(fit)[, c("term", "role", "measure", "reliability", "n_people", "n_defined")]
#>   term    role      measure reliability n_people n_defined
#> 1 Days outcome     marginal   0.8485308       18        18
#> 2 Days outcome residualized   0.8275788       18        18
diagnostics <- repair_diagnostics(fit)
diagnostics$weight_summary
#>   n_people max_fraction weight_equivalent_n
#> 1       18   0.05555556                  18
diagnostics$conditioning
#>          step dimension min_eigenvalue max_eigenvalue condition_number
#> 1 preliminary         2     4.10807536     2122.89192         516.7607
#> 2       final         2     0.09612987       49.67614         516.7607
```

``` r
person_reliability <- repair_reliability(fit, level = "person")
plot(person_reliability)
plot(person_reliability, type = "observations")
```

Posterior reliability is distinct from the likelihood-only scoring-error
covariance used by the estimator. Undefined values retain reasons and
explicit denominators; diagnostics do not impose automatic adequacy
cutoffs. See the [reliability and diagnostics
guide](vignettes/reliability-diagnostics.Rmd) for marginal versus
intercept-residualized reliability, cohort selection, and interpretation
of numerical and weight diagnostics.

## Interpretation and limits

- Standard errors are model-based and conditional on estimated
  measurement-model quantities. They use the modified final predictor
  moment, consistent with the manuscript and primary simulation
  estimator. Measurement-model parameter uncertainty is not propagated;
  intervals and tests use a normal reference.
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
- Multi-group lme4 scoring remains experimental; multi-group structural
  fits and multilevel nlme extraction are rejected. Shared nuisance
  effects can create cross-person score errors that the structural
  estimator cannot represent.

## Development and reproducibility

The [design plan](dev/package_plan.Rmd) describes the implemented
contracts and extension sequence. The [audit and
decisions](dev/design-audit.md) trace findings to corrections and tests.
The test suite includes frozen REPAIR simulation and independent
algebra-audit fixtures for both R and C++ backends, plus three
single-model nlme and two dual-process reproduction datasets; see
[fixture provenance](tests/testthat/fixtures/README.md).

Pull requests and pushes to main run package checks on Linux, macOS, and Windows,
with additional Linux checks for previous, development, and minimum-supported R.
Tests, examples, and all vignettes run with suggested dependencies installed;
errors and warnings fail CI. See the [CI guide](dev/ci-checks.md) for the matrix,
check artifacts, local reproduction, and the `Package checks` branch-protection target.

``` r
# From a source checkout:
testthat::test_local()
```
