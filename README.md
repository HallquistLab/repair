
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

## Package website

Browse the [package website](https://zacharyvig.github.io/repair/) for a
guided introduction, searchable function reference, executable articles,
and release notes. The site is rebuilt from `main`; it documents the
development version.

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
#> People: 60 analyzed / 60 requested; 0 excluded
#> Measurement model: growth (lme4; group person; 0 IDs outside cohort)
#> Cross-model scoring errors: not applicable (one model)
#> Scores: likelihood-only random deviations
#>   model  group term    role
#>  growth person  day outcome
#>
#> Posterior reliability (analyzed cohort; equal person weights):
#>   model term    role      measure reliability n_people n_defined
#>  growth  day outcome     marginal       0.943       60        60
#>  growth  day outcome residualized       0.897       60        60
#>
#>             estimate      se  zstat      pval  ci_low ci_high
#> (Intercept)  -0.7895 0.11473 -6.882 5.905e-12 -1.0144 -0.5647
#> practice      0.1755 0.02142  8.193 2.548e-16  0.1335  0.2175
#>
#> Regularization (preliminary / final): minimum / minimum
#> Error-moment subtraction multipliers: 0.95 / 0.95
#> Largest analysis-weight fraction: 0.0179 ; weight-equivalent count: 59.8
#> Moment conditioning and individual weights: repair_diagnostics(fit).
#> Conditional on fitted measurement-model parameters; model-based normal SEs and intervals.
```

The marker `re("day", model = "growth")` selects the random slope for
`day`. It represents a deviation from the fitted fixed slope. Use
`re("(Intercept)", model = "growth")` for the random intercept. When
there is only one matching effect, `model` may be omitted.

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
head(repair_scores(analysis))
#>     id  model  group term                       label    role      score
#> 1 P001 growth person  day re("day", model = "growth") outcome  0.4451909
#> 2 P002 growth person  day re("day", model = "growth") outcome  0.4105728
#> 3 P003 growth person  day re("day", model = "growth") outcome -0.4013259
#> 4 P004 growth person  day re("day", model = "growth") outcome  0.6289102
#> 5 P005 growth person  day re("day", model = "growth") outcome  0.6052267
#> 6 P006 growth person  day re("day", model = "growth") outcome -1.2120879
#>   std.error error.variance   cohort analyzed status reason
#> 1 0.1401763     0.01964941 analyzed     TRUE     ok
#> 2 0.1401763     0.01964941 analyzed     TRUE     ok
#> 3 0.2268540     0.05146273 analyzed     TRUE     ok
#> 4 0.1618617     0.02619921 analyzed     TRUE     ok
#> 5 0.1401763     0.01964941 analyzed     TRUE     ok
#> 6 0.2268540     0.05146273 analyzed     TRUE     ok
head(repair_score_covariance(analysis))
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
The person-level data frame defines the requested cohort.
Missing/nonfinite values, unmatched IDs, and failed scores stop by
default. Use `na_action = "omit"` to explicitly exclude them with a
warning and an ID/reason audit trail.

Use public table extractors for reporting, joins, and export:

``` r
repair_coefficients(fit, conf.level = .95)
#>          term   estimate std.error statistic      p.value   conf.low  conf.high
#> 1 (Intercept) -0.7895289 0.1147252 -6.881913 5.905407e-12 -1.0143862 -0.5646716
#> 2    practice  0.1754854 0.0214190  8.192979 2.548385e-16  0.1335049  0.2174659
#>   conf.level
#> 1       0.95
#> 2       0.95
head(repair_scores(fit))
#>     id  model  group term                       label    role      score
#> 1 P001 growth person  day re("day", model = "growth") outcome  0.4451909
#> 2 P002 growth person  day re("day", model = "growth") outcome  0.4105728
#> 3 P003 growth person  day re("day", model = "growth") outcome -0.4013259
#> 4 P004 growth person  day re("day", model = "growth") outcome  0.6289102
#> 5 P005 growth person  day re("day", model = "growth") outcome  0.6052267
#> 6 P006 growth person  day re("day", model = "growth") outcome -1.2120879
#>   std.error error.variance   cohort analyzed status reason
#> 1 0.1401763     0.01964941 analyzed     TRUE     ok
#> 2 0.1401763     0.01964941 analyzed     TRUE     ok
#> 3 0.2268540     0.05146273 analyzed     TRUE     ok
#> 4 0.1618617     0.02619921 analyzed     TRUE     ok
#> 5 0.1401763     0.01964941 analyzed     TRUE     ok
#> 6 0.2268540     0.05146273 analyzed     TRUE     ok
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
#> 1  day outcome     marginal   0.9433813       60        60
#> 2  day outcome residualized   0.8967932       60        60
diagnostics <- repair_diagnostics(fit)
diagnostics$weight_summary
#>   n_people max_fraction weight_equivalent_n
#> 1       60   0.01788085            59.82168
diagnostics$conditioning
#>          step dimension min_eigenvalue max_eigenvalue condition_number
#> 1 preliminary         2       17.52838       1747.322         99.68531
#> 2       final         2       74.13875       7551.600        101.85768
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

## Support, troubleshooting, and citation

The [support guide](vignettes/support-and-troubleshooting.Rmd)
consolidates the supported model/inference matrix, a complete
`repair_example` analysis, common diagnostic messages, and saved-fit
guidance. The teaching data include both repeated observations and
meaningful synthetic person-level variables.

Run `citation("repair")` for separate software and method references.
Zachary A. Vig and Michael N. Hallquist are authors and co-maintainers;
Zach is the primary R package maintainer contact. This remains a
development version, with no tagged release claimed by these
instructions.

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
