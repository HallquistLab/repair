
# repair

`repair` implements REPAIR: structural regression involving random
coefficients from a Gaussian mixed model, with a modified
method-of-moments correction for sampling error in the estimated
coefficients. It uses **likelihood-only scores** and their
sampling-error covariance, rather than treating shrunken BLUPs as
error-free observed variables.

This is a development package. The supported structural workflow is an
additive regression with an intercept, using random coefficients from
**one fitted `lme4::lmer()` model with one grouping factor**. Random
coefficients may be predictors, the outcome, or both. Numeric, logical,
and factor observed predictors are supported. See the [workflow
guide](vignettes/repair-workflow.Rmd) for matching, omission,
interpretation, and current limits.

## Install

``` r
# install.packages("remotes")
remotes::install_github("zacharyvig/repair")
```

Installation from source requires R \>= 4.3 and a C++20 compiler with
the R build toolchain. On Windows, use the Rtools version appropriate
for your R installation.

## Fit, inspect, and summarize

Fit the measurement model using the usual `lmer()` interface. Its
weights, offsets, contrasts, estimation choices, and effective
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
#> Scores: likelihood-only random deviations
#>   model   group term    role
#>  growth Subject Days outcome
#> 
#>             estimate    se   zstat   pval  ci_low ci_high
#> (Intercept)  -2.9721 3.215 -0.9245 0.3552 -9.2729  3.3287
#> observed      0.3129 0.297  1.0534 0.2922 -0.2692  0.8949
#> 
#> Regularization (preliminary / final): minimum / minimum 
#> Error-moment subtraction multipliers: 0.8333 / 0.8333 
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
head(analysis@ids)
#> [1] "308" "309" "310" "330" "331" "332"
head(analysis@predictors)
#>     observed
#> 308        1
#> 309        2
#> 310        3
#> 330        4
#> 331        5
#> 332        6
analysis@diagnostics$excluded
#> [1] id     reason
#> <0 rows> (or 0-length row.names)
```

IDs determine alignment; data-frame row order never determines a match.
The person-level data frame defines the requested cohort.
Missing/nonfinite values, unmatched IDs, and failed scores stop by
default. Use `na_action = "omit"` to explicitly exclude them with a
warning and an ID/reason audit trail.

`fit@analysis_data` retains the matched scores, error covariance,
selected effect sources, and exclusions. `fit@struc_fit` retains M-MOM
steps and numerical diagnostics. For custom score/error inputs, use the
lower-level `fit_mmom()` API.

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
  coefficients, models without an intercept, joint multiple-equation
  fitting, and cross-model error covariances are not implemented.
  Unsupported formulas fail explicitly.
- The lower-level scorer has experimental multi-group and `nlme`
  support. Shared nuisance effects can create cross-person score errors;
  those cases are not admitted by the structural workflow.

## Development and reproducibility

The [design plan](dev/package_plan.Rmd) describes the implemented
contracts and extension sequence. The [audit and
decisions](dev/design-audit.md) trace findings to corrections and tests.
The test suite includes frozen REPAIR simulation and independent
algebra-audit fixtures for both R and C++ backends; see [fixture
provenance](tests/testthat/fixtures/README.md).

``` r
# From a source checkout:
testthat::test_local()
```
