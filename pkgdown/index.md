# repair: regression with estimated random coefficients

**Estimate relationships involving random intercepts and slopes while accounting
for their sampling error.** `repair` implements the REPAIR modified
method-of-moments estimator for structural regressions using fitted Gaussian
mixed models from **lme4** or **nlme**.

REPAIR uses likelihood-only scores and their sampling-error covariance. A
regression that treats estimated BLUPs as error-free data does not make this
correction.

[**Get started with a complete analysis →**](articles/repair-workflow.html)

## Install

```r
# install.packages("remotes")
remotes::install_github("HallquistLab/repair")
```

This is a development package, currently distributed through GitHub. Installation
requires **R 4.3 or later** and a C++20 compiler with the R build toolchain
(including the matching Rtools installation on Windows).

## Fit a first model

Fit your measurement model as usual, then provide a separate table with one row
per person. The included data are entirely simulated; `practice` represents prior practice
hours per week. See `?repair_example` for provenance and generating parameters.

```r
library(repair)

data(repair_example)
growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
people <- repair_example$people

fit <- repair(
  re("day", model = "growth") ~ practice,
  models = list(growth = growth),
  data = people,
  id = "person"
)
repair_coefficients(fit)
```

`re("day", model = "growth")` selects the random slope for `day`, expressed as
a deviation from the fitted fixed slope. IDs align people across the measurement
model and structural data; row order never determines a match. The
[support guide](articles/support-and-troubleshooting.html) executes this example.
The [getting-started guide](articles/repair-workflow.html) explains the broader
workflow, omissions, regularization, and interpretation.

## Know what each step does

| Task | Interface | What you get |
|---|---|---|
| Discover available random coefficients | `repair_effects()` | Model, grouping factor, and coefficient names to use in `re()` |
| Inspect the analysis before fitting | `repair_data()` | Aligned scores, error covariances, and an audit of the analysis cohort |
| Fit the structural regression | `repair()` | Corrected coefficients and conditional model-based inference |
| Report or export results | `repair_coefficients()`, `repair_scores()`, `repair_score_covariance()` | Ordinary data frames with stable column contracts and explicit IDs where applicable |
| Assess information and numerical behavior | `repair_reliability()`, `repair_diagnostics()` | Reliability, cohort accounting, and numerical diagnostics |

See the [function reference](reference/index.html) for arguments, return values,
and executable examples. Familiar `summary()`, `coef()`, `vcov()`, `confint()`, and
`nobs()` methods are also available.

## Choose a guide

| If you want to… | Read… |
|---|---|
| Check supported models or diagnose a problem | [Supported analyses and troubleshooting](articles/support-and-troubleshooting.html) |
| Refit or revise an analysis with reference verification | [Validated refitting](articles/refitting.html) |
| Run and understand a complete analysis | [A structural regression with REPAIR](articles/repair-workflow.html) |
| Create reporting tables, join by ID, or export results | [Extracting tables for reporting and reuse](articles/extracting-results.html) |
| Assess marginal and intercept-residualized reliability | [Reliability and analysis diagnostics](articles/reliability-diagnostics.html) |
| Relate random coefficients from separate measurement models | [REPAIR across multiple measurement models](articles/multiple-models.html) |
| Include day or occasion effects nested within each person | [Nested nuisance effects within people](articles/nested-nuisance.html) |
| Use correlated or heterogeneous residuals in nlme | [Validated nlme scoring and reproduction](articles/nlme-workflow.html) |

## Interpret the correction carefully

- Structural equations are additive, include an intercept, and can use random
  coefficients as predictors, the outcome, or both. Each measurement model must
  select one target grouping factor. lmer nuisance groups may be contained within
  target people; shared nuisance levels remain rejected. Observed structural
  variables are treated as measured without error.
- Standard errors condition on estimated measurement-model quantities. They do
  not propagate measurement-model parameter uncertainty; tests and intervals
  use a normal reference.
- Multiple measurement models require an explicit cross-model scoring-error
  covariance policy. Fitting them separately does not establish independent
  scoring errors.
- Posterior reliability and likelihood-only scoring-error covariance describe
  different quantities. `vcov(fit)` describes structural coefficient uncertainty;
  `repair_score_covariance(fit)` describes the assembled scoring errors.

The [workflow guide](articles/repair-workflow.html) documents supported formulas
and model structures. The [changelog](news/index.html) records implementation
changes. For reproducibility and design decisions, see the repository's
[design audit](https://github.com/HallquistLab/repair/blob/main/dev/design-audit.md)
and [fixture provenance](https://github.com/HallquistLab/repair/blob/main/tests/testthat/fixtures/README.md).
