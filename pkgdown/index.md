# repair: regression with estimated random coefficients

`repair` implements **Random-Effect Prediction Adjustment for Inference in
Regression (REPAIR)**, a two-stage method for estimating relationships involving
person- or cluster-specific processes measured with repeated observations.

## Why use REPAIR?

Do people who improve more rapidly during treatment have better subsequent
outcomes? Does a person's response to stress predict their symptoms? Are rates
of change measured in two different tasks related? These questions concern
individual differences in processes, such as change over time or within-unit
associations, and how those differences relate to other variables.

Multilevel (mixed) models describe these differences using random intercepts
and slopes. A common next step is to extract the predicted random effects and
use them as predictors or outcomes in another regression. However, these
predictions are estimates, with precision that can vary across people because
of their observation counts, measurement schedules, and predictor distributions.
Empirical Bayes predictions, often called BLUPs, also shrink uncertain individual
values toward the population mean. Shrinkage changes their variances and
covariances, and treating the predictions as perfectly measured can distort
regression coefficients and their standard errors. The scientific goal is to
learn about relationships among the underlying processes, rather than just
relationships among their predicted values.

REPAIR carries information about measurement precision into that second
regression. It first obtains **likelihood-only predictions** and their unit-specific
sampling-error covariances from fitted Gaussian mixed models. These predictions
remove shrinkage but remain noisy estimates. It then uses a **modified
method-of-moments (M-MOM)** regression to account for their sampling error when
random coefficients serve as predictors, the outcome, or both.

This workflow is useful when you want to connect processes measured across
different tasks, settings, or timescales while keeping the MLMs
separate from the regression of interest. The same principle applies to other
clusters, such as students nested within classrooms. The motivation follows
Vig and Hallquist's *Regression among random effects* manuscript; run
`citation("repair")` for the method and software references.

[**Get started with a complete analysis →**](articles/repair-workflow.html)

## Install

```r
# install.packages("remotes")
remotes::install_github("HallquistLab/repair")
```

Version 0.5.0 is currently distributed through GitHub. Installation
requires **R 4.3 or later** and a C++20 compiler with the R build toolchain
(including the matching Rtools installation on Windows).

## Fit a first model

Fit your MLM as usual, then provide a separate table with one row
per unit. The included data are entirely simulated; `practice` represents prior practice
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
| Inspect the analysis before fitting | `repair_data()` | Aligned predictions, error covariances, and an audit of the analysis scope |
| Fit the structural regression | `repair()` | Corrected coefficients and conditional model-based inference |
| Report or export results | `repair_coefficients()`, `repair_predictions()`, `repair_prediction_covariance()` | Ordinary data frames with stable column contracts and explicit IDs where applicable |
| Assess information and numerical behavior | `repair_reliability()`, `repair_diagnostics()` | Reliability, scope accounting, and numerical diagnostics |

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
| Relate random coefficients from separate MLMs | [REPAIR across multiple MLMs](articles/multiple-models.html) |
| Include day or occasion effects nested within each person | [Nested nuisance effects within people](articles/nested-nuisance.html) |
| Use correlated or heterogeneous residuals in nlme | [Validated nlme prediction extraction and reproduction](articles/nlme-workflow.html) |

## Interpret the correction carefully

- Structural equations are additive, include an intercept, and can use random
  coefficients as predictors, the outcome, or both. Each MLM must
  select one target grouping factor. lmer nuisance groups may be contained within
  target people; shared nuisance levels remain rejected. Observed structural
  variables are treated as measured without error.
- Standard errors condition on estimated MLM quantities. They do
  not propagate MLM parameter uncertainty; tests and intervals
  use a normal reference.
- Multiple MLMs require an explicit cross-model prediction-error
  covariance policy. Fitting them separately does not establish independent
  prediction errors.
- Posterior reliability and likelihood-only prediction-error covariance describe
  different quantities. `vcov(fit)` describes structural coefficient uncertainty;
  `repair_prediction_covariance(fit)` describes the assembled prediction errors.

The [workflow guide](articles/repair-workflow.html) documents supported formulas
and model structures. The [changelog](news/index.html) records implementation
changes. For reproducibility and design decisions, see the repository's
[design audit](https://github.com/HallquistLab/repair/blob/main/dev/design-audit.md)
and [fixture provenance](https://github.com/HallquistLab/repair/blob/main/tests/testthat/fixtures/README.md).
