# Reliability and diagnostics: audit, contracts, and validation

Baseline: merged multiple-model workflow, `a77b0ac`. This increment implements
measurement reliability and structural diagnostic interfaces without changing
scoring, M-MOM estimation, covariance assumptions, regularization, or SEs.

## Findings and corrections

| Finding | Practical consequence | Implementation and evidence |
|---|---|---|
| The manuscript recommends empirical posterior reliability, but the package had only a deferred helper placeholder | Analysts could fit the model without computing its central information diagnostic | `repair_reliability()` implements the manuscript appendix and simulation definitions; frozen references are generated independently of repair |
| Posterior uncertainty and likelihood-only score error were easy to conflate | A user could substitute posterior variance for the estimator's error input | Separate API, named variance columns, explicit documentation, and tests establishing the quantities differ |
| A formula based on inverse G fails for singular random-effect covariance | Boundary models could crash summaries or receive fabricated zero reliability | Use the marginal observation covariance formula; retain valid marginals and label undefined prior-variance ratios |
| Unequal observation counts and multiple model cohorts leave aggregation ambiguous | Reliability averages could silently use different people or observation weights | Default to analyzed IDs, equal person weights, explicit cohort and defined/undefined counts; measurement-cohort summaries use each model's effective levels |
| Intermediate structural quantities were exposed only through nested internal slots | Users had difficulty tracing sample exclusions and numerical stability | `repair_diagnostics()` provides labeled tables for cohort, scoring, variances, regularization, corrected moments, and weight concentration |
| Numerical failure discards intermediate structural state | A diagnostic accessor could present unavailable results as zeros | Empty structural tables plus an explicit failure step/message; measurement and cohort diagnostics remain available |

## Reliability contract

The public input is a `repair()` result, including a numerical structural failure.
Only random coefficients selected in its formula are summarized. The default
`type = "both"` reports marginal reliability and, where applicable, slope
reliability residualized on the intercept. `level = "person"` returns ID-aligned
rows with effective observation count, prior/posterior variance, reliability,
status, and reason; the default summary returns an ordinary data frame with
class `repair_reliability` for plotting dispatch.

For each single-group measurement model, normalize its effective rows, full
random-effect design, fitted residual covariance R_i, and fitted random-effect
covariance G through the existing validated backend. Compute

    V_i = G - G Z_i' (R_i + Z_i G Z_i')^-1 Z_i G.

Cholesky solves avoid explicitly inverting either G or the observation covariance.
The formula permits positive-semidefinite G and rank-deficient Z_i. This is the
posterior covariance conditional on fitted parameters, not the covariance of
likelihood-only scoring errors. Fixed-effect/variance-parameter uncertainty is
not propagated.

Marginal reliability is `1 - V_i[k,k] / G[k,k]`. For a selected slope in a model
with exactly two random coefficients, one an intercept, the residualized contrast
has intercept coefficient `-G[intercept,slope] / G[intercept,intercept]` and slope
coefficient one. Its reliability is `1 - a' V_i a / (a' G a)`. This residualizes on
the intercept within that measurement model only. More general partial contrasts,
structural covariate adjustment, and cross-model posterior borrowing are not
implemented. Explicitly requested inapplicable residualized effects carry a
not-applicable status.

Default summaries use analyzed IDs, equal weights, and the same person cohort
across models. `cohort = "measurement"` uses all effective grouping levels within
each model; this includes structural exclusions and may have different denominators
across models. A failed likelihood-only score does not imply an undefined posterior
reliability. This distinction is tested with a rank-deficient within-person design.

Zero or numerically degenerate prior contrast variance yields NA with a reason.
Degeneracy uses a cancellation tolerance of 100 machine epsilons times the sum
of absolute contrast-variance contributions, in the same units as the contrast
variance. It is not a substantive reliability threshold. Computed reliabilities
within 1e-8 of the theoretical boundaries are rounded to [0, 1]; larger violations
and posterior calculation failures are undefined with reasons. If any requested
person is undefined, the overall reliability is NA. `mean_defined` separately
reports the mean of defined values, with explicit counts. There is no silent
`na.rm = TRUE` replacement of the requested-cohort estimand.

## Diagnostic and UI contract

`repair_diagnostics()` returns tables, preserving readable column names and IDs:

- `status`, `cohort`, `excluded`, `measurement`, and `score_failures` describe fit
  completion, inclusion, model sources, and failed scores (even outside the cohort).
- `error_variances` contains structural input-scale likelihood-only variances,
  including zeros for observed variables under the model's assumptions.
- `regularization` identifies the two moment policies, alpha, and error-subtraction
  multipliers. A multiplier of one still subtracts measurement error in full.
- `conditioning` contains eigenvalue extrema and the spectral condition number
  for each corrected predictor moment, including the structural intercept. Values
  depend on predictor units. Nonpositive minimum eigenvalues give Inf.
- `weights` contains individual inverse case-variance weights and normalized
  fractions. `weight_summary` contains maximum fraction and reciprocal squared
  fraction sum, labeled weight-equivalent count rather than inferential sample size.
- `reliability` contains the analyzed-cohort reliability table; `cross_model_error`
  retains the existing policy and supplied blocks.

`summary(fit)` now includes a compact role-labeled reliability table and weight
summary. `plot(repair_reliability(fit, level = "person"))` draws distributions;
`type = "observations"` draws reliability against effective observation count.
Panels identify model/effect/measure and count defined versus total people.
Undefined rows are not drawn; no adequacy cutoffs are imposed. Graphics state is
restored, including cex/mex, which base R resets when changing panel layout.

## Independent evidence and design workflow

1. Read the manuscript's empirical posterior reliability appendix and the
   companion `R/reliability_calibration.R` functions before implementing an API.
2. Define marginal/residualized contrasts, cohort defaults, aggregation, and
   undefined behavior; keep posterior uncertainty separate from score error.
3. Implement posterior calculation through normalized backend inputs and add
   read-only diagnostic tables without changing the estimator.
4. Generate `inst/extdata/reliability-reference.rds` using
   `dev/generate_reliability_fixtures.R`, which never loads repair. The fixture
   records source commit/hash, input-fixture hash, and R/nlme versions.
5. Compare against simulation helpers on balanced, unbalanced, and AR(1)/varIdent
   datasets, and compare lmer posterior variances to `ranef(condVar = TRUE)`.
   Test singular priors analytically and through controlled covariance blocks,
   rank-deficient clusters, cohort permutation, undefined denominators, multiple
   models, structural failures, and independent moment reconstruction.
6. Execute README/help/vignette examples, inspect plots, run the full regression
   suite and source-package check, and record results in the PR.

Validation on R 4.5.2 / macOS arm64: both source and installed-package tests
pass 1,098 assertions without failures, errors, warnings, or skips. README
rendering and all four vignette builds/rebuilds passed. The new plot output was
visually inspected for readable labels, cohort counts, and panel layout.
`R CMD check --no-manual` completed with zero errors, zero notes, and one
pre-existing installation warning from R's header referring to the unsupported
Apple-clang warning group `-Wfixed-enum-extension`. Restricted-network dependency
index lookups logged connectivity messages; installed dependency checks passed.
`git diff --check` passed after removing generated README output trailing spaces.
The existing coefficient/SE and multiple-model fixtures are unchanged. These checks establish implementation agreement, not
confidence intervals for reliability or new inferential coverage guarantees.

## Boundaries and next work

No bootstrap, automatic reliability cutoffs, measurement-fitting wrappers,
additional measurement-model structures, general contrast API, or cross-person
covariance is introduced. Reliability uses plug-in model quantities, so a
well-defined value does not certify model adequacy. Stable broader result
extractors and cross-platform CI are suitable subsequent engineering work;
bootstrap inference and multigroup support need separate methodological validation.
