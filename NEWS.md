# repair (development version)

* Add explicit experimental `repair_bayes_plugin()` snapshots for narrowly
  supported sampled Gaussian identity-link brms/rstanarm measurement models.
  Reuse the existing likelihood-only GLS and M-MOM workflow with posterior means
  of shared parameters, full covariance matrices and E(sigma^2).
* Record priors, native-engine version, sampler diagnostics and inference limits.
  Reject unsupported likelihood/design features and Bayesian refitting/bootstrap
  requests. Structural intervals remain conditional model-based confidence
  intervals; full refitted-workflow coverage is unvalidated.
* Add a Bayesian plug-in guide, focused software tests, opt-in native sampling
  tests and a prospective full refitted-workflow simulation plan in `dev/`.

# repair 1.0.0 (2026-11-04)

First release candidate for Random-Effect Prediction Adjustment for Inference
in Regression (REPAIR). This version provides the core estimator, adaptive
step-up regularization, and experimental whole-unit bootstrap inference.

## Estimation and supported models

* Add optional `repair_mlm()` specifications for fitting first-stage Gaussian
  lmer or linear lme models inside `repair()`, `repair_data()`, and
  `repair_effects()`. Retain fitted-model inputs and allow named lists mixing
  specifications and fitted models. Forward native fitting choices and preserve
  data expressions and resolved controls for supported reference replay.
  Structural estimation, regularization, inference and bootstrap failure rules
  are unchanged.

* Add `repair()`, `re()`, `repair_effects()`, and `repair_data()` for additive
  structural regressions with an intercept. Random coefficients can be
  predictors, the outcome, or both; observed numeric, logical, and factor
  predictors are supported.
* Extract likelihood-only random-coefficient predictions and their full
  sampling-error covariances from Gaussian `lme4::lmer()` and supported linear
  `nlme::lme()` models. Support contained lmer nuisance groups and supported
  correlated or heterogeneous nlme residual structures. Reject shared nuisance
  groups and unsupported model structures explicitly.
* Support coefficients from multiple measurement models with an explicit
  independent-error assumption or supplied cross-model covariance blocks.
  Align observations by identifier and retain exclusions and their reasons.
* Provide the modified method-of-moments estimator with minimum, none, custom,
  and sequential step-up (M-MOMS) regularization, using R or Rcpp computation.
  Record selected alphas, conditioning checks, evaluated candidates, and explicit
  minimum-alpha fallbacks. Refits and bootstrap attempts reselect adaptive alphas.
* Accelerate step-up with early stopping and dataset-local caching of invariant
  moments, determinant roots, weights, and covariance calculations. Verify exact
  agreement with the pre-optimization implementation for selected alphas,
  estimates, SEs, covariance, bootstrap results, and failure decisions.

## Refitting and experimental bootstrap

* Add `repair_refit()` for identity-verified replay and replacement-data fitting
  of one supported measurement model and its structural analysis. Preserve
  factor coding, model specification, effective rows, and fitting recipes.
* Add `repair_bootstrap()` for serial whole-unit sampling with separate saved RNG
  streams, fresh target/nuisance identities, and complete attempt accounting.
  Provide empirical covariance, marginal percentile intervals, and result,
  diagnostic, and draw-table extractors. Withhold inference after failed attempts
  unless the user explicitly requests inference conditional on success.
* Add opt-in `repair_lmer_refitter()` for deterministic, likelihood-checked
  BOBYQA fitting with initialization and bounded retries. Retain every internal
  attempt, convergence warning, likelihood, and singularity diagnostic.
* Keep bootstrap inference experimental. The paired six-profile robust-refitter
  screen improves interval availability but does not establish general coverage
  or validate adaptive M-MOMS. Ordinary structural SEs remain conditional on
  fitted measurement quantities. Influence-function, leverage, and delete-one
  inference are excluded from this release.

## Results, documentation, and release preparation

* Add standard summary, print, coefficient, covariance, confidence-interval,
  and observation-count methods; stable plain-table extractors; posterior
  reliability summaries; and sample, conditioning, and weight diagnostics.
* Include the entirely simulated `repair_example` teaching data, executable help
  and workflow guides, support/troubleshooting documentation, software and
  manuscript citations, and a searchable package website.
* Verify independent numerical references, R/Rcpp agreement, bootstrap/RNG
  behavior, and fresh-session persistence. Check source packages across Linux,
  macOS, and Windows, including development and minimum-supported R versions.
* Prepare version 1.0.0 metadata and release instructions, identify both authors
  as copyright holders, and cite Fuller's measurement-error methods in the
  package description. Keep statistical criteria and defaults unchanged.
* Declare tested dependency floors for S7, lme4, Rcpp, and RcppArmadillo. Include
  HTML-manual validation alongside PDF and incoming checks in release CI.
