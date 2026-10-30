# repair 0.5.0

* Introduce Random-Effect Prediction Adjustment for Inference in Regression
  (REPAIR) by its full name in the landing page, workflow guide, and method
  reference documentation. Explain the applied motivation, shrinkage, and
  sampling error before introducing the technical workflow.

# repair 0.0.0.9000

* Add `repair_refit()` for identity-verified replay of one supported lmer/lme
  MLM and its structural analysis, with explicit reference and
  replacement inputs, factor/row/schema checks, and recorded fitting recipes.
* Add deterministic whole-unit draws with fresh target/nuisance identities,
  audited warnings and controls, and explicit scope restrictions. This is
  refitting groundwork; bootstrap inference is not yet implemented.
* Record structural factor coding for reproducible replay and add a refit audit
  to summaries. Add tests against manual refits and fresh-process replay, and an
  executable refitting guide.

* Support Gaussian lmer nuisance groups contained within each selected target
  person, including nested random slopes and multiple nuisance groups. Verify
  containment on all effective fitted rows; retain rejection of shared levels.
* Add optional `group` qualification to `re()` for ambiguous random coefficient
  names. Each referenced model selects one target group; existing single-group
  labels and stable table schemas are preserved.
* Include nuisance covariance in both prediction extraction and posterior reliability; report
  verified nuisance groups in structural summaries and measurement metadata.
* Add independent nested GLS, posterior, and production/algebra structural
  references, integration/rejection/persistence tests, and an executable guide.

* Add a documented synthetic teaching dataset, `repair_example`, with repeated
  observations, unit-level predictors, and explicit generating provenance.
* Complete executable examples for all exported objects and add a support and
  troubleshooting guide covering model/inference boundaries and saved results.
* Verify saved single- and multiple-model fits in fresh R sessions, including
  standard generics, stable tables, supplied covariance, and diagnostics.
* Add separate software and method citations, identify the C++20 requirement,
  and index all installed reference datasets. Credit Michael N. Hallquist as an
  author and co-maintainer, with Zachary A. Vig remaining primary maintainer.
* Document the release-preparation workflow and outstanding publication checks;
  the package remains a development version.

* Add a pkgdown documentation website with a guided home page, task-oriented
  reference, all five executable guides, search, and release notes. Build and
  validate it on pull requests and deploy main to GitHub Pages.

* Add automated source-package checks on Linux, macOS, and Windows, with Linux
  coverage for previous, development, and minimum-supported R. Include full
  dependencies, tests, examples, vignettes, and a checked-package README run.
* Preserve check artifacts and a stable `Package checks` aggregate status; document
  local reproduction, warning policy, and troubleshooting.

* Add `repair_coefficients()` for plain structural coefficient tables with
  documented columns and configurable normal confidence intervals.
* Add `repair_predictions()` for ID-keyed likelihood-only score tables and
  `repair_prediction_covariance()` for assembled within- and cross-model error
  covariance. Both work on fitted results and pre-fit `repair_data()` assemblies.
* Add exact named-table selection to `repair_diagnostics()`, preserving the
  existing list interface and avoiding unrequested reliability computation.
* Document table schemas, types, keys, ordering, cohorts, failure behavior, and
  reporting/join/export workflows; test the interfaces against existing references.

* Add `repair_reliability()` with manuscript-defined marginal and intercept-
  residualized posterior reliability, unit-level results, explicit scope
  selection, equally weighted summaries, and reasons for undefined values.
* Add `repair_diagnostics()` tables for sample accounting, prediction extraction failures,
  input error variances, regularization, corrected-moment conditioning, and
  analysis-weight concentration. Preserve useful diagnostics on structural failure.
* Include compact reliability and weight summaries in `summary(fit)`; add plots
  of unit-specific reliability distributions and observation-count relationships.
* Add independent simulation-derived reliability fixtures, boundary/diagnostic
  tests, and an executable reliability guide distinguishing posterior uncertainty
  from the likelihood-only prediction-error covariance used by REPAIR.

* Support structural regressions using coefficients from multiple fitted Gaussian
  lmer or linear lme MLMs, each with one grouping factor.
* Require an explicit cross-model prediction-error policy: independent errors or
  named per-ID covariance blocks for every model pair. Validate complete model,
  coefficient, and ID coverage and joint positive semidefiniteness.
* Preserve within-model covariance, align each model to the requested scope,
  and report model-specific exclusions and MLM diagnostics. Record the
  covariance policy and supplied blocks in analysis diagnostics and summaries.
* Add two independently generated dual-process reproduction fixtures, tests for
  both computation backends, and an executable multiple-model vignette covering
  covariance assumptions and the simulation's population-SD reporting scale.

* Add validated single-group linear `nlme::lme()` support to prediction extraction, effect
  discovery, structural-data assembly, and fitting. Support independent/AR(1)/
  continuous-time AR(1) correlation with constant/varIdent/varPower residual
  variance, using stored effective rows, fitted contrasts, and fitted residual SDs.
* Reject unsupported nlme structures and missing stored data explicitly. Report
  the MLM backend and available variance-parameter uncertainty messages.
* Add three end-to-end reproduction datasets from the paper's production prediction extraction
  and M-MOM functions, independent GLS tests, and an executable nlme vignette.

* Add the exported `repair()` workflow for additive structural regressions from
  one fitted, single-group Gaussian lmer model. `re()` references random
  coefficients on either side of the equation; `repair_effects()` lists their
  model, group, and scalar coefficient names.
* Add `repair_data()` for explicit ID matching, covariance alignment, and
  inspectable sample accounting. Missing values and failed predictions stop by
  default; explicit omission warns and records excluded IDs and reasons.
* Add summary, print, coef, vcov, confint, and nobs methods. Summaries report
  conditional normal inference, regularization, sample counts, effect sources,
  and MLM diagnostics.
* Support independent random-effect blocks for a shared grouping factor,
  including lme4 `||` syntax, with level/term/covariance alignment.
* Remove the shadowed legacy prediction extraction function and correct the documented S7
  prediction extraction interface. Reject non-Gaussian prediction extraction and incomplete fixed-effect
  matches; reject the unimplemented step-up regularization option immediately.
* Document estimator assumptions, supported formula grammar, experimental
  backend boundaries, and design decisions in an executable README, vignette,
  audit, and implementation plan. Preserve REPAIR numerical reference fixtures.
