# repair 0.0.0.9000

* Add experimental `repair_leverage()` for a specified full-system HC3-style
  correction, including all measurement/M-MOM nuisance equations and cross-stage
  covariance. Normalize common weights to reproduce OLS HC3; retain the original
  full-influence gate and require every deletion sensitivity to be usable.
* Add `repair_jackknife()` for validated whole-person N-1 refitting, separate
  coefficient/interval extraction, strict complete-vector failure handling,
  optimizer diagnostics and diagnostic bias estimates. Original coefficients
  remain unchanged. Normal intervals are primary; N-1 t intervals are optional
  experimental sensitivity analyses.
* Document a fresh independent-stream 600-dataset comparison. HC3-style
  adjustment improves selected-case coverage over plain IF, with substantial
  width tails; delete-one refits remain more widely available but their Wald
  intervals often under-cover. Retain the frozen protocol, person-shift arrays,
  endpoint-based verifier and full availability/precision tradeoffs under `dev/`.

* Add experimental `repair_sandwich()` and separate coefficient/covariance/interval
  methods for full person-level influence inference. Propagate measurement
  estimation, preliminary/final M-MOM moments, estimated weights and regularization,
  including signed cross-stage covariance. Retain the original point estimates
  and ordinary conditional inference. Report frozen-measurement influence
  separately and withhold unstable/boundary full inference explicitly.
* Implement exact Gaussian ML/REML scores with the REML determinant's
  case-frequency derivative, validated numerical sensitivities, ID-aligned
  influences and covariance decomposition. Initial support is one Gaussian lmer,
  including contained nuisance groups, weights and offsets; nlme/multiple-model
  derivatives remain unsupported.
* Document a paired 600-dataset comparison with person bootstrap and the existing
  research OLS-outcome stacked HC0–HC3 SEs. Full influence results were available
  for 1,351/2,400 structural analyses; difficult noisy-predictor coverage remained
  below nominal and did not improve on bootstrap. Preserve endpoint-based review
  evidence and matched Monte Carlo comparisons under `dev/influence-simulation`.

* Add opt-in `repair_lmer_refitter()` for deterministic, likelihood-checked
  nloptwrap-initialized BOBYQA fitting and bounded retries. Retain every attempt,
  convergence warning, likelihood and singularity diagnostic. Persist the recipe
  for automatic original/reference/replacement replay and saved-result auditing.
  Do not disable convergence checks or change the bootstrap failed-attempt gate.
* Document a paired six-profile simulation screen: default ordinary-bootstrap
  interval availability improved from 1--53% to 95--100%, with condition-dependent
  coverage and persistent small/short-panel bias. Include compact evidence and
  an independent interval-summary verifier under `dev/`; broader validation
  remains outstanding.


* Add experimental `repair_bootstrap()` for serial whole-person sampling through
  `repair_refit()`, with separate saved RNG streams, complete attempt accounting,
  fresh target/nuisance IDs, and unchanged original conditional inference.
* Add a separate bootstrap result class, coefficient and attempt/draw table
  extractors, empirical covariance and percentile intervals. Withhold inference
  after failed attempts unless explicitly requested conditional on success.
* Add independent manual-loop, RNG, failure, nested/nlme and saved-result tests,
  structured refit failure stages, an executable user guide, and a design audit.
  Scientific coverage validation and multiple-model sampling remain future work.

* Add `repair_refit()` for identity-verified replay of one supported lmer/lme
  measurement model and its structural analysis, with explicit reference and
  replacement inputs, factor/row/schema checks, and recorded fitting recipes.
* Add deterministic whole-person draws with fresh target/nuisance identities,
  audited warnings and controls, and explicit cohort restrictions. This is
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
* Include nuisance covariance in both scoring and posterior reliability; report
  verified nuisance groups in structural summaries and measurement metadata.
* Add independent nested GLS, posterior, and production/algebra structural
  references, integration/rejection/persistence tests, and an executable guide.

* Add a documented synthetic teaching dataset, `repair_example`, with repeated
  observations, person-level predictors, and explicit generating provenance.
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
* Add `repair_scores()` for ID-keyed likelihood-only score tables and
  `repair_score_covariance()` for assembled within- and cross-model error
  covariance. Both work on fitted results and pre-fit `repair_data()` assemblies.
* Add exact named-table selection to `repair_diagnostics()`, preserving the
  existing list interface and avoiding unrequested reliability computation.
* Document table schemas, types, keys, ordering, cohorts, failure behavior, and
  reporting/join/export workflows; test the interfaces against existing references.

* Add `repair_reliability()` with manuscript-defined marginal and intercept-
  residualized posterior reliability, person-level results, explicit cohort
  selection, equally weighted summaries, and reasons for undefined values.
* Add `repair_diagnostics()` tables for sample accounting, scoring failures,
  input error variances, regularization, corrected-moment conditioning, and
  analysis-weight concentration. Preserve useful diagnostics on structural failure.
* Include compact reliability and weight summaries in `summary(fit)`; add plots
  of person-specific reliability distributions and observation-count relationships.
* Add independent simulation-derived reliability fixtures, boundary/diagnostic
  tests, and an executable reliability guide distinguishing posterior uncertainty
  from the likelihood-only score-error covariance used by REPAIR.

* Support structural regressions using coefficients from multiple fitted Gaussian
  lmer or linear lme measurement models, each with one grouping factor.
* Require an explicit cross-model scoring-error policy: independent errors or
  named per-ID covariance blocks for every model pair. Validate complete model,
  coefficient, and ID coverage and joint positive semidefiniteness.
* Preserve within-model covariance, align each model to the requested cohort,
  and report model-specific exclusions and measurement diagnostics. Record the
  covariance policy and supplied blocks in analysis diagnostics and summaries.
* Add two independently generated dual-process reproduction fixtures, tests for
  both computation backends, and an executable multiple-model vignette covering
  covariance assumptions and the simulation's population-SD reporting scale.

* Add validated single-group linear `nlme::lme()` support to scoring, effect
  discovery, structural-data assembly, and fitting. Support independent/AR(1)/
  continuous-time AR(1) correlation with constant/varIdent/varPower residual
  variance, using stored effective rows, fitted contrasts, and fitted residual SDs.
* Reject unsupported nlme structures and missing stored data explicitly. Report
  the measurement backend and available variance-parameter uncertainty messages.
* Add three end-to-end reproduction datasets from the paper's production scoring
  and M-MOM functions, independent GLS tests, and an executable nlme vignette.

* Add the exported `repair()` workflow for additive structural regressions from
  one fitted, single-group Gaussian lmer model. `re()` references random
  coefficients on either side of the equation; `repair_effects()` lists their
  model, group, and scalar coefficient names.
* Add `repair_data()` for explicit ID matching, covariance alignment, and
  inspectable sample accounting. Missing values and failed scores stop by
  default; explicit omission warns and records excluded IDs and reasons.
* Add summary, print, coef, vcov, confint, and nobs methods. Summaries report
  conditional normal inference, regularization, sample counts, effect sources,
  and measurement-model diagnostics.
* Support independent random-effect blocks for a shared grouping factor,
  including lme4 `||` syntax, with level/term/covariance alignment.
* Remove the shadowed legacy scoring function and correct the documented S7
  scoring interface. Reject non-Gaussian scoring and incomplete fixed-effect
  matches; reject the unimplemented step-up regularization option immediately.
* Document estimator assumptions, supported formula grammar, experimental
  backend boundaries, and design decisions in an executable README, vignette,
  audit, and implementation plan. Preserve REPAIR numerical reference fixtures.
