# repair 0.0.0.9000

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
