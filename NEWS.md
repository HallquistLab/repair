# repair 0.0.0.9000

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
