# Design audit and implementation decisions

Audit baseline: `47382f5` (2026-10-29 working review), following the October 28
review of `f7ab6bf`. The latter predated the new target-group scorer. This PR is
built in an isolated checkout; uncommitted edits in the author's original
`R/mmom.R` are not included or overwritten.

## Findings, corrections, and evidence

| Finding | Consequence | Applied correction | Verification |
|---|---|---|---|
| `repair()`, formula parsing, measurement fitting, structural fitting, and assembly were empty or chained empty stubs | A plausible call could yield a model containing NULL fits | Export a complete fitted-model workflow; replace unused fitting stubs with model normalization and effect discovery | Executable README/help/vignette; `test-workflow.R` |
| The outline combined fitting choices and structural specification without ID or covariance contracts | Users would have to infer joins and how selected effects map to error blocks | Explicit `id`, one-row-per-ID structural table, inspectable `repair_data()` object | Permuted/subset IDs, independent GLS and joint-error assembly tests |
| `re("x")` lacked model qualification and a scalar-coefficient definition | Ambiguous effects and factor slopes could select the wrong quantity | `repair_effects()` catalog, qualified literal marker, canonical intercept name, ambiguity errors | Unknown/ambiguous names, aliases, both-side references |
| Formula scope was broader than the M-MOM engine | No-intercept or transformed score formulas could be interpreted incorrectly | Additive grammar with intercept; reject interactions, nonlinear expressions, offsets, `.`, and unsupported models | Grammar rejection tests |
| Missingness and rank failure had no structural policy | Silent sample changes or invalid matching | Fail by default; explicit omission warns and retains every excluded ID/reason; cohort counts and unused score IDs | Missing/nonfinite/unmatched IDs, failed score blocks, duplicate IDs, minimum sample and rank tests |
| The score contract had changed from a list to S7, but a second old function definition and stale help remained | Users saw nonexistent arguments and two contradictory implementations | Remove shadowed legacy function; document the actual group-keyed S7 contract and arguments | Regenerated help, current scoring suite, documentation execution |
| Multi-group scoring supplied within-cluster blocks but no cross-cluster covariance | Shared nuisance effects can violate the structural estimator's independent-row assumption | Keep scorer available as experimental; reject multi-group models in structural assembly and document the extension gate | Explicit multi-group structural rejection; target-group memo addendum |
| S4 extraction dispatch accepted `merMod`, including GLMMs | Non-Gaussian fits could enter Gaussian GLS formulas | Guard public scoring for Gaussian lmer/lme classes | GLMM rejection regression test |
| Adding fixed effects silently skipped missing coefficient matches | Mixed deviations/total coefficients could appear in one result | Require all matching fixed coefficients when `include_fixed = TRUE`; structural workflow consistently uses deviations | Complete/incomplete matching tests |
| Extraction assumed one random-effect block per grouping factor | Ordinary `||` models failed despite a single grouping factor | Merge independent blocks in level-major order and align their block-diagonal covariance | Independent per-person GLS and OLS tests for `||` models |
| No summaries or standard accessors existed | Analysts had to navigate internal S7 properties for routine results | Summary/print, coef, vcov, confint, nobs; explicit failed-fit behavior | Accessor values, confidence levels, printed sample counts, failure tests |
| `step_up` was accepted but never implemented | A control choice produced a later numerical failure | Reject it immediately with an implementation message | Control validation test |
| Placeholder metadata and README; stale outline PDF | Neither supported use nor inference limits were communicated | Rewrite metadata, executable README, packaged vignette, help, this audit, and the implementation plan; remove obsolete PDF | README/vignette execution and package check |
| Conditional standard-error convention was no longer recorded in diagnostics | Users could miss the plug-in inference limitation | Restore explicit diagnostics and help; retain the corrected final-moment covariance; use stable normal-tail p-values | Existing simulation/algebra fixtures on both backends plus workflow checks |

## Why this scope is coherent

The PR completes one scientific path: fitted measurement model → effect discovery
→ ID/covariance assembly → M-MOM → interpretable output. It retains the existing
S7 result and scorer contracts, rather than inventing a parallel fit hierarchy.
The assembly class adds explicit dimensions and identity; M-MOM validates the
joint covariance numerically. The engine continues to use the paper's estimator.

The structural workflow supports one referenced Gaussian lmer model, one grouping
factor, and one additive equation with an intercept. Scores may be on either or
both sides; observed factor contrast columns are permitted and carry zero error.
Multiple fits may be cataloged, but cross-model equations fail instead of
assuming zero score-error covariance. Fitted-model input preserves users' fitting
controls and lets the same measurement model support several structural analyses.

Rejected scope is explicit: joint equations, implicit long-data aggregation,
nonlinear score-error propagation, no-intercept M-MOM, cross-model covariance,
structural fitting with shared nuisance groups, and automatic reliability
summaries. The unused reliability generic is removed; its statistical definition
must precede implementation. Automatic measurement fitting is deferred until the
fitted-model workflow is stable. See `package_plan.Rmd` for extension gates.

## Scientific reference and tests

The existing frozen fixtures originate in the companion `blup_unshrinking`
repository at commit `b1f09f430d2214d5c4c8815b1386b62cb4ae159e`:

- manuscript: `documentation/paper_drafts/9oct2026/VigHallquist2026.tex`;
- matrix estimator: `R/fuller_mlm.R`, `fit_fuller_mlm()`;
- primary simulation estimator: `R/stage2_estimators.R`, `fit_fuller_dual_core()`;
- independent audit: `vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R`,
  `fit_fuller_audit(preliminary = "current_cross", se_bread = "current_final")`.

`test-repair-reference.R` compares coefficients, covariance, confidence intervals,
preliminary estimates, weights, roots, and corrections on six scenarios, through
both computation backends. Fixtures contain source hashes and estimator settings;
normal tests require no sibling repository. `dev/generate_repair_fixtures.R` is the
explicit maintenance path, not part of ordinary tests.

New integration tests calculate within-person GLS independently, assemble all
outcome/predictor error terms, and compare the public workflow to explicit M-MOM.
The equal-precision, observed-predictor case is checked against OLS. Tests cover
alignment, factors, omissions, failure diagnostics, grammar, and normal intervals.
These establish implementation agreement; they do not establish new statistical
coverage properties or validate experimental backend extensions.

Run `testthat::test_local()`, render `README.Rmd`, build the package with its
vignette, and run `R CMD check --no-manual` on the tarball. The PR records the
actual results and environment-specific warnings. Keep acceptance tests and
support claims aligned when extending the design.

## Verification recorded for this PR (2026-10-29)

R 4.5.2 on macOS arm64: 399 assertions passed with zero failures, errors,
warnings, or skips. The installed-package test run also reports 399 passes.
README rendering, examples, vignette creation/rebuild, help consistency, and
method registration checks passed. Rendered output was inspected to verify that
`summary(fit)` displays the compact summary rather than an unformatted list;
a regression test checks dispatch from outside the package namespace.

`R CMD check --no-manual` completed with zero errors, zero notes, and one
installation warning from R's header suppressing the unsupported Apple-clang
warning group `-Wfixed-enum-extension`. This compiler-header warning also occurred
before this design work. `git diff --check` passed. The fixture RDS is unchanged.
