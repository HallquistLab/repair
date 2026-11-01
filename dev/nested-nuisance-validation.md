# Nested nuisance effects: audit, implementation, and validation

Date: 2026-10-30. Baseline: main `c409bb9` after release-hardening PR #10.
This increment implements the contained-nuisance special case of
[target-specific covariance](target_grouping_vars.Rmd). It admits Gaussian lmer
models while preserving the independent-person structural covariance contract.

## Findings and corrections

| Finding | Correction | Verification |
|---|---|---|
| Structural assembly rejected every model with multiple grouping factors | Permit one target group and nuisance levels contained within that target | Full workflow from repeated observations, two frozen references, both estimation backends |
| An intercept term can belong to both target and nuisance groups | Add optional literal `group` argument to `re()`; require one selected target group per model | Ambiguity, named/positional arguments, invalid references, and mixed-target rejection tests |
| Nuisance group names do not establish containment | Check actual membership using every effective fitted row, before structural selection/omission | Globally unique visit IDs admitted; shared local day/school IDs and reverse target direction rejected |
| Posterior reliability used only ordinary residual covariance | Share target-specific covariance construction between GLS prediction extraction and reliability | Independent full posterior matrices, marginal/residualized contrasts, and a regression showing that ignoring nuisance covariance overstates information |
| New model structure needs interpretable saved results | Print verified nuisance groups and record them in MLM metadata | Stable table identities, scope diagnostics, joint covariance, and fresh-process persistence |

## Statistical contract

For every nuisance group g and fitted nuisance level h, all effective observations
with level h must have the same target value. The rule applies to the complete
MLM fit, including people not requested in structural data. A zero
estimated nuisance variance or `na_action = "omit"` does not waive the rule.
Only excluded measurement observations are absent from the membership check.

Under the lmer model's independent random-effect blocks and independent ordinary
residuals, containment makes the conditional prediction-error covariance zero
between target people. Within each person, the covariance used by prediction extraction is

`R_target = R_error + sum(Z_nuisance G_nuisance Z_nuisance')`.

The complete target intercept/slope vector is estimated jointly; its G is
excluded from R_target. Posterior reliability then uses that R_target with the
target G included in the posterior calculation. Multiple nuisance groups,
correlated nuisance coefficients, independent target blocks, and nuisance groups
that cross *within* people are allowed when every membership check passes.

The structural M-MOM equations, regularization, backend arithmetic, and inference
convention are unchanged. SEs/normal intervals condition on fitted fixed effects,
variance components, and any supplied cross-model covariance. This does not
propagate measurement-parameter estimation uncertainty or add a bootstrap.
Shared nuisance levels can induce cross-person error covariance and are rejected.
Multilevel nlme remains unsupported.

## Interface and compatibility

`re(term, model = NULL, group = NULL)` preserves the positional model argument.
The group name comes from `repair_effects()`. Terms, models, and groups are parsed
as literal strings; no arbitrary expression evaluation is introduced. Every
selected effect in one model must use the same target group. Separate models can
have different target-group names but must identify the same analysis people.
They still require an explicit cross-model prediction-error policy.

For a model with multiple groups, canonical labels include model and group.
Existing single-group labels and public table column schemas are preserved.
`summary(fit)$mlm_models[[model]]$nuisance_groups` records verified names,
and printed summaries identify them. No new diagnostic table or column is needed.
Saved nested fits are checked in a new R process, including recomputed reliability.
The existing same-version persistence boundary remains in force.

## Independent references

`dev/generate_nested_fixtures.R` does not load repair. It generates two synthetic
48-person datasets and fits MLMs with the same documented optimizer
controls used by the refitting tests. It calculates GLS and posterior covariance
from fitted quantities directly. In particular, same-day indicator matrices and
within-day covariates form nuisance covariance; this does not call repair's
extraction/prediction extraction helpers.

The structural inputs include a random slope outcome, random intercept predictor,
an observed practice predictor, and nonzero outcome/predictor error covariance.
The companion `fit_fuller_mlm()` gives full coefficient covariance; its focal
estimate/SE/interval is required to agree with `fit_fuller_dual_core()` and
`fit_fuller_audit()` before the reference is saved. Source commit and hashes,
generator hash, R/lme4 versions, and inference definitions are embedded in
`inst/extdata/nested-workflow-reference.rds`.

Cases cover weighted/offset models, unequal numbers of days and observations,
shuffled rows, an excluded response, unused factor levels, day random intercepts,
and multiple local nuisance groups including day random slopes. Tests refit from
observations and compare joint predictions, nuisance/error/posterior covariance,
reliability, all structural coefficients, SEs, full coefficient covariance, and
intervals. Both R and C++ backends are exercised. Refit comparisons use tolerance
1e-5 for optimizer/platform variation; same-fit analytical checks use 1e-9.
Existing reference fixtures and tolerances are unchanged. These are per-dataset
implementation checks, not a full coverage or performance simulation.

## Repository and website location

The canonical repository is now `HallquistLab/repair`, and the configured Pages
URL is `https://hallquistlab.github.io/repair/`. Current install, citation, support,
metadata, and site configuration URLs are updated so the new guide points to the
current project. Version/release selection remains the separate metadata PR.

## Workflow and acceptance

1. Inspect the target-group memo, actual scorer, parser, assembly, and reliability.
2. Define containment using effective fitted membership, with explicit rejection
   of sharing and of different targets selected from the same model.
3. Factor the target covariance into a shared prediction extraction/reliability helper.
4. Freeze independent GLS/posterior/research-estimator expectations without
   loading repair; retain their provenance and existing fixture maintenance rules.
5. Exercise the public workflow, labels, scope permutation/omission, failed
   designs, shared/zero-variance effects, multiple models, and saved-fit dispatch.
6. Add the executable nested-nuisance guide; reconcile current README/help/site/
   support matrix and append dated follow-ups to historical design audits.
7. Run the full source-package check, README, all examples/vignettes, pkgdown
   content/link/search validation, and the unchanged strict six-environment CI.

Record final counts, tested revision, platform results, and any local toolchain
limitations in the PR. Reproduction commands follow `dev/ci-checks.md`.
