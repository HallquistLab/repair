# nlme extraction audit and validation

This follow-up to PR #2 completes the next design milestone: single-group linear
nlme extraction → ID-aligned structural assembly → paper-code reproduction.
The baseline is merge commit `94bdaef`; the structural M-MOM equations and prior
reference fixtures are unchanged.

## Confirmed problems and corrections

1. **Effective rows:** `nlme::getData()` may retain exclusions or re-evaluate
   `subset`; `getResponse()` may pad excluded observations. The previous code
   combined those with unpadded fitted groups. Extraction now joins retained
   data to the model's effective row names and uses the saved fitted/residual
   response. It verifies all design/response/group rows agree.
2. **Contrasts:** reconstructing model matrices with current defaults can differ
   from fitting-time contrasts. Stored contrasts are applied to fixed and random
   designs; reconstructed fixed and conditional fitted values must agree with
   the fit. Missing retained data and inconsistent reconstructed designs fail
   with actionable messages.
3. **Residual covariance:** covariance blocks now combine fitted residual SDs in
   effective-row order with fitted correlation matrices whose group ordering is
   checked explicitly. They retain heterogeneous variance and correlated errors.
   `varPower(~ fitted(.))` uses fitted weights: final fitted values, and even the
   saved variance covariate/power expression, can differ from the last retained
   weights by convergence tolerance. No variance function is re-estimated.
4. **Covariance assembly:** the nlme block assignment in `get_ranef_base_covariance`
   incorrectly included `drop = FALSE` on the assignment side. It now uses valid
   two-dimensional replacement indexing.
5. **Backend assumptions:** effect discovery, grouping validation, and measurement
   diagnostics no longer call lme4 methods on lme objects. Summaries name the
   backend and retain nlme variance-parameter uncertainty messages. nlme
   singularity is explicitly not assessed using the lme4 criterion.
6. **Unsupported scope:** the old extractor attempted multilevel nlme extraction
   despite `getVarCov.lme()` rejecting that case. Supported correlation classes
   are independent, AR(1), and continuous AR(1); supported variance classes are
   constant, varIdent, and varPower. AR(1) internally converted to corARMA(1,0)
   is accepted. Finer correlation grouping, other covariance classes, nonlinear
   fits, missing stored data, and multilevel nlme extraction now fail explicitly.

## Validation design

Same-fit analytical tests construct covariance from residual scale, AR
coefficients/time distances, and variance multipliers (or the authoritative
fitted weights for fitted-value varPower). GLS scores and covariance are computed
with base matrix solves, independently of package extraction. Tests include
unequal cluster sizes, shuffled observations, unused factor levels, subset,
na.omit/na.exclude, fitted contrasts, changed caller data, all admitted residual
classes, intercept-only scoring, rank failure and structural omission. Finer
correlation groups and unsupported structures are rejection tests.

Three bundled fixtures start from raw synthetic repeated observations and refit
nlme in ordinary tests. Expectations are generated without repair by the paper's
production functions in the companion repository at
`b1f09f430d2214d5c4c8815b1386b62cb4ae159e`:

- `R/blup_helpers.R`: `get_closed_form_corrected_scores()`;
- `R/fuller_mlm.R`: `fit_fuller_mlm()` for complete matrices;
- `R/stage2_estimators.R`: `fit_fuller_dual_core()` for focal inference;
- `vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R`: independent
  `fit_fuller_audit()` required to agree before saving.

Cases cover balanced predictor error, unbalanced outcome error, and same-model
correlated outcome/predictor error with AR(1) and heterogeneous residuals. The
last is a covariance-propagation extension, not a cross-model analysis. Tests
compare score values, error covariances, coefficients, SEs, CIs, and coefficient
covariance through both numerical backends. Stored provenance contains source
hashes, estimator settings, and R/nlme versions. Refitted comparisons use 1e-5
tolerance for platform/optimizer variation; same-fit GLS checks are stricter.

The production comparison uses modified preliminary/final moments and modified
variance bread with `skip_internal_scaling = TRUE`, matching package input-scale
inference. The predictor fixture also verifies the external population-slope-SD
reporting convention used by `vh_study2.R`; the known generating SD is 0.8.
The scale multiplies estimates, SEs, and CI endpoints together. No automatic
standardization or alpha-search implementation is added.

## Developer workflow and remaining limits

Run `testthat::test_local()`, render the README, build both executable vignettes,
and run `R CMD check --no-manual`. For intentional fixture maintenance, run
`Rscript dev/generate_nlme_fixtures.R /path/to/blup_unshrinking` and review the
new outputs/provenance. Ordinary tests and examples use the bundled fixtures
without any sibling-repository dependency.

These fixtures establish per-dataset agreement, not reproduction of every
published Monte Carlo table or new inference coverage. First-stage parameter
uncertainty remains unpropagated. At this milestone structural fitting required one model
and one grouping factor; the later [multiple-model milestone](multimodel-validation.md)
relaxes the model count while retaining one grouping factor per model. At that stage, shared nuisance groups, cross-model covariance, other
residual families, and arbitrary nonlinear score transformations require separate
statistical and extraction validation before their interfaces can be enabled.

## Verification (2026-10-29)

On R 4.5.2 / nlme 3.1.170 / macOS arm64, all 804 assertions pass with zero
failures, errors, warnings, or skips, including the existing manuscript fixtures.
The built package's examples and both vignettes execute successfully. Rendered
reference tables show matching input-scale and population-SD-scale inference.
`R CMD check --no-manual` reports zero errors, zero notes, and one pre-existing
compiler-header warning: Apple's clang does not recognize R's suppression of
`-Wfixed-enum-extension`. `git diff --check` passes.
