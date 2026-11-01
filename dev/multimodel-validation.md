# Multiple MLMs: audit and design decisions

This milestone extends the validated single-group workflow to a structural
regression using coefficients from several fitted MLMs. It does
not change the REPAIR estimator, prediction extraction equations, regularization, or SE bread.

## Audit findings and corrections

| Finding | Consequence | Correction and evidence |
|---|---|---|
| The parser cataloged several models but rejected equations using more than one | The paper's dual-process workflow required manual arrays | Retain all referenced model identities and score each model once; end-to-end dual-process fixtures exercise public `repair()` |
| Assembly selected one model, group, and prediction object | No coherent cross-model ID or coefficient alignment | Construct a joint selected-effect covariance by named model and match each model to the requested scope; tests reverse model, coefficient, block, and person order |
| Separate fits do not determine between-model prediction-error covariance | Implicit zero blocks could encode an unjustified scientific assumption | Require `cross_model_error = "independent"` or supplied blocks for every pair; never infer independence |
| Pairwise covariance checks are insufficient for three models | Individually plausible pairs may form an impossible joint covariance | Check the entire per-person selected-effect covariance for positive semidefiniteness; reject a three-model counterexample |
| Exclusion and MLM diagnostics identified a single model | Missingness and prediction extraction failures across models would be hard to trace | Accumulate reasons naming the responsible models; expose per-model backend, group, unused IDs, messages, and singularity |
| Single-model references could not verify dual-process inference or its reporting scale | A correct API could still misassemble predictor/outcome covariance or report the wrong units | Freeze two independently cross-checked observation-to-inference cases; test predictions, all estimator inputs, coefficients, vcov, CIs, and scaled focal inference on both backends |

## Public contract

The new argument is appended to `repair()` and `repair_data()` so existing
positional calls retain their meaning. A single referenced model keeps its
existing defaults and prediction object. Passing a cross-model policy for a single
model is rejected to catch stale specifications after a formula change. Unused
models in the input catalog do not affect covariance requirements.

For several referenced models:

1. Each fit must have one grouping factor. Group column names may differ; their
   values must refer to the same units. No row-order matching or ID crosswalk is inferred.
2. `"independent"` explicitly sets between-model prediction-error covariance to zero.
   It does not imply independent latent processes or erase within-model covariance.
3. Supplied blocks have `models = c("a", "b")` and `covariance`, a numeric 3D
   array named by selected coefficients of a, selected coefficients of b, and IDs.
   All pairs are required, including explicit zero blocks. Pair orientation is
   arbitrary but must agree with the array axes. Duplicate pairs are errors.
4. Selected term sets must match exactly. Every retained ID must be supplied;
   extra IDs are allowed. Arrays must be finite and fully named even for IDs
   outside the scope. Invalid specifications do not become omission candidates.
5. Within-model blocks always come from the scorer. The full selected-effect
   covariance is checked at relative eigenvalue tolerance 1e-8. No supplied
   covariance is silently repaired or estimated.
6. `@predictions` becomes a named list only for multiple models. Diagnostics always
   include a named `mlm_models` list and `cross_model_error`, recording
   the policy and original supplied blocks. Final aligned error inputs remain
   in `@err_x`, `@err_y`, and `@err_xy`. The summary prints the policy and each
   MLM's diagnostics.

The existing M-MOM input contract already represents arbitrary within-unit
predictor and outcome errors. The extension therefore belongs in parsing,
assembly, and diagnostics. Multiple equations, cross-person errors, automatic
cross-model covariance estimation, reliability summaries, and uncertainty from
measurement/covariance estimation are outside this milestone.

## Reproduction evidence

`dev/generate_multimodel_fixtures.R` does not load `repair`. It sources the
companion repository's production prediction extraction helper, full M-MOM estimator, primary
`fit_fuller_dual_core()` implementation, and independent algebra audit. It checks
agreement before writing `inst/extdata/multimodel-workflow-reference.rds`.
The fixture records source commit and hashes, package versions, and estimator
settings; ordinary tests need neither the companion repository nor its dependencies.

Two deterministic synthetic datasets have 70 people, unequal cluster sizes,
different time covariates and grouping column names, correlated latent slopes,
and either independent residual errors or residual correlation 0.25. For the
second case, cross-model prediction-error covariance is derived as L_a R_ab L_b'
using the known correlation and fitted marginal residual scales. This is an
illustrative externally supplied covariance model, not an estimate inferred
from the two separate fits. Both cases preserve the predictor model's full
intercept/slope error covariance.

Tests start with repeated observations, refit both nlme models, and compare both
package backends to frozen production results, including all assembled error
inputs, coefficients, full coefficient covariance, SEs, and confidence intervals.
The dual-process population-SD reporting ratio is applied to estimate, SE, and
interval endpoints together. The reference tolerances are 1e-5 to accommodate
optimizer/platform differences; order-invariance tests use 1e-9. These are
per-dataset reproducibility tests, not claims to reproduce published Monte Carlo
coverage or every simulation scenario.

## Design and verification workflow

1. Inspect the merged parser, scorer, assembly, inference, and simulation usage.
2. Define an explicit covariance/identity contract before relaxing the parser.
3. Generalize assembly while preserving the existing single-model API.
4. Generate independent production/algebra fixtures; add alignment, omission,
   invalid-input, mixed-backend, and joint-covariance tests.
5. Regenerate public help; update the README, workflow guide, and design plan;
   add an executable multiple-model vignette explaining statistical assumptions.
6. Run the full tests, render the README, and build/check the source package with
   all vignettes. Record results and any remaining warnings in the PR.

Validation on R 4.5.2 / macOS arm64: 926 assertions passed in both the source
and installed-package test runs, with no failures, errors, warnings, or skips.
README rendering and all three vignette builds/rebuilds passed. Rendered summary
output was inspected for model identities, scope accounting, and covariance
policy. `R CMD check --no-manual` completed with zero errors, zero notes, and
one pre-existing installation warning: R's header refers to the unsupported
Apple-clang warning group `-Wfixed-enum-extension`. Network-restricted dependency
index lookups logged connectivity messages but dependency checks passed with
installed packages. `git diff --check` passed. No estimator or existing reference
fixture was modified.


## Nested nuisance follow-up (2026-10-30)

The original one-group-per-fit restriction above is now relaxed for lmer nuisance
levels verified to be contained within the selected target group. The cross-model
covariance contract is unchanged; shared nuisance levels remain rejected. See
[nested-nuisance-validation.md](nested-nuisance-validation.md).
