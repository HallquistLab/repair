# Stable table extractors: audit, contracts, and verification

Baseline: merged reliability/diagnostics workflow, `58b10b9`. This increment
provides public reporting and inspection tables without changing prediction extraction,
measurement assumptions, M-MOM estimation, reliability, or inference.

## Audit findings and changes

| Finding | Consequence | Correction |
|---|---|---|
| Structural coefficient tables lived inside S7 slots or a larger summary | Reporting scripts depended on internal property names and unrelated summary calculations | Export `repair_coefficients()` with a documented plain-data-frame schema and configurable normal confidence level |
| The saved score container differs for one versus multiple models | Every user's extraction code needed model-count-specific branches and manual ID alignment | Export `repair_predictions()` for fitted and assembled objects, returning explicit model/group/term/ID keys |
| Marginal score SDs alone lose within- and cross-model dependence | A convenient score export could encourage invalid independent-error reuse | Export `repair_prediction_covariance()` from the actual assembled error inputs, retaining covariance source and selected-effect identities |
| Diagnostic access returned an entire list and always calculated reliability | A user needing only exclusions or weights had to compute and navigate unrelated results | Add exact `table` selection while preserving the old list API; skip reliability unless requested |
| Sample exclusions could be mistaken for failed MLM predictions | An analyst might discard valid measurement information or lose failed cases | Separate scope, analyzed membership, score status, and reason; preserve MLM-scope failed rows |

## Public contract decisions

### Structural coefficient table

`repair_coefficients(fit, conf.level = .95)` returns exactly the current documented
columns in order: character `term`, followed by numeric `estimate`, `std.error`,
`statistic`, `p.value`, `conf.low`, `conf.high`, and `conf.level`. Terms follow the
standard coefficient vector: structural intercept, then model-matrix columns,
including named factor contrasts and model-qualified random-effect labels.
Intervals use the existing `confint()` calculation; changing their level changes
no fitted values or p-values. Failed structural fits error rather than presenting
partial coefficients. No standardization or extra uncertainty is implied.

### Score table

`repair_predictions(object, scope = "structural")` accepts a `repair()` fit or
`repair_data()` assembly. Columns are character `id`, `model`, `group`, `term`,
`label`, `role`; numeric `score`, `std.error`, `prediction.error.variance`; character `scope`;
logical `analyzed`; character `status`, `reason`. Selected effects follow formula
selection order, with IDs inside each effect. The key is model/group/term/ID.
Model and coefficient identity are explicit columns, not parsed from display labels.

Default IDs follow retained structural-data order. The MLM scope uses
each model's saved effective grouping levels, including out-of-scope and failed
predictions. Only selected effects are returned. Excluded and failed are distinct:
status is `ok` or `failed`, analyzed is a separate membership indicator, and failed
rows retain saved missing values and explanatory reasons. Absent MLM IDs
have no prediction row and remain discoverable in the exclusion table. Scores are
likelihood-only deviations and their SDs describe conditional sampling error.
They are neither BLUPs nor posterior summaries.

### Prediction extraction-error covariance table

`repair_prediction_covariance(object)` accepts the same fit/assembly types. It returns
character `id`, two sets of character model/group/term/label/role identity columns
(suffixed 1 and 2), numeric `covariance`, and character `source`. Each unordered
pair including the diagonal appears once per analyzed person. Pair order is the
lower triangle in selected-effect order; ID order follows the analysis. The
source distinguishes within-model prediction extraction, explicit assumed independence, and
supplied cross-model covariance, even when a supplied covariance is zero.

Values come from aligned structural error inputs, not a new covariance estimate.
The table reconstructs the full selected-random-coefficient error matrix by
mirroring off-diagonal entries. Observed error-free columns are omitted. There is
no MLM-scope covariance option because supplied cross-model blocks need
only cover analyzed IDs. `vcov(fit)` continues to represent structural coefficient
covariance, a different uncertainty object.

### Diagnostic selection and compatibility

`repair_diagnostics(fit, table = NULL)` preserves the previous list, including
its names, order, metadata, and reliability class. An exact table name returns a
plain data frame with automatic row names. Unsupported/partial names error.
Reliability is calculated only for the complete list or the reliability table.
The covariance policy is list metadata, not a selectable table. Existing typed
zero-row structural tables remain the failure behavior, with status explaining
the missing results. Diagnostic table columns are cataloged in the user guide.

All new tables use explicit identifier columns and automatic row names. No
extractor refits or rescores a model. Coefficient/score/covariance extraction
avoids posterior reliability calculation. Existing S7 objects and standard
`coef`, `vcov`, `confint`, `nobs`, and summary methods remain available. The new
schema contracts are the supported routine reporting interface; any future
schema change should be documented and accompanied by deliberate contract tests.

## Development and verification workflow

1. Inspect the merged object classes, numerical outputs, saved prediction containers,
   scope handling, and diagnostic list before choosing new names or shapes.
2. Define table columns, scalar types, unique keys, row order, scope meaning,
   uncertainty meaning, and failure behavior before implementation.
3. Add a shared fit/assembly resolver and thin extractors reading saved results;
   extend the existing diagnostic function without changing its default interface.
4. Test exact schemas, factor labels, fit/assembly equivalence, custom confidence
   levels, scope permutation, failed predictions, absent IDs, numerical failures,
   covariance reconstruction, source provenance, and portable table serialization.
5. Compare extracted coefficients, SEs, and intervals against existing independent
   dual-process fixtures on both R and C++ backends. Reconstruct all selected-effect
   covariances from exported tables and compare to frozen reference arrays.
6. Write executable reporting/join/export examples, regenerate help/README, and
   run all tests plus source-package build/check with all vignettes.

Validation on R 4.5.2 / macOS arm64: source and installed-package tests pass
1,278 assertions without failures, errors, warnings, or skips. README rendering
and all five vignette builds/rebuilds passed. `R CMD check --no-manual` completed
with zero errors, zero notes, and one pre-existing installation warning from R's
header referring to the unsupported Apple-clang warning group
`-Wfixed-enum-extension`. Restricted-network dependency index lookups logged
connectivity messages; installed dependency checks passed. `git diff --check`
passed. Existing numerical and reliability fixtures are unchanged. The extractor
tests establish interface and numerical agreement, not new coverage properties.

## Deliberate boundaries

No broom dependency or duplicate tidy/glance/augment methods, automatic CSV writes,
new fitted-value or prediction estimands, arbitrary score transformations, or
standardization is introduced. The user guide demonstrates ordinary data-frame
joins and export, with explicit identity keys and the prediction-error/posterior/
structural-uncertainty distinction. Bootstrap inference and automated cross-platform
checks remain separate subsequent work.
