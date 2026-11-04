# Optional first-stage MLM specifications

Date: November 4, 2026. The original October 20 outline proposed fitting both
stages in one call. `repair_mlm()` now supplies the optional first-stage fitting
specification, while the outer formula continues to describe the structural
regression. The helper name follows the existing MLM terminology and the
maintainer's explicit naming choice.

## Scope and implementation

- `repair_mlm(formula, data, engine = c("lmer", "lme"), ...)` creates an unfitted
  S3 specification. Native named options are validated against the selected
  fitter's formal arguments; implicit partial matching in `...` is rejected.
- Repeated data and evaluated non-data arguments are retained at construction.
  lmer weights/offset/subset and lme subset remain native data expressions.
  Native defaults are applied when fitting. Each public consuming call fits
  each specification once; specifications are not a mutable fitting cache.
- `normalize_repair_models()` accepts bare specifications and uniquely named
  lists mixing specifications and fitted inputs. Existing fitted objects are
  returned unchanged. Ordinary model objects then enter the existing extraction,
  alignment and structural pipeline.
- Native calls use a retained data/formula environment. Stored call metadata
  resolves helper-local argument symbols to their actual values, preserving
  controls, random/covariance objects and data expressions for supported native
  replay. Replay evaluation now also accepts evaluated literal lists, preserving
  their attributes; symbolic external dependencies retain their existing limits.
- No changes to M-MOM calculations, adaptive selection, inference, optimizer
  retry policies, or bootstrap failure rules accompany this API addition.
  PRs 21 and 22 remain excluded.

## User contracts

`repair()`, `repair_data()` and `repair_effects()` consume specifications. A bare
input is named `mlm`. Reuse a result's `@mlm_fit` for repeated structural analyses
or discovery without another native fit. Model-level data are not aggregated into
structural rows; callers still supply the ID-keyed structural table and any
cross-model prediction-error covariance policy explicitly.

The native fitter's warnings remain visible. Errors identify the named model
and engine. Only fitted Gaussian lmer/linear lme results are supported. For lme,
an explicit random-effects specification and `keep.data = TRUE` are required.
All existing grouping, residual-covariance and dependence restrictions remain.

Refits/bootstrap still require explicit original data and pass reference identity,
factor/schema, convergence and effective-row checks. Data-only native recipes can
replay; external vectors, custom transformations or unsupported fitting choices
can require an explicit callback. This is not a general mechanism for serializing
arbitrary global R dependencies. Custom fitting strategies such as the robust
lmer refitter remain available through fitted-model inputs.

## Validation

`tests/testthat/test-mlm-spec.R` covers deferred fitting and per-call counts,
unchanged fitted inputs, named mixtures and explicit cross-model policy,
R/Rcpp agreement with manual native fits, complete structural calculations,
ML/REML, frozen control objects and data, factor contrasts, weights, offsets,
subsetting and missing observations, AR(1)/varIdent residual and pdDiag random-effect lme structures, native replay,
paired public bootstrap draws/results, and fresh-process reconstruction/replay
for both engines. Invalid arguments and a real native rank-deficiency warning
exercise error/warning boundaries. Existing independent estimator fixtures are
unchanged.

Executable examples are included in installed help, README, the main workflow
and nlme guide. Candidate source-archive checks, updated counts, the artifact
hash and hosted results are recorded in `cran-comments.md` and the local ignored
`check/cran/` bundle. Earlier candidate/Win-builder results identify their own
archive and do not validate a subsequently rebuilt archive. No CRAN submission
or release tag is created by this implementation.

## November 4 candidate evidence

Source archive: `repair_1.0.0.tar.gz`, 674,611 bytes. SHA-256:
`6abb6caf09c3bff5a0be7c9d77252c49b2a380879884061aae37174bcdee78dc`.

Full `R CMD check --as-cran` on this archive reports zero errors, zero warnings,
and two NOTEs: New submission and unavailable external current-time verification.
All 4,029 assertions pass with zero failures, warnings or skips. Incoming checks,
examples, vignette rebuilds and PDF/HTML manual validation run. A separate archive
check recompiles against the declared older dependency stack and passes the same
4,029 assertions with no errors or warnings and only the clock-verification NOTE.
The installed README executes; pkgdown validates 48 HTML pages, internal links,
assets and search targets.

Win-builder R-devel received this exact archive on November 4 at 14:34 UTC
(09:34 America/New_York), returning HTTP 200 with the matching filename and byte
count. The result URL is pending delivery to the designated maintainer email.
This receipt does not establish a passing Win-builder check. The earlier source
archive and receipt remain under `check/cran-1de2d65/`; current logs and the new
receipt are in `check/cran/`. No CRAN submission or release tag is created.
