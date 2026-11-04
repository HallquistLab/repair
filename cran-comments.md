## Candidate

repair 1.0.0, first CRAN submission candidate, checked November 4, 2026.

Scope: core REPAIR, optional first-stage lmer/lme fitting specifications through
repair_mlm(), conditional model-based inference, sequential step-up
(M-MOMS) regularization, and experimental whole-unit bootstrap inference with
an optional deterministic robust lmer refitter. Full influence-function,
leverage, and delete-one inference are excluded.

## Test environment and results

* macOS 27.2, arm64; R 4.5.2; GCC 16.2.0 (C++20), TeX Live 2026,
  HTML Tidy 5.8.0.
* R CMD build followed by R CMD check --as-cran on the exact source tarball.
  CRAN incoming checks, dependency/URL checks, PDF manual, HTML manual,
  examples, tests, and vignette rebuilds were enabled.
* 0 errors | 0 warnings | 2 notes

The notes are "New submission" and "unable to verify current time". The latter
is an external time-verification issue; it is retained in the check log rather
than suppressed. PDF/HTML manual validation, tests and vignettes pass.

Installed-package tests: 4,076 passing assertions, with no test failures,
warnings, or skips. Existing independent numerical fixtures are unchanged.
New step-up references compare both computation backends with the separate
research implementation, including sequential alpha selection and fallback.
Computation-only changes also match the frozen pre-optimization estimator
exactly for alphas, coefficients, SEs, covariance, intermediate calculations
and failure decisions. Matched public bootstrap runs cover robust/nested lmer,
correlated-residual lme, and natural failed draws, with identical draws,
results, inference eligibility and RNG preservation. Candidate histories omit
unevaluated larger grid points; statistical criteria and fallback are unchanged.

A second archive check recompiles against the declared older dependency stack
(S7 0.2.0, lme4 1.1-35.5, Rcpp 1.0.12, RcppArmadillo 0.12.8.4.0)
and passes all 4,076 assertions with no errors or warnings.

The new optional first-stage interface is compared with manual native fitting
for both engines, including full structural calculations, controls, contrasts,
weights/offsets, effective rows, lme covariance structures, native replay, paired
bootstrap results and fresh-session persistence. Fitted-model inputs and the
existing estimator remain unchanged.

Bare and quoted random-coefficient references produce identical structural
fits, selected constants, prediction/covariance tables, replay and bootstrap
results. Names are parsed literally; model/group ambiguity remains an error.

The README executes against the checked installation. The documentation site
build and internal-link validation cover 48 HTML pages.

## Source artifact

File: repair_1.0.0.tar.gz
SHA-256: 3cf98f3d57b5df3e64838b1183a7ef9523b01362ed799ea68b27c16450fcb171

## Win-builder

The most recent Win-builder R-devel upload (November 4, 2026, 14:34 UTC)
applies to the previous candidate archive, SHA-256
6abb6caf09c3bff5a0be7c9d77252c49b2a380879884061aae37174bcdee78dc.
Its receipt is retained with that candidate under check/cran-bb0ed01/; the result
URL is sent to the designated maintainer email and has not been reviewed here.
That upload does not validate the current archive above. Recheck the exact final
archive on Win-builder before CRAN submission.

No CRAN submission or release tag has been created during preparation.
