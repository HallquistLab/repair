## Candidate

repair 1.0.0, first CRAN submission candidate, checked November 3, 2026.

Scope: core REPAIR, conditional model-based inference, sequential step-up
(M-MOMS) regularization, and experimental whole-unit bootstrap inference with
an optional deterministic robust lmer refitter. Full influence-function,
leverage, and delete-one inference are excluded.

## Test environment and results

* macOS 27.2, arm64; R 4.5.2; GCC 16.2.0 (C++20), TeX Live 2026,
  HTML Tidy 5.8.0.
* R CMD build followed by R CMD check --as-cran on the exact source tarball.
  CRAN incoming checks, dependency/URL checks, PDF manual, HTML manual,
  examples, tests, and vignette rebuilds were enabled.
* 0 errors | 0 warnings | 1 note

The note is "New submission".

Installed-package tests: 3,875 passing assertions, with no test failures,
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
and passes all 3,875 assertions with no errors or warnings.

The README executes against the checked installation. The documentation site
build and internal-link validation cover 47 HTML pages.

## Source artifact

File: repair_1.0.0.tar.gz
SHA-256: b55a0256bb626094e7a2349c3c1e2995d89fb2dc5a4a147192a2f69200f54929

## Win-builder

The exact archive above was uploaded to Win-builder R-devel on November 3, 2026
(2026-11-04 02:18 UTC). The upload page returned a receipt for the matching
filename and 664,027-byte size. Results are pending and are sent to the
designated maintainer email in DESCRIPTION; a receipt is not a passing check.

No CRAN submission or release tag has been created during preparation.
