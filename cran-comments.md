## Candidate

repair 0.5.0, first CRAN submission candidate, checked November 3, 2026.

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

The README executes against the checked installation. The documentation site
build and internal-link validation cover 47 HTML pages.

## Source artifact

File: repair_0.5.0.tar.gz
SHA-256: 5aadd0c101efc5f0ae9705bb5b0768a1468ac03422b84061185ab8e988bbfda5

No CRAN submission or release tag has been created during preparation.
