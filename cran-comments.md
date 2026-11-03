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

Installed-package tests: 2,308 passing assertions, with no test failures,
warnings, or skips. Existing independent numerical fixtures are unchanged.
New step-up references compare both computation backends with the separate
research implementation, including sequential alpha selection and fallback.

The README executes against the checked installation. The documentation site
build and internal-link validation cover 47 HTML pages.

## Source artifact

File: repair_0.5.0.tar.gz
SHA-256: de59565ac833f315fd17021d78ecdc6ee61a68888d12b45db1621f84bf9eba91

No CRAN submission or release tag has been created during preparation.
