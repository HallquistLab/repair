# Version 1.0.0 release hardening

Prepared November 3, 2026. This is the first release candidate, not a CRAN
submission or a published release tag. Core REPAIR, step-up and experimental
whole-unit bootstrap remain in scope. PRs 21 and 22, influence-function,
leverage, delete-one and all statistical rule changes remain excluded.

## Metadata and documentation

- Consolidate unreleased development changes into NEWS for 1.0.0 and update the
  README/site version without implying CRAN availability.
- Retain Zachary A. Vig as designated maintainer. Name Zachary A. Vig and
  Michael N. Hallquist as copyright holders in both license files and Authors@R;
  the MIT license is unchanged.
- Cite the estimator's measurement-error foundation, Fuller (1987), in
  DESCRIPTION. The title/year/DOI are verified against the
  [publisher's record](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470316665).
  This reference does not describe the package's adaptations as Fuller's exact
  unmodified estimator. Retain a separate unpublished REPAIR manuscript citation,
  updated against the companion repository's November 3, 2026 draft. No journal,
  publication status or manuscript DOI is inferred.
- Update the refitting guide to link the implemented bootstrap and remove
  obsolete manual-skipping advice from the release checklist.
- Install HTML Tidy in the full release CI job so PDF, HTML and incoming checks
  are all exercised, rather than recording an unavailable-validator note.

## Dependency review

The supported older stack is tested together in an isolated temporary library:

| Dependency | Declared/tested floor |
|---|---|
| S7 | 0.2.0 |
| lme4 | 1.1-35.5 |
| Rcpp | 1.0.12 |
| RcppArmadillo | 0.12.8.4.0 |

These are conservative supported floors, not claims that every earlier version
is incompatible. S7 0.1.0 and 0.1.1 emit method-registration warnings when loading
the package source; 0.2.0 removes these warnings. The floor stack retains the
installed versions of the other dependencies (including Matrix, RcppEigen,
rlang and cli), recorded with the logs. This is not an exhaustive Cartesian
dependency-version study or a claim about every old dependency on every R version.

The complete source suite passes on this stack. The archive check additionally
recompiles the package against the older Rcpp/Armadillo headers. Reproduce with
`Rscript dev/check_dependency_floors.R /path/to/repair_1.0.0.tar.gz`; compiler
and system prerequisites must already be installed. The script uses a new
temporary library and retains its results separately from the user's library.
Local validation uses R 4.5.2/macOS arm64 and GCC 16.2.0, without suppressing
compiler diagnostics. The normal hosted matrix exercises current dependencies
on Linux, macOS, Windows, R-devel and R 4.3.

Downloaded CRAN archive SHA-256 values:

```
b8675a7fac7a396e524b21cd353ef0823d2acf76088b5f229d2a55a182a4d49b  S7_0.2.0.tar.gz
7d6664db7ea96429562efe1058da58985d779d6fe79ec6f4e86ba68047135170  lme4_1.1-35.5.tar.gz
0c7359cc43beee02761aa3df2baccede1182d29d28c9cd49964b609305062bd0  Rcpp_1.0.12.tar.gz
cabf6865073204c184ddad4dec1f30aae1b8020b4c9d6e74d1a6b0baaef03046  RcppArmadillo_0.12.8.4.0.tar.gz
```

## Source archive and external checks

Archive: `repair_1.0.0.tar.gz` (664,027 bytes). SHA-256:
`b55a0256bb626094e7a2349c3c1e2995d89fb2dc5a4a147192a2f69200f54929`.

The exact archive passes full `R CMD check --as-cran` with current installed
dependencies: zero errors, zero warnings, only the New submission NOTE. PDF and
HTML manuals, incoming/dependency/URL checks, examples, 3,875 assertions (zero
failures, warnings or skips), and vignette rebuilds pass. The README runs against
the checked installation. pkgdown builds and validates 47 HTML pages, internal
links/assets and search targets.

The second archive check compiles against the older headers and passes the same
3,875 assertions, examples and vignette rebuilds with no errors or warnings. It
uses `--no-manual` because the complete manual check is covered above. Its two
environment NOTEs are unavailable current-time verification under the restricted
sandbox and an Apple `xcrun_db` temporary cache; neither occurs in the unrestricted
full check. Complete logs, compiler header paths and dependency versions are
retained in the ignored `check/cran/` bundle.

Win-builder R-devel received this exact archive on November 3, 2026
(2026-11-04 02:18:28 UTC). The upload page returned the matching filename and
664,027-byte size (HTTP 200). The response and upload metadata are retained in
`check/cran/win-builder-receipt.txt` and `win-builder-upload.json`. The result is
pending; the receipt does not establish that the package passed. Win-builder sends
its result URL to the maintainer address in DESCRIPTION; uploading there requests
a Windows test build, not a CRAN submission. A receipt is not a passing check.
