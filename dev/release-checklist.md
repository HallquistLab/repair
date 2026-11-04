# Preparing the first release

Candidate: `1.0.0`, prepared November 3, 2026. Scope: core REPAIR, conditional
inference, sequential step-up regularization, whole-unit bootstrap sampling and
its deterministic robust lmer refitter. PRs 21 and 22 (full influence, leverage,
and delete-one inference) are excluded. Submission and release publication have
not been performed by this checklist.

## Freeze the intended scope and identity

- Review the support/troubleshooting guide and conditional-inference boundaries.
- Retain Zachary A. Vig as primary maintainer and both Zachary A. Vig and Michael
  N. Hallquist as authors/co-maintainers, per the confirmed maintainer decision.
  Retain the designated maintainer email and identify both authors as copyright
  holders in DESCRIPTION and LICENSE.
- Reconcile `inst/CITATION` with the final manuscript; update its status/date/DOI
  only when verified. Check `citation("repair")` in the candidate installation.
- Set the agreed version in DESCRIPTION, consolidate/date NEWS, regenerate help
  as needed, and render README against that candidate. Add version-pinned install
  instructions once a real tag is selected; do not imply it already exists.

## Validate the candidate

- Run the complete package and website workflows at the candidate revision.
  Require all six package environments and the documentation build to pass;
  review notes, examples, all articles, README, and check artifacts.
- Build a source artifact (`R CMD build .`) and check that exact tarball in a clean
  library (`R CMD check --as-cran repair_<version>.tar.gz`), including the PDF
  manual and HTML validation with a suitable Tidy installation. Record its
  SHA-256, Git revision, R/platform/dependency versions, and check logs. A source
  checkout passing tests does not replace checking the release artifact.
- Verify installed teaching-data discovery, `citation("repair")`, author metadata,
  all help examples, and fresh-session persistence. Verify independent fixtures
  have not been regenerated from repair or loosened to achieve passing checks.
- Inspect site home/help/support/citation pages and search at desktop and narrow
  widths. Retain GitHub Actions as the repository's Pages source and confirm a
  successful main deployment and reachable public URLs before advertising them.
- Keep historical audit counts as dated evidence; record new results rather than
  silently rewriting those audits. Resolve genuine warnings. A local toolchain
  warning and clean hosted checks should be reported separately.

## Publish the reviewed artifact

Version 1.0.0 and the scope above are selected. After the maintainers authorize
publication, create its tag and release notes for the tested commit/artifact. The site tracks main,
so say explicitly whether readers are viewing development or release docs.
Release preparation does not tag, publish a release, or submit to CRAN.

Before submitting, review the dated results and candidate hash in
`cran-comments.md`. The Ubuntu/release CI job now includes the PDF manual and
CRAN incoming checks; retain all six platform/R checks and the documentation
workflow on the final revision. Review any new warning or significant note
against current CRAN policy. Submit the tested source tarball through the CRAN
form only after the designated maintainer reviews the candidate.

## November 3, 2026 candidate evidence

The bootstrap/refitter branch is integrated without PRs 21 and 22. The merge
retains the current prediction extractors, `mlm_data` refit API, and unit-copy
schema. Step-up is available through `repair()` and `fit_mmom()`; the default
minimum estimator and existing reference fixtures are unchanged. Each bootstrap
attempt reselects adaptive constants and retains the chosen alphas and guard/
fallback indicators. See `vignette("step-up")` and the independent generator
`dev/generate_step_up_fixtures.R` for algorithm settings and provenance.

The final tarball passed full `R CMD check --as-cran` on R 4.5.2/macOS arm64
with GCC 16.2.0, PDF and HTML manual checks, and network/incoming checks enabled:
zero errors, zero warnings, one new-submission note; 2,308 assertions passed
without failures, warnings or skips. Local validation used a temporary TeX font
tree, HTML Tidy 5.8.0, and a temporary Makevars selecting GCC. No package compiler
flags or diagnostics were suppressed. The earlier Apple-clang/R-header warning
is toolchain-specific; it is not hidden by the package.

The README executes against the installed candidate and pkgdown validates all
47 HTML pages, internal links/assets, and search targets. Candidate hashes and
dated check counts are in `cran-comments.md`; the tarball and logs are retained
locally under ignored `check/cran/`. The automated matrix additionally checks
Linux release/oldrel/devel/R 4.3, macOS, and Windows on the integrated revision.
The release CI now includes PDF-manual and incoming checks on Linux release.

## Computation-only step-up update, November 3, 2026

This update supersedes the preceding local candidate: full source-tarball
`R CMD check --as-cran` passes with zero errors, zero warnings, and the single
new-submission note, including PDF/HTML manuals. All 3,875 assertions pass with
no failures, warnings or skips. Exact old/new structural and public-bootstrap
comparisons, deterministic operation counts and environment-specific timing
results are documented in `dev/README.md` and `dev/step-up-efficiency.csv`.
The frozen oracle is revision 60c46ec. Statistical guard/fallback changes remain
deferred, and PRs 21 and 22 remain excluded. The updated exact artifact hash is
in `cran-comments.md`; the previous local bundle is retained in ignored
`check/cran-60c46ec/` and the current candidate in `check/cran/`.

## Version 1.0.0 candidate, November 3, 2026

Version 1.0.0 supersedes the dated 0.5.0 candidate evidence above. NEWS is
consolidated for the first release; README and the website identify it as a
candidate, without claiming CRAN availability. Both authors are explicitly
identified as copyright holders under the unchanged MIT license; the method
citation remains unpublished and reflects the November 3 manuscript draft.

The exact source archive passes full local `--as-cran` checks with PDF/HTML
manuals and incoming checks: zero errors, zero warnings and only New submission.
Both current and declared older dependency stacks pass all 3,875 assertions.
README execution and validation of all 47 site pages pass. See
[version 1.0 evidence](release-1.0-validation.md) for dependency bounds and logs;
`cran-comments.md` records the tested archive hash. The previous local bundle
is preserved under `check/cran-914e613/`; the 1.0.0 bundle is `check/cran/`.

Win-builder R-devel upload has a receipt for the exact checked archive; the
check result is pending delivery to the designated maintainer email. Preserve
and review that result before treating Win-builder as passed. No CRAN submission,
release tag or release publication has been performed. Statistical changes and
PRs 21 and 22 remain excluded.
