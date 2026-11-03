# Preparing the first release

Candidate: `0.5.0`, prepared November 3, 2026. Scope: core REPAIR, conditional
inference, sequential step-up regularization, whole-unit bootstrap sampling and
its deterministic robust lmer refitter. PRs 21 and 22 (full influence, leverage,
and delete-one inference) are excluded. Submission and release publication have
not been performed by this checklist.

## Freeze the intended scope and identity

- Review the support/troubleshooting guide and conditional-inference boundaries.
- Retain Zachary A. Vig as primary maintainer and both Zachary A. Vig and Michael
  N. Hallquist as authors/co-maintainers, per the confirmed maintainer decision.
  Confirm the current primary email and license holder before publication.
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
  library (`R CMD check --as-cran --no-manual repair_<version>.tar.gz`). Record its
  SHA-256, Git revision, R/platform/dependency versions, and check logs. A source
  checkout passing tests does not replace checking the release artifact.
- Verify installed teaching-data discovery, `citation("repair")`, author metadata,
  all help examples, and fresh-session persistence. Verify independent fixtures
  have not been regenerated from repair or loosened to achieve passing checks.
- Inspect site home/help/support/citation pages and search at desktop and narrow
  widths. Enable GitHub Actions as the repository's Pages source and confirm a
  successful main deployment and reachable public URLs before advertising them.
- Keep historical audit counts as dated evidence; record new results rather than
  silently rewriting those audits. Resolve genuine warnings. A local toolchain
  warning and clean hosted checks should be reported separately.

## Publish the reviewed artifact

After the maintainers choose the release version and publication scope, create
its tag and release notes for the tested commit/artifact. The site tracks main,
so say explicitly whether readers are viewing development or release docs.
This hardening PR does not tag, publish a release, or submit to CRAN.

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
