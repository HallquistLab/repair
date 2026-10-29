# Preparing the first release

This checklist is preparation guidance, not a claim that a release or CRAN
submission exists. Current development version: `0.0.0.9000`. Proposed first
GitHub release: `0.1.0`, subject to the maintainers' release decision.

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

For a later CRAN submission, perform the additional incoming, PDF manual,
URL/spelling, package-name, dependency/toolchain, and external-platform checks
against current CRAN policy. Replace `cran-comments.md` with dated results from
that actual submission candidate. Current development CI alone is not a complete
CRAN submission rehearsal.
