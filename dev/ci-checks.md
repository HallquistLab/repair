# Automated package checks

`.github/workflows/R-CMD-check.yaml` runs on every pull request, pushes to `main`,
and manual dispatch. There are no path filters: documentation, tests, source,
fixtures, and workflow changes all receive the same checks. The workflow uses
pinned revisions of the official [r-lib actions](https://github.com/r-lib/actions)
and [checkout action](https://github.com/actions/checkout).

## Check matrix and failure policy

| Runner | R version | Purpose |
|---|---|---|
| Ubuntu | release | Full CRAN check including PDF manual, incoming checks, and executable README |
| macOS | release | Apple compiler, BLAS, and platform behavior |
| Windows | release | Rtools/C++20 and Windows package behavior |
| Ubuntu | oldrel-1 | Previous R release compatibility |
| Ubuntu | devel | Detect changes in upcoming R |
| Ubuntu | 4.3 | Test the minimum minor version declared in DESCRIPTION |

All six jobs are required by the final **Package checks** job. That stable job
name can be selected as a required status check in branch protection. This PR
does not modify repository branch-protection settings. R-devel failures are
visible failures, not silently accepted experimental results. The minimum-R job
uses the latest patch release within 4.3; it does not test every patch release.

Each job installs the full declared dependency set, including Suggests and
LinkingTo, and Pandoc. Dependency caching is handled by setup-r-dependencies.
`_R_CHECK_FORCE_SUGGESTS_=true` and an explicit prerequisite script prevent missing
nlme/test/vignette dependencies from quietly removing coverage. The setup action
prints package versions and session information for reproducing failures.

The check builds a source package and runs:

```r
rcmdcheck::rcmdcheck(
  args = c("--no-manual", "--as-cran"),
  build_args = "--no-manual",
  error_on = "warning",
  check_dir = "check"
)
```

Examples, installed-package tests, and all vignettes are enabled. Existing tests
explicitly exercise both R and Rcpp estimator backends; the workflow does not
substitute an R-only test run or regenerate frozen expectations. Ubuntu/release
installs TinyTeX and manual fonts, omits `--no-manual` from both build and check,
and explicitly enables `_R_CHECK_CRAN_INCOMING_=true`. Other matrix members
retain the faster `--no-manual` check above. Errors and
warnings fail the job. Notes remain visible but are not automatic failures;
new-package incoming checks can emit contextual notes that require review.
The final submission must still be made from the exact reviewed tarball, with
its dated check results recorded in `cran-comments.md`.

The Ubuntu/release job additionally executes README.Rmd against the package
installed in `check/repair.Rcheck`, rather than a development namespace or an
unrelated installed version. Rendering occurs in a temporary directory and does
not update committed generated documentation. Other executable guides are
already covered by the package vignette checks on every matrix member.

## Logs and troubleshooting

Every matrix member uploads its check directory on success or failure under a
unique `check-<runner>-<R-version>` artifact name, when check output exists.
The directory includes `00check.log`, installation output, installed test output,
and example/vignette diagnostics. Dependency/setup failures occur before that
directory exists; consult the corresponding workflow step's log and session
information. README failures appear in their own step log.

1. Open the failed matrix job and find the first failed step.
2. For package failures, download its check artifact and read `00check.log`,
   `00install.out`, and `tests/testthat.Rout*` in the package's `.Rcheck` directory.
3. Reproduce with the same R/dependency versions and the commands below.
4. Resolve extraction, compiler, numerical, or documentation failures from their
   evidence. Do not regenerate fixtures to match a failure, broadly relax numeric
   tolerances, or disable warnings merely to make a runner green.

Earlier local macOS checks emitted an installation warning from R's header about
`-Wfixed-enum-extension` on the local Apple compiler. The workflow deliberately
does not suppress it or globally downgrade warnings: hosted-runner results must
be checked independently. Any actual CI exception requires a separate, documented
cause and narrowly scoped resolution.

## Local reproduction

Use a checkout with an R build toolchain (C++20, BLAS/LAPACK/Fortran; Rtools on
Windows) and Pandoc. To install check dependencies with pak:

```r
install.packages("pak")
pak::pkg_install(c("deps::.", "any::rcmdcheck"), dependencies = TRUE)
```

Then, from the package root:

```sh
Rscript .github/scripts/check-prerequisites.R
```

```r
Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = "true", "_R_CHECK_CRAN_INCOMING_" = "false")
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"),
  build_args = "--no-manual", error_on = "warning", check_dir = "check")
```

After a successful check:

```sh
Rscript .github/scripts/check-readme.R
```

`check/` is ignored by Git and excluded from source builds, so rerunning a check
cannot accidentally package its previous output. The workflow itself can be
validated with [actionlint](https://github.com/rhysd/actionlint):

```sh
actionlint .github/workflows/R-CMD-check.yaml
```

When updating pinned actions, verify the upstream revision, run actionlint, and
let the full matrix complete. `release`, `oldrel-1`, and `devel` follow upstream
R versions; R 4.3 intentionally stays fixed until the package minimum changes.
All jobs have time limits, older runs on the same PR/ref are canceled when
superseded, and matrix failures do not cancel the other diagnostic runs.

Ubuntu/release and Ubuntu/devel have 90-minute job limits; the other package
matrix members retain 45 minutes. The first Bayesian plug-in run (PR 23,
Actions run 41856471343, revision 51f3376) exceeded the old 45-minute limits:
release spent 38m54s downloading system dependencies (70.5 MB at 30.2 kB/s),
and devel spent 37m20s installing dependencies, including 19m41s compiling
rstanarm. Devel then passed all 4,136 assertions but timed out during vignette
rebuilding; release timed out during tests. Native sampling was not reached.
The larger cold-start budget preserves every dependency, package/manual check,
README execution, and the native Bayesian integration step. A timeout is not a
test assertion failure, and the aggregate correctly fails if any member is
canceled. Hosted results must still demonstrate that these jobs complete.

Vignette links to development plans must resolve before a PR is merged. Use a
published immutable commit for a plan not yet on main, rather than a main-branch
URL that returns 404 during incoming checks.

The follow-up run at revision 59435cd completed all six package checks. The
Ubuntu/release native sampling step then failed before fitting brms because
RStan could not find Boost. Binary dependency resolution can omit LinkingTo-only
packages, so the native job explicitly installs BH and RcppEigen and verifies
their actual Boost/Eigen header files and RStan's configured include directories
before starting checks. The brms and rstanarm native tests are separate, so a
failure in one engine cannot prevent exercising the other. Both remain required.

The same follow-up's documentation job spent about 28 minutes downloading system
dependencies and reached its 30-minute limit just as website construction began.
Its budget is now 90 minutes as well, preserving examples, articles, link checks
and preview artifacts. These changes address compiler prerequisites and cold
setup time; they do not disable sampling or relax numerical assertions.

## Scope and audit

The package previously relied on manually executed macOS checks. The new workflow
addresses missing platform/R-version coverage, missing-dependency skips,
inaccessible failure logs, and the lack of a consistent merge-check name. It
uses read-only repository permissions, ordinary pull-request events, and no
publishing credentials. It does not deploy packages, change branch protection,
create recurring app tasks, or modify estimator behavior.

Implementation validation combines actionlint, the prerequisite/README scripts,
a local source-package check, and the first hosted matrix run. Results and any
platform findings are recorded in the PR. Independent numerical reference
fixtures remain the acceptance evidence for scientific behavior.


## Website checks

The separate pkgdown workflow builds the package website on PRs and validates
its generated pages, internal links, and assets. Main-branch builds deploy to
GitHub Pages. See [website.md](website.md) for source editing, HTML previews,
publishing prerequisites, and local verification.
