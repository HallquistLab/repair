# First-release documentation and reproducibility hardening

Date: 2026-10-29. Based on main `a80bebe` after website PR #8. This implements the
bounded next step in the [readiness audit PR #9](https://github.com/zacharyvig/repair/pull/9).
The development version remains `0.0.0.9000`; a release/tag is a separate step.

## Findings and corrections

| Finding | Correction | Evidence |
|---|---|---|
| Ten of twelve exported-object help topics lacked executable examples | Add a compact example to every missing topic; explain the advanced S7 class and its return value | Generated Rd; R CMD check examples |
| Beginner examples invented a sleepstudy predictor with no shared teaching dataset | Add `repair_example`, with repeated learning scores, practice hours, and documented generating truth | `R/data.R`, `dev/generate_example_data.R`, README and support-guide execution |
| Supported model and inference limits were scattered across guides | Add one support/troubleshooting matrix, interpretation guidance, and a diagnostic-to-action table | `vignettes/support-and-troubleshooting.Rmd`, linked from help, README, plan, and site |
| Only plain extracted tables had explicit persistence tests | Save single lmer and paired nlme fits, read in fresh processes, compare summaries/generics/all public result tables | `tests/testthat/test-persistence.R`; includes an omitted person and supplied cross-model covariance |
| Software and manuscript had no explicit citation file | Add separate software and unpublished-method references | `inst/CITATION`; software version comes from installed metadata |
| Maintainer credits needed an explicit decision | Credit Zachary A. Vig and Michael N. Hallquist as authors and co-maintainers, with Zach retaining the primary contact | User's explicit instruction; DESCRIPTION, package help, README, site authors page |
| Source installation used C++20 without declaring it | Add `SystemRequirements: C++20` | Existing Makevars toolchain requirement; no untested package-version minima invented |
| Installed fixture README described only nlme references | Index nlme, multiple-model, and reliability fixtures and distinguish teaching data | `inst/extdata/README.md`; frozen expectations remain unchanged |
| Generic CRAN comments resembled completed submission evidence | Replace placeholder counts with explicit development status and a release checklist | `cran-comments.md`, `dev/release-checklist.md` |

## Design decisions and workflow

1. Preserve the estimator, tolerances, and independent reference fixtures. This
   increment strengthens access, documentation, and reproducibility; it does not
   extend inferential scope or implement the proposed bootstrap.
2. Reuse one teaching example across installed help, README, and the support
   guide. The 60 people have 8–12 observations each, so repeated and structural
   sample sizes are visibly different. Practice predicts slope deviation in the
   generating model. Truth is documented for teaching, never passed into a fit.
3. Preserve the existing advanced exports. Explain how they differ from the
   routine `repair()` and stable table workflow instead of breaking callers.
4. Use `callr` to cross an actual process boundary. The collector has a base
   environment and receives only a saved file, not the original test bindings.
   Installed checks load the checked package; `test_local()` loads that source
   checkout with `pkgload` so an older installed version cannot mask regressions.
   Both packages are declared test dependencies. Comparisons cover values,
   identities, cohorts, covariance sources, printed summaries, and diagnostics;
   the paired-model coefficients/SEs also retain independent-reference agreement.
5. Regenerate help from roxygen sources, install the package, execute README and
   examples/articles, build and link-check pkgdown, then run the unchanged strict
   six-environment package-check matrix and website workflow on the PR.

R permits one primary maintainer contact (`cre`). Zach retains that role and his
existing address; both authors' metadata comments identify co-maintenance.
The manuscript citation was checked against the title page of the companion
repository's `documentation/paper_drafts/27oct2026/` draft: *Regression among
random effects*, Zachary A. Vig and Michael N. Hallquist. It is explicitly an
unpublished manuscript dated October 27, 2026, without an invented DOI or journal
publication claim. Reconcile it with the final manuscript before a release.

## Acceptance and remaining boundaries

The new tests check same-version persistence, including recomputed diagnostics
from stored measurement fits. They do not promise serialization compatibility
across future S7/package/dependency versions. Keep analysis code and session
information and export plain tables for durable reporting.

Fixture agreement is implementation evidence, not a full rerun of the paper's
Monte Carlo study or a coverage guarantee. Structural inference still conditions
on fitted measurement quantities and any supplied error-covariance inputs.
Dependency lower bounds beyond the existing R minimum have not been systematically
tested; adding speculative minima would not establish compatibility.

Website publication remains distinct from a successful build. The first main
[deployment after PR #8](https://github.com/zacharyvig/repair/actions/runs/40930760883)
failed with HTTP 404 and an instruction to enable Pages; the documentation build
succeeded. An authorized administrator must select **Settings → Pages → Source →
GitHub Actions**, then verify deployment and the public site. Do not treat a PR
preview as publication. See [website operations](website.md).

The PR validation record should identify the tested commit, exact results, local
toolchain limitations, and hosted run links. The [release checklist](release-checklist.md)
requires a fresh check of the eventual source artifact, not reuse of old counts.
