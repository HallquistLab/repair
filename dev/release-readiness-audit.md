# First-release readiness audit

Date: 2026-10-29. Status: planning audit, not a release approval or new statistical
validation. Suggested initial target: a tagged GitHub `0.1.0` release. CRAN
submission is a separate decision and has additional checks below.

## Scope and evidence

This audit inspected the public API, documentation, fixture provenance,
representative tests, metadata, workflows, tags/releases, and Pages state.
It uses the already completed numerical and hosted checks; it does not claim a
new exhaustive code review, coverage simulation, performance benchmark, or
lower-bound dependency compatibility study.

Snapshot:

- Merged main: `150e56fa2afe2fe9c1ab75e8bc7b1ce9d1822abf`, through PR #7.
- Website PR #8: `4a93e7c173fcd4835bdc18145c15c76e374200a7`, open at audit time.
  Its documentation build is verified, but its changes are not yet on main.
- DESCRIPTION version remains `0.0.0.9000`; no remote tags or GitHub releases were
  returned by the repository API at audit time.
- GitHub Pages is not configured in the inspected repository metadata. The Pages
  endpoint/setup request returns HTTP 404 with the connected account, which has
  push access but no administrator/maintainer role. Public deployment has not
  been verified; enabling the publishing source requires an authorized owner.

The [main package-check run](https://github.com/zacharyvig/repair/actions/runs/40910349364)
completed successfully. On the website PR, all six
[package-check environments](https://github.com/zacharyvig/repair/actions/runs/40911936623)
reported `Status: OK`, each with 1,278 passing assertions and zero failures,
warnings, or skips. The [website workflow](https://github.com/zacharyvig/repair/actions/runs/40911936762)
built and validated 28 HTML pages, links, assets, and search targets; deployment
was correctly skipped for the PR. Local browser checks covered home/reference/
article rendering, search navigation, and home/guide layout at 390px width.

These counts are dated evidence for the revisions above, not evergreen package
claims. Refresh the checks on the actual release candidate.

## What has been accomplished

| Increment | Outcome | Evidence and boundaries |
|---|---|---|
| [PR #1](https://github.com/zacharyvig/repair/pull/1): scoring and inference correction | Likelihood-only scoring and the M-MOM covariance calculation were corrected and reconciled with independent references | Six frozen structural-input cases, both R/Rcpp backends, intermediate matrices/weights, and the weak-information variance regression test in `tests/testthat/test-repair-reference.R` |
| [PR #2](https://github.com/zacharyvig/repair/pull/2): coherent workflow | `repair_effects()`, `re()`, `repair_data()`, and `repair()` form a usable discovery/assembly/fitting path with ID-based alignment and explicit omissions | [Design audit](design-audit.md), formula/cohort tests, generic methods, and the main workflow vignette |
| [PR #3](https://github.com/zacharyvig/repair/pull/3): nlme support | Validated single-group linear lme extraction, including supported correlated/heterogeneous residual structures | [nlme audit](nlme-validation.md); three observation-to-inference reference datasets and analytical GLS checks |
| [PR #4](https://github.com/zacharyvig/repair/pull/4): multiple models | Model-qualified effects and explicit cross-model scoring-error policies, with full within-model covariance retained | [Multiple-model audit](multimodel-validation.md); independent and supplied cross-covariance cases, paired-model ID/covariance validation |
| [PR #5](https://github.com/zacharyvig/repair/pull/5): reliability and diagnostics | Manuscript-defined posterior reliability, cohort accounting, numerical/weight diagnostics, plots, and summary integration | [Reliability audit](reliability-validation.md); independent person-level and aggregated references |
| [PR #6](https://github.com/zacharyvig/repair/pull/6): stable extractors | Plain coefficient, ID-keyed score, covariance, and diagnostic tables with documented schemas | [Extractor contracts](extractor-contracts.md); ordering, failure, covariance reconstruction, and table serialization tests |
| [PR #7](https://github.com/zacharyvig/repair/pull/7): package automation | Six OS/R configurations, installed tests/examples/vignettes, README execution, artifacts, and a stable aggregate check | [CI guide](ci-checks.md); strict warning policy passed on hosted runners |
| [PR #8](https://github.com/zacharyvig/repair/pull/8): website, still open | pkgdown home page, task-based reference, five executable guides, search, NEWS, preview artifacts, and main-branch deployment workflow | Build verified; merge, Pages configuration, and first public deployment remain outstanding |

The main numerical achievement is implementation agreement with independently
produced reference calculations, including the production simulation estimator
and algebra audit. Tests cover predictor, outcome, joint/correlated,
heterogeneous-error, and weak-information inputs. End-to-end cases refit
measurement models from repeated observations. Fixtures are bundled and normal
tests do not require the companion repository.

This does not mean the package reproduces the entire published Monte Carlo study
or has newly established coverage in every supported setting. The fixture
provenance explicitly distinguishes per-dataset reproduction from simulation
performance. Conditional inference remains conditional; bootstrap refitting is
only a [future-version design](bootstrap-refit-design.md).

## Overall assessment

The package has enough coherent core functionality for a bounded first release.
The next highest-value work is release hardening and improving the path from
installed help to a successful analysis. A major estimator extension is not
needed before releasing the existing, accurately documented conditional-inference
workflow.

The audit found concrete release/documentation gaps below. It did not identify a
new demonstrated numerical defect in the inspected paths. That statement is not
an exhaustive numerical correctness guarantee. Retain the existing independent
fixtures and cross-platform acceptance tests during preparation.

## Recommended release prerequisites

### 1. Establish a versioned, citable release identity

Evidence: DESCRIPTION and NEWS still identify `0.0.0.9000`; remote tags/releases
are empty; there is no `inst/CITATION`. `Authors@R` currently lists Zach Vig alone
as author/maintainer, and LICENSE identifies the generic `repair authors` holder.
The audit does not establish that those credits are wrong; they need an explicit
maintainer decision before publication.

Recommended work:

- Confirm release channel and scope; use `0.1.0` as the proposed first GitHub
  release rather than implying a CRAN release has occurred.
- Confirm author/contributor roles, maintainer address, and copyright holder with
  the maintainers. Do not infer author order or silently add names.
- Supply a reviewed `inst/CITATION` distinguishing the software from the method
  manuscript. Cite the actual manuscript status/version and DOI only if verified;
  do not invent a publication reference.
- Date and consolidate release NEWS, record the release commit and reference
  provenance, and test the exact candidate source tarball before creating a tag.
- Provide version-pinned installation instructions for the release alongside the
  development installation. A website tracking main is not a frozen release manual.

Acceptance: `packageVersion()` and `citation()` communicate the intended release;
installation of the tagged/source artifact works in a clean library; release
notes state supported scope and inference limitations. Tagging/publishing remains
a future authorized release action, not an action performed by this audit.

### 2. Complete publication and verify the candidate on hosted CI

Evidence: website PR #8 is open; its build passes but Pages publication has not
occurred. Its DESCRIPTION website/issue URLs are also still pending merge.

Merge the website changes, have an authorized repository administrator select
Settings -> Pages -> Source -> GitHub Actions, and verify the first main-branch
Build documentation/Deploy documentation run and public links. Include navigation,
search, code outputs, and a narrow-width check on the actual public site. If the
site is advertised in the release, it should be reachable before announcing it.

Run the full package/site workflows on the exact release candidate. `Package
checks` and `Build documentation` can be made required checks; their existence
does not itself establish that branch protection is configured. Settings were not
changed by this audit. Preserve warning failures and fixture tolerances.

The hosted macOS checks do not reproduce the earlier local Apple-clang/R-header
warning. Keep those two environments distinct instead of waiving warnings or
calling the local warning a demonstrated package-wide defect.

### 3. Finish installed help and define the public compatibility boundary

Evidence: NAMESPACE exports 12 named objects. A scan of their corresponding Rd
files found examples only for `repair()` and `repair_coefficients()`; the other
10 export topics lack an examples section. The `class_repair_mod` topic also lacks
a Value section. This is a discoverability/documentation gap despite the five
executable guides and does not mean those functions lack implementation tests.

Add small executable examples where users seek help: `repair_effects()`,
`repair_data()`, `repair_scores()`, `repair_score_covariance()`,
`repair_reliability()`, and `repair_diagnostics()` first. Reuse one compact setup;
show what the returned object means, not only how to print it. Also cover `re()`,
`fit_mmom()`, and `repair_ranef_gls()` appropriately. Update roxygen sources, then
regenerate Rd rather than editing generated files directly.

Decide whether the exported S7 class is a supported constructor or an advanced
implementation detail. Document it accordingly without removing an export
casually. Clearly distinguish stable high-level functions/tables, advanced
estimation interfaces, internal slots, and experimental multi-group scoring.
Document error/failure behavior and future compatibility promises for each tier.

Acceptance: a new user can complete fit/inspect/report operations from installed
help without reading source or navigating S7 internals. Standard summaries and
stable extractor schemas retain their tested contracts.

### 4. Publish a compact support and inference matrix

Evidence: boundaries are present in README, vignettes, and audits, but distributed
across them. The current main workflow guide still contains the old sentence
rejecting separate-model equations; website PR #8 corrects it. Historical audit
sections should remain historical, with follow-up links, rather than being
mistaken for current support statements.

Add one authoritative user-facing matrix covering Gaussian single-group lmer,
supported single-group linear lme residual structures, multiple models and their
covariance policy, observed predictors, effect identities/scales, omissions,
regularization, and inference. Link to it from README/help/site. Reconcile it with
tests and the final manuscript reference before freezing the release.

Keep these distinctions prominent:

- Likelihood-only score error, posterior reliability, and structural coefficient
  uncertainty are different quantities.
- Scores are random deviations; observed structural predictors are assumed
  measured without error; reported coefficients retain input scale.
- Conditional normal intervals do not propagate measurement-model estimation or
  uncertainty in estimated supplied cross-model blocks.
- Minimum regularization stabilizes an estimator; it is not a guarantee against
  weak information or poor coverage. `none` still corrects measurement error.
- Single-group model support is not support for dependent/crossed analysis people.

Acceptance: all entry points agree about supported models, inference, and
experimental boundaries; reproducing fixture outputs is not described as proving
nominal coverage.

## Valuable small additions before or just after the first release

| Priority | Addition | Why it helps | Bounded acceptance criterion |
|---|---|---|---|
| High | One documented synthetic example with repeated observations and meaningful person-level variables | The beginner sleepstudy example deliberately invents a predictor; existing bundled RDS files are reference fixtures rather than an obvious teaching dataset | A simple discoverable dataset/list and a short fit-inspect-report example, with known simulation provenance and no suggestion that invented variables are real measurements |
| High | A troubleshooting guide tied to existing diagnostics | Users will encounter missing IDs, rank-deficient person designs, singular measurement fits, weak corrected moments, and concentrated weights | Each common message points to a diagnostic/action; never recommend arbitrary exclusion, variance deletion, or covariance-policy changes just to make a fit succeed |
| Medium | Clean-session fitted-object persistence check | Existing extractor tests round-trip plain tables through RDS, but there is no corresponding explicit saved-fit/fresh-process test in the inspected suite | Save/read a fit in a new R process, load repair, and verify summary/generics/table outputs; document that cross-version S7 object compatibility is not automatically promised |
| Medium | Installation metadata/dependency review | Source installation uses C++20 and BLAS/LAPACK/Fortran via Makevars; DESCRIPTION has no SystemRequirements field and mostly unversioned dependency requirements | Review appropriate toolchain metadata and genuinely required package API minima; test any claimed minima instead of inventing version constraints |
| Medium | Clear bundled-data provenance index | `inst/extdata/README.md` describes nlme data, whereas the directory also contains multimodel and reliability references; the test fixture README has fuller coverage | Installed-data documentation identifies every fixture, its role, source revision, and regeneration policy |
| Optional | A small, versioned benchmark | Refitting/diagnostics and covariance-table exports can become costly as people/effects grow; no benchmark evidence was reviewed | Record time/memory over representative n/effect counts before claiming performance; avoid optimizing without a measured problem |

The first two additions have the strongest direct user value. They can share the
same documented example and support matrix rather than introducing another
parallel reporting API. A new broom layer, broad plotting framework, or automatic
model-selection interface is not needed to make the current workflow usable.

## What to defer

- Person bootstrap/refit inference until its replay, cohort, failure, and
  scientific-validation contracts are implemented: see the separate design draft.
- Multi-group/crossed structural dependence and cross-person covariance support.
- Automatic estimation of cross-model scoring-error covariance.
- Random-coefficient interactions, nonlinear transformations, joint equations,
  and automated standardization without explicit error-propagation contracts.
- Prediction/fitted-value APIs until the target is specified (latent outcome,
  observed score, existing person, or new person are not interchangeable).
- Additional model families, residual structures, and backend expansion without
  independent references and tests.

These are substantive methodological or API projects, not small omissions to
fill immediately before a first release. A clearly bounded release is more useful
than an expanded support claim without the corresponding evidence.

## If the first release is intended for CRAN

The current CI is strong development evidence but intentionally skips PDF manual
creation and uses r-lib's default disabling of CRAN incoming checks. It is not a
complete submission rehearsal. `cran-comments.md` currently contains the generic
`0 errors | 0 warnings | 1 note` text and no dated environment/results record;
regenerate it from the actual submission candidate rather than treating it as
verified submission evidence.

For CRAN, additionally review the current repository policy/new-package guidance,
build the PDF manual, exercise the final package with appropriate incoming checks,
inspect URLs/spelling/docs and portable dependencies, and use suitable external
R check services. Confirm package-name availability and maintainer/citation/license
metadata. Evaluate dependency minima and check time/size from the candidate.
No CRAN name reservation, submission, or complete compliance audit has been done.

Useful release guidance: [R Packages: Releasing to CRAN](https://r-pkgs.org/release.html).
Recheck the [official CRAN policy](https://cran.r-project.org/web/packages/policies.html)
when preparing an actual submission; this document does not reproduce or certify
its full requirements.

## Proposed next implementation work

The next coherent PR should be **first-release documentation and reproducibility
hardening**: installed-help examples, the shared support/troubleshooting guide,
reviewed citation/metadata, a simple documented example, and a clean-session
saved-fit test. Keep estimator calculations and independent expectations fixed.
Authorship and manuscript citation details require maintainer decisions; prepare
the concrete metadata for review rather than guessing.

Then prepare the release candidate: merge/publish the website, reconcile the
final paper/reference version, set the agreed version, validate the exact source
artifact and hosted checks, and create a reviewed tag/release. Bootstrap work can
proceed afterward using the separate design and its explicit acceptance gates.
