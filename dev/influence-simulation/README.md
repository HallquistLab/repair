# Influence-function / bootstrap / legacy stacked SE review bundle

Start with [the implementation audit](../influence-sandwich-design.md) and
[the complete report](report.html), also available as [R Markdown](report.Rmd).
This is a paired exploratory extension of the prior robust-refit screen, with
600 measurement datasets, two structural arms and two truths. It is not new
independent confirmatory coverage evidence. All files are synthetic summaries.
The existing `^dev$` rule excludes this bundle from the installed package.

## Verify the results

From a repair checkout, base R alone suffices:

```sh
Rscript dev/influence-simulation/verify-results.R
```

The verifier independently computes coverage/zero exclusion from the 17,200
retained interval endpoints, verifies all 172 performance cells and matched
decision contrasts, checks unavailable-as-NA behavior, numerical evidence and
bundle checksums. It does not claim to rerun raw native optimizations or derive
the likelihood gradients independently; those checks are in package tests and
the full research validation archive.

- `dataset-intervals.rds`: profile/rep/arm/truth/method keyed original point
  estimates, SE/endpoints, availability/reasons and decisions. Missing intervals
  are not valid non-rejections. Legacy OLS points can differ in unequal designs.
- `performance.csv`: denominators, coverage/rejection Wilson MC intervals,
  interval delivery measures, bias, SD, SE calibration and widths.
- `paired_comparisons.csv`: common-eligible method contrasts and paired MC
  uncertainty. No observed discordance is not proof of equivalence.
- `analytic_diagnostics.csv`, `covariance_decomposition.csv`: singularity,
  derivative/branch reasons and focal signed covariance contributions.
- `legacy_normal_diagnostic.csv`: secondary diagnostic changing only the legacy
  critical value from t to normal; primary comparisons retain the original t.
- `balanced_ols_validation.csv`, `public_bootstrap_validation.csv`,
  `guarded_api_validation.csv`, `validation.txt`: executed identity/API evidence.
- `bootstrap_sensitivity.csv`, `derivative_sensitivity.csv`: prespecified
  first-three/profile B499/B1999 and half-step comparisons.
- `summary.rds`, PNG/PDF figures: report inputs and exportable plots.
- `00check.log`, `00install.out`, `testthat.Rout`: source-package validation;
  2,099 passing assertions, examples/vignettes passed, zero errors/notes and
  one external R-header compiler warning (retained, not suppressed).
- `drivers/`: executed comparison/postprocessing/validation sources. These source
  external research helpers; the GPL legacy implementation is not copied here.
- `protocol.md`, `study-README.md`, `source-manifest.csv`, `raw-manifest.csv`,
  `provenance.rds`, `sessionInfo.txt`: frozen design and numerical environment.
  Paths in bundled manifests are normalized; raw archives remain external.
- `bundle-manifest.csv`: checksums of all other bundle files.

## Scope and reproducibility

Executed numerical code: repair 558bdcead1c15cedd1e332aed6f01e86d49f3721, on
robust-refit head 7770874, R 4.5.2 and lme4 2.0.6. Later pointwise-recipe scope
guards were checked against 72 preselected public calls with identical results.
Do not overwrite the frozen installed numerical library when rebuilding the PR.

The full source/data/native-bootstrap archive is in the research workspace:

```
blup_unshrinking/vig_hallquist_2026/analysis/repair_bootstrap_pilot/
  robust_refit/          # preserved original native fits/draw maps
  influence_comparison/ # executed drivers, 600+18 raw results and HTML report
```

The local archive's README records reproduction commands and dependencies. The
research checkout baseline was b1f09f430d2214d5c4c8815b1386b62cb4ae159e; the study
drivers/raw outputs are local artifacts, not asserted to be at that public
commit. Full reruns require that archive and its unchanged
R/stacked_sandwich_helpers.R and R/derivative_backends.R. Copying the bundled
drivers alone into an arbitrary checkout is not a complete raw-data reproduction.
Cross-platform bitwise agreement is not promised.

Arm A retains noisy-predictor M-MOM; arm B constructs an observed predictor with
known regression on the recovered slope outcome. The old stacked OLS helper is
applicable only to B. Its ML fitting, HC0–HC3 conventions and t_(N−2) intervals
are preserved; REPAIR uses REML and primary IF normal intervals. Under unequal
information M-MOM estimates precision weights but legacy OLS is unweighted.
These differences are explicitly reported rather than described as pure SE gains.

Arm-B focal bootstrap arrays reuse actual previous native fits, their residual
SDs, draws and convergence decisions; every structural fit is recomputed.
Common fixed-slope centering cancels for this focal coefficient. This does not
reconstruct or claim intercept bootstrap inference. Twenty-four complete public
B99 bootstraps verify the shortcut (2,376 attempts, max focal gap 6.66e-16).

The full sandwich is available for 1,351/2,400 analyses, versus 2,400/2,400
stage2 results; 1,048 refusals were singular native fits and one an unstable
derivative. Difficult noisy-predictor nonzero coverage was 73–83% for available
full IF intervals and 79–96% for bootstrap intervals. Bootstrap retains failures
and small/short-panel bias; neither a finite interval nor propagated covariance
establishes calibrated inference. See the report for exact denominators and
matched uncertainty rather than pooling unlike profiles/estimands.
