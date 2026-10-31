# Fresh leverage / delete-one / bootstrap review evidence

Read [the method design and audit](../leverage-jackknife-design.md) and the
complete [HTML report](report.html), also available as [R Markdown](report.Rmd).
This bundle covers six profiles ×100 **fresh** measurement datasets, with paired
null/effect structural outcomes. It is an exploratory validation screen, not a
claim of uniform validity or a substitute for more precise coverage studies.
All person-level records here are synthetic. Existing `^dev$` build exclusion
keeps this review material out of the installed package.

## Verify independently

From a repair checkout, using base R alone:

```sh
Rscript dev/leverage-simulation/verify-results.R
```

The verifier recalculates coverage/zero exclusion from all 12,200 endpoints,
checks 122 performance cells and matched decision contrasts, strict missing-
interval handling, actual centered delete-one covariance and uncentered HC3
covariance from retained person vectors, signed contributions, public-API checks,
fresh/old input-hash disjointness, sensitivity record counts and bundle checksums.
It does not claim to rerun native optimization or independently derive the
leverage equations; those are validated by package tests and the full archive.

- `dataset-intervals.rds`, `performance.csv`: original points, endpoints,
  availability/reasons, coverage/rejection Wilson MC intervals, bias/SD, SEs,
  widths and delivery per planned data. Missing intervals remain NA.
- `person-shifts.rds`: all actual delete-one coefficient vectors/usability flags
  and available analytic HC3 shifts, aligned by person ID, for both truths.
  Jackknife variance centers vectors; HC3 does not center shifts. All methods
  retain the original point estimate; no jackknife bias correction is applied.
- `paired_comparisons.csv`: common-eligible method contrasts with paired MC
  uncertainty. No discordance is not evidence of equivalence.
- `analytic_diagnostics.csv`, `leverage_refusals.csv`: boundary/derivative/
  conditioning gates and explicit failure reasons, including augmented HC0
  agreement with the existing full influence calculation.
- `covariance_decomposition.csv`, `leverage_exact_comparison.csv`: signed direct,
  measurement and cross terms, and analytic versus actual deletion shifts.
- `branch_performance.csv`, `deletion_diagnostics.rds`, `deletion_statistics.csv`:
  original branch strata, actual native singularity, floor/correction changes
  and diagnostic jackknife bias. Subset rates are descriptive, not causal.
- `public_api_validation.csv`: fixed-subset complete public API checks, including
  24 B99 bootstrap calls and 12 whole-person jackknives, all coefficient vectors.
- `bootstrap_sensitivity.csv`, `derivative_sensitivity.csv`: preselected
  first-three/profile B499/B1999 and half-step comparisons.
- `posthoc_extreme_audit.csv`: separately labeled numerical audit of the 12
  largest HC3/HC0 variance ratios at h/2 and h/4. It never replaces main results
  or tunes thresholds; all 24 remain available, max relative change 1.06e-4.
- `freshness.csv`, `previous-input-hashes.csv`: 600 distinct new inputs and no
  identical old main-study input. Seed labels alone are not the only check.
- `summary.rds`, PNG/PDF figures: report inputs and exportable scientific plots.
- `bias_plot_*`, `matched_plot_cohorts.csv`, `matched_coverage_plot_data.csv`,
  `width_tail_plot_data.csv`, `width_exceedance_plot_data.csv`, `studentized_*`:
  auditable inputs to the bias, 50–100% coverage, and untrimmed tail figures.
  `plot-validation.txt` and `drivers/check-plots.R` verify the plotted summaries,
  common eligibility cohorts, all observed maxima and exact exceedance counts.
- `protocol.md`, `study-README.md`, `drivers/`, source/raw manifests,
  `provenance.rds`, `sessionInfo.txt`, `validation.txt`: frozen design, commands,
  numerical environment and executed checks. The drivers source external
  research helpers; their GPL implementations are not copied into this package.
- `initial-00check.log`, `initial-00install.out`, `testthat.Rout`,
  `ci-validation.json`: 2,164 passing assertions and implementation-head CI.
  Initial local logs retain an external R-header warning and the tail() namespace
  note; the qualification fix passed all six remote R-platform checks. These
  labeled initial logs are not presented as a later zero-note local rerun.

## Numerical provenance and interpretation

Numerical baseline is repair `cbffe31`, on influence PR #21 head `6bb14dc`.
The later `utils::tail` qualification fixes a namespace check note without
changing the equations. The frozen installed simulation library remains
unchanged. R 4.5.2 and lme4 2.0.6 were used; cross-platform bitwise reproduction
is not promised.

Full native optimization histories and generated data remain in:

```
blup_unshrinking/vig_hallquist_2026/analysis/repair_bootstrap_pilot/leverage_fresh/
```

The local archive's README records complete reproduction commands and external
sources. The compact bundle permits verification without hundreds of thousands
of model fits; copying its drivers alone into an arbitrary checkout does not
reproduce the full optimization archive. Earlier robust/influence studies are
preserved unchanged. Master seed 202610319 and bootstrap seeds 340000000+rep
are distinct from the earlier experiment. Development replication 801 is excluded.

The HC3-style system representation includes measurement equations, raw/weighted
moments, a mean-weight normalization and final coefficient equations. It retains
cross-stage dependence and exact deleted REML terms. Original N/policy are fixed
in the analytic correction; actual delete-one refits N-1 and updates finite-N
terms. Normal Wald inference is primary; N-1 t intervals are secondary heuristics,
not derived Satterthwaite tests. No failed vector/sensitivity is omitted, no
positive-only variance addition is made and no correction is chosen by coverage.

The OLS HC3 and augmented HC0 identities establish numerical anchors. They do
not establish finite-sample unbiasedness of this representation-dependent
nonlinear correction. The report evaluates coverage together with eligibility,
bias, width, detection and numerical stability.
