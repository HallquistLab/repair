# Fresh leverage / delete-one / bootstrap comparison

The frozen [protocol](design.md) specifies six profiles × 100 **new** measurement
datasets, two paired structural truths, B499 main / B1999 first-three/profile,
and the exact HC3-style system representation. This is an independent-stream
exploratory validation screen, not a replacement for the earlier experiments.

## Reproduce

Numerical package baseline: `cbffe31`, built on full-influence PR #21 head
`6bb14dc`. Frozen installed library: `/private/tmp/repair-leverage-lib`.
R 4.5.2 / lme4 2.0.6. Do not overwrite that library while validating manifests.
Source, reference and binary hashes are checked on resume. Later documentation
and evidence commits do not change this numerical implementation.

From this directory:

```sh
REPAIR_STUDY_WORKERS=2 Rscript run.R benchmark
REPAIR_STUDY_WORKERS=12 Rscript run.R main
REPAIR_STUDY_WORKERS=4 Rscript run.R sensitivity
Rscript summarize.R
Rscript check-plots.R
Rscript validate.R
Rscript audit_extremes.R
Rscript -e 'rmarkdown::render("report.Rmd")'
```

Set `REPAIR_STUDY_LIB` to a matching library when rerunning elsewhere. The
`audit_extremes.R` command is separately labeled post hoc: it checks the 12
largest reported HC3/HC0 variance ratios at smaller steps, without replacing
main estimates or changing thresholds. The main sensitivity subset remains
the prespecified first three datasets/profile.

`summarize.R` calls `plots.R`; the latter can also run on existing summaries
without recomputing models. It exports bias, null/nonzero and common-eligible
coverage, interval-width quantile/exceedance, and studentized-tail figures as
PNG/PDF files, with plotted numerical inputs in CSV. Coverage views use 50–100%
and include every Monte Carlo bound. `check-plots.R` verifies bias summaries,
eligibility intersections, coverage, tail maxima and exact exceedance counts
from the retained endpoints. It accounts for CSV rounding at empirical jumps.

The
generator/calibration and native bootstrap driver are sourced unchanged from
`../robust_refit/engine.R` and the research repository. Full reruns therefore
need this research archive; the compact PR bundle is a review/verification
archive, not a claim to contain all native bootstrap histories or external GPL
helper implementations. Cross-platform bitwise equivalence is not promised.

The new L'Ecuyer master seed is 202610319, distinct from the old 202610301;
bootstrap seeds are 340000000+replication. Streams are paired across profiles
and null/effect outcomes. Development pilots use replication 801 and B9 and
are excluded. Main/sensitivity input hashes are compared against all 600 old
main inputs to verify that no dataset was reused. No failed dataset is replaced.

## Archive and computational route

`results/main/native` contains generated synthetic data, actual native bootstrap
draw maps, complete coefficient arrays, sigma/optimizer/convergence histories,
frozen and unequal-profile stratified comparisons, and public API checks.
`results/main` contains analytic results, actual delete-one arrays and native
histories, interval records and API checks. The sensitivity phase mirrors this
layout at B1999 with half-step diagnostics. Both retain source/library provenance
and session information. The excluded pilot archive remains separate.

For speed the exact iid matching-design score calculation caches immutable
individual OLS/design-inverse ingredients. Every native bootstrap/delete-one
measurement fit is genuinely optimized; its fixed coefficients/sigma update and
the entire structural estimator is re-estimated. All complete coefficient
vectors and usability flags are checked against the public APIs on fixed
subsets: 24 B99 bootstrap calls (2,376 attempts) and 12 whole-person jackknives
(760 deletions). This is not a closed-form converged-measurement substitute or
a scalar focal-only eligibility gate.

HC3-style covariance uses uncentered analytic deletion-shift outer products,
with original N fixed in structural policy terms. Actual jackknife refits N-1
people and uses centered deletion cross-products scaled by (N-1)/N. Their point
estimate remains original REPAIR. Normal Wald intervals are primary; N-1 t
references are prespecified heuristics. Bootstrap-SD Wald intervals use the
ordinary resample SD, **not** a studentized percentile-t bootstrap.

Coverage/zero-exclusion denominators use available intervals. Missing intervals
are NA and are never treated as valid non-rejections. `results/summary` contains
independently verifiable endpoints/rates, matched MC contrasts, branch/boundary
strata, leverage/exact-deletion comparisons, signed covariance contributions,
API and freshness validation, sensitivity tables, plots and manifests. The
full report distinguishes numerical availability, bias and interval calibration.
