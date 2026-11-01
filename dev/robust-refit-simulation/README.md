# Review bundle: robust bootstrap refitting

Start with [the audit and results](../robust-refit-audit.md). This compact bundle
preserves the numerical evidence for the 2026-10-30 screen. It is excluded from
the installed package by the existing `^dev$` build rule.

## Verify the published results

From a repair checkout, run (base R only):

```sh
Rscript dev/robust-refit-simulation/verify-results.R
```

`dataset-intervals.rds` contains 6,400 paired old/new dataset-level interval
summaries, keyed by profile, outer replication, focal truth and method. It stores
original estimates, SEs/endpoints, availability, draw/failure counts and decisions,
including the labeled successes-only pathways. It contains synthetic-study
summaries, not raw person records or fitted models. The verifier independently
recomputes coverage and zero exclusion **from the retained endpoints**, checks
all 64 method/profile/effect summaries and before/after availability, and verifies
recovery totals and the prespecified B sensitivity decisions.

- `performance.csv`: robust results, including denominators, Monte Carlo
  uncertainty, bias/SD/SE/width and delivery measures.
- `old_robust_performance.csv`: old versus robust rates, returned-interval counts,
  rescued/lost datasets and original-estimate changes.
- `paired_comparisons.csv`: matched-eligible robust method contrasts and MC
  uncertainty; no observed discordance is not treated as proof of equivalence.
- `optimizer_recovery.csv`: distinct measurement calls/internal optimizer
  attempts, retries, recovered/retained old failures and new rejections.
- `unresolved_measurement_fits.csv`: remaining measurement rejection histories.
- `tail_sensitivity.csv`: the first-three/profile B499/B1999 comparison.
- `profiles.csv`: sample sizes, schedules, covariance/reliability calibration.
- `api_validation.csv`, `stratified_api_validation.csv`, `validation.txt`:
  executed public-API checks and their numerical agreement.
- `sessionInfo.txt`, `source-manifest.csv`, `bundle-manifest.csv`: tested
  numerical environment, original executed source digests, and this bundle's
  checksums. Original source paths are normalized to the study directory.

## Full-run provenance and design

Executed baseline: repair `4cc48f6ebb585f019c55dfaf1b7d920f14c801cd`.
Executed robust fitting code: `3c4f6016de4cc2d034af8f66acb52f5081ff48e1`;
`cb7ed20` subsequently clarified audit prose without changing fitting code.
The PR's later commits add documentation/evidence without changing that numerical
strategy. R 4.5.2 and lme4 2.0.6 were used. Cross-platform bitwise equivalence
is not promised.

The full executed study, raw per-resample histories, draw maps, source snapshots,
HTML report, source package and check logs are retained in the research workspace:

```
blup_unshrinking/vig_hallquist_2026/analysis/repair_bootstrap_pilot/
  paper_like/     # original six-profile screen
  robust_refit/   # frozen numerical rerun and complete reproduction README
```

The research checkout was at `b1f09f430d2214d5c4c8815b1386b62cb4ae159e`.
The study drivers are local research artifacts, not asserted to be present at
that public commit. This PR preserves summary verification, rather than claiming
to include the full raw-simulation archive. A complete optimization rerun requires
that research archive and the matching installed package libraries. Its recorded
commands are `run.R benchmark`, `run.R main`, `run.R sensitivity`, followed by
`summarize.R`, `compare.R`, `validate.R`, and rendering `report.Rmd`. The research
README documents checkpoint/hash matching and `REPAIR_STUDY_LIB` /
`REPAIR_STUDY_WORKERS`. The original and rerun outputs are retained separately.

The shared DGP uses time centered/sample-SD standardized within person,
`Y = .5*time + b0 + b1*time + eps`, random-coefficient SDs .9 and correlation
+/−.5. Gaussian iid residual SD is calibrated using full posterior slope
reliability (profile CSV). The observed structural outcome is
`Z = 1.5 + .4*b0 + (gamma/.9)*b1 + eta`, gamma=0/.40, with eta variance adjusted
so Var(Z)=1. Both random coefficients enter REPAIR. Structural errors are
homoskedastic Gaussian. Reported focal coefficients/SEs/intervals multiply by
the population SD .9, not an estimated SD.

Six profiles have 100 outer replications each; null/alternative share measurement
data and paired structural noise. Across profiles, common random numbers use
L'Ecuyer-CMRG master seed 202610301, successive outer streams; bootstrap seed
240000000 + replication. B499 is primary; the preselected first three/profile
are extended to B1999. Unequal information fixes half at n=3 and half at n=17,
allocated independently of latent values/errors. Its additional stratified
comparator preserves those counts and uses bootstrap seed +700000000. This is
not a package strata option. Excluded timing pilots use replication 801.

Every measurement resample uses actual native optimization. The executed driver
caches only immutable OLS/design-inverse ingredients for exact iid matching-design
scores; fixed effects and residual variance update on every refit and structural
fits use public `fit_mmom()`. Public `repair()`/`repair_refit()` identity checks,
24 complete B99 public bootstrap checks, 20 stratified public refits and independent
raw-rate recalculation verify this route. It is not a closed-form converged-fit
substitute. Original/reference replay also uses the robust recipe.

See [the protocol clarification](protocol-clarification.md): explicit nloptwrap
1e-8 tolerances were already the installed defaults, so its development candidate
increased the budget rather than tightening stopping thresholds. The default
nloptwrap algorithm and minqa BOBYQA are different implementations of the same
optimizer family. The screen reuses the original data and is not independent
confirmatory coverage evidence. General coverage, adaptive M-MOMS, multiple
measurement models and random outcomes remain outside this screen.
