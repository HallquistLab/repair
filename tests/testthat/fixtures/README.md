# REPAIR numerical reference fixtures

`repair-reference.rds` contains six synthetic inputs and expected results
produced independently of this R package. Metadata in the RDS records source
paths, hashes, source-repository commit, estimator settings, and R version.

The companion `blup_unshrinking` repository at commit
`b1f09f430d2214d5c4c8815b1386b62cb4ae159e` supplied:

- `R/fuller_mlm.R`: full M-MOM matrices and intermediate quantities;
- `R/stage2_estimators.R`: primary simulation coefficient, SE, and interval;
- `vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R`: independent
  algebra-audit coefficient, SE, and interval.

Cases cover predictor-only, outcome-only, both-sided, heteroskedastic, correlated
error, and weak-information inputs. The weak-information case distinguishes the
manuscript's modified covariance bread from the former unregularized bread.

Ordinary tests read the fixture without the sibling repository. Regeneration is
an intentional maintenance action using `dev/generate_repair_fixtures.R`; inspect
that script's usage and review any resulting fixture change. Do not regenerate
expected outputs from the package being tested.
