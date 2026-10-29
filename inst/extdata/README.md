# Bundled reference data

These RDS files are independent implementation references, not teaching data or
published Monte Carlo tables. For a discoverable teaching dataset use
`data(repair_example)` and its help page; that dataset is stored under `data/`.

| File | Role | Generation script |
|---|---|---|
| `nlme-workflow-reference.rds` | Three observation-to-inference nlme cases | `dev/generate_nlme_fixtures.R` |
| `multimodel-workflow-reference.rds` | Two paired-process cases with independent or supplied cross-model error covariance | `dev/generate_multimodel_fixtures.R` |
| `reliability-reference.rds` | Posterior covariance, marginal and intercept-residualized reliability references | `dev/generate_reliability_fixtures.R` |

Each reference contains provenance with source revisions/hashes and environment
information. Generators use the companion research code rather than repair to
produce expected values. Do not regenerate expectations from the implementation
being tested. See `tests/testthat/fixtures/README.md` in the source repository for
the complete reference inventory and maintenance policy. Normal installed
examples/tests require no companion repository.

## End-to-end nlme REPAIR reference datasets

`nlme-workflow-reference.rds` contains three synthetic repeated-observation
and person-level datasets. It includes model specifications, expected
likelihood-only scores, joint score-error covariances, and structural estimates,
standard errors, confidence intervals, and covariance matrices.

Expectations were generated WITHOUT loading repair, using production code in
`blup_unshrinking` at commit `b1f09f430d2214d5c4c8815b1386b62cb4ae159e`:
`R/blup_helpers.R`, `R/stage2_estimators.R`, and `R/fuller_mlm.R`. The independent
algebra audit in `vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R`
was required to agree. Source hashes, R/nlme versions, estimator controls, and
scaling settings are embedded in `provenance`.

These are exact per-dataset implementation checks, not the paper's published
Monte Carlo table values or a new coverage study. Cases cover balanced predictor
errors, unbalanced outcome errors, and correlated errors on both sides with
AR(1) and heterogeneous residual variances. The last case is a same-model
correlated-error extension, not a cross-model structural fit.

The package reports coefficients on the input scale. The production comparison
uses `skip_internal_scaling = TRUE`. The predictor case additionally records
results multiplied by the generating population slope SD (0.8), matching the
external reporting-scale convention in `vh_study2.R`. All reported estimates,
SEs, and interval endpoints are scaled together. The other cases use scale 1.

Regenerate intentionally with `dev/generate_nlme_fixtures.R`; never replace these
expectations with package-generated outputs. Normal tests and the vignette read
this bundled file and require no sibling repository. Tests refit measurement
models from the raw data, with tolerance 1e-5 to accommodate optimizer/platform
variation; same-fit analytical GLS tests use stricter numerical tolerances.
