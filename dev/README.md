# Development plans and audits

These documents describe design decisions, evidence, and future work. Historical
milestone audits retain their original scope; later follow-ups supersede earlier
support boundaries. README and executable vignettes describe current user-facing
behavior. A design draft does not make its proposed functions available.

## Current planning

* [Bayesian plug-in full refitted-workflow simulation plan](bayesian-plugin-simulation-plan.md):
  prospective inference validation and software/native sampling checks; no
  unconditional coverage or Bayesian structural posterior is claimed.

- [First-release readiness audit](release-readiness-audit.md): accomplished work,
  dated CI evidence, merged versus pending website status, concrete release gaps,
  short-term user improvements, and a separate CRAN preparation path.
- [Person-level bootstrap and refitting design](bootstrap-refit-design.md): a
  future-version proposal covering resampling identity, replay/scope contracts,
  cross-model covariance, failure handling, inference, and validation. Not a
  first-release commitment and not an implemented API.
- [Package design plan](package_plan.Rmd): implemented interfaces, statistical
  boundaries, and the development/acceptance workflow.

## Implementation evidence

- [Version 1.0.0 hardening](release-1.0-validation.md): release identity, attribution,
  tested dependency floors, source-archive checks and Win-builder status.
- [Design audit](design-audit.md): original workflow corrections and follow-ups.
- [nlme validation](nlme-validation.md): effective rows, residual covariance,
  independent prediction extraction and end-to-end reference checks.
- [Multiple-model validation](multimodel-validation.md): effect identities,
  cross-model covariance contracts, paired processes, and reproduction.
- [Reliability validation](reliability-validation.md): manuscript definitions,
  scope/undefined-value behavior, and independent references.
- [Extractor contracts](extractor-contracts.md): public table schemas, ordering,
  keys, covariance reconstruction, and failure behavior.
- [CI checks](ci-checks.md): platform/R-version matrix, artifacts, local
  reproduction, and warning policy.
- [Fixture provenance](../tests/testthat/fixtures/README.md): independent
  numerical/reliability reference sources and deliberate regeneration policy.

The website work is tracked in [PR #8](https://github.com/zacharyvig/repair/pull/8).
Its build is verified but publication is still pending as of the release audit;
consult that audit's timestamp rather than treating this as a live status page.

## Step-up computational equivalence and timing

The computation-only optimization is compared with the frozen estimator at
`60c46ece17071f7401717a2eefe1a169e5a598d6`; see the test fixture provenance.
The comparison requires identical selected alphas, complete coefficient tables
(including SEs and intervals), covariance matrices, weights, intermediate steps,
status/failure decisions, and guard measurements at every evaluated point. It
covers both R/Rcpp backends, all independent step-up fixtures, mixed policies,
custom grids, zero errors, failed fits, determinant-root branches, thresholds at
floating-point boundaries, and heterogeneous joint errors with 1/3/4 predictors.

Matched public bootstrap runs pair 18 attempted samples across robust lmer,
nested lmer, AR(1)/varIdent nlme, and natural rare-factor failures. An observer
captures every structural fit during both versions' public bootstrap calls,
including reference replays and every replacement, so each replicate's SEs,
covariance and intermediate calculations are compared as well as the saved
coefficient vectors. Draw maps, streams/caller RNG, selected constants, warnings,
failure decisions, percentile intervals and covariance extraction must match
exactly. Elapsed times and unevaluated candidate histories are excluded.

Run `dev/benchmark_step_up.R` against the installed candidate to reproduce
[step-up-efficiency.csv](step-up-efficiency.csv); the accompanying session file
records the environment. The same-process microbenchmark uses five timing
repetitions, both backends, independent N=80 fixtures and an N=800 workload made
by repeating synthetic inputs (not additional independent statistical evidence).
Representative operation counts:

| Case | Complete fits, old -> optimized | Determinant-root solves, old -> optimized |
|---|---:|---:|
| Minimum already passes | 17 -> 1 | 34 -> 2 |
| Both stages refine | 39 -> 30 | 78 -> 16 |
| No passing candidate | 17 -> 17 | 34 -> 11 |

In this environment the median structural-fit speed ratios were about 8--12x,
1.6x and 1.1--1.3x, respectively. The larger repeated-input workload improved
about 1.6--1.9x. These timings are descriptive and environment-specific; they
are not end-to-end bootstrap speedups. Measurement-model fitting is unchanged.

The cache is local to one normalized dataset and one search. Preliminary-control
changes invalidate all weight-dependent quantities; final candidates reuse only
quantities invariant to final alpha. Every evaluated candidate retains the full
coefficient/SE eligibility checks and their original failure ordering. Ascending
early stopping preserves the first passing coarse point and its preceding
failure bracket, without assuming global monotonicity of candidate availability.
Histories report only points actually evaluated. Statistical guard changes,
centering/scaling, sample-size-dependent thresholds and fallback changes are
explicitly deferred.
