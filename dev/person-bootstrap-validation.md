# Person bootstrap: implementation audit and verification workflow

Date: 2026-10-30. Baseline: validated refitting PR #13, `e263b9a`.
Implements step 2 of the [bootstrap design](bootstrap-refit-design.md), within
the refitter's single-model scope. Scientific coverage validation is outstanding;
the public inference label is experimental.

## Findings and decisions

| Finding | Decision | Verification |
|---|---|---|
| Fixed scores/weights omit measurement-estimation variability | Call public `repair_refit()` for every attempt, including full reference and replacement fitting | Independent manual lmer, nested-slope and correlated/heterogeneous lme loops |
| Duplicate people can merge nuisance effects | Use validated deterministic draw/re-keying; sample uniform source indices with replacement | Full effective histories, unequal counts, nested effects and source maps |
| Sampling to a success target conceals unstable samples | B counts attempts; retain all failures and maps, no replacement/retry/trim | Natural rare-level loss and explicit measurement failures |
| Failure selection changes the distribution | Withhold covariance/intervals by default; explicit accessor-level `allow_failed = TRUE` labels success-conditional inference | Common-vector covariance, table/matrix metadata and inspectable unavailable summaries |
| Bootstrap output has different inference semantics | Separate S3 result and dedicated table extractors; no bootstrap z-statistics/p-values | Generic dispatch, schemas and fresh-process persistence |
| Fitter RNG use/failure can affect later samples | Save independent L'Ecuyer-CMRG streams for preflight and attempts; restore caller state | Prefix stability, RNG-consuming recipe, error restoration and absent seed |
| Hidden normal-generator state is not fully saved in .Random.seed | Reject Box-Muller and user-supplied caller RNGs before mutation | Explicit state-preserving rejection test |
| Errors need stage and warning provenance | Add private condition fields to refitter failures without changing existing messages | Failed measurement/structural stages and the existing refit suite |
| Percentile definitions differ between implementations | Use R quantile type 7 and sample covariance on the same complete vectors for all terms | Independent quantile/covariance checks and coefficient selection |

## Actual API and contracts

```r
boot <- repair_bootstrap(fit, measurement_data, data, id = NULL,
  refitter = NULL, B = 1999L, seed = 1L, tolerance = 1e-6)
coef(boot)       # original estimate
vcov(boot)       # empirical covariance, not added to conditional variance
confint(boot)    # marginal percentile intervals
repair_bootstrap_coefficients(boot)
repair_bootstrap_diagnostics(boot, "replicates")
repair_bootstrap_diagnostics(boot, "draws")
```

The coefficient table gives original estimates, bootstrap SDs, interval endpoints,
level/method, attempt/usable/failed counts, success-conditional status and an
experimental label. Ordinary `repair_coefficients()` is unchanged. Covariance and
interval matrices carry count/status/method attributes. All terms use the same
complete finite coefficient vectors, with at least two required for inference.
This is a numerical minimum, not adequate precision or a validity guarantee.

The result retains original estimates/conditional comparison, the B-by-p matrix
(failed rows all NA), integer source-draw indices, character source cohort,
attempt diagnostics, reference audit, seed streams, versions and call. It retains
neither every replicate fit nor copied raw datasets. Saved extraction needs no
original fitting bindings; source maps are derived without parsing generated IDs.

Scope matches `repair_refit()`: one Gaussian lmer with contained nuisance effects
or supported single-group linear lme. Measurement/analyzed cohorts must coincide;
extra or structurally omitted measurement people require a revised reference fit.
Preflight verifies identity and draw eligibility before sampling. Rows/order,
coding/schema, ML/REML, recorded recipes, regularization policy/custom alpha and
computation backend follow the refitter. Sample-dependent corrections/weights
are recomputed. Sampling conditions on the chosen cohort and observed missingness.

The public refitter repeats reference validation for each attempt. This costs an
extra reference fit but avoids an unchecked private fast path. A future prepared
recipe optimization must retain the same checks and independent test agreement.

## Failures, warnings and RNG

Measurement nonconvergence, score/schema failure and structural failure invalidate
the whole vector. Singular-but-usable lmer fits and informational warnings remain
eligible. nlme singularity remains unknown. Diagnostics identify attempt, stage,
reason, intended row/person counts, analyzed count when available, singularity,
control provenance, realized correction multipliers, elapsed seconds, and character
warning/message list columns. Unavailable diagnostics remain NA. Extreme finite
estimates remain; there are no automatic fallbacks or within-attempt omissions.
Summary stays inspectable when inference is withheld. Opting into successes does
not make selection ignorable. Interruptions stop and restore RNG; no recovery file
or unexpected raw-data log is written.

The seed argument defaults to 1. Preflight uses stream 0; attempt b uses the b-th
`parallel::nextRNGStream()` successor, with Inversion normals/Rejection sampling.
Each attempt samples before fitting. RNG consumption/failure cannot change later
draws. Kind and seed, including absence, are restored on success/error. Hidden
Box-Muller/user-supplied states are explicitly unsupported. Reproducibility is
within the same numerical environment and explicit deterministic recipe, not
cross-version bitwise optimizer equivalence.

## Verification workflow

1. Read the existing design/refit contract and official R RNG/quantile semantics.
2. Compare with separately written sampling, record copying, lmer/lme fitting,
   REPAIR estimation and quantile/covariance calculations.
3. Exercise natural factor-loss failures, explicit fitting failures, stage/warning
   records, gates, RNG restoration, nested identities and saved results.
4. Keep manuscript/simulation fixtures and tolerances unchanged; run the complete
   suite and source checks with all examples and executable vignettes.
5. Install the candidate and execute README, then build/link-check pkgdown
   sequentially with package builds to avoid generated-file races.
6. Require the existing strict six-environment matrix/documentation build; record
   exact results and local toolchain limitations in the PR. No CI policy is relaxed.

## Scientific gate and deferred functionality

Software tests establish implementation agreement, not nominal coverage. Before
promoting this out of experimental status, prespecify simulations spanning sample
size, reliability, irregular histories, correlated random coefficients, weak
moments, variance boundaries and regularization transitions. Compare bias,
empirical SD versus estimated SE, widths, coverage, failure rates and whole-run
interval availability, with Monte Carlo uncertainty. Report the probability that
an attempted analysis both produces and covers with an interval, avoiding selected
successful subsets. Increasing B cannot repair invalid assumptions or selection.

Multiple models/supplied covariance, strata, higher-level clusters, parallelism,
checkpoints/resume, replicate model retention, BCa, studentized and simultaneous
inference remain outside this PR.

Sources: [R RNG streams](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/RngStream.html),
[R quantiles](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html),
and [refitting audit](refitting-validation.md).
