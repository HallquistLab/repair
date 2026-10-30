# Validated refitting: design, audit, and acceptance

Date: 2026-10-30. Baseline: merged nested-nuisance PR #12, main `a84c353`.
This implements the single-model replay foundation proposed in the
[bootstrap/refit design](bootstrap-refit-design.md). Random bootstrap inference
remains a later increment.

## Findings and decisions

| Finding | Decision | Evidence |
|---|---|---|
| An original fitted call can refer to missing or mutable external objects | Require explicit reference inputs; reconstruct native data-only recipes or accept an explicit fitting callback | Missing original bindings and fresh-process replay tests; nonliteral recipe rejection |
| Replacement data alone cannot establish identity with the original estimator | Separate reference `measurement_data`/`data` from `new_measurement_data`/`new_data`; verify the former before fitting the latter | Deliberately changed reference observations, structural values, and cohort order fail |
| Original subsets and NA padding are not a safe replay recipe | Select original effective measurement rows by fitted row names; require all replacement rows to be used | Excluded/shuffled observation cases; callback row-drop and replacement NA rejection |
| Contrasts can change with global options or lost factor levels | Record structural factor matrices at assembly, freeze native fixed-factor coding, and require stable coefficient schemas | Changed contrast defaults, sum coding, absent/new levels, and random-effect identity tests |
| Covariance/control choices are not all retained as complete backend metadata | Compare recorded recipes and observable configuration; distinguish retained optimizer checks from explicitly declared callback recipes | Optimizer/ML/REML/covariance/contrasts checks; symbolic control and callback tests |
| Duplicate source people can merge target or nuisance levels | Add deterministic whole-person copying with fresh target/nuisance IDs and an explicit source map | Person:day, nuisance slopes, multiple nuisance groups, and standalone visit-ID tests; manual-refit agreement |
| Successful refits need provenance and actionable failure handling | Return an ordinary REPAIR object with a summary audit; stop on unusable stages and retain informational warnings/messages | Failure-code, callback error, singular-but-usable, warning and persistence tests |

## Concrete interface

`repair_refit(object, measurement_data, data, new_measurement_data = NULL,
new_data = NULL, id = NULL, refitter = NULL, draw = NULL, tolerance = 1e-6)`.

The two unprefixed data arguments always describe the original reference fit.
With no replacements, return the refitted identity replay. Replacements default
to the reference inputs for whichever side is unchanged. A deterministic `draw`
is an alternative to explicit replacements, not an additional source of data.
This avoids storing another full reference dataset in the structural object or
inferring a historical raw dataset from a call. Native refitted measurement
models retain their explicit fitting data in a private formula environment.

Initial scope is one referenced measurement model: Gaussian lmer (including
contained nuisance groups) or supported single-group linear lme. Multiple-model
joint sampling and supplied cross-model covariance replay remain deferred.
Regularization policy/custom alpha, structural formula, and actual computation
backend come from the fitted result, not from evaluation of the original call.

## Identity, settings, and data contract

Reference validation compares effective rows and grouping identities; fixed/random
designs; responses, weights and offsets; covariance; measurement-cohort scores;
structural inputs and exclusions; coefficients and their covariance. Numeric
agreement uses a documented positive tolerance. Cohort/effect identity and
specification/schema checks are exact. Reference fitting precedes replacements.

Native adapters interpret a small grammar of explicit columns, deterministic
formula operators, literal control/covariance constructors, and atomic settings.
They never evaluate original data/subset expressions or arbitrary global bindings.
ML/REML and native fixed-factor matrices come from fitted metadata. Stateful or
sample-dependent transformations, random-only/transformed categorical coding,
and symbolic fitting objects require an explicit deterministic callback.

Recorded control recipes are normalized and compared exactly. lmer additionally
retains optimizer names/control settings in its fitted object. Some nlme control
choices exist only in the original call: when those refer to unavailable symbolic
objects, a callback must declare the original recipe. The audit labels this
`callback_declared` rather than implying that unrecoverable settings were checked.
Identity agreement alone cannot prove which historical global defaults were used.
Callbacks own unrecorded choices and must not conceal convergence failures.

All replacement measurement rows must be consumed in their supplied order.
Replacement structural assembly is strict (`na_action = "fail"`), even when the
reference deliberately omitted people. Original exclusions are retained in the
reference audit. Every fixed factor level and coefficient meaning must remain
estimable; no silent term dropping, level renaming, or within-replicate rescue
omission is introduced.

## Deterministic person draws

`draw` makes no RNG calls or inferential summaries. It copies a source person's
complete effective records plus its structural row for every occurrence, with a
fresh target ID and an explicit source-to-copy mapping. Within-person observation
order is retained. Grouping expressions containing the target ID (person:day)
become distinct automatically; standalone eligible nuisance ID columns are also
re-keyed. Nuisance IDs used as covariates, computed target IDs, and unsupported
composite nuisance expressions require explicit replacement data.

The measurement and analyzed source cohorts must coincide for this convenience
interface. It cannot silently discard measurement-only or originally omitted
people. A deliberately restricted reference fit is required first. The ordinary
shared-nuisance validator remains active. This is a deterministic building block,
not evidence of valid bootstrap uncertainty; stochastic replicate management,
percentile intervals, and scientific coverage validation remain future work.

## Results, failures, and compatibility

Results use the existing REPAIR class and generics/table schemas. A public
`summary(result)$refit` audit records mode, adapter, stage counts, control-check
status, warnings/messages, optional draw mapping, excluded reference observations,
and software versions. Native recipe internals stay outside this public summary.
Default ordinary-fit diagnostics and table columns are preserved.

Measurement errors identify their stage. Available optimizer failure codes and
known nonconvergence warnings stop; singularity alone does not. Informational
warnings/messages are captured for inspection. Structural numerical failure and
changed schemas stop, making future replicate engines able to record invalid
attempts rather than silently replacing them. Current SEs remain the conditional
model-based SEs of the refitted analysis. No additional variance is added.

## Verification workflow

1. Inspect fitted-model metadata and backend APIs before defining replay scope.
2. Test native identity replay and independent manual replacement refits.
3. Exercise nuisance duplication, fresh ID maps, fixed coding, ML/REML, weights,
   offsets, covariance structures, controls, exclusions, and negative cases.
4. Verify native replay with original call bindings unavailable in a new process,
   and round-trip completed refits with their public audit and result tables.
5. Keep all existing independent numerical fixtures/tolerances unchanged.
6. Regenerate help, install and execute README, run source-package checks with all
   examples/vignettes, then build/link-check the site sequentially.
7. Require the strict six-environment package CI and documentation build on the
   PR revision; record exact counts, run links, and local toolchain limitations.

The existing frozen references continue to validate REPAIR calculations. Replay
comparisons verify the refitting procedure; they do not establish interval
coverage or validate future bootstrap inference.
