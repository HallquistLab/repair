# Deterministic robust lmer refitting

## Trigger and diagnosis

The paper-inspired bootstrap screen (2026-10-30) used correlated random
intercept/slope predictors, six measurement profiles, 100 outer replications,
paired focal null/alternative outcomes and B=499. Of 348,302 distinct native
measurement resamples, 2,363 (0.678%) were rejected by the existing convergence
policy. Default interval availability was only 1–53% across profiles, because
one failed coefficient vector withholds inference. Singular-but-usable fits were
already retained. The failures were gradient/Hessian/identifiability warnings,
not structural moment failures or primarily exhausted evaluation budgets.

Development used the first 20 failures/profile, in fixed replication/mode/draw
order. Cold BOBYQA resolved 119/120 but sometimes gave worse likelihoods.
Tighter nloptwrap left all 120 unresolved. BOBYQA initialized at the nloptwrap
theta resolved 119/120 and matched/improved the original objective in all cases.
The remaining case had a near-unidentifiable Hessian. These are diagnostic
results on selected failures, not estimated full-grid recovery probabilities.

## Implementation and design workflow

1. Preserve the original package commit, simulation outputs, source hashes and
   failure histories; work on a separate `codex/robust-lmer-refits` branch.
2. Diagnose causes before changing the fitting or eligibility policy. Compare
   deterministic optimizer candidates on a fixed development subset.
3. Add opt-in `repair_lmer_refitter()` rather than changing the meaning of
   existing original fits. Require rebuilding the original measurement and
   REPAIR analysis with the new recipe.
4. Use nloptwrap for initialization, then tighter BOBYQA. Retry cold and from a
   deterministic interior perturbation only if necessary. Select a converged
   BOBYQA result only when its deviance is within 1e-6 of the best finite attempt.
   Keep unresolved fits rejected, including converged candidates that fail that
   objective check. Preserve model structure and ordinary convergence checks.
5. Persist the recipe in each fitted model and apply it automatically in native
   reference/replacement replay. Compare the declared policy as part of the
   configuration contract; callbacks cannot silently change it. Keep each
   attempt, selected index, objective, warnings/messages and singularity in the
   refit and bootstrap audit. An unresolved finite model remains inspectable with
   a warning; reference/replacement inference rejects it.
6. Test real recorded failures against an independent manual optimizer fit,
   deterministic recovery/retry branches, unresolved and inferior-fit rejection,
   weights/offsets/ML/contrasts/nested groups, RNG isolation, identity replay,
   bootstrap attempt records and saved-model replay. Retain existing tests.
7. Freeze the numerical recipe and rerun all original datasets/draws. Compare
   default availability and coverage with explicit denominators and matched
   datasets, retaining successes-only as a labeled diagnostic. Assess B1999 on
   the same preselected first-three/profile subset. No replacement resamples.

The complete study scripts/results are in the blup_unshrinking repository at
`vig_hallquist_2026/analysis/repair_bootstrap_pilot/robust_refit/`. This numerical
rerun reuses the screen's data and is not independent confirmatory evidence.

## Scope and residual limitations

The helper covers supported Gaussian lmer formulas. nlme and existing optimizer
recipes keep their previous behavior. Factors need supported coding and explicit
complete measurement rows; missing observations are not silently removed.
Improved numerical convergence does not establish percentile-interval validity,
resolve weak measurement/structural identification, or repair M-MOM point bias.
The default any-failure gate and explicit successes-only opt-in are unchanged.

The strategy follows the restart/tolerance/cross-optimizer recommendations in
[lme4's convergence guide](https://lme4.github.io/lme4/reference/convergence.html),
while retaining the package's stricter replay and likelihood checks.
