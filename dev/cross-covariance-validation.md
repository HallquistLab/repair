# Explicit cross-model covariance propagation and sensitivity

November 4, 2026. These tools support explicit assumptions without changing
the REPAIR estimator, estimating new covariance parameters, or jointly fitting
measurement models. PRs 21/22 remain excluded. This work is a separate PR;
it is not a CRAN submission.

## Contracts

`repair_cross_covariance()` derives likelihood-only extraction operators using
the scorer's target design and fitted effective error covariance. It propagates
observation cross-covariance with `L_a R_ab L_b'`, then selects exactly the
formula's coefficients. Covariances of latent target effects never enter R_ab.
Effective nuisance covariance does enter within-model R and must be represented
in supplied R_ab when cross-correlated. Fitted fixed effects and variance
parameters remain conditional quantities.

Occasion keys are explicit character vectors named by original observation row
identities. Effective fitted rows are matched by identity after native subset/NA
handling; target IDs come from the fitted grouping variables. Occasion matching
is within person, never positional or nearest-time. Explicit matrix dimensions
must cover the exact effective occasion keys and all requested IDs. Every model
pair is required. Both full joint observation covariance and propagated selected
prediction covariance must be PSD within relative tolerance 1e-8. No projection
to a nearby covariance or silent omission is performed.

The optional scalar-rho shortcut is limited to two models with homoscedastic iid
errors and no nuisance groups. It uses fitted marginal error scales and equal
occasion keys; all unequal/unmatched occasion pairs have zero cross-covariance.
The chosen rho is an assumption. More general residual structures require
explicit compatible matrices. This is not a residual-correlation estimator.

`repair_sensitivity()` normalizes/fits measurement models once, validates common
scope and settings, and calls the existing `repair()` for each named covariance
assumption. It records invalid inputs, numerical structural failures, warnings,
selected alphas, and fallback flags. Successful fits supply the usual conditional
coefficient tables. Full fits and original assumptions are retained. Errors while
constructing blocks occur before this call; the vignette preserves those in a
separate construction table. No scenarios are automatically chosen or averaged.
There is no multiple-model bootstrap or propagation of covariance-parameter
uncertainty added here.

## Validation

`test-cross-covariance.R` compares the iid shortcut and explicit propagation
against the independently generated dual-process fixture (including coefficients,
SE/covariance and intervals for both R and C++ estimators). Zero cross-covariance
produces the identical structural fit as `"independent"`. Tests also cover
reversed model pairs, reordered matrix axes/IDs/input rows, one-coefficient and
one-ID dimensions, unmatched occasions, negative/boundary correlations, mixed
lmer/nlme backends, correlated heterogeneous nlme residuals with omitted rows,
weighted lmer residuals, effective nested nuisance covariance, and singular
extraction failure. Observation-level PSD rejection is tested in directions that
prediction extraction would remove, and with three pairwise-valid but jointly
invalid processes. Identity, coverage and numerical-input errors are explicit.

`test-sensitivity.R` verifies exact fit, coefficient-table and adaptive-alpha
equivalence to direct calls in both computation backends; usable minimum-alpha
fallbacks; errors beside successful scenarios and all-error output schemas;
numerical failure/warning accounting; common validation errors; and fitting a
specification only once while preserving pre-fitted models.

The executable multiple-model guide includes design choices, known-simulation
rho caveats, explicit observation matrices, propagation equivalence checks, a
sensitivity grid retaining construction errors, and conditional-inference limits.
Automatic covariance estimation is deliberately left for separately validated
methodological work.

## Repeated-sampling follow-up

The [basic simulation screen](covariance-screen/README.md) adds 500 replicates
in each of 48 iid two-process conditions, with known and fitted marginal
measurement quantities, three assumed correlations and minimum/default step-up
regularization. It records all 288,000 structural fits, Monte Carlo uncertainty,
prediction-error covariance checks, and public-interface/backend audits. Correct
assumptions give approximately nominal coverage in this screen; wrong but
admissible cross-covariances can severely bias inference. Every step-up selection
remained at alpha 3, so this does not validate activated adaptive selection.
