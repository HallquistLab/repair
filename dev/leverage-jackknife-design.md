# Full-system leverage and delete-one inference: design and audit

## Scope and acceptance workflow

Build on PR #21 without changing its HC0 default or ordinary conditional SEs.
Use one isolated branch, freeze numerical code before the main comparison, and
validate on teaching data / excluded replication 801. Use fresh independent
random-number streams for the six-profile noisy-predictor comparison. Preserve
all previous studies and report missing intervals as missing, with actual
denominators and matched comparisons. The study protocol is retained verbatim
in the review bundle; no correction is selected by main-grid coverage.

The analytic method supports the same single-Gaussian-lmer context as
repair_sandwich(). The exact jackknife supports validated single-model person
draws through repair_refit(), including contained nuisance effects and supported
nlme structures. Neither implements shared nuisance effects or cross-model
derivatives/replay. Original point estimates remain unchanged.

## Explicit augmented system

For person i, form z_i=(y_i,1,x_i')', M_i=z_i z_i', and E_i the augmented score-
error covariance. Structurally zero error entries, including the error-free
intercept and observed variables, are not free nuisance coordinates. Retain
only the noisy submatrix of E; all symmetric observed moments are retained.

The parameter vector contains native transformed measurement parameters psi,
unweighted mean moments M/E, weighted mean moments W/F, log mean weight log(a),
and final structural coefficients b. Preliminary b0 is the canonical corrected
solve using M/E. With c=(1,-b0')', residual variance is the canonical forced-zero
or max(0, N/(N-p)*c'Mc-c'Ec) rule. Individual inverse case variances w_i use
that residual variance, b0 and each reconstructed error covariance. Both
preliminary/final multipliers retain the user's policy/alpha and generalized
eigenvalue rule, including their active branches.

The per-person structural equations are:

    vech(M_i)-vech(M)
    vech(E_i)-vech(E)                # only allowed error coordinates
    w_i*vech(M_i)-vech(W)
    w_i*vech(E_i)-vech(F)
    w_i/a - 1
    (w_i/a) * [corrected_Xy_i - corrected_XX_i*b]

The final correction multiplier in corrected_Xy/XX is computed from W/F.
The log mean-weight equation cancels irrelevant common estimated weights in the
OLS special case. Without normalization a common nuisance scale can create a
spurious HC3 adjustment. This was caught by the pre-main OLS acceptance test;
the final specified representation reproduces ordinary OLS HC3.

The measurement block uses the exact marginal Gaussian score, including the
global REML determinant. Unlike ML, it is not falsely represented as independent
additive per-person likelihood scores. Its case-frequency derivatives provide
the HC0 meat; its exact deleted score/determinant provides deletion sensitivity.

## Correction and decomposition

Write eta=(psi,moments,log(a),b) and phi_i for the structural equations. Define
the deleted system at the original N (normalization by N is immaterial to roots):

    G_-i(eta) = ( S_-i(psi)/N, sum_{j!=i} phi_j(eta)/N ).

Subtract the small original native score residual as a constant anchor from
the measurement block. At eta_hat, the residual is

    ((S_-i(psi_hat)-S_full(psi_hat))/N, -phi_i(eta_hat)/N).

The deleted Jacobian J_-i uses the exact REML H-H_i and D-D_i determinant terms,
and derivatives of all retained structural equations. Solve
delta_i=-J_-i^{-1} G_-i(eta_hat); let d_i select its b coordinates. The candidate
covariance is sum_i d_i d_i'. These are **uncentered** HC3-style outer products,
not the centered exact jackknife. In OLS, d_i is the exact deletion coefficient
shift, giving ordinary HC3. Original N and policy are fixed in the analytic
structural equations; exact deletion below updates N to N-1.

Split the deleted residual into structural-only and measurement-only parts.
Their linear solutions give direct and propagated shifts; covariance decomposes
into direct, measurement and both signed cross-products. No independence
assumption or positive-only variance addition is imposed. The result retains
original coefficients, person shifts, conditioning and numerical diagnostics.

This is a newly specified HC3-style extension, not a proof of finite-sample
unbiasedness or a CR2/Satterthwaite implementation. Leverage corrections in
nonlinear systems are representation dependent. The chosen representation and
normalization are part of the method. Mancl–DeRouen (2001) and MacKinnon–White
(1985) motivate deletion/leverage corrections; their guarantees do not transfer
automatically to regularized M-MOM.

## Numerical checks and failure gates

Use validated central/half-step/one-sided derivatives on unchanged active
branches. Collect all deleted measurement scores in one evaluation to avoid
rebuilding every marginal covariance independently for each deletion. Scale
both parameter columns and equation rows before conditioning/solving; do not
ridge, truncate or pseudoinvert weak systems. The augmented HC0 influence must
agree with the canonical full influence (relative tolerance .001); much tighter
teaching/OLS checks test the derivation. Every deleted system must be usable,
and the original full influence must be available, before corrected inference
is returned. Unavailable corrected covariance fields are NA, never plain HC0.

Tests cover exact native ML/REML deletion scores; canonical moment reconstruction
under minimum/none/custom policies; augmented HC0 agreement; OLS HC3; signed
decomposition; saved extraction; ID/RNG/original-result preservation; conditioning
and explicit unavailability. Scope guards from PR #21 remain in force.

## Actual whole-person deletion

repair_jackknife() validates reference replay, then calls repair_refit() with each
N-1 person draw. All observations and contained nuisance IDs follow the person.
The measurement model is optimized and the complete M-MOM estimator recalculated,
including N-1 finite-sample multipliers. Store all N complete vectors or failed
rows, optimizer/warning/message histories, singularity and corrections. Covariance
is (N-1)/N times centered deletion cross-products. Bias=(N-1)*(mean deletes-original)
is diagnostic only; original estimates remain the interval centers. Any failed
vector invalidates covariance/interval extraction, with no successes-only option.
Standard RNG state is restored even when a callback uses RNG or fails.

Normal Wald intervals are primary. N-1 t intervals are explicitly labeled
secondary heuristics, not derived degrees of freedom. Neither symmetric method
is presented as a universal correction for skewness, bias or active boundaries.
