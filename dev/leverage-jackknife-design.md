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

## Implementation verification

All 2,164 assertions passed with no test failures, warnings or skips; examples
and vignette rebuilds passed. The initial local source check recorded the known
external R_ext/Boolean.h Clang warning and a missing namespace qualification for
tail(). The latter is fixed with utils::tail; the six-platform remote check
matrix passed on implementation head 9dac903 (Linux release/devel/oldrel/R4.3,
Windows release and macOS release). Website CI also passed, and local validation
checked 69 HTML pages, links/assets and search targets. Initial logs are labeled
as such in the review bundle; no compiler-warning suppression was introduced.

The numerical study library remains frozen at cbffe31. Qualifying a utility
function and adding documentation do not alter the estimating equations. Fresh
native fits, analytic outputs, source/binary hashes, input hashes and endpoint-
based verification are retained separately from the earlier studies.

## Fresh simulation findings (2026-10-31)

Read [the report and review bundle](leverage-simulation/README.md). The executed
protocol is preserved verbatim, including normalization, original-N analytic
policy, normal primary references and strict failure gates. New master seed
202610319 produces 600 inputs distinct from all 600 old main inputs. Main B499,
first-three/profile B1999 and half-step checks use frozen code cbffe31. No
dataset, failed coefficient vector or extreme interval was omitted or replaced.

There are 1,200 paired null/effect structural analyses and 12,200 method-specific
interval records (122 cells, each 100 planned datasets). Full IF and HC3 each
return 756/1,200; all 444 refusals are singular original measurement fits, with
no extra HC3 conditioning/branch failures. Exact delete-one returns 1,192/1,200
(eight failed outcome-specific deletion vectors, four measurement datasets);
person percentile returns 1,158/1,200. Availability is not calibration.

Available-case nonzero-effect coverage:

| Profile | Full IF (n) | HC3-style (same n) | Delete-one (n) | Person percentile (n) |
| --- | ---: | ---: | ---: | ---: |
| Low, balanced | 82.1% (56) | 92.9% (56) | 81.6% (98) | 94.8% (97) |
| Negative correlation | 85.4% (48) | 95.8% (48) | 91.9% (99) | 94.9% (98) |
| Small, short | 71.1% (45) | 86.7% (45) | 81.8% (99) | 81.5% (92) |
| Low, larger N | 85.9% (64) | 89.1% (64) | 82.0% (100) | 96.9% (96) |
| High reliability | 89.0% (100) | 92.0% (100) | 90.0% (100) | 89.0% (100) |
| Unequal information | 83.1% (65) | 93.8% (65) | 88.0% (100) | 94.8% (96) |

On identical IF/HC3 datasets, four hard profiles gain about 10–16 percentage
points of coverage (unadjusted paired MC intervals exclude zero). The larger-N
gain is only 3.1 points with MC bounds including zero. Common-eligible HC3 minus
bootstrap coverage is −3.6 points balanced (n=55), +2.1 negative (n=48), +4.9
small/short (n=41), −6.5 larger N (n=62), +3.0 high control (n=100), and −1.6
unequal (n=61). Only the larger-N contrast excludes zero in this exploratory
unadjusted comparison. These do not establish general superiority.

The cost is substantial: median matched HC3/bootstrap width ratios are 1.62–1.85
in hard profiles. HC3's 95th width percentiles at truth .40 range 8.7–28.0,
versus roughly 1.0–1.64 for returned bootstrap intervals. Corrected/HC0 variance
ratios range .957–2252 across all available analyses. Mean SEs hide the heavy
tails; effect detection often drops. There is no post hoc leverage cap or
variance/interval trimming. Very large interval inflation is not declared a
successful calibrated method merely because coverage improves.

Exact deletion does not solve symmetric Wald inference: effect coverage remains
81.6–91.9% in the hard profiles despite 98–100 returned intervals/profile. N-1
t references add little. Point bias remains −.178 for small/short and −.156
for negative correlation; neither method subtracts diagnostic jackknife bias.
Bootstrap-SD Wald coverage is lower than percentile coverage in hard profiles,
confirming that interval shape/location matters beyond estimating a variance.

The analytic/actual deletion comparison is informative. At truth .40 the median
focal relative shift RMSE is .012 in high control. Among hard-profile datasets
without actual singularity/floor/final-branch transitions it is .027–.201; among
those with a transition it is .55–1.33. These selected diagnostic strata support
investigation of local approximation limits, not a causal claim about a single
mechanism. Full details retain per-deletion native and moment diagnostics.

All endpoints and paired rates were independently checked. Complete public
checks comprise 24 B99 bootstraps (2,376 attempts) and 12 jackknives (760 deletions),
with maximum coefficient discrepancy 2.93e-14. Signed decomposition error is at
most 6.31e-16; augmented HC0 relative gap is at most 6.16e-6. Half steps change no
availability decisions on the preselected subset and change matched covariance
by at most 4.86e-7. These validate numerical implementation, not coverage theory.

In the B sensitivity subset, all 78 percentile intervals are available at B499;
six are withheld at B1999 under the unchanged all-vector policy. The 72 common
intervals have unchanged coverage decisions (all cover). Prefixes of native
draws, full/frozen coefficients and usability flags agree exactly. Lost
intervals are retained in the sensitivity evidence, not counted as non-rejections.

After observing the wide tails, a separately labeled post hoc numerical audit
checked the 12 largest HC3/HC0 variance ratios at h/2 and h/4. All 24 calculations
remain available; maximum relative covariance change is 1.06e-4 (about .011%).
No main estimate, threshold, step or interval was replaced. This rules out gross
differencing error as an explanation of those extremes without turning the audit
into a prespecified coverage comparison or a guarantee of exact precision.

Decision: keep both new APIs experimental and separate from ordinary inference.
Full-system leverage can expose sensitivity and improve selected-case coverage,
but its extreme widths and boundary gate prevent a general default recommendation.
Delete-one refitting is a useful complete-person influence/branch diagnostic;
these results do not validate its Wald intervals as a universal bootstrap substitute.
