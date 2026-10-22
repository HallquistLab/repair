# Run from the repair package root:
# Rscript dev/generate_repair_fixtures.R /path/to/blup_unshrinking
# This is a maintenance script, not part of the test suite. The frozen fixture
# keeps package tests independent of the simulation repository and its packages.
# Reference results come from the simulation's PRIMARY stabilized variant, not
# from repair. The independent algebra audit must agree before anything is saved.
args <- commandArgs(trailingOnly = TRUE)
reference_root <- normalizePath(if (length(args)) args[[1L]] else '../blup_unshrinking')
reference <- new.env(parent = globalenv())
paths <- c(
  'R/stage2_estimators.R', 'R/fuller_mlm.R',
  'vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R'
)
for (path in paths) sys.source(file.path(reference_root, path), reference)

make_case <- function(seed, scenario, single = FALSE) {
  set.seed(seed)
  sim <- reference$simulate_audit_replication(80L, 'regular', scenario)
  if (single) {
    keep <- c(1L, 3L)
    sim$x <- sim$x[, keep, drop = FALSE]
    sim$omega_x <- sim$omega_x[, keep, keep, drop = FALSE]
    sim$omega_xy <- sim$omega_xy[, keep, drop = FALSE]
  }
  sim
}
cases <- list(
  predictor_homoskedastic = make_case(202610221L, 'predictor_only', TRUE),
  outcome_only = make_case(202610223L, 'outcome_heteroskedastic'),
  both_sides = make_case(202610224L, 'both_heteroskedastic')
)
# Add heterogeneous predictor error to a replication with an observed outcome.
set.seed(202610222L)
sim <- reference$simulate_audit_replication(80L, 'regular', 'no_measurement_error')
measurement <- reference$measurement_parameters(
  reference$geometry_parameters('regular')$g, 'both_heteroskedastic', 80L
)
for (i in seq_len(80L)) {
  sim$x[i, -1L] <- sim$x[i, -1L] +
    drop(rnorm(2L) %*% chol(measurement$omega_x[i, -1L, -1L]))
}
sim$omega_x <- measurement$omega_x
cases$predictor_heteroskedastic <- sim

# Correlated sampling errors must enter the cross-moment, feasible weights and
# variance meat. Generate a coherent joint error distribution with nonzero
# predictor-predictor and predictor-outcome covariances.
set.seed(202610225L)
n <- 80L
h <- rep(c(.3, 1, 2), length.out = n)
joint <- diag(c(.12, .16, .16)) + .025
errors <- matrix(rnorm(n * 3L), n) %*% chol(joint) * sqrt(h)
latent <- matrix(rnorm(n * 2L), n)
sim <- list(
  y = .3 + drop(latent %*% c(.6, -.25)) + rnorm(n, sd = .7) + errors[, 1L],
  x = cbind(1, latent + errors[, -1L]),
  omega_x = array(0, c(n, 3L, 3L)),
  omega_y = joint[1L, 1L] * h,
  omega_xy = cbind(0, outer(h, joint[1L, -1L]))
)
for (i in seq_len(n)) sim$omega_x[i, -1L, -1L] <- h[i] * joint[-1L, -1L]
cases$both_sides_correlated <- sim

# This low-information case sharply distinguishes the modified covariance bread
# from the old package calculation while retaining finite primary REPAIR SEs.
set.seed(42)
n <- 60L
latent <- rnorm(n)
cases$weak_predictor <- list(
  x = cbind(1, latent + rnorm(n, sd = 3)),
  y = .5 + .8 * latent + rnorm(n, sd = .6),
  omega_x = array(rep(c(0, 0, 0, 9), n), c(2L, 2L, n)),
  omega_y = rep(0, n), omega_xy = matrix(0, n, 2L)
)
cases$weak_predictor$omega_x <- aperm(cases$weak_predictor$omega_x, c(3, 1, 2))

assert_close <- function(a, b) {
  stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-9, check.attributes = FALSE)))
}
fixtures <- lapply(cases, function(sim) {
  q <- ncol(sim$x) - 1L
  colnames(sim$x) <- c('(Intercept)', paste0('x', seq_len(q)))
  inputs <- list(
    outcome = sim$y, predictors = sim$x[, -1L, drop = FALSE],
    err_x = aperm(sim$omega_x[, -1L, -1L, drop = FALSE], c(2, 3, 1)),
    err_y = sim$omega_y, err_xy = sim$omega_xy[, -1L, drop = FALSE]
  )
  full <- reference$fit_fuller_mlm(
    inputs$outcome, inputs$predictors, inputs$err_x,
    inputs$err_y, inputs$err_xy, backend = 'r'
  )
  stopifnot(full$converged)
  dat <- data.frame(
    y = sim$y, x1 = inputs$predictors[, q], s22 = inputs$err_x[q, q, ],
    syy = inputs$err_y, c1 = inputs$err_xy[, q]
  )
  production_args <- list(
    stage2_df = dat, outcome = 'y', predictor_u1 = 'x1', meas22 = 's22',
    outcome_meas_var = 'syy', predictor_outcome_meas_cov_u1 = 'c1',
    skip_internal_scaling = TRUE
  )
  if (q == 2L) {
    production_args$stage2_df$x0 <- inputs$predictors[, 1L]
    production_args$stage2_df$s11 <- inputs$err_x[1L, 1L, ]
    production_args$stage2_df$s12 <- inputs$err_x[1L, 2L, ]
    production_args$stage2_df$c0 <- inputs$err_xy[, 1L]
    production_args <- c(production_args, list(
      predictor_u0 = 'x0', meas11 = 's11', meas12 = 's12',
      predictor_outcome_meas_cov_u0 = 'c0'
    ))
  }
  production <- do.call(reference$fit_fuller_dual_core, production_args)
  audit <- reference$fit_fuller_audit(
    sim$y, sim$x, sim$omega_x, sim$omega_y, sim$omega_xy,
    preliminary = 'current_cross', se_bread = 'current_final'
  )
  stopifnot(production$status_code == 0L, audit$status == 'ok')
  focal <- unlist(full$coefficients[q + 1L, c('estimate', 'se', 'ci_low', 'ci_high')])
  assert_close(focal, unlist(production[1L, names(focal)]))
  assert_close(focal, unlist(as.data.frame(audit)[1L, names(focal)]))
  assert_close(full$diagnostics$sigma2, production$fuller_sigma2)
  assert_close(full$diagnostics$lambda1, audit$lambda1)
  assert_close(full$diagnostics$lambda2, audit$lambda2)
  list(
    inputs = inputs, coefficients = full$coefficients, vcov = full$vcov,
    preliminary_coef = full$preliminary_coefficients,
    weights = full$analysis_weights,
    diagnostics = full$diagnostics[c('lambda1', 'lambda2', 'correction1', 'correction3', 'sigma2')],
    simulation_focal = as.data.frame(production[, c('estimate', 'se', 'ci_low', 'ci_high')]),
    independent_audit_focal = as.data.frame(audit)[, c('estimate', 'se', 'ci_low', 'ci_high')]
  )
})
fixture <- list(
  provenance = list(
    simulation_commit = system2('git', c('-C', shQuote(reference_root), 'rev-parse', 'HEAD'), stdout = TRUE),
    source_md5 = setNames(unname(tools::md5sum(file.path(reference_root, paths))), paths),
    preliminary_moment = 'modified', variance_bread = 'modified',
    skip_internal_scaling = TRUE, R_version = R.version.string,
    note = 'Exact per-dataset checks of primary REPAIR, not Monte Carlo coverage validation.'
  ),
  cases = fixtures
)
saveRDS(fixture, 'tests/testthat/fixtures/repair-reference.rds', version = 2)
cat('Saved', length(fixtures), 'independently cross-checked REPAIR fixtures.\n')
