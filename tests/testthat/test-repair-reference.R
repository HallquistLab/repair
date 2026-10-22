# Frozen results were produced without repair, checked against both the primary
# simulation estimator and its independent algebra audit. Regenerate only with
# dev/generate_repair_fixtures.R; normal tests need no sibling repository.
reference <- readRDS(test_path('fixtures', 'repair-reference.rds'))

for (case_name in names(reference$cases)) {
  local({
    label <- case_name
    fixture <- reference$cases[[label]]
    test_that(paste('REPAIR reproduces simulation estimates and inference:', label), {
      for (backend in c('r', 'rcpp')) {
        fit <- fit_mmom(fixture$inputs, backend = backend)
        expect_true(fit@converged, info = paste(label, backend, fit@message))
        expect_equal(fit@coefficients, fixture$coefficients, tolerance = 1e-8)
        expect_equal(fit@vcov, fixture$vcov, tolerance = 1e-8)
        expect_equal(fit@preliminary_coef, fixture$preliminary_coef, tolerance = 1e-8)
        expect_equal(fit@analysis_weights, fixture$weights, tolerance = 1e-8)
        expect_equal(fit@diagnostics$lambda_prelim, fixture$diagnostics$lambda1, tolerance = 1e-8)
        expect_equal(fit@diagnostics$lambda_final, fixture$diagnostics$lambda2, tolerance = 1e-8)
        expect_equal(fit@diagnostics$correction_prelim, fixture$diagnostics$correction1, tolerance = 1e-8)
        expect_equal(fit@diagnostics$correction_final, fixture$diagnostics$correction3, tolerance = 1e-8)
        expect_equal(fit@diagnostics$resid_var, fixture$diagnostics$sigma2, tolerance = 1e-8)
        focal <- tail(fit@coefficients, 1L)[, c('estimate', 'se', 'ci_low', 'ci_high')]
        rownames(focal) <- NULL
        expect_equal(focal, fixture$simulation_focal, tolerance = 1e-8)
        expect_equal(focal, fixture$independent_audit_focal, tolerance = 1e-8)
        expect_identical(fit@diagnostics$variance_bread, 'modified')
        expect_true(fit@diagnostics$conditional_on_stage1)
      }
    })
  })
}

test_that('the weak-information fixture detects the old variance-bread error', {
  fit <- fit_mmom(reference$cases$weak_predictor$inputs, backend = 'r')
  # An independently saved simulation result is ~0.695, whereas the former
  # unregularized bread gives ~5.347 on precisely these data.
  expect_equal(tail(fit@coefficients$se, 1L), 0.69516986258826, tolerance = 1e-9)
  expect_lt(tail(fit@coefficients$se, 1L), 1)
})

test_that('the manuscript covariance uses final moments with custom or no correction', {
  input <- reference$cases$both_sides_correlated$inputs
  for (correction in c('minimum', 'custom', 'none')) {
    control <- list(correction = correction, alpha = 8)
    for (backend in c('r', 'rcpp')) {
      fit <- fit_mmom(input, final_moment_control = control, backend = backend)
      expect_true(fit@converged)
      # Direct row-wise evaluation of the manuscript covariance, independent of
      # all package covariance kernels and of its normalization by n.
      x <- cbind(1, input$predictors)
      w <- fit@analysis_weights
      gamma0 <- fit@preliminary_coef[-1L]
      meat <- matrix(0, ncol(x), ncol(x))
      for (j in seq_len(nrow(x))) {
        tilde <- c(0, input$err_xy[j, ] - input$err_x[, , j] %*% gamma0)
        meat <- meat + w[j] * tcrossprod(x[j, ]) + w[j]^2 * tcrossprod(tilde)
      }
      inverse <- solve(fit@corrected_design_moment)
      expected <- inverse %*% meat %*% inverse
      expect_equal(unname(fit@vcov), unname(expected), tolerance = 1e-10)
    }
  }
})
