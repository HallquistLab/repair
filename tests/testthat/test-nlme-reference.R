# These expectations were produced without loading repair. Every test starts
# with repeated observations, refits nlme, and invokes the public structural API.
test_that("the public nlme workflow reproduces production prediction extraction and REPAIR inference", {
  skip_if_not_installed("nlme")
  reference <- readRDS(system.file("extdata", "nlme-workflow-reference.rds", package = "repair"))
  for (case in reference$cases) {
    correlated <- case$kind == "both_ar1_ident"
    model <- nlme::lme(y ~ z, random = ~ z | id, data = case$data,
      correlation = if (correlated) nlme::corAR1(form = ~ time | id) else NULL,
      weights = if (correlated) nlme::varIdent(form = ~ 1 | phase) else NULL,
      control = nlme::lmeControl(msMaxIter = 300))
    expect_equal(nlme::fixef(model), case$fitted_beta, tolerance = 1e-5)
    predictions <- repair_ranef_gls(model)
    expected <- case$scores[match(predictions@levels$id, case$scores$id), ]
    expect_equal(unname(predictions@predictions$id[, "z"]), expected$corrected_slope_full, tolerance = 1e-5)
    expect_equal(unname(predictions@predictions$id[, "(Intercept)"]), expected$corrected_intercept_full, tolerance = 1e-5)
    expect_equal(unname(predictions@prediction_error_covariance$id["z", "z", ]), expected$gls_var22, tolerance = 1e-5)
    expect_equal(unname(predictions@prediction_error_covariance$id["z", "(Intercept)", ]), expected$gls_var12, tolerance = 1e-5)
    for (backend in c("r", "rcpp")) {
      fit <- repair(as.formula(case$formula), list(growth = model), case$people, "id", backend = backend)
      expect_true(fit@struc_fit@converged, info = paste(case$kind, backend))
      expect_equal(unname(coef(fit)), case$coefficients$estimate, tolerance = 1e-5)
      expect_equal(unname(vcov(fit)), unname(case$vcov), tolerance = 1e-5)
      expect_equal(unname(confint(fit)), unname(as.matrix(case$coefficients[, c("ci_low", "ci_high")])), tolerance = 1e-5)
      a <- fit@analysis_data
      expect_equal(unname(a@outcome), unname(case$input$outcome), tolerance = 1e-5)
      expect_equal(unname(a@err_x), unname(case$input$err_x), tolerance = 1e-5)
      expect_equal(unname(a@err_y), unname(case$input$err_y), tolerance = 1e-5)
      expect_equal(unname(a@err_xy), unname(case$input$err_xy), tolerance = 1e-5)
      # VH Study 2 applies the generating population slope SD after estimation.
      # The package estimates on the input scale; scale SE and CI consistently.
      focal <- tail(fit@struc_fit@coefficients, 1)[, c("estimate", "se", "ci_low", "ci_high")]
      rownames(focal) <- NULL
      expect_equal(focal * case$reporting_scale, case$reported_focal, tolerance = 1e-5)
    }
  }
})
