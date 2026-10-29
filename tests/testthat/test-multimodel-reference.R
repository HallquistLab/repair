test_that("dual-process fits reproduce independent production and algebra references", {
  skip_if_not_installed("nlme")
  for (kind in c("independent", "supplied")) {
    ex <- multimodel_example(kind); case <- ex$case
    expect_equal(nlme::fixef(ex$models$growth), case$fitted_beta$a, tolerance = 1e-5)
    expect_equal(nlme::fixef(ex$models$response), case$fitted_beta$b, tolerance = 1e-5)
    for (backend in c("r", "rcpp")) {
      fit <- repair(ex$formula, ex$models, case$people, "person", backend = backend,
        cross_model_error = case$cross_model_error)
      expect_true(fit@struc_fit@converged)
      expect_equal(unname(coef(fit)), case$coefficients$estimate, tolerance = 1e-5)
      expect_equal(unname(vcov(fit)), unname(case$vcov), tolerance = 1e-5)
      expect_equal(unname(confint(fit)), unname(as.matrix(case$coefficients[, c("ci_low", "ci_high")])), tolerance = 1e-5)
      a <- fit@analysis_data
      for (field in names(case$input)) {
        expect_equal(unname(S7::prop(a, field)), unname(case$input[[field]]), tolerance = 1e-5)
      }
      expect_equal(unname(a@scores$growth@predictions$id[, "z"]), case$scores$a$corrected_slope_full, tolerance = 1e-5)
      expect_equal(unname(a@scores$response@predictions$subject[, "w"]), case$scores$b$corrected_slope_full, tolerance = 1e-5)
      focal <- tail(fit@struc_fit@coefficients, 1)[, c("estimate", "se", "ci_low", "ci_high")]
      rownames(focal) <- NULL
      expect_equal(focal * case$reporting_scale, case$reported_focal, tolerance = 1e-5)
      expect_identical(summary(fit)$cross_model_error$policy, kind)
      expect_setequal(names(fit@mlm_fit), c("growth", "response"))
      expect_output(print(fit), if (kind == "independent") "assumed independent" else "user-supplied")
      expect_equal(summary(fit)$measurement_models$growth$group, "id")
      expect_equal(summary(fit)$measurement_models$response$group, "subject")
    }
  }
})
