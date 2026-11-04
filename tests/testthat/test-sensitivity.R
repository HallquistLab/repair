test_that("sensitivity results equal direct fits, including adaptive alphas and fallback", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example()
  scenarios <- list(independent = "independent", supplied = ex$case$cross_model_error)
  for (backend in c("r", "rcpp")) {
    result <- repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", scenarios,
      backend = backend, regularization = "step_up", conf.level = .9,
      step_up_control = list(min_relative_eigenvalue = .99999))
    expect_identical(result$scenarios$scenario, names(scenarios))
    expect_identical(result$assumptions, scenarios)
    expect_true(all(result$scenarios$status == "ok"))
    expect_identical(result$conf.level, .9)
    for (j in seq_along(scenarios)) {
      direct <- repair(ex$formula, ex$models, ex$case$people, "person",
        cross_model_error = scenarios[[j]], backend = backend, regularization = "step_up",
        step_up_control = list(min_relative_eigenvalue = .99999))
      expect_identical(result$fits[[j]]@struc_fit, direct@struc_fit)
      table <- result$coefficients[result$coefficients$scenario == names(scenarios)[j], -1]
      rownames(table) <- NULL
      expect_identical(table, repair_coefficients(direct, .9))
      d <- direct@struc_fit@diagnostics
      expect_identical(result$scenarios$alpha_preliminary[j], d$alpha_prelim)
      expect_identical(result$scenarios$alpha_final[j], d$alpha_final)
      expect_identical(result$scenarios$fallback_preliminary[j], d$step_up$preliminary$fallback)
      expect_identical(result$scenarios$fallback_final[j], d$step_up$final$fallback)
    }
  }
  fallback <- repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", scenarios,
    regularization = "step_up", step_up_control = list(min_eigenvalue = 1e10))
  expect_true(all(fallback$scenarios$fallback_preliminary))
  expect_true(all(fallback$scenarios$fallback_final))
  expect_true(all(fallback$scenarios$status == "ok"))
})

test_that("invalid scenarios are retained beside successful fits and have no coefficient rows", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); bad <- ex$case$cross_model_error
  bad[[1]]$covariance[] <- 1e6
  scenarios <- list(independent = "independent", missing_policy = NULL,
    impossible = bad, incomplete = list())
  result <- repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", scenarios)
  expect_identical(result$scenarios$status, c("ok", "error", "error", "error"))
  expect_match(result$scenarios$message[2], "explicit")
  expect_match(result$scenarios$message[3], "positive semidefinite")
  expect_identical(names(result$fits), names(scenarios))
  expect_null(result$fits$impossible)
  expect_identical(unique(result$coefficients$scenario), "independent")
  expect_true(all(is.na(result$scenarios$nobs[-1])))
  empty <- repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", scenarios[-1])
  expect_identical(nrow(empty$coefficients), 0L)
  expect_identical(names(empty$coefficients), names(result$coefficients))
})

test_that("numerical failures and warnings are recorded rather than discarded", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example()
  # Exercise a numerical failure through the existing estimator, without making
  # input validation or covariance PSD invalid.
  local_mocked_bindings(run_struc = function(...) class_mmom_fit(
    converged = FALSE, message = "test numerical failure", diagnostics = list(step = "step2")))
  result <- repair_sensitivity(ex$formula, ex$models, ex$case$people, "person",
    list(independent = "independent"))
  expect_identical(result$scenarios$status, "failed")
  expect_match(result$scenarios$message, "test numerical failure")
  expect_match(result$scenarios$warnings, "Structural fit failed")
  expect_equal(result$scenarios$nobs, 70)
  expect_identical(nrow(result$coefficients), 0L)
  expect_false(result$fits$independent@struc_fit@converged)
})

test_that("shared input and control errors stop outside the scenario handler", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); scenarios <- list(independent = "independent")
  for (bad in list(list(), list("independent"), setNames(list(1, 2), c("x", "x")))) {
    expect_error(repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", bad), "uniquely named")
  }
  expect_error(repair_sensitivity(ex$formula, ex$models, ex$case$people, "wrong", scenarios), "id")
  expect_error(repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", scenarios,
    conf.level = 1), "conf.level")
  expect_error(repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", scenarios,
    regularization = "custom", alpha = 1), "npredictors")
  expect_error(repair_sensitivity(ex$formula, ex$models, ex$case$people, "person", scenarios,
    regularization = "step_up", step_up_control = list(unknown = 1)), "step_up_control")
})

test_that("specifications fit once per sensitivity call and preserve fitted inputs", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); count <- 0L
  native <- repair_fit_mlm_spec
  local_mocked_bindings(repair_fit_mlm_spec = function(...) {
    count <<- count + 1L; native(...)
  })
  models <- list(growth = repair_mlm(y ~ z, ex$case$a, engine = "lme", random = ~ z | id),
    response = ex$models$response)
  result <- repair_sensitivity(ex$formula, models, ex$case$people, "person",
    list(zero1 = "independent", zero2 = "independent"))
  expect_identical(count, 1L)
  expect_identical(result$fits$zero1@mlm_fit$response, ex$models$response)
  expect_identical(result$fits$zero1@struc_fit, result$fits$zero2@struc_fit)
})
