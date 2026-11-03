step_up_reference <- readRDS(test_path("fixtures", "step-up-reference.rds"))

for (case_name in names(step_up_reference$cases)) local({
  name <- case_name
  case <- step_up_reference$cases[[name]]
  test_that(paste("step-up matches the independent research implementation:", name), {
    control <- list(correction = "step_up", step_up = case$settings)
    expected <- case$expected
    for (backend in c("r", "rcpp")) {
      fit <- fit_mmom(case$inputs, control, control, backend = backend)
      expect_true(fit@converged)
      expect_equal(tail(fit@coefficients$estimate, 1), expected$estimate, tolerance = 1e-9)
      expect_equal(tail(fit@coefficients$se, 1), expected$se, tolerance = 1e-9)
      expect_equal(tail(fit@coefficients$ci_low, 1), expected$ci_low, tolerance = 1e-9)
      expect_equal(tail(fit@coefficients$ci_high, 1), expected$ci_high, tolerance = 1e-9)
      expect_equal(fit@diagnostics$alpha_prelim, expected$fuller_alpha_step1_used, tolerance = 1e-10)
      expect_equal(fit@diagnostics$alpha_final, expected$fuller_alpha_step3_used, tolerance = 1e-10)
      expect_equal(fit@diagnostics$lambda_prelim, expected$fuller_lambda1, tolerance = 1e-10)
      expect_equal(fit@diagnostics$lambda_final, expected$fuller_lambda2, tolerance = 1e-10)
      expect_equal(fit@diagnostics$resid_var, expected$fuller_sigma2, tolerance = 1e-10)
      expect_equal(fit@diagnostics$step_up$preliminary$selected$relative_min_eigenvalue,
        expected$fuller_sx1_star_relative_min_eigen, tolerance = 1e-10)
      expect_equal(fit@diagnostics$step_up$final$selected$relative_min_eigenvalue,
        expected$fuller_sx_star_relative_min_eigen, tolerance = 1e-10)
    }
  })
})

test_that("selected constants reproduce canonical moments, weights and covariance", {
  case <- step_up_reference$cases$strict_correlated
  control <- list(correction = "step_up", step_up = case$settings)
  fit <- fit_mmom(case$inputs, control, control)
  expect_gt(fit@diagnostics$alpha_prelim, 4)
  expect_gt(fit@diagnostics$alpha_final, 4)
  fixed <- fit_mmom(case$inputs,
    list(correction = "custom", alpha = fit@diagnostics$alpha_prelim),
    list(correction = "custom", alpha = fit@diagnostics$alpha_final))
  expect_equal(fit@coefficients, fixed@coefficients, tolerance = 1e-12)
  expect_equal(fit@vcov, fixed@vcov, tolerance = 1e-12)
  expect_equal(fit@analysis_weights, fixed@analysis_weights, tolerance = 1e-12)
  expect_equal(fit@preliminary_coef, fixed@preliminary_coef, tolerance = 1e-12)
  expect_equal(fit@corrected_design_moment, fixed@corrected_design_moment, tolerance = 1e-12)
  for (stage in names(fit@diagnostics$step_up)) {
    search <- fit@diagnostics$step_up[[stage]]
    expect_true(search$met_criteria)
    expect_false(search$fallback)
    expect_equal(search$alpha, min(search$candidates$alpha[search$candidates$guard_pass]))
    expect_lte(search$alpha, search$upper)
    expect_equal(search$settings$search_tolerance, .005)
    # Independently reconstruct the guarded predictor block; the intercept is
    # excluded and preliminary moments alone are divided by N.
    x <- cbind(1, case$inputs$predictors)
    w <- if (stage == "preliminary") rep(1, nrow(x)) else fit@analysis_weights
    e <- matrix(0, ncol(x), ncol(x))
    for (i in seq_len(nrow(x))) e[-1, -1] <- e[-1, -1] + w[i] * case$inputs$err_x[, , i]
    multiplier <- if (stage == "preliminary") fit@diagnostics$correction_prelim else fit@diagnostics$correction_final
    scale <- if (stage == "preliminary") nrow(x) else 1
    moment <- (crossprod(x * sqrt(w)) - multiplier * e)[-1, -1, drop = FALSE] / scale
    ev <- eigen(moment, symmetric = TRUE)$values
    observed <- eigen(crossprod(x * sqrt(w))[-1, -1, drop = FALSE] / scale, symmetric = TRUE)$values
    expect_equal(search$selected$min_eigenvalue, min(ev), tolerance = 1e-12)
    expect_equal(search$selected$relative_min_eigenvalue, min(ev) / max(observed), tolerance = 1e-12)
    expect_equal(search$selected$condition_number, max(ev) / min(ev), tolerance = 1e-12)
  }
})

test_that("no passing alpha is an explicit minimum-alpha fallback", {
  case <- step_up_reference$cases$impossible_guard
  control <- list(correction = "step_up", step_up = case$settings)
  fit <- fit_mmom(case$inputs, control, control)
  baseline <- fit_mmom(case$inputs)
  expect_equal(fit@coefficients, baseline@coefficients)
  for (search in fit@diagnostics$step_up) {
    expect_true(search$fallback)
    expect_false(search$met_criteria)
    expect_equal(search$alpha, 4)
    expect_equal(search$reason, "no_passing_candidate")
    expect_false(any(search$candidates$guard_pass))
  }
})

test_that("step-up handles zero errors, restricted grids, mixed policies and failures", {
  input <- step_up_reference$cases$strict_correlated$inputs
  grid <- list(correction = "step_up", step_up = list(candidate_alphas = c(80, 4, 4), max_refinements = 0L,
    min_relative_eigenvalue = .65))
  fit <- fit_mmom(input, grid, list(correction = "none"))
  expect_identical(fit@diagnostics$final_moment_control$correction, "none")
  expect_identical(names(fit@diagnostics$step_up), "preliminary")
  expect_true(all(fit@diagnostics$step_up$preliminary$candidates$alpha %in% c(4, 80)))
  input$err_x[] <- 0; input$err_y[] <- 0; input$err_xy[] <- 0
  control <- list(correction = "step_up")
  fit <- fit_mmom(input, control, control)
  ols <- lm(input$outcome ~ input$predictors)
  expect_equal(unname(fit@coefficients$estimate), unname(coef(ols)), tolerance = 1e-12)
  expect_equal(unname(fit@vcov), unname(vcov(ols)), tolerance = 1e-12)
  input$predictors[, 2] <- input$predictors[, 1]
  failed <- fit_mmom(input, control, control)
  expect_false(failed@converged)
  expect_true(failed@diagnostics$step_up$preliminary$fallback)
  expect_false(failed@diagnostics$step_up$preliminary$selected$usable)
})

test_that("invalid step-up settings fail before fitting", {
  for (bad in list(list(unknown = 1), list(1), list(coarse_grid_size = 2),
    list(max_refinements = 1.5), list(search_tolerance = 0), list(min_eigenvalue = NA_real_),
    list(max_condition_number = .5), list(candidate_alphas = numeric()), list(candidate_alphas = Inf)))
    expect_error(mmom_step_up_control(bad), "Step-up|step_up_control")
})

test_that("release exports core, step-up controls and bootstrap without analytic extensions", {
  exports <- getNamespaceExports("repair")
  expect_true(all(c("repair", "fit_mmom", "repair_bootstrap", "repair_lmer_refitter") %in% exports))
  expect_false(any(c("repair_sandwich", "repair_sandwich_coefficients", "repair_leverage",
    "repair_jackknife", "repair_jackknife_coefficients") %in% exports))
})

test_that("public refits and bootstrap rerun adaptive selection and preserve saved results", {
  data(repair_example)
  model <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
  settings <- list(min_relative_eigenvalue = .999)
  fit <- repair(practice ~ re("day"), model,
    repair_example$people, "person", regularization = "step_up", step_up_control = settings)
  expect_identical(fit@struc_fit@diagnostics$preliminary_moment_control,
    fit@struc_fit@diagnostics$final_moment_control)
  diagnostics <- repair_diagnostics(fit, table = "regularization")
  expect_equal(diagnostics$alpha, c(fit@struc_fit@diagnostics$alpha_prelim, fit@struc_fit@diagnostics$alpha_final))
  replay <- repair_refit(fit, repair_example$repeated, repair_example$people)
  expect_equal(coef(replay), coef(fit), tolerance = 1e-8)
  expect_output(print(fit), "Step-up preliminary")
  boot <- repair_bootstrap(fit, mlm_data = repair_example$repeated, data = repair_example$people, B = 3, seed = 41)
  expect_true(all(boot$diagnostics$usable))
  expect_true(any(abs(boot$diagnostics$preliminary_alpha - fit@struc_fit@diagnostics$alpha_prelim) > .01))
  for (i in seq_len(3)) {
    draw <- boot$cohort[boot$draws[i, ]]
    manual <- repair_refit(fit, repair_example$repeated, repair_example$people, draw = draw)
    expect_equal(unname(boot$replicates[i, ]), unname(coef(manual)), tolerance = 1e-9)
    expect_identical(manual@struc_fit@diagnostics$preliminary_moment_control$correction, "step_up")
    # Direct fitting on copied data independently verifies replay reselects alpha.
    copied <- repair_refit_draw(fit, repair_example$repeated, repair_example$people, "person", draw)
    m <- lme4::lmer(score ~ day + (day | person), copied$mlm)
    direct <- repair(practice ~ re("day"), m, copied$units, "person",
      regularization = "step_up", step_up_control = settings)
    expect_equal(manual@struc_fit@diagnostics$alpha_prelim, direct@struc_fit@diagnostics$alpha_prelim, tolerance = 1e-8)
    expect_equal(manual@struc_fit@diagnostics$alpha_final, direct@struc_fit@diagnostics$alpha_final, tolerance = 1e-8)
    expect_equal(boot$diagnostics$preliminary_alpha[i], direct@struc_fit@diagnostics$alpha_prelim, tolerance = 1e-8)
    expect_equal(boot$diagnostics$final_alpha[i], direct@struc_fit@diagnostics$alpha_final, tolerance = 1e-8)
  }
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  saved <- readRDS(path)
  expect_equal(summary(saved)$diagnostics$step_up, summary(fit)$diagnostics$step_up)
  expect_error(repair(re("day") ~ practice, model, repair_example$people, "person",
    step_up_control = settings), "only with step-up")
})
