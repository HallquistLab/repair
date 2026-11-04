# Exercise the same public refitting/bootstrap code with the frozen structural
# estimator substituted only for the legacy run. Measurement models, draw maps,
# convergence policy, and failure accounting follow the ordinary public APIs.
step_up_bootstrap_pair <- function(fit, repeated, people, id, B, seed) {
  optimized <- mmom_step_up
  run <- function(legacy) {
    implementation <- if (legacy) step_up_legacy$mmom_step_up else optimized
    structural <- list()
    local_mocked_bindings(mmom_step_up = function(inputs) {
      value <- implementation(inputs)
      # Capture every structural fit, including reference replays and every
      # replacement, to compare SEs/weights/moments not stored in boot vectors.
      structural[[length(structural) + 1L]] <<- list(converged = value@converged,
        status_code = value@status_code, message = value@message,
        coefficients = value@coefficients, vcov = value@vcov, steps = value@steps,
        alpha_prelim = value@diagnostics$alpha_prelim, alpha_final = value@diagnostics$alpha_final)
      value
    }, .package = "repair")
    warnings <- character()
    value <- withCallingHandlers(repair_bootstrap(fit, repeated, people, id = id, B = B, seed = seed),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(value = value, warnings = warnings, structural = structural)
  }
  before <- .Random.seed; kind <- RNGkind()
  actual <- run(FALSE)
  expect_identical(.Random.seed, before); expect_identical(RNGkind(), kind)
  expected <- run(TRUE)
  expect_identical(.Random.seed, before); expect_identical(RNGkind(), kind)
  expect_identical(actual$warnings, expected$warnings)
  expect_identical(actual$structural, expected$structural)
  new <- actual$value; old <- expected$value
  expect_identical(names(new), names(old))
  for (field in setdiff(names(new), c("diagnostics", "call")))
    expect_identical(new[[field]], old[[field]], info = field)
  fields <- setdiff(names(new$diagnostics), "elapsed")
  expect_identical(new$diagnostics[fields], old$diagnostics[fields])
  if (sum(new$diagnostics$usable) >= 2L) {
    expect_identical(vcov(new, allow_failed = TRUE), vcov(old, allow_failed = TRUE))
    expect_identical(confint(new, allow_failed = TRUE), confint(old, allow_failed = TRUE))
    expect_identical(repair_bootstrap_coefficients(new, allow_failed = TRUE),
      repair_bootstrap_coefficients(old, allow_failed = TRUE))
  }
  if (any(!new$diagnostics$usable)) {
    expect_error(vcov(new), "withheld"); expect_error(vcov(old), "withheld")
  }
  new
}

test_that("adaptive robust lmer bootstrap retains all draws, alphas, estimates and intervals", {
  set.seed(901)
  example <- repair_example
  model <- repair_lmer_refitter(score ~ day + (day | person))(example$repeated)
  fit <- repair(practice ~ re("day"), model, example$people, "person",
    regularization = "step_up", step_up_control = list(min_relative_eigenvalue = .999))
  boot <- step_up_bootstrap_pair(fit, example$repeated, example$people, "person", 5, 41)
  expect_true(all(boot$diagnostics$usable))
  expect_true(any(abs(boot$diagnostics$preliminary_alpha - fit@struc_fit@diagnostics$alpha_prelim) > .01))
})

test_that("nested and correlated-residual bootstrap outputs are unchanged", {
  set.seed(902)
  example <- nested_example("nested_slope")
  fit <- repair(example$formula, list(growth = example$model), example$case$people, "person",
    regularization = "step_up", step_up_control = list(min_relative_eigenvalue = .1))
  nested <- step_up_bootstrap_pair(fit, example$case$data, example$case$people, "person", 2, 304)
  expect_true(all(nested$diagnostics$usable))
  case <- readRDS(system.file("extdata", "nlme-workflow-reference.rds", package = "repair"))$cases$both_ar1_ident
  model <- nlme::lme(y ~ z, random = ~ z | id, data = case$data,
    correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varIdent(form = ~ 1 | phase),
    control = nlme::lmeControl(msMaxIter = 300))
  fit <- repair(as.formula(case$formula), list(growth = model), case$people, "id",
    regularization = "step_up", step_up_control = list(min_relative_eigenvalue = .1))
  correlated <- step_up_bootstrap_pair(fit, case$data, case$people, "id", 2, 143)
  expect_true(all(correlated$diagnostics$usable))
})

test_that("natural failed draws retain the strict bootstrap gate and success-conditional inference", {
  set.seed(903)
  example <- repair_example
  people <- example$people
  people$category <- factor(c("rare", rep("common", nrow(people) - 1L)))
  model <- lme4::lmer(score ~ day + (day | person), example$repeated)
  fit <- repair(practice ~ re("day") + category, model, people, "person", regularization = "step_up")
  boot <- step_up_bootstrap_pair(fit, example$repeated, people, "person", 9, 41)
  expect_true(any(!boot$diagnostics$usable))
  expect_gte(sum(boot$diagnostics$usable), 2L)
})
