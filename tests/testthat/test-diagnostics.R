test_that("structural diagnostics reconstruct moments and identify analysis weights", {
  fit <- reliability_example(); a <- fit@analysis_data; d <- repair_diagnostics(fit)
  expect_equal(d$cohort$n_analyzed, nobs(fit))
  expect_equal(d$measurement$backend, "lme4")
  expect_equal(d$error_variances$variance[d$error_variances$role == "outcome"], unname(a@err_y))
  expect_true(all(d$error_variances$variance[d$error_variances$term == "x"] == 0))
  # Independent weighted predictor moments, full unscaled error subtraction.
  design <- cbind(1, a@predictors)
  for (step in c("preliminary", "final")) {
    w <- if (step == "preliminary") rep(1, nrow(design)) else d$weights$weight
    error <- matrix(0, ncol(design), ncol(design))
    for (i in seq_len(nrow(design))) error[-1, -1] <- error[-1, -1] + w[i] * a@err_x[, , i]
    multiplier <- d$regularization$subtraction_multiplier[d$regularization$step == step]
    moment <- crossprod(design * sqrt(w)) - multiplier * error
    ev <- eigen(moment, symmetric = TRUE)$values
    row <- d$conditioning[d$conditioning$step == step, ]
    expect_equal(row$min_eigenvalue, min(ev), tolerance = 1e-10)
    expect_equal(row$max_eigenvalue, max(ev), tolerance = 1e-10)
    expect_equal(row$condition_number, max(ev) / min(ev), tolerance = 1e-10)
  }
  expect_identical(d$weights$id, a@ids)
  expect_equal(sum(d$weights$fraction), 1)
  expect_equal(d$weight_summary$weight_equivalent_n, sum(d$weights$weight)^2 / sum(d$weights$weight^2))
  expect_equal(d$weight_summary$max_fraction, max(d$weights$fraction))
  expect_equal(d$reliability, repair_reliability(fit))
})

test_that("disabled regularization is distinct from disabled measurement correction", {
  ex <- workflow_example()
  fit <- repair(re("Days") ~ x, ex$models, ex$data, "person", regularization = "none")
  d <- repair_diagnostics(fit)
  expect_equal(d$regularization$subtraction_multiplier, c(1, 1))
  expect_false(any(d$regularization$regularized))
  expect_true(any(d$error_variances$variance > 0))
  expect_equal(d$weight_summary$weight_equivalent_n, 18, tolerance = 1e-10)
})

test_that("failed structural fits still expose measurement and cohort diagnostics", {
  ex <- workflow_example(); ex$data$x <- 0
  expect_warning(fit <- repair(x ~ re("Days"), ex$models, ex$data, "person"), "Structural fit failed")
  d <- repair_diagnostics(fit)
  expect_false(d$status$converged); expect_equal(d$status$step, "case_weights")
  for (table in c("weights", "weight_summary", "conditioning", "regularization")) expect_equal(nrow(d[[table]]), 0L)
  expect_true(all(is.finite(d$reliability$reliability)))
  expect_equal(d$cohort$n_analyzed, 18L)
  expect_output(print(fit), "Structural fit FAILED")
  expect_error(repair_diagnostics(list()), "repair")
})

test_that("unequal analysis weights follow IDs through cohort permutation", {
  skip_if_not_installed("nlme")
  case <- readRDS(system.file("extdata", "nlme-workflow-reference.rds", package = "repair"))$cases$outcome_unbalanced
  model <- nlme::lme(y ~ z, random = ~ z | id, data = case$data)
  fit <- repair(re("z") ~ observed, model, case$people, "id")
  reverse <- repair(re("z") ~ observed, model, case$people[60:1, ], "id")
  a <- repair_diagnostics(fit); b <- repair_diagnostics(reverse)
  expect_equal(b$weights$id, rev(a$weights$id))
  expect_equal(b$weights$fraction, rev(a$weights$fraction), tolerance = 1e-10)
  expect_gt(diff(range(a$weights$fraction)), 1e-4)
  expect_lt(a$weight_summary$weight_equivalent_n, nobs(fit))
  expect_equal(a$weight_summary, b$weight_summary, tolerance = 1e-10)
})
