test_that("Bayesian snapshots preserve IDs and reproduce independent per-unit GLS", {
  p <- bayes_test_parts()
  snapshot <- bayes_test_snapshot(p)
  pred <- repair_ranef_gls(snapshot)
  expect_true(pred@success)
  expect_identical(pred@levels$id, p$levels)
  residual <- p$y - drop(p$x %*% p$beta)
  for (id in p$levels) {
    rows <- which(p$grouping == id)
    z <- p$z[rows, , drop = FALSE]
    expect_equal(unname(pred@predictions$id[id, ]), unname(lm.fit(z, residual[rows])$coefficients), tolerance = 1e-10)
    expect_equal(unname(pred@prediction_error_covariance$id[, , id]),
      unname(p$sigma2 * solve(crossprod(z))), tolerance = 1e-10)
  }
  totals <- repair_ranef_gls(snapshot, include_fixed = TRUE)
  expect_equal(totals@predictions$id, sweep(pred@predictions$id, 2L, p$beta[p$terms], "+"))
  changed <- p; changed$g <- diag(c(2, 3))
  expect_equal(repair_ranef_gls(bayes_test_snapshot(changed))@predictions, pred@predictions)
  expect_equal(repair_ranef_gls(bayes_test_snapshot(changed))@prediction_error_covariance,
    pred@prediction_error_covariance)
})

test_that("snapshots work through structural fitting, covariance, diagnostics and persistence", {
  p <- bayes_test_parts()
  snapshot <- bayes_test_snapshot(p)
  people <- data.frame(id = rev(p$levels), x = seq_along(p$levels))
  catalog <- repair_effects(snapshot)
  expect_identical(catalog$term, p$terms)
  for (engine in c("r", "rcpp")) {
    fit <- repair(re("z") ~ re("(Intercept)") + x, snapshot, people, "id", backend = engine)
    expect_true(fit@struc_fit@converged)
    expect_identical(fit@analysis_data@ids, people$id)
    expect_identical(repair_diagnostics(fit, table = "mlm")$backend, "brms")
    expect_match(repair_diagnostics(fit, table = "mlm")$messages, "Experimental")
    expect_identical(summary(fit)$mlm_models$mlm$stage1_uncertainty, "not_propagated")
    expect_output(print(summary(fit)), "Experimental Bayesian plug-in")
    covariance <- repair_prediction_covariance(fit)
    expect_true(any(covariance$term1 != covariance$term2 & abs(covariance$covariance) > 0))
    expect_true(all(is.finite(repair_reliability(fit)$reliability)))
    expect_error(repair_refit(fit, data.frame(), people), "not yet supported")
    expect_error(repair_bootstrap(fit, data.frame(), people, B = 1), "not yet supported")
    path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
    saveRDS(fit, path)
    expect_equal(repair_coefficients(readRDS(path)), repair_coefficients(fit))
  }
  expect_output(print(snapshot), "Experimental REPAIR Bayesian plug-in snapshot")
})

test_that("Gaussian snapshots match lmer scoring at identical shared parameters", {
  data <- expand.grid(z = seq(-1, 1, length.out = 6), id = factor(1:20))
  set.seed(3)
  data$y <- .4 + .7 * data$z + rep(rnorm(20), each = 6) +
    rep(rnorm(20, sd = .5), each = 6) * data$z + rnorm(nrow(data), sd = .3)
  model <- lme4::lmer(y ~ z + (z | id), data)
  e <- repair_get_mlm_elements(model)
  p <- list(spec = list(group = "id"), x = e@fixed_design, beta = e@fixed_coefs,
    y = e@response, z = model.matrix(~ z, data), grouping = data$id,
    levels = e@group_levels$id, terms = e@random_terms$id,
    g = as.matrix(lme4::VarCorr(model)$id), sigma2 = sigma(model)^2)
  snapshot <- bayes_test_snapshot(p)
  expect_equal(repair_ranef_gls(snapshot)@predictions, repair_ranef_gls(model)@predictions)
  expect_equal(repair_ranef_gls(snapshot)@prediction_error_covariance,
    repair_ranef_gls(model)@prediction_error_covariance)
  people <- data.frame(id = levels(data$id), x = seq_len(20))
  for (f in list(re("z") ~ x, x ~ re("z"), re("z") ~ re("(Intercept)"))) {
    frequentist <- repair(f, model, people, "id")
    plugin <- repair(f, snapshot, people, "id")
    expect_equal(repair_coefficients(plugin), repair_coefficients(frequentist))
    expect_equal(repair_reliability(plugin)$reliability, repair_reliability(frequentist)$reliability)
  }
})

test_that("proper pooling cannot rescue a rank-deficient likelihood design", {
  p <- bayes_test_parts()
  bad <- p$levels[1L]
  p$z[p$grouping == bad, 2L] <- 1
  snapshot <- bayes_test_snapshot(p)
  expect_warning(pred <- repair_ranef_gls(snapshot), "predictions failed for 1 cluster")
  expect_false(pred@success)
  expect_identical(pred@diagnostics$id$failed_levels, bad)
  people <- data.frame(id = p$levels, x = seq_along(p$levels))
  expect_warning(expect_error(repair_data(re("z") ~ x, snapshot, people, "id"), "incomplete"), "predictions failed")
  expect_warning(assembled <- suppressWarnings(repair_data(re("z") ~ x, snapshot, people, "id", na_action = "omit")), NA)
  expect_identical(assembled@diagnostics$excluded$id, bad)
})

test_that("the experimental measurement formula grammar is deliberately narrow", {
  data <- data.frame(y = 1:12, z = 1:12, w = 12:1, id = rep(1:3, each = 4), category = factor(rep(1:2, 6)))
  expect_identical(repair_bayes_formula(y ~ z + w + (1 + z | id), data)$group, "id")
  expect_identical(repair_bayes_formula(y ~ 0 + z + (0 + z | id), data)$group, "id")
  for (f in list(log(y) ~ z + (z | id), y ~ log(z) + (z | id),
    y ~ z * w + (z | id), y ~ z + offset(w) + (z | id),
    y ~ category + (z | id), y ~ z + (z || id), y ~ z + (1 | id) + (0 + z | id),
    y ~ z + (1 | id:w), y ~ id + (z | id), y ~ z + (id | id), y ~ z)) {
    expect_error(repair_bayes_formula(f, data), "requires|requires|support requires")
  }
  expect_error(repair_bayes_formula(y ~ z + (z | id), data["y"]), "retained")
  bad <- bayes_test_parts(); bad$beta <- rev(bad$beta)
  expect_error(bayes_test_snapshot(bad), "invalid")
  bad <- bayes_test_parts(); bad$g <- -diag(2)
  expect_error(bayes_test_snapshot(bad), "positive-semidefinite")
})

test_that("native fits require explicit opt-in and unsupported engines fail", {
  expect_error(normalize_repair_models(structure(list(), class = "brmsfit")), "explicit opt-in")
  expect_error(normalize_repair_models(structure(list(), class = "stanreg")), "explicit opt-in")
  skip_if_not_installed("posterior")
  expect_error(repair_bayes_plugin(lm(mpg ~ wt, mtcars)), "supported sampled")
})

test_that("Bayesian snapshots coexist with native MLM specifications", {
  snapshot <- bayes_test_snapshot()
  p <- bayes_test_parts()
  data <- data.frame(y = p$y, z = p$z[, "z"], id = p$grouping)
  models <- normalize_repair_models(list(bayes = snapshot,
    frequentist = repair_mlm(y ~ z + (z | id), data)))
  expect_identical(models$bayes, snapshot)
  expect_s4_class(models$frequentist, "lmerMod")
  expect_setequal(repair_effects(models)$model, c("bayes", "frequentist"))
  people <- data.frame(id = p$levels, x = seq_along(p$levels))
  mixed <- repair(re("z", "bayes") ~ x, models, people, "id")
  single <- repair(re("z") ~ x, snapshot, people, "id")
  expect_equal(repair_coefficients(mixed), repair_coefficients(single))
})
