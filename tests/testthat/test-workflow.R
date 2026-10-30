test_that("the public workflow reproduces explicit correlated-error M-MOM inputs", {
  ex <- workflow_example()
  f <- re("Days", model = "growth") ~ re("(Intercept)", model = "growth") + x
  input <- repair_data(f, ex$models, ex$data, "person")
  # Independent per-person least squares of fixed-effect-adjusted responses.
  # The intercept/slope covariance must be retained, including across y and x.
  model <- ex$models$growth
  d <- model.frame(model)
  residual <- model.response(d) - drop(lme4::getME(model, "X") %*% lme4::fixef(model))
  expected_y <- numeric(18)
  expected_x <- matrix(0, 18, 2)
  expected_ex <- array(0, c(2, 2, 18))
  expected_ey <- numeric(18)
  expected_exy <- matrix(0, 18, 2)
  for (i in seq_len(18)) {
    rows <- as.character(d$Subject) == ex$data$person[i]
    z <- cbind(1, d$Days[rows])
    coefficients <- lm.fit(z, residual[rows])$coefficients
    covariance <- sigma(model)^2 * solve(crossprod(z))
    expected_y[i] <- coefficients[2]
    expected_x[i, ] <- c(coefficients[1], ex$data$x[i])
    expected_ex[1, 1, i] <- covariance[1, 1]
    expected_ey[i] <- covariance[2, 2]
    expected_exy[i, 1] <- covariance[1, 2]
  }
  expect_equal(unname(input@outcome), expected_y, tolerance = 1e-9)
  expect_equal(unname(input@predictors), expected_x, tolerance = 1e-9)
  expect_equal(unname(input@err_x), expected_ex, tolerance = 1e-9)
  expect_equal(unname(input@err_y), expected_ey, tolerance = 1e-9)
  expect_equal(unname(input@err_xy), expected_exy, tolerance = 1e-9)
  colnames(expected_x) <- colnames(input@predictors)
  for (backend in c("r", "rcpp")) {
    manual <- fit_mmom(list(outcome = expected_y, predictors = expected_x,
      err_x = expected_ex, err_y = expected_ey, err_xy = expected_exy), backend = backend)
    fit <- repair(f, ex$models, ex$data, "person", backend = backend)
    expect_true(fit@struc_fit@converged)
    expect_equal(unname(coef(fit)), manual@coefficients$estimate, tolerance = 1e-9)
    expect_equal(vcov(fit), manual@vcov, tolerance = 1e-9)
    expect_equal(unname(confint(fit)), as.matrix(manual@coefficients[, c("ci_low", "ci_high")]),
      ignore_attr = TRUE, tolerance = 1e-9)
    expect_identical(fit@struc_fit@diagnostics$variance_bread, "modified")
    expect_true(fit@struc_fit@diagnostics$conditional_on_stage1)
  }
})

test_that("IDs define alignment and cohort independently of row order", {
  ex <- workflow_example()
  f <- re("Days") ~ x
  full <- repair(f, ex$models, ex$data, "person", backend = "r")
  reversed <- repair(f, ex$models, ex$data[18:1, ], "person", backend = "r")
  expect_equal(coef(full), coef(reversed), tolerance = 1e-10)
  expect_equal(vcov(full), vcov(reversed), tolerance = 1e-10)
  expect_identical(reversed@analysis_data@ids, ex$data$person[18:1])
  expect_equal(unname(reversed@analysis_data@err_y), rev(unname(full@analysis_data@err_y)))
  subset <- repair_data(f, ex$models, ex$data[3:18, ], "person")
  expect_identical(subset@diagnostics$unused_score_ids, ex$data$person[1:2])
  expect_identical(subset@ids, ex$data$person[3:18])
  ordinary <- lm(unname(full@analysis_data@outcome) ~ ex$data$x)
  expect_equal(unname(coef(full)), unname(coef(ordinary)), tolerance = 1e-9)
  expect_equal(unname(vcov(full)), unname(vcov(ordinary)), tolerance = 1e-9)
})

test_that("observed outcomes and factor predictors carry zero error", {
  ex <- workflow_example()
  a <- repair_data(x ~ re("Days") + category, ex$models, ex$data, "person")
  expect_equal(ncol(a@predictors), 3L)
  expect_true(all(a@err_y == 0))
  expect_true(all(a@err_xy == 0))
  expect_true(all(a@err_x[2:3, , ] == 0))
  expect_equal(colnames(a@predictors)[2:3], c("categoryB", "categoryC"))
  expect_equal(a@predictors[, 2:3], model.matrix(~ category, ex$data)[, -1], ignore_attr = TRUE)
})

test_that("missingness and failed IDs require explicit omission and remain inspectable", {
  ex <- workflow_example()
  dat <- ex$data
  dat$person[1] <- "not-measured"
  dat$x[2] <- NA_real_
  dat$x[3] <- Inf
  expect_error(repair_data(re("Days") ~ x, ex$models, dat, "person"), "3 analysis row")
  expect_warning(a <- repair_data(re("Days") ~ x, ex$models, dat, "person", na_action = "omit"), "Omitting 3")
  expect_identical(a@ids, ex$data$person[4:18])
  expect_identical(a@diagnostics$excluded$id, dat$person[1:3])
  expect_match(a@diagnostics$excluded$reason[1], "absent")
  expect_match(a@diagnostics$excluded$reason[2], "observed variable x")
  expect_error(repair_data(re("Days") ~ x, ex$models, rbind(ex$data, ex$data[1, ]), "person"), "unique")
  dat$person[1] <- NA_character_
  expect_error(repair_data(re("Days") ~ x, ex$models, dat, "person", na_action = "omit"), "nonmissing")
  dat <- ex$data; dat$person[1] <- " "
  expect_error(repair_data(re("Days") ~ x, ex$models, dat, "person"), "nonempty")
  dat <- ex$data; dat$x <- NA_real_
  expect_error(repair_data(re("Days") ~ x, ex$models, dat, "person", "omit"), "No complete")
})

test_that("rank-deficient measurement clusters propagate to exclusion diagnostics", {
  d <- lme4::sleepstudy
  d$Days[d$Subject == "308"] <- 0
  model <- suppressWarnings(lme4::lmer(Reaction ~ Days + (Days | Subject), d))
  dat <- workflow_example()$data
  expect_warning(expect_error(repair_data(re("Days") ~ x, model, dat, "person"), "scoring failed"), "predictions failed")
  expect_warning(expect_warning(a <- repair_data(re("Days") ~ x, model, dat, "person", "omit"), "Omitting 1"), "predictions failed")
  expect_identical(a@diagnostics$excluded$id, "308")
  expect_identical(a@scores@diagnostics$Subject$failed_levels, "308")
})

test_that("effect names and formula grammar fail clearly before numerical fitting", {
  ex <- workflow_example()
  catalog <- repair_effects(ex$models)
  expect_identical(catalog$term, c("(Intercept)", "Days"))
  expect_identical(catalog$model, rep("growth", 2))
  expect_error(re("Days"), "formula marker")
  invalid <- list(~ x, x ~ 1, x ~ missing, x ~ re("typo"),
    x ~ re("Days", model = "typo"), x ~ re(Days), x ~ I(re("Days")^2),
    x ~ re("Days") * category, x ~ 0 + re("Days"), x ~ re("Days") - 1,
    x ~ ., x ~ offset(re("Days")), cbind(x, x) ~ re("Days"),
    re("Days") ~ re("Days"), re("Days") ~ x + person)
  for (f in invalid) expect_error(repair_data(f, ex$models, ex$data, "person"), class = "repair_error")
  models <- list(a = ex$models$growth, b = ex$models$growth)
  expect_error(repair_data(x ~ re("Days"), models, ex$data, "person"), "ambiguous")
  expect_error(repair_data(re("Days", "a") ~ re("(Intercept)", "b"), models, ex$data, "person"), "Cross-model")
  expect_equal(repair_data(x ~ re("1"), ex$models, ex$data, "person")@err_x,
    repair_data(x ~ re("(Intercept)"), ex$models, ex$data, "person")@err_x)
  expect_error(repair_data(x ~ re("Days"), list(ex$models$growth), ex$data, "person"), "uniquely named")
  expect_error(repair_data(x ~ re("Days"), list(a = lm(Reaction ~ Days, lme4::sleepstudy)), ex$data, "person"), "Gaussian")
  expect_error(repair_data(x ~ re("Days") + absent, ex$models, ex$data, "person"), "missing from")
  dat <- ex$data; dat$constant <- 1
  expect_error(repair_data(x ~ re("Days") + constant, ex$models, dat, "person"), "rank-deficient")
  expect_error(repair_data(x ~ re("Days"), ex$models, ex$data[1:7, ], "person"), "at least")
})

test_that("summaries and accessors expose inference and failure without silent NULLs", {
  ex <- workflow_example()
  fit <- repair(re("Days") ~ x, ex$models, ex$data, "person")
  s <- summary(fit)
  expect_s3_class(s, "summary_repair")
  expect_equal(s$n_analyzed, nobs(fit))
  expect_equal(nobs(fit), 18L)
  expect_match(s$inference, "Conditional")
  expect_output(print(fit), "18 analyzed / 18 requested")
  expect_output(print(fit@analysis_data), "18 people")
  expect_equal(confint(fit, parm = 2, level = .90), confint(fit, parm = "x", level = .90))
  expect_equal(unname(confint(fit, "x", .90)), matrix(coef(fit)["x"] + c(-1, 1) * qnorm(.95) * sqrt(vcov(fit)["x", "x"]), 1), tolerance = 1e-10)
  expect_error(confint(fit, level = 1), "level")
  expect_error(confint(fit, parm = "absent"), "Unknown")
  expect_error(confint(fit, parm = 0), "positions")
  expect_error(repair(re("Days") ~ x, ex$models, ex$data, "person", regularization = "custom"), "Supply.*alpha")
  expect_error(repair(re("Days") ~ x, ex$models, ex$data, "person", alpha = 5), "only with custom")
  for (mode in c("none", "custom")) {
    a <- if (mode == "custom") 8 else NULL
    fit2 <- repair(re("Days") ~ x, ex$models, ex$data, "person", regularization = mode, alpha = a)
    expect_equal(summary(fit2)$regularization$final_moment_control$correction, mode)
  }
  failed <- fit
  failed@struc_fit <- mmom_failure("case_weights", "No positive case variance")
  expect_false(summary(failed)$converged)
  expect_output(print(failed), "FAILED")
  expect_error(coef(failed), "did not converge")
  expect_error(vcov(failed), "did not converge")
  expect_error(confint(failed), "did not converge")
  expect_error(mmom_control(list(correction = "step_up"), 1, "control"), "not implemented")
})

test_that("Gaussian and independent-group boundaries are enforced", {
  d <- lme4::sleepstudy
  d$site <- factor(rep(1:3, length.out = nrow(d)))
  multi <- suppressMessages(lme4::lmer(Reaction ~ Days + (Days | Subject) + (1 | site), d))
  expect_error(repair_data(re("Days") ~ x, multi, workflow_example()$data, "person"), "shared across target")
  binomial <- suppressMessages(lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
    family = binomial, data = lme4::cbpp))
  expect_error(repair_ranef_gls(binomial), "GLMMs are not supported")
  ex <- workflow_example()
  deviations <- repair_ranef_gls(ex$models$growth)
  totals <- repair_ranef_gls(ex$models$growth, include_fixed = TRUE)
  expect_equal(totals@predictions$Subject,
    sweep(deviations@predictions$Subject, 2, lme4::fixef(ex$models$growth), "+"))
  expect_equal(totals@error_covariance, deviations@error_covariance)
  unmatched <- lme4::lmer(Reaction ~ 1 + (Days | Subject), d)
  expect_error(repair_ranef_gls(unmatched, include_fixed = TRUE), "matching fixed coefficients")
})

test_that("separate random blocks for one grouping factor preserve score ordering", {
  d <- lme4::sleepstudy
  independent <- lme4::lmer(Reaction ~ Days + (Days || Subject), d)
  score <- repair_ranef_gls(independent)
  expect_identical(score@terms$Subject, c("(Intercept)", "Days"))
  residual <- d$Reaction - drop(lme4::getME(independent, "X") %*% lme4::fixef(independent))
  for (id in score@levels$Subject) {
    rows <- as.character(d$Subject) == id
    z <- cbind(1, d$Days[rows])
    expect_equal(unname(score@predictions$Subject[id, ]), unname(lm.fit(z, residual[rows])$coefficients), tolerance = 1e-9)
    expect_equal(unname(score@error_covariance$Subject[, , id]), sigma(independent)^2 * solve(crossprod(z)), tolerance = 1e-9)
  }
  fit <- repair(re("Days") ~ x, independent, workflow_example()$data, "person")
  ordinary <- lm(score@predictions$Subject[, "Days"] ~ workflow_example()$data$x)
  expect_equal(unname(coef(fit)), unname(coef(ordinary)), tolerance = 1e-9)
})

test_that("marker argument matching follows ordinary R named argument order", {
  ex <- workflow_example()
  expected <- repair_data(x ~ re("Days", "growth"), ex$models, ex$data, "person")
  actual <- repair_data(x ~ re(model = "growth", term = "Days"), ex$models, ex$data, "person")
  expect_equal(actual@predictors, expected@predictors)
  expect_error(repair_data(x ~ re("Days", unknown = "growth"), ex$models, ex$data, "person"), "literal strings")
})

test_that("a numerical M-MOM failure returns a diagnosed fit with a warning", {
  ex <- workflow_example()
  ex$data$x <- 0
  expect_warning(fit <- repair(x ~ re("Days"), ex$models, ex$data, "person"), "Structural fit failed")
  expect_false(fit@struc_fit@converged)
  expect_identical(summary(fit)$diagnostics$step, "case_weights")
  expect_equal(nobs(fit), 18L)
  expect_error(coef(fit), "did not converge")
})

test_that("observed factors retain user contrasts when all levels are present", {
  ex <- workflow_example()
  contrasts(ex$data$category) <- contr.sum(3)
  assembled <- repair_data(x ~ re("Days") + category, ex$models, ex$data, "person")
  expected <- model.matrix(~ category, ex$data)[, -1, drop = FALSE]
  expect_equal(unname(assembled@predictors[, 2:3]), unname(expected))
  expect_true(all(assembled@err_x[2:3, , ] == 0))
})

test_that("summary printing is registered outside the package namespace", {
  ex <- workflow_example()
  fit <- repair(re("Days") ~ x, ex$models, ex$data, "person")
  # Simulate a renderer that cannot find unexported methods by caller lookup.
  external <- list2env(list(x = summary(fit)), parent = baseenv())
  expect_output(evalq(print(x), external), "REPAIR structural regression")
})
