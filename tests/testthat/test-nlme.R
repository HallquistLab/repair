test_that("nlme effective rows and fitted residual structures reproduce independent GLS", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  for (structure in c("iid", "ar1", "ident", "ar1_ident", "car1_power", "fitted_power")) {
    correlation <- switch(structure,
      ar1 = nlme::corAR1(form = ~ time | id),
      ar1_ident = nlme::corAR1(form = ~ time | id),
      car1_power = nlme::corCAR1(form = ~ z | id), NULL)
    weights <- switch(structure,
      ident = nlme::varIdent(form = ~ 1 | phase),
      ar1_ident = nlme::varIdent(form = ~ 1 | phase),
      car1_power = nlme::varPower(form = ~ scale_var),
      fitted_power = nlme::varPower(form = ~ fitted(.)), NULL)
    model <- nlme::lme(y ~ z + phase, random = ~ z | id, data = dat,
      subset = time > 0, na.action = na.exclude, correlation = correlation,
      weights = weights, contrasts = list(phase = contr.sum(2)),
      control = nlme::lmeControl(msMaxIter = 300))
    expected <- independent_nlme_gls(model, dat, structure)
    scores <- repair_ranef_gls(model)
    expect_true(scores@success, info = structure)
    expect_equal(scores@predictions$id, expected$estimates, tolerance = 1e-8, info = structure)
    expect_equal(scores@error_covariance$id, expected$errors, tolerance = 1e-8, info = structure)
    extracted <- repair_get_mlm_elements(model)
    expect_identical(extracted@diagnostics$effective_row_ids, rownames(model$groups))
    expect_equal(extracted@response, dat$y[match(rownames(model$groups), rownames(dat))])
    for (id in names(expected$residual_blocks)) {
      expect_equal(unname(extracted@residual_covariance$blocks[[id]]),
        unname(expected$residual_blocks[[id]]), tolerance = 1e-9, info = paste(structure, id))
    }
    people <- data.frame(id = rev(rownames(nlme::ranef(model))), x = seq_len(35))
    for (backend in c("r", "rcpp")) {
      fit <- repair(re("z", "growth") ~ re("(Intercept)", "growth") + x,
        list(growth = model), people, "id", backend = backend)
      a <- fit@analysis_data
      rows <- match(people$id, rownames(expected$estimates))
      expect_equal(unname(a@outcome), expected$estimates[rows, "z"], ignore_attr = TRUE)
      expect_equal(unname(a@err_xy[, 1L]), expected$errors["z", "(Intercept)", rows], ignore_attr = TRUE)
      input <- list(outcome = expected$estimates[rows, "z"],
        predictors = cbind(expected$estimates[rows, "(Intercept)"], people$x),
        err_y = expected$errors["z", "z", rows], err_x = array(0, c(2, 2, 35)),
        err_xy = cbind(expected$errors["z", "(Intercept)", rows], 0))
      colnames(input$predictors) <- colnames(a@predictors)
      input$err_x[1, 1, ] <- expected$errors["(Intercept)", "(Intercept)", rows]
      reference <- fit_mmom(input, backend = backend)
      expect_true(fit@struc_fit@converged)
      expect_equal(unname(coef(fit)), reference@coefficients$estimate, tolerance = 1e-8)
      expect_equal(vcov(fit), reference@vcov, tolerance = 1e-8)
      expect_identical(a@diagnostics$measurement_backend, "nlme")
    }
  }
})

test_that("nlme omission, subset and changed caller data cannot change fitted rows", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  omitted <- nlme::lme(y ~ z + phase, random = ~ z | id, data = dat,
    subset = time > 0, na.action = na.omit)
  excluded <- nlme::lme(y ~ z + phase, random = ~ z | id, data = dat,
    subset = time > 0, na.action = na.exclude)
  expect_equal(repair_ranef_gls(omitted)@predictions, repair_ranef_gls(excluded)@predictions)
  before <- repair_ranef_gls(excluded)
  dat$y <- dat$z <- 999
  expect_equal(repair_ranef_gls(excluded)@predictions, before@predictions)
  no_data <- excluded; no_data$data <- NULL
  expect_error(repair_ranef_gls(no_data), "keep.data = TRUE", fixed = TRUE)
  mutated <- excluded; mutated$data$z <- 999
  expect_error(repair_ranef_gls(mutated), "do not reproduce fitted values")
  broken <- excluded; rownames(broken$groups)[1L] <- "missing-row"
  expect_error(repair_ranef_gls(broken), "observation IDs")
  expect_identical(repair_effects(list(growth = excluded))$term, c("(Intercept)", "z"))
  # Generic backend messages do not call lme4 diagnostics on an lme object.
  bad_apvar <- excluded; bad_apvar$apVar <- "Non-positive definite approximate variance-covariance"
  expect_match(repair_measurement_diagnostics(bad_apvar)$measurement_messages, "variance-parameter uncertainty")
})

test_that("unsupported nlme residual structures and grouping fail explicitly", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  model <- nlme::lme(y ~ z + phase, random = ~ z | id, data = dat, na.action = na.omit)
  nested <- model; nested$groups$nested <- nested$groups$id
  expect_error(repair_ranef_gls(nested), "exactly one grouping")
  unsupported <- model; unsupported$modelStruct$varStruct <- nlme::varExp()
  expect_error(repair_ranef_gls(unsupported), "variance structures")
  unsupported <- model; unsupported$modelStruct$corStruct <- nlme::corCompSymm()
  expect_error(repair_ranef_gls(unsupported), "Supported nlme correlations")
  nonlinear <- model; class(nonlinear) <- c("nlme", "lme")
  expect_error(repair_ranef_gls(nonlinear), "Nonlinear")
  expect_error(repair_data(re("z") ~ x, nonlinear, data.frame(id = 1:35, x = 1:35), "id"), "Gaussian")
})

test_that("nlme rank failures are exposed to structural omission", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  dat$z[dat$id == "p01"] <- 0
  model <- nlme::lme(y ~ z + phase, random = ~ z | id, data = dat, na.action = na.omit)
  people <- data.frame(id = rownames(nlme::ranef(model)), x = seq_len(35))
  expect_warning(expect_error(repair_data(re("z") ~ x, model, people, "id"), "scoring failed"), "predictions failed")
  expect_warning(expect_warning(fit <- repair(re("z") ~ x, model, people, "id", na_action = "omit"), "Omitting 1"), "predictions failed")
  expect_identical(fit@analysis_data@diagnostics$excluded$id, "p01")
  expect_equal(nobs(fit), 34L)
})

test_that("single random-intercept scoring reduces to residualized per-person means", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  model <- nlme::lme(y ~ z + phase, random = ~ 1 | id, data = dat, na.action = na.exclude)
  score <- repair_ranef_gls(model)
  d <- dat[match(rownames(model$groups), rownames(dat)), ]
  residual <- d$y - drop(model.matrix(~ z + phase, d) %*% nlme::fixef(model))
  expected <- tapply(residual, d$id, mean)[score@levels$id]
  expect_equal(unname(score@predictions$id[, 1]), as.numeric(expected))
  expect_equal(unname(score@error_covariance$id[1, 1, ]), as.numeric(model$sigma^2 / table(d$id)[score@levels$id]))
})

test_that("finer correlation groups are not mistaken for independent analysis people", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  model <- nlme::lme(y ~ z + phase, random = ~ z | id, data = dat, na.action = na.omit,
    correlation = nlme::corAR1(form = ~ time | id/phase), control = nlme::lmeControl(msMaxIter = 300))
  expect_error(repair_ranef_gls(model), "correlation groups must match")
})

test_that("random-factor coefficients retain their fitted contrasts", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  set.seed(101)
  effect <- stats::setNames(rnorm(35, sd = 1.5), sprintf("p%02d", 1:35))
  dat$y <- dat$y + effect[as.character(dat$id)] * as.numeric(dat$phase == "1")
  model <- nlme::lme(y ~ z + phase, random = ~ phase | id, data = dat,
    na.action = na.omit, contrasts = list(phase = contr.sum(2)),
    control = nlme::lmeControl(msMaxIter = 300))
  d <- dat[match(rownames(model$groups), rownames(dat)), ]
  x <- model.matrix(~ z + phase, d, contrasts.arg = list(phase = contr.sum(2)))
  z <- model.matrix(~ phase, d, contrasts.arg = list(phase = contr.sum(2)))
  residual <- d$y - drop(x %*% nlme::fixef(model))
  ids <- rownames(nlme::ranef(model))
  values <- matrix(0, length(ids), 2)
  errors <- array(0, c(2, 2, length(ids)))
  for (i in seq_along(ids)) {
    rows <- as.character(d$id) == ids[i]
    values[i, ] <- lm.fit(z[rows, ], residual[rows])$coefficients
    errors[, , i] <- model$sigma^2 * solve(crossprod(z[rows, ]))
  }
  scores <- repair_ranef_gls(model)
  expect_equal(unname(scores@predictions$id), values, tolerance = 1e-9)
  expect_equal(unname(scores@error_covariance$id), errors, tolerance = 1e-9)
  expect_identical(repair_effects(model)$term, c("(Intercept)", "phase1"))
})
