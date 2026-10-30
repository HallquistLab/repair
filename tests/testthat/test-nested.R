test_that("nested nuisance fits reproduce independent GLS and structural references", {
  for (kind in c("nested_intercept", "nested_slope")) {
    ex <- nested_example(kind); case <- ex$case
    expect_equal(lme4::fixef(ex$model), case$fitted_beta, tolerance = 1e-5)
    for (backend in c("r", "rcpp")) {
      fit <- repair(ex$formula, list(growth = ex$model), case$people, "person", backend = backend)
      a <- fit@analysis_data
      expect_true(fit@struc_fit@converged)
      expect_equal(unname(coef(fit)), case$coefficients$estimate, tolerance = 1e-5)
      expect_equal(unname(vcov(fit)), unname(case$vcov), tolerance = 1e-5)
      expect_equal(repair_coefficients(fit)$std.error, case$coefficients$se, tolerance = 1e-5)
      expect_equal(unname(confint(fit)), unname(as.matrix(case$coefficients[c("ci_low", "ci_high")])), tolerance = 1e-5)
      for (field in names(case$input)) {
        expect_equal(unname(S7::prop(a, field)), unname(case$input[[field]]), tolerance = 1e-5)
      }
      score <- a@predictions
      expect_true(score@success)
      expect_equal(unname(score@predictions$person[case$people$person, ]),
        unname(case$scores), tolerance = 1e-5)
      expect_equal(unname(score@prediction_error_covariance$person[, , case$people$person]),
        unname(case$error), tolerance = 1e-5)
      for (id in case$people$person) {
        expect_equal(unname(score@diagnostics$person$per_level[[id]]$covariance),
          unname(case$covariance[[id]]), tolerance = 1e-5)
      }
      covariance <- repair_prediction_covariance(fit)
      expect_equal(nrow(covariance), 3L * nrow(case$people))
      expect_true(all(covariance$group1 == "person" & covariance$group2 == "person"))
      expect_true(all(covariance$source == "within_model"))
      off <- covariance[covariance$term1 != covariance$term2, ]
      expect_equal(off$covariance, unname(case$error[1, 2, ]), tolerance = 1e-5)
      expect_identical(unique(repair_predictions(fit)$group), "person")
      expect_true(all(grepl('group = "person"', a@effects$label, fixed = TRUE)))
      expect_setequal(summary(fit)$mlm_models$growth$nuisance_groups,
        if (kind == "nested_slope") c("person:day", "person:block") else "person:day")
      expect_output(print(summary(fit)), "Nested nuisance groups:")
    }
  }
})

test_that("posterior reliability includes every nested nuisance covariance contribution", {
  for (kind in c("nested_intercept", "nested_slope")) {
    ex <- nested_example(kind); case <- ex$case
    fit <- repair(ex$formula, list(growth = ex$model), case$people, "person")
    person <- repair_reliability(fit, detail = "unit")
    g <- case$target_g
    for (term in c("(Intercept)", "x")) {
      k <- match(term, colnames(case$scores))
      for (measure in if (term == "x") c("marginal", "residualized") else "marginal") {
        a <- c(0, 0); a[k] <- 1
        if (measure == "residualized") a[1] <- -g[1, 2] / g[1, 1]
        prior <- drop(crossprod(a, g %*% a))
        post <- vapply(seq_len(nrow(case$people)), function(i)
          drop(crossprod(a, case$posterior[, , i] %*% a)), numeric(1))
        selected <- person[person$term == term & person$measure == measure, ]
        expect_identical(selected$id, case$people$person)
        expect_identical(selected$n_observations, case$n_observations)
        expect_true(all(selected$status == "ok"))
        expect_equal(selected$prior_variance, rep(prior, nrow(selected)), tolerance = 1e-5)
        expect_equal(selected$posterior_variance, post, tolerance = 1e-5)
        expect_equal(selected$reliability, 1 - post / prior, tolerance = 1e-5)
      }
    }
    # Deliberately ignoring nuisance variance would overstate information.
    e <- repair_get_mlm_elements(ex$model); id <- case$people$person[1]
    obs <- e@group_rows$person[[id]]
    naive <- repair_posterior_covariance(g, get_ranef_design(e, "person", obs, id),
      get_ranef_base_covariance(e, obs))
    expect_gt(max(abs(naive - case$posterior[, , 1])), 1e-3)
  }
})

test_that("target qualification resolves ambiguity and preserves single-group labels", {
  ex <- nested_example(); models <- list(growth = ex$model); people <- ex$case$people
  expect_error(repair_data(re("x") ~ re("(Intercept)"), models, people, "person"), "unknown or ambiguous")
  a <- repair_data(re("x") ~ re("1", group = "person"), models, people, "person")
  b <- repair_data(re(model = "growth", group = "person", term = "x") ~
    re("(Intercept)", "growth", "person"), models, people, "person")
  expect_identical(a@effects, b@effects)
  expect_equal(a@predictors, b@predictors)
  expect_equal(a@effects$label, c('re("x", model = "growth", group = "person")',
    're("(Intercept)", model = "growth", group = "person")'))
  expect_error(repair_data(re("x", group = "absent") ~ practice, models, people, "person"), "unknown or ambiguous")
  for (group in list(NA_character_, "", 1, c("person", "person:day"), as.name("target"))) {
    f <- as.formula(call("~", call("re", "x", group = group), as.name("practice")))
    expect_error(repair_data(f, models, people, "person"), "literal strings")
  }
  expect_error(repair_data(re("x", group = "person") ~
    re("(Intercept)", group = "person:day"), models, people, "person"), "one target grouping factor")
  single <- workflow_example()
  old <- repair_data(x ~ re("Days", "growth"), single$models, single$data, "person")
  qualified <- repair_data(x ~ re("Days", "growth", "Subject"), single$models, single$data, "person")
  expect_identical(qualified@effects, old@effects)
  expect_equal(qualified@predictors, old@predictors)
})

test_that("nested scope alignment and omissions preserve predictions and joint covariance", {
  ex <- nested_example(); case <- ex$case
  people <- case$people[48:1, ]; people$practice[2] <- NA_real_
  expect_warning(fit <- repair(ex$formula, list(growth = ex$model), people, "person",
    na_action = "omit"), "Omitting 1")
  ids <- people$person[-2]; a <- fit@analysis_data
  expect_identical(a@ids, ids)
  expect_equal(unname(a@outcome), unname(case$scores[ids, "x"]), tolerance = 1e-5)
  expect_equal(unname(a@err_xy[, 1]), unname(case$error[2, 1, ids]), tolerance = 1e-5)
  expect_equal(nrow(repair_predictions(fit, "all")), 96L)
  expect_identical(repair_diagnostics(fit, "excluded")$id, people$person[2])
  expect_equal(repair_diagnostics(fit, "mlm")$n_levels, 48L)
  expect_equal(repair_diagnostics(fit, "scope")$n_analyzed, 47L)
  expect_equal(repair_reliability(fit)$n_units, rep(47L, 3))
})

test_that("containment follows fitted membership rather than group names or unused rows", {
  ex <- nested_example(); d <- ex$case$data
  # Globally unique visit IDs are genuinely nested even without a colon name.
  d$visit <- interaction(d$person, d$day, drop = TRUE)
  # A conflicting membership on an excluded response must not enter the fit.
  extra <- d[which(!is.na(d$y))[1], ]; extra$person <- "P002"
  extra$visit <- levels(d$visit)[1]; extra$y <- NA_real_
  d <- rbind(d, extra)
  model <- lme4::lmer(y ~ x + offset(off) + (x | person) + (1 | visit), d,
    weights = prior, na.action = na.exclude,
    control = lme4::lmerControl(optimizer = "bobyqa"))
  fit <- repair(ex$formula, list(growth = model), ex$case$people, "person")
  expect_equal(unname(coef(fit)), ex$case$coefficients$estimate, tolerance = 1e-5)
  expect_identical(summary(fit)$mlm_models$growth$nuisance_groups, "visit")
})

test_that("shared nuisance levels are rejected across the entire fitted scope", {
  ex <- nested_example(); d <- ex$case$data
  # Day numbers are reused across people; day alone represents a shared effect.
  shared_day <- suppressMessages(lme4::lmer(y ~ x + (x | person) + (1 | day), d))
  expect_error(repair(ex$formula, list(growth = shared_day), ex$case$people, "person"),
    "shared across target group")
  d$school <- factor(ceiling(as.integer(d$person) / 4))
  school <- suppressMessages(lme4::lmer(y ~ x + (x | person) + (1 | school), d))
  # One target person per school does not erase sharing in the MLM fit.
  people <- ex$case$people[seq(1, 48, 4), ]
  for (action in c("fail", "omit")) {
    expect_error(repair_data(ex$formula, list(growth = school), people, "person", na_action = action),
      "shared across target group")
  }
  # A colon in a name does not prove nesting in the requested direction.
  expect_error(repair_data(re("(Intercept)", group = "person:day") ~ practice,
    list(growth = ex$model), ex$case$people, "person"), "shared across target group")
})

test_that("multiple MLMs retain explicit cross-model covariance policies", {
  a <- nested_example(); b <- nested_example("nested_slope")
  models <- list(growth = a$model, response = b$model)
  f <- re("x", "growth", "person") ~ re("x", "response", "person") + practice
  expect_error(repair(f, models, a$case$people, "person"), "cross_model_error")
  fit <- repair(f, models, a$case$people, "person", cross_model_error = "independent")
  expect_true(fit@struc_fit@converged)
  expect_equal(fit@analysis_data@err_xy, matrix(0, 48, 2), ignore_attr = TRUE)
  cov <- repair_prediction_covariance(fit)
  cross <- cov[cov$model1 != cov$model2, ]
  expect_true(all(cross$source == "assumed_independent" & cross$covariance == 0))
  expect_equal(repair_diagnostics(fit, "mlm")$n_levels, c(48L, 48L))
  expect_equal(repair_reliability(fit)$n_units, rep(48L, 4))
})

test_that("independent target blocks retain nested nuisance covariance", {
  ex <- nested_example(); d <- ex$case$data
  model <- lme4::lmer(y ~ x + offset(off) + (x || person) + (1 | person:day), d,
    weights = prior, na.action = na.exclude)
  fit <- repair(ex$formula, list(growth = model), ex$case$people, "person")
  mf <- model.frame(model); id <- ex$case$people$person[1]
  rows <- which(as.character(mf$person) == id)
  z <- cbind(1, mf$x[rows])
  r <- diag(sigma(model)^2 / model.weights(mf)[rows]) +
    as.numeric(lme4::VarCorr(model)[[match("person:day", names(lme4::getME(model, "cnms")))]]) *
      outer(mf$day[rows], mf$day[rows], "==")
  expected <- solve(crossprod(z, solve(r, z)))
  expect_true(fit@struc_fit@converged)
  expect_equal(unname(fit@analysis_data@predictions@prediction_error_covariance$person[c("(Intercept)", "x"), c("(Intercept)", "x"), id]),
    unname(expected), tolerance = 1e-9)
})

test_that("failed target designs retain explicit omission and score-failure audits", {
  ex <- nested_example(); d <- ex$case$data
  d$x[d$person == "P001"] <- 0
  model <- lme4::lmer(y ~ x + offset(off) + (x | person) + (1 | person:day), d,
    weights = prior, na.action = na.exclude)
  expect_warning(expect_warning(fit <- repair(ex$formula, list(growth = model),
    ex$case$people, "person", na_action = "omit"), "Omitting 1"), "predictions failed")
  expect_equal(nobs(fit), 47L)
  bad <- repair_predictions(fit, scope = "all")
  expect_true(all(bad$status[bad$id == "P001"] == "failed"))
  expect_true(all(!bad$analyzed[bad$id == "P001"]))
  expect_match(repair_diagnostics(fit, "prediction_failures")$reason, "rank-deficient")
  expect_identical(repair_diagnostics(fit, "excluded")$id, "P001")
})

test_that("zero estimated shared variance does not waive the membership rule", {
  set.seed(20261030)
  d <- expand.grid(trial = 1:4, day = 1:4, person = factor(1:24))
  d$x <- (d$day - 2.5) / 2 + (d$trial - 2.5) / 3
  b0 <- rnorm(24); b1 <- rnorm(24, sd = .4)
  d$y <- 2 + .7 * d$x + b0[d$person] + b1[d$person] * d$x + c(.2, -.2, -.2, .2)[d$trial]
  model <- suppressMessages(lme4::lmer(y ~ x + (x | person) + (1 | day), d,
    control = lme4::lmerControl(optimizer = "bobyqa")))
  expect_lt(as.numeric(lme4::VarCorr(model)$day), 1e-8)
  people <- data.frame(person = levels(d$person), observed = seq_len(24))
  expect_error(repair_data(re("x", group = "person") ~ observed, model, people, "person"),
    "shared across target group")
})
