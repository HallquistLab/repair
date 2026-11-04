# These tests sample real models. Run explicitly; ordinary package checks exercise
# the Gaussian contract without a compiler or MCMC in their test loop.
bayes_native_test <- function(backend) {
  skip_if(Sys.getenv("REPAIR_RUN_STAN_TESTS") != "true", "set REPAIR_RUN_STAN_TESTS=true for native sampling tests")
  skip_if_not_installed("brms")
  skip_if_not_installed("rstanarm")
  skip_if_not_installed("posterior")
  set.seed(803)
  data <- expand.grid(t = seq(-.4, 1.4, length.out = 8), id = factor(paste0("person", 1:30)))
  data$w <- rep(rnorm(30), each = 8)
  data$y <- .5 + .4 * data$t + .2 * data$w + rep(rnorm(30, sd = .8), each = 8) +
    rep(rnorm(30, sd = .6), each = 8) * data$t + rnorm(nrow(data), sd = .5)
  # Native row handling, unused levels, shuffled IDs, and omission must survive.
  data <- data[sample(nrow(data)), ]
  data$id <- factor(data$id, levels = c(rev(levels(data$id)), "unused"))
  data$y[1L] <- NA_real_
  formula <- y ~ t + w + (1 + t | id)
  # stan_lmer() in rstanarm 2.32.1 evaluates an unqualified stan_glmer call
  # in its caller, including when stan_lmer itself was namespace-qualified.
  stan_glmer <- rstanarm::stan_glmer
  model <- switch(backend,
    brms = suppressWarnings(brms::brm(formula, data = data, family = gaussian(),
      backend = Sys.getenv("REPAIR_BRMS_BACKEND", "rstan"), seed = 811,
      stan_model_args = if (Sys.getenv("REPAIR_BRMS_BACKEND", "rstan") == "cmdstanr")
        list(cpp_options = list(PRECOMPILED_HEADERS = "false")) else list(),
      chains = 2, cores = 2, iter = 1200, warmup = 600, refresh = 0,
      prior = c(brms::prior(normal(0, 2), class = "b"),
        brms::prior(normal(0, 2), class = "Intercept"),
        brms::prior(exponential(1), class = "sd"), brms::prior(exponential(1), class = "sigma")),
      control = list(adapt_delta = .95, max_treedepth = 12))),
    rstanarm = suppressWarnings(rstanarm::stan_lmer(formula, data = data,
      seed = 812, chains = 2, cores = 2, iter = 1200, warmup = 600,
      refresh = 0, adapt_delta = .95, control = list(max_treedepth = 12)))
  )
  effective <- data[complete.cases(data[c("y", "t", "w", "id")]), ]
  people <- data.frame(id = rev(unique(as.character(effective$id))), x = seq_len(30))
  expect_warning(snapshot <- repair_bayes_plugin(model), "Experimental Bayesian plug-in")
  expect_identical(snapshot$diagnostics$mlm_backend, backend)
  expect_identical(snapshot$diagnostics$sampler$chains, 2L)
  expect_true(is.finite(snapshot$diagnostics$sampler$max_observed_treedepth))
  expect_equal(snapshot$diagnostics$sampler$max_treedepth, 12)
  expect_true(length(snapshot$diagnostics$priors) > 0L)
  expect_output(print(snapshot), "Experimental REPAIR Bayesian plug-in snapshot")
  e <- snapshot$elements
  raw <- if (backend == "brms") brms::as_draws_df(model)$sigma else as.matrix(model)[, "sigma"]
  expect_equal(e@diagnostics$residual_variance, mean(raw^2), tolerance = 1e-12)
  expect_gt(abs(e@diagnostics$residual_variance - mean(raw)^2), 0)
  pred <- repair_ranef_gls(snapshot)
  expect_true(pred@success)
  expect_false("unused" %in% pred@levels$id)
  for (id in pred@levels$id) {
    rows <- which(as.character(effective$id) == id)
    observed <- effective[rows, ]
    residual <- observed$y - drop(model.matrix(~ t + w, observed) %*% e@fixed_coefs)
    z <- model.matrix(~ t, observed)
    expect_equal(unname(pred@predictions$id[id, ]), unname(lm.fit(z, residual)$coefficients), tolerance = 1e-8)
    expect_equal(unname(pred@prediction_error_covariance$id[, , id]),
      unname(mean(raw^2) * solve(crossprod(z))), tolerance = 1e-8)
  }
  fit <- repair(re("t") ~ re("(Intercept)") + x, snapshot, people, "id")
  expect_true(fit@struc_fit@converged)
  expect_true(all(is.finite(repair_coefficients(fit)$std.error)))
  expect_identical(summary(fit)$mlm_models$mlm$coverage, "full_refitted_workflow_unvalidated")
  expect_match(repair_diagnostics(fit, table = "mlm")$messages, "not propagated")
  expect_true(all(is.finite(repair_reliability(fit)$reliability)))
  expect_error(repair_refit(fit, effective, people), "not yet supported")
  if (backend == "brms") {
    bad <- model; bad$family$link <- "log"
    expect_error(repair_bayes_plugin(bad), "Gaussian identity-link")
    bad <- model; bad$algorithm <- "meanfield"
    expect_error(repair_bayes_plugin(bad), "sampled")
    bad <- model; bad$formula$pforms <- list(sigma = ~ t)
    expect_error(repair_bayes_plugin(bad), "distributional/nonlinear")
    bad <- model; bad$formula$formula <- y | weights(w) ~ t + (t | id)
    expect_error(repair_bayes_plugin(bad), "support requires")
  } else {
    bad <- model; bad$weights <- rep(1, nrow(effective))
    expect_error(repair_bayes_plugin(bad), "Weights and offsets")
    bad <- model; bad$call$prior_PD <- TRUE
    expect_error(repair_bayes_plugin(bad), "Prior-predictive")
    bad <- model; bad$algorithm <- "meanfield"
    expect_error(repair_bayes_plugin(bad), "sampled")
  }
}

test_that("actual brms samples feed the complete plug-in workflow", {
  bayes_native_test("brms")
})

test_that("actual rstanarm samples feed the complete plug-in workflow", {
  bayes_native_test("rstanarm")
})
