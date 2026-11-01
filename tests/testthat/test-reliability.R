test_that("posterior reliability reproduces independent manuscript/simulation fixtures", {
  skip_if_not_installed("nlme")
  reference <- readRDS(system.file("extdata", "reliability-reference.rds", package = "repair"))
  for (case in reference$cases) {
    correlated <- case$kind == "both_ar1_ident"
    model <- nlme::lme(y ~ z, random = ~ z | id, data = case$data,
      correlation = if (correlated) nlme::corAR1(form = ~ time | id) else NULL,
      weights = if (correlated) nlme::varIdent(form = ~ 1 | phase) else NULL,
      control = nlme::lmeControl(msMaxIter = 300))
    fit <- repair(re("z") ~ re("(Intercept)"), model,
      data.frame(person = rev(case$per_person$id)), "person")
    person <-     repair_reliability(fit, "unit")
    summary <- repair_reliability(fit)
    cases <- list(intercept = list(term = "(Intercept)", measure = "marginal"),
      slope = list(term = "z", measure = "marginal"),
      residualized = list(term = "z", measure = "residualized"))
    for (name in names(cases)) {
      spec <- cases[[name]]
      p <- person[person$term == spec$term & person$measure == spec$measure, ]
      s <- summary[summary$term == spec$term & summary$measure == spec$measure, ]
      order <- match(p$id, case$per_person$id)
      expect_equal(p$reliability, case$per_person[[name]][order], tolerance = 1e-5)
      expect_identical(p$n_observations, case$per_person$n_observations[order])
      expect_equal(s$reliability, unname(case$summary[name]), tolerance = 1e-5)
      expect_equal(s$reliability, mean(p$reliability), tolerance = 1e-12)
      expect_equal(s$n_units, nrow(case$per_person)); expect_equal(s$n_undefined, 0L)
      expect_true(all(p$status == "ok"))
    }
    # Cohort summaries use people, not number of observations, as weights.
    if (case$kind != "predictor_balanced") {
      p <- person[person$term == "z" & person$measure == "marginal", ]
      expect_gt(abs(mean(p$reliability) - weighted.mean(p$reliability, p$n_observations)), 1e-5)
    }
  }
})

test_that("lmer posterior variance agrees with conditional ranef and differs from GLS error", {
  fit <- reliability_example()
  model <- fit@mlm_fit$growth
  post <- attr(lme4::ranef(model, condVar = TRUE)$Subject, "postVar")
  g <- as.matrix(lme4::VarCorr(model)$Subject)
  p <-   repair_reliability(fit, "unit", type = "marginal")
  for (j in 1:2) {
    rows <- p[p$term == colnames(g)[j], ]
    expect_equal(rows$posterior_variance, as.numeric(post[j, j, ]), tolerance = 1e-8)
    expect_equal(rows$reliability, 1 - as.numeric(post[j, j, ]) / g[j, j], tolerance = 1e-8)
  }
  slope <- p[p$term == "Days", ]
  expect_gt(max(abs(slope$posterior_variance - fit@analysis_data@err_y)), 1e-3)
  expect_equal(summary(fit)$reliability, repair_reliability(fit))
  expect_output(print(fit), "Posterior reliability.*structural scope")
})

test_that("scopes preserve analyzed IDs and expose MLM-only units", {
  ex <- workflow_example()
  dat <- ex$data[18:3, ]; dat$x[1] <- NA_real_
  expect_warning(fit <- repair(re("Days") ~ x, ex$models, dat, "person", "omit"), "Omitting")
  p <-   repair_reliability(fit, "unit", type = "marginal")
  expect_identical(p$id, dat$person[-1])
  expect_equal(repair_reliability(fit)$n_units, c(15L, 15L))
  all <- repair_reliability(fit, "unit", "all", "marginal")
  expect_setequal(all$id, ex$data$person)
  expect_equal(repair_reliability(fit, scope = "all")$n_units, c(18L, 18L))
  expect_error(repair_reliability(ex$models), "repair")
  expect_error(repair_reliability(fit, scope = "requested"), "arg")
})

test_that("singular prior covariances retain defined marginals and explain undefined contrasts", {
  fit <- reliability_example()
  elements <- repair_get_mlm_elements(fit@mlm_fit$growth)
  original <- elements
  # Controlled singular covariance blocks isolate boundary definitions from
  # optimizer behavior; fitted design and residual covariance are unchanged.
  local_mocked_bindings(repair_get_mlm_elements = function(model) elements)
  for (g in list(diag(c(1, 0)), matrix(1, 2, 2), diag(c(0, 1)), matrix(0, 2, 2))) {
    elements <- original
    elements@random_covariance <- list(Subject = kronecker(diag(18), g))
    p <-     repair_reliability(fit, "unit")
    marg <- p[p$measure == "marginal", ]
    for (j in 1:2) {
      rows <- marg[marg$term == c("(Intercept)", "Days")[j], ]
      if (g[j, j] == 0) {
        expect_true(all(is.na(rows$reliability)))
        expect_true(all(rows$status == "undefined_variance"))
      } else expect_true(all(is.finite(rows$reliability)))
    }
    r <- p[p$measure == "residualized", ]
    expect_true(all(is.na(r$reliability)))
    expect_true(all(nzchar(r$reason)))
    s <- repair_reliability(fit)
    expect_true(all(is.na(s$reliability[s$measure == "residualized"])))
    expect_true(all(s$n_defined[s$measure == "residualized"] == 0))
  }
})

test_that("posterior covariance handles singular priors and unidentifiable likelihood designs", {
  z <- cbind(1, 0:4); r <- diag(5); a <- c(1, .5); g <- tcrossprod(a)
  expected <- g / (1 + drop(crossprod(z %*% a)))
  expect_equal(repair_posterior_covariance(g, z, r), expected, tolerance = 1e-12)
  expect_equal(repair_posterior_covariance(matrix(0, 2, 2), z, r), matrix(0, 2, 2))
  d <- lme4::sleepstudy; d$Days[d$Subject == "308"] <- 0
  model <- suppressWarnings(lme4::lmer(Reaction ~ Days + (Days | Subject), d))
  people <- workflow_example()$data
  expect_warning(expect_warning(fit <- repair(re("Days") ~ x, model, people, "person", "omit"), "Omitting"), "predictions failed")
  p <-   repair_reliability(fit, "unit", "all")
  expect_true(all(is.finite(p$reliability[p$id == "308"])))
  expect_false("308" %in% repair_reliability(fit, "unit")$id)
  diag <- repair_diagnostics(fit)
  expect_identical(diag$prediction_failures$id, "308")
  expect_false(diag$prediction_failures$analyzed)
  expect_match(diag$prediction_failures$reason, "rank-deficient")
})

test_that("undefined people are not silently removed from summary denominators", {
  p <-   repair_reliability(reliability_example(), "unit", type = "marginal")
  p$reliability[1] <- NA_real_; p$status[1] <- "posterior_failure"; p$reason[1] <- "test numerical failure"
  s <- repair_summarize_reliability(p)
  bad <- s[s$term == p$term[1], ]
  expect_true(is.na(bad$reliability)); expect_equal(bad$n_units, 18L)
  expect_equal(bad$n_defined, 17L); expect_equal(bad$n_undefined, 1L)
  expect_equal(bad$mean_defined, mean(p$reliability[p$term == p$term[1]], na.rm = TRUE))
  expect_match(bad$reason, "test numerical failure")
})

test_that("multiple-model reliability is model-specific and uses a common structural scope", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example()
  fit <- repair(ex$formula, ex$models, ex$case$people[1:60, , drop = FALSE], "person",
    cross_model_error = ex$case$cross_model_error)
  p <-   repair_reliability(fit, "unit")
  expect_setequal(unique(p$model), c("growth", "response"))
  expect_setequal(unique(p$id), ex$case$people$person[1:60])
  expect_true(all(repair_reliability(fit)$n_units == 60L))
  expect_true(all(repair_reliability(fit, scope = "all")$n_units == 70L))
  independent <- repair(ex$formula, ex$models, ex$case$people[1:60, , drop = FALSE], "person", cross_model_error = "independent")
  expect_equal(repair_reliability(independent), repair_reliability(fit))
  expect_equal(repair_diagnostics(fit)$cross_model_error, fit@analysis_data@diagnostics$cross_model_error)
})

test_that("residualization is explicit for effects outside the supported two-term case", {
  ex <- workflow_example()
  fit <- repair(re("(Intercept)") ~ x, ex$models, ex$data, "person")
  r <- repair_reliability(fit, type = "residualized")
  expect_true(is.na(r$reliability)); expect_match(r$reason, "two-coefficient")
  expect_identical(repair_reliability(fit)$measure, "marginal")
})

test_that("reliability plots dispatch outside the namespace and restore graphics state", {
  p <-   repair_reliability(reliability_example(), "unit")
  path <- tempfile(fileext = ".pdf"); grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  graphics::par(cex = .8)
  old <- graphics::par(c("mfrow", "cex", "mar"))
  external <- list2env(list(x = p), parent = asNamespace("graphics"))
  expect_equal(evalq(plot(x), external), p)
  expect_equal(graphics::par(c("mfrow", "cex", "mar")), old)
  expect_equal(plot(p, type = "observations", pch = 19), p)
  p$reliability[] <- NA_real_
  expect_equal(plot(p), p)
  expect_error(plot(repair_reliability(reliability_example())), "unit-level")
  expect_error(plot(p[FALSE, ]), "No reliability rows")
})
