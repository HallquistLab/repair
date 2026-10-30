coefficient_columns <- c("term", "estimate", "std.error", "statistic", "p.value",
  "conf.low", "conf.high", "conf.level")
prediction_columns <- c("id", "model", "group", "term", "label", "role", "prediction", "std.error",
  "prediction.error.variance", "scope", "analyzed", "status", "reason")
covariance_columns <- c("id", "model1", "group1", "term1", "label1", "role1",
  "model2", "group2", "term2", "label2", "role2", "covariance", "source")

test_that("coefficient tables have a stable schema and reproduce standard inference", {
  ex <- workflow_example()
  fit <- repair(x ~ re("Days", "growth") + category, ex$models, ex$data, "person")
  snapshot <- function(x) list(coefficients = x@struc_fit@coefficients, vcov = x@struc_fit@vcov,
    ids = x@analysis_data@ids, err_x = x@analysis_data@err_x, err_y = x@analysis_data@err_y,
    err_xy = x@analysis_data@err_xy)
  original <- snapshot(fit)
  table <- repair_coefficients(fit)
  expect_identical(class(table), "data.frame")
  expect_identical(names(table), coefficient_columns)
  expect_identical(vapply(table, typeof, character(1)),
    setNames(c("character", rep("double", 7)), coefficient_columns))
  expect_identical(table$term, names(coef(fit)))
  expect_equal(table$estimate, unname(coef(fit)))
  expect_equal(table$std.error, unname(sqrt(diag(vcov(fit)))))
  expect_equal(as.matrix(table[c("conf.low", "conf.high")]), unname(confint(fit)), ignore_attr = TRUE)
  expect_equal(table$statistic, table$estimate / table$std.error)
  expect_equal(table$p.value, 2 * pnorm(abs(table$statistic), lower.tail = FALSE))
  expect_identical(tail(table$term, 2), c("categoryB", "categoryC"))
  alternative <- repair_coefficients(fit, conf.level = .9)
  expect_identical(table[1:5], alternative[1:5])
  expect_equal(alternative$conf.low, table$estimate - qnorm(.95) * table$std.error)
  expect_equal(alternative$conf.high, table$estimate + qnorm(.95) * table$std.error)
  expect_true(all(alternative$conf.level == .9))
  expect_identical(snapshot(fit), original)
  for (level in list(NA_real_, 0, 1, Inf, c(.8, .9), ".9")) {
    expect_error(repair_coefficients(fit, level), "level")
  }
  expect_error(repair_coefficients(fit@analysis_data), "repair")
})

test_that("score tables reproduce saved predictions without model-count-dependent structure", {
  ex <- workflow_example()
  f <- re("Days", "growth") ~ re("(Intercept)", "growth") + x
  fit <- repair(f, ex$models, ex$data[18:1, ], "person")
  input <- repair_data(f, ex$models, ex$data[18:1, ], "person")
  predictions <- repair_predictions(fit)
  expect_identical(class(predictions), "data.frame")
  expect_identical(names(predictions), prediction_columns)
  expect_identical(vapply(predictions, typeof, character(1)), setNames(
    c(rep("character", 6), rep("double", 3), "character", "logical", "character", "character"), prediction_columns))
  expect_identical(predictions, repair_predictions(input))
  slope <- predictions[predictions$role == "outcome", ]
  intercept <- predictions[predictions$role == "predictor", ]
  expect_identical(slope$id, input@ids)
  expect_equal(slope$prediction, unname(input@outcome))
  expect_equal(slope$prediction.error.variance, unname(input@err_y))
  expect_equal(intercept$prediction, unname(input@predictors[, 1]))
  expect_equal(intercept$prediction.error.variance, unname(input@err_x[1, 1, ]))
  expect_equal(predictions$std.error^2, predictions$prediction.error.variance)
  expect_true(all(predictions$analyzed)); expect_true(all(predictions$status == "ok"))
  expect_true(all(predictions$reason == "")); expect_equal(nrow(predictions), 36L)
  expect_equal(anyDuplicated(predictions[c("model", "group", "term", "id")]), 0L)
  expect_error(repair_predictions(list()), "repair_data")
  expect_error(repair_predictions(fit, "requested"), "arg")
})

test_that("MLM scopes distinguish failed predictions, exclusions, and unused people", {
  ex <- workflow_example(); d <- lme4::sleepstudy
  d$Days[d$Subject == "308"] <- 0
  model <- suppressWarnings(lme4::lmer(Reaction ~ Days + (Days | Subject), d))
  people <- ex$data[-18, ]; people$x[2] <- NA_real_; people$person[3] <- "absent"
  expect_warning(expect_warning(fit <- repair(re("Days") ~ x, model, people, "person", "omit"), "Omitting 3"), "predictions failed")
  p <- repair_predictions(fit, "all")
  expect_identical(p$id, fit@analysis_data@predictions@levels$Subject)
  expect_true(is.na(p$prediction[p$id == "308"]))
  expect_true(is.na(p$std.error[p$id == "308"]))
  expect_equal(p$status[p$id == "308"], "failed")
  expect_match(p$reason[p$id == "308"], "rank-deficient")
  expect_false(p$analyzed[p$id == people$person[2]])
  expect_equal(p$status[p$id == people$person[2]], "ok")
  expect_false(p$analyzed[p$id == ex$data$person[18]])
  expect_false("absent" %in% p$id)
  expect_true("absent" %in% repair_diagnostics(fit, "excluded")$id)
  expect_identical(repair_predictions(fit)$id, fit@analysis_data@ids)
})

test_that("covariance tables reconstruct within-model and cross-model error inputs", {
  skip_if_not_installed("nlme")
  for (kind in c("independent", "supplied")) {
    ex <- multimodel_example(kind); case <- ex$case
    for (backend in c("r", "rcpp")) {
      fit <- repair(ex$formula, ex$models, case$people, "person", backend = backend,
        cross_model_error = case$cross_model_error)
      b <- repair_coefficients(fit)
      expect_equal(b$estimate, case$coefficients$estimate, tolerance = 1e-5)
      expect_equal(b$std.error, case$coefficients$se, tolerance = 1e-5)
      expect_equal(b$conf.low, case$coefficients$ci_low, tolerance = 1e-5)
      expect_equal(b$conf.high, case$coefficients$ci_high, tolerance = 1e-5)
      s <- repair_predictions(fit); cov <- repair_prediction_covariance(fit)
      expect_identical(names(s), prediction_columns)
      expect_identical(class(cov), "data.frame")
      expect_identical(names(cov), covariance_columns)
      expect_identical(vapply(cov, typeof, character(1)),
        setNames(c(rep("character", 11), "double", "character"), covariance_columns))
      expect_equal(nrow(cov), 70L * 6L)
      expect_equal(anyDuplicated(cov[c("id", "model1", "group1", "term1", "model2", "group2", "term2")]), 0L)
      expect_identical(cov, repair_prediction_covariance(fit@analysis_data))
      # Reconstruct the full joint matrix using only exported table identities.
      labels <- unique(s$label)
      expected <- array(0, c(3, 3, 70))
      expected[1, 1, ] <- case$input$err_y
      expected[2:3, 2:3, ] <- case$input$err_x
      expected[1, 2:3, ] <- expected[2:3, 1, ] <- t(case$input$err_xy)
      actual <- array(0, c(3, 3, 70))
      for (row in seq_len(nrow(cov))) {
        i <- match(cov$id[row], case$people$person)
        a <- match(cov$label1[row], labels); b <- match(cov$label2[row], labels)
        actual[a, b, i] <- actual[b, a, i] <- cov$covariance[row]
      }
      expect_equal(actual, expected, tolerance = 1e-5)
      expect_true(all(cov$source[cov$model1 == cov$model2] == "within_model"))
      cross <- cov[cov$model1 != cov$model2, ]
      expect_true(all(cross$source == if (kind == "independent") "assumed_independent" else "supplied"))
      diagonal <- cov[cov$label1 == cov$label2, ]
      expect_equal(diagonal$covariance, s$prediction.error.variance)
      expect_identical(diagonal$id, s$id)
    }
  }
})

test_that("covariance extraction handles observed outcomes and single selected effects", {
  ex <- workflow_example()
  single <- repair_data(x ~ re("Days"), ex$models, ex$data, "person")
  cov <- repair_prediction_covariance(single)
  expect_equal(nrow(cov), 18L)
  expect_true(all(cov$term1 == "Days" & cov$term2 == "Days"))
  expect_equal(cov$covariance, unname(single@err_x[1, 1, ]))
  f <- x ~ re("Days") + re("(Intercept)") + category
  both <- repair_data(f, ex$models, ex$data[18:1, ], "person")
  cov <- repair_prediction_covariance(both)
  expect_equal(nrow(cov), 54L)
  expect_identical(cov$id[1:18], ex$data$person[18:1])
  expect_equal(cov$covariance[19:36], unname(both@err_x[2, 1, ]))
  expect_false(any(grepl("category", cov$term1)))
  expect_error(repair_prediction_covariance(list()), "repair_data")
})

test_that("diagnostic selectors preserve tables and avoid unrequested reliability computation", {
  fit <- reliability_example()
  all <- repair_diagnostics(fit)
  expect_identical(names(all), c("status", "scope", "excluded", "mlm", "prediction_failures",
    "prediction_error_variances", "regularization", "conditioning", "weights", "weight_summary", "reliability", "cross_model_error"))
  for (name in setdiff(names(all), "cross_model_error")) {
    table <- repair_diagnostics(fit, table = name)
    expect_identical(class(table), "data.frame")
    expect_equal(table, as.data.frame(all[[name]]), ignore_attr = TRUE)
    expect_identical(rownames(table), as.character(seq_len(nrow(table))))
  }
  for (bad in list("weight", "cross_model_error", NA_character_, 1, c("weights", "status"), "")) {
    expect_error(repair_diagnostics(fit, bad), "exact diagnostic table name")
  }
  local_mocked_bindings(repair_reliability = function(...) stop("unexpected reliability computation"))
  expect_equal(repair_diagnostics(fit, "scope"), all$scope)
  expect_equal(repair_coefficients(fit)$estimate, unname(coef(fit)))
  expect_equal(nrow(repair_predictions(fit)), 36L)
  expect_equal(nrow(repair_prediction_covariance(fit)), 54L)
})

test_that("measurement extractors remain usable after structural failure", {
  ex <- workflow_example(); ex$data$x <- 0
  expect_warning(fit <- repair(x ~ re("Days"), ex$models, ex$data, "person"), "Structural fit failed")
  expect_error(repair_coefficients(fit), "did not converge")
  expect_equal(nrow(repair_predictions(fit)), 18L)
  expect_true(all(repair_predictions(fit)$status == "ok"))
  expect_equal(nrow(repair_prediction_covariance(fit)), 18L)
  expect_false(repair_diagnostics(fit, "status")$converged)
  empty <- repair_diagnostics(fit, "weights")
  expect_identical(names(empty), c("id", "weight", "fraction"))
  expect_identical(vapply(empty, typeof, character(1)), c(id = "character", weight = "double", fraction = "double"))
  expect_equal(nrow(empty), 0L)
})

test_that("extractor identities survive unusual model names and portable table serialization", {
  ex <- workflow_example()
  models <- setNames(ex$models, "growth process")
  f <- re("Days", model = "growth process") ~ x
  fit <- repair(f, models, ex$data, "person")
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  tables <- list(predictions = repair_predictions(fit), covariance = repair_prediction_covariance(fit),
                 coefficients = repair_coefficients(fit))
  saveRDS(tables, path); restored <- readRDS(path)
  expect_identical(restored, tables)
  expect_identical(unique(restored$predictions$model), "growth process")
  expect_identical(unique(restored$predictions$label), 're("Days", model = "growth process")')
})
