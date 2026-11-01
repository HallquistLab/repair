test_that("cross-model covariance aligns by model, coefficient, and person identity", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); case <- ex$case; cross <- case$cross_model_error
  baseline <- repair(ex$formula, ex$models, case$people, "person", cross_model_error = cross)
  # Reverse every independently orderable dimension, and reverse block orientation.
  swapped <- list(list(models = rev(cross[[1]]$models),
    covariance = aperm(cross[[1]]$covariance[2:1, , 70:1, drop = FALSE], c(2, 1, 3))))
  f <- re("w", "response") ~ re("z", "growth") + re("(Intercept)", "growth")
  changed <- repair(f, rev(ex$models), case$people[70:1, , drop = FALSE], "person", cross_model_error = swapped)
  order <- c(1, 3, 2)
  expect_equal(coef(changed), coef(baseline)[order], tolerance = 1e-9)
  expect_equal(vcov(changed), vcov(baseline)[order, order], tolerance = 1e-9)
  expect_equal(changed@analysis_data@err_xy, baseline@analysis_data@err_xy[70:1, 2:1], tolerance = 1e-9)
  expect_equal(summary(changed)$cross_model_error$blocks, swapped)
  zero <- cross; zero[[1]]$covariance[] <- 0
  independent <- repair(ex$formula, ex$models, case$people, "person", cross_model_error = "independent")
  explicit <- repair(ex$formula, ex$models, case$people, "person", cross_model_error = zero)
  expect_equal(coef(explicit), coef(independent))
  expect_equal(vcov(explicit), vcov(independent))
  expect_gt(max(abs(coef(baseline) - coef(independent))), 1e-4)
  expect_gt(max(abs(vcov(baseline) - vcov(independent))), 1e-4)
  # An observed outcome moves all selected cross-model covariance into err_x.
  people <- case$people; people$observed <- seq_len(nrow(people))
  a <- repair_data(observed ~ re("w", "response") + re("z", "growth") + re("(Intercept)", "growth"),
    ex$models, people, "person", cross_model_error = cross)
  expect_equal(unname(a@err_x[1, 2:3, ]), unname(t(baseline@analysis_data@err_xy[, 2:1])))
  expect_true(all(a@err_y == 0)); expect_true(all(a@err_xy == 0))
})

test_that("all referenced models contribute ID and failed-score exclusions", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); case <- ex$case
  # Missing from different model cohorts; different grouping column names.
  ex$models$growth <- nlme::lme(y ~ z, random = ~ z | id, data = case$a[case$a$id != "p01", ])
  ex$models$response <- nlme::lme(q ~ w, random = ~ w | subject, data = case$b[case$b$subject != "p02", ])
  expect_error(repair_data(ex$formula, ex$models, case$people, "person", cross_model_error = "independent"), "2 analysis row")
  expect_warning(a <- repair_data(ex$formula, ex$models, case$people, "person", "omit", "independent"), "Omitting 2")
  expect_identical(a@ids, case$people$person[-(1:2)])
  expect_match(a@diagnostics$excluded$reason[1], "MLM growth")
  expect_match(a@diagnostics$excluded$reason[2], "MLM response")
  # Supplied blocks need only cover retained IDs; omissions cannot hide bad blocks.
  cross <- case$cross_model_error
  cross[[1]]$covariance <- cross[[1]]$covariance[, , -(1:2), drop = FALSE]
  expect_warning(b <- repair_data(ex$formula, ex$models, case$people, "person", "omit", cross), "Omitting 2")
  expect_identical(b@ids, a@ids)
  cross[[1]]$covariance <- cross[[1]]$covariance[, , -1, drop = FALSE]
  expect_error(repair_data(ex$formula, ex$models, case$people, "person", "omit", cross), "missing retained IDs.*p03")
  subset <- repair_data(ex$formula, ex$models, case$people[3:65, , drop = FALSE], "person", cross_model_error = "independent")
  expect_identical(subset@diagnostics$mlm_models$growth$unused_prediction_ids, c("p02", sprintf("p%02d", 66:70)))
  # A cluster present in the model but without an identifiable random slope.
  d <- lme4::sleepstudy; d$Days[d$Subject == "308"] <- 0
  bad <- suppressWarnings(lme4::lmer(Reaction ~ Days + (Days | Subject), d))
  good <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  people <- data.frame(person = levels(d$Subject))
  expect_warning(expect_warning(failed <- repair_data(re("Days", "bad") ~ re("(Intercept)", "good"),
    list(bad = bad, good = good), people, "person", "omit", "independent"), "Omitting 1"), "predictions failed")
  expect_match(failed@diagnostics$excluded$reason, "prediction extraction failed in model bad")
})

test_that("invalid or implicit cross-model covariance never falls back to independence", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); case <- ex$case
  assemble <- function(value) repair_data(ex$formula, ex$models, case$people, "person", cross_model_error = value)
  for (value in list(NULL, "auto", TRUE, list())) expect_error(assemble(value), "explicit")
  cross <- case$cross_model_error
  expect_error(assemble(c(cross, cross)), "Duplicate")
  bad <- cross; bad[[1]]$models <- c("growth", "growth")
  expect_error(assemble(bad), "two distinct")
  bad <- cross; bad[[1]]$models <- c("growth", "absent")
  expect_error(assemble(bad), "two distinct")
  bad <- cross; bad[[1]]$extra <- 1
  expect_error(assemble(bad), "exactly")
  bad <- cross; bad[[1]]$covariance[1] <- NA_real_
  expect_error(assemble(bad), "finite numeric")
  bad <- cross; dimnames(bad[[1]]$covariance) <- NULL
  expect_error(assemble(bad), "dimnames")
  bad <- cross; dimnames(bad[[1]]$covariance)[[3]][2] <- "p01"
  expect_error(assemble(bad), "unique")
  bad <- cross; dimnames(bad[[1]]$covariance)[[1]][1] <- "wrong"
  expect_error(assemble(bad), "term names")
  bad <- cross; bad[[1]]$covariance[, , 1] <- 100
  expect_error(assemble(bad), "positive semidefinite.*p01")
  expect_error(repair_data(re("z", "growth") ~ observed, ex$models,
    transform(case$people, observed = 1:70), "person", cross_model_error = "independent"), "only when.*multiple")
  # Unused catalog entries do not create spurious requirements.
  one <- repair_data(re("z", "growth") ~ observed, ex$models, transform(case$people, observed = 1:70), "person")
  expect_identical(one@diagnostics$cross_model_error$policy, "not_applicable")
})

test_that("three-model assembly requires all pairs and checks joint positive semidefiniteness", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); case <- ex$case
  models <- list(a = ex$models$growth, b = ex$models$response, c = ex$models$response)
  f <- re("w", "b") ~ re("z", "a") + re("(Intercept)", "c")
  independent <- repair_data(f, models, case$people, "person", cross_model_error = "independent")
  var <- list(a = independent@err_x[1, 1, ], b = independent@err_y, c = independent@err_x[2, 2, ])
  terms <- c(a = "z", b = "w", c = "(Intercept)")
  pair <- function(a, b, rho) list(models = c(a, b), covariance = array(rho * sqrt(var[[a]] * var[[b]]),
    c(1, 1, 70), dimnames = list(terms[[a]], terms[[b]], case$people$person)))
  blocks <- list(pair("a", "b", .9), pair("a", "c", .9), pair("b", "c", -.9))
  expect_error(repair_data(f, models, case$people, "person", cross_model_error = blocks[-3]), "every pair")
  expect_error(repair_data(f, models, case$people, "person", cross_model_error = blocks), "positive semidefinite")
  blocks[[3]] <- pair("b", "c", .9)
  valid <- repair_data(f, models, case$people, "person", cross_model_error = blocks)
  expect_equal(unname(valid@err_xy[, 1]), unname(.9 * sqrt(var$a * var$b)))
  expect_equal(unname(valid@err_xy[, 2]), unname(.9 * sqrt(var$c * var$b)))
  expect_equal(unname(valid@err_x[1, 2, ]), unname(.9 * sqrt(var$a * var$c)))
})

test_that("mixed lmer and nlme models preserve their own prediction covariances", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example("independent"); case <- ex$case
  ex$models$growth <- lme4::lmer(y ~ z + (z | id), case$a)
  fit <- repair(ex$formula, ex$models, case$people, "person", cross_model_error = "independent")
  a <- fit@analysis_data
  expect_true(fit@struc_fit@converged)
  expect_equal(summary(fit)$mlm_models$growth$mlm_backend, "lme4")
  expect_equal(summary(fit)$mlm_models$response$mlm_backend, "nlme")
  direct <- repair_ranef_gls(ex$models$growth)
  expect_equal(unname(a@err_x), unname(direct@prediction_error_covariance$id[, , a@ids]))
  expect_true(all(a@err_xy == 0))
  expect_output(print(fit), "growth \\(lme4")
  expect_output(print(fit), "response \\(nlme")
})
