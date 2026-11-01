scoring_data <- function() {
  set.seed(20261022)
  sizes <- rep(4:9, 4)
  dat <- do.call(rbind, lapply(seq_along(sizes), function(i) {
    z <- seq(-1, 1, length.out = sizes[i]) + i / 50
    data.frame(
      id = paste0("person", i), z = z,
      y = .4 + .7 * z + rnorm(1) + rnorm(1, sd = .6) * z + rnorm(length(z), sd = .3),
      prior = seq(.5, 2, length.out = sizes[i]),
      off = sin(z)
    )
  }))
  dat$id <- factor(dat$id, levels = c(rev(unique(dat$id)), "unused"))
  dat <- dat[sample(nrow(dat)), ]
  dat$y[1L] <- NA_real_
  dat
}

test_that("direct predictions and covariances match per-cluster GLS after row alignment", {
  skip_if_not_installed("lme4")
  dat <- scoring_data()
  fit <- lme4::lmer(y ~ z + offset(off) + (1 + z | id), dat,
    weights = prior, na.action = na.omit
  )
  predictions <- repair_ranef_gls(fit)
  mf <- model.frame(fit)
  residual <- model.response(mf) - model.offset(mf) -
    drop(lme4::getME(fit, "X") %*% lme4::fixef(fit))
  expect_true(predictions@success)
  expect_true(predictions@diagnostics[["id"]]$success)
  expect_length(predictions@diagnostics[["id"]]$failed_levels, 0L)
  expect_true(all(is.finite(predictions@predictions[["id"]])))
  expect_true(all(is.finite(predictions@prediction_error_covariance[["id"]])))
  expect_identical(predictions@levels[["id"]], rownames(lme4::ranef(fit)[[1L]]))
  for (id in predictions@levels[["id"]]) {
    rows <- which(as.character(mf$id) == id)
    z <- cbind(1, mf$z[rows])
    weights <- model.weights(mf)[rows]
    expected <- lm.wfit(z, residual[rows], weights)$coefficients
    expected_cov <- sigma(fit)^2 * solve(crossprod(z * sqrt(weights)))
    expect_equal(unname(predictions@predictions[["id"]][id, ]), unname(expected), tolerance = 1e-9)
    expect_equal(unname(predictions@prediction_error_covariance[["id"]][, , id]), unname(expected_cov), tolerance = 1e-9)
  }
})

test_that("rank-deficient clusters are reported without discarding valid clusters", {
  skip_if_not_installed("lme4")
  dat <- scoring_data()
  bad_id <- "person1"
  dat$z[dat$id == bad_id] <- 1
  fit <- lme4::lmer(y ~ z + (1 + z | id), dat, na.action = na.omit)
  expect_warning(predictions <- repair_ranef_gls(fit), "predictions failed for 1 cluster")
  expect_false(predictions@success)
  expect_identical(predictions@diagnostics[["id"]]$failed_levels, bad_id)
  expect_true(all(is.na(predictions@predictions[["id"]][bad_id, ])))
  expect_true(all(is.na(predictions@prediction_error_covariance[["id"]][, , bad_id])))
  expect_match(predictions@diagnostics[["id"]]$per_level[[bad_id]]$error, "rank-deficient")
  expect_equal(predictions@diagnostics[["id"]]$per_level[[bad_id]]$information_solver$rank, 1L)
  expect_true(all(is.finite(predictions@predictions[["id"]][setdiff(predictions@levels[["id"]], bad_id), ])))
})

test_that("a singular random-effect covariance does not prevent direct prediction extraction", {
  skip_if_not_installed("lme4")
  dat <- expand.grid(z = seq(-1, 1, length.out = 6), id = factor(seq_len(12)))
  dat$y <- .4 + .7 * dat$z + rep(c(.1, -.2, .1, .1, -.2, .1), 12)
  fit <- suppressMessages(lme4::lmer(y ~ z + (1 + z | id), dat))
  expect_true(lme4::isSingular(fit))
  predictions <- repair_ranef_gls(fit)
  expect_true(predictions@success)
  expect_true(all(is.finite(predictions@predictions[["id"]])))
  expect_true(all(is.finite(predictions@prediction_error_covariance[["id"]])))
})

test_that("prediction extraction and M-MOM together recover the OLS no-predictor-error limit", {
  skip_if_not_installed("lme4")
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  predictions <- repair_ranef_gls(fit)
  # The balanced design gives equal outcome precision. With an observed
  # predictor, REPAIR must reduce to OLS of the likelihood-only slope predictions.
  x <- seq_along(predictions@levels[["Subject"]])
  result <- fit_mmom(list(
    outcome = predictions@predictions[["Subject"]][, "Days"], predictors = cbind(observed = x),
    err_x = matrix(0, 1, 1), err_y = predictions@prediction_error_covariance[["Subject"]]["Days", "Days", ]
  ))
  ordinary <- lm(predictions@predictions[["Subject"]][, "Days"] ~ x)
  expect_true(result@converged)
  expect_equal(result@coefficients$estimate, unname(coef(ordinary)), tolerance = 1e-9)
  expect_equal(unname(result@vcov), unname(vcov(ordinary)), tolerance = 1e-9)
})
