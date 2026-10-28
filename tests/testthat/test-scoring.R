scoring_data <- function() {
  set.seed(20261022)
  sizes <- rep(4:9, 4)
  dat <- do.call(rbind, lapply(seq_along(sizes), function(i) {
    z <- seq(-1, 1, length.out = sizes[i]) + i / 50
    data.frame(id = paste0('person', i), z = z,
      y = .4 + .7 * z + rnorm(1) + rnorm(1, sd = .6) * z + rnorm(length(z), sd = .3),
      prior = seq(.5, 2, length.out = sizes[i]),
      off = sin(z))
  }))
  dat$id <- factor(dat$id, levels = c(rev(unique(dat$id)), 'unused'))
  dat <- dat[sample(nrow(dat)), ]
  dat$y[1L] <- NA_real_
  dat
}

test_that('direct scores and covariances match per-cluster GLS after row alignment', {
  skip_if_not_installed('lme4')
  dat <- scoring_data()
  fit <- lme4::lmer(y ~ z + offset(off) + (1 + z | id), dat,
    weights = prior, na.action = na.omit)
  scores <- ranef_like(fit)
  mf <- model.frame(fit)
  residual <- model.response(mf) - model.offset(mf) -
    drop(lme4::getME(fit, 'X') %*% lme4::fixef(fit))
  expect_true(scores$success)
  expect_length(scores$diagnostics$failed_clusters, 0L)
  expect_true(all(is.finite(scores$estimates)))
  expect_true(all(is.finite(scores$err_cov)))
  expect_identical(scores$levels, rownames(lme4::ranef(fit)[[1L]]))
  expect_identical(scores$score_cov, scores$err_cov)
  for (id in scores$levels) {
    rows <- which(as.character(mf$id) == id)
    z <- cbind(1, mf$z[rows])
    weights <- model.weights(mf)[rows]
    expected <- lm.wfit(z, residual[rows], weights)$coefficients
    expected_cov <- sigma(fit)^2 * solve(crossprod(z * sqrt(weights)))
    expect_equal(unname(scores$estimates[id, ]), unname(expected), tolerance = 1e-9)
    expect_equal(unname(scores$err_cov[, , id]), unname(expected_cov), tolerance = 1e-9)
  }
})

test_that('rank-deficient clusters are reported without discarding valid clusters', {
  skip_if_not_installed('lme4')
  dat <- scoring_data()
  bad_id <- 'person1'
  dat$z[dat$id == bad_id] <- 1
  fit <- lme4::lmer(y ~ z + (1 + z | id), dat, na.action = na.omit)
  expect_warning(scores <- ranef_like(fit), 'scoring failed for 1 cluster')
  expect_false(scores$success)
  expect_identical(scores$diagnostics$failed_clusters, bad_id)
  expect_true(all(is.na(scores$estimates[bad_id, ])))
  expect_true(all(is.na(scores$err_cov[, , bad_id])))
  expect_match(scores$diagnostics$clusters[[bad_id]]$error, 'rank-deficient')
  expect_equal(scores$diagnostics$clusters[[bad_id]]$rank, 1L)
  expect_true(all(is.finite(scores$estimates[setdiff(scores$levels, bad_id), ])))
})

test_that('a singular random-effect covariance does not prevent direct scoring', {
  skip_if_not_installed('lme4')
  dat <- expand.grid(z = seq(-1, 1, length.out = 6), id = factor(seq_len(12)))
  dat$y <- .4 + .7 * dat$z + rep(c(.1, -.2, .1, .1, -.2, .1), 12)
  fit <- suppressMessages(lme4::lmer(y ~ z + (1 + z | id), dat))
  expect_true(lme4::isSingular(fit))
  scores <- ranef_like(fit)
  expect_true(scores$success)
  expect_true(all(is.finite(scores$estimates)))
  expect_true(all(is.finite(scores$err_cov)))
})

test_that('unsupported model families and multiple grouping factors fail explicitly', {
  skip_if_not_installed('lme4')
  glmm <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
    lme4::cbpp, family = binomial)
  expect_error(ranef_like(glmm), 'Gaussian')
  dat <- scoring_data()
  dat$batch <- factor(rep(1:3, length.out = nrow(dat)))
  fit <- suppressMessages(lme4::lmer(y ~ z + (1 | id) + (1 | batch), dat))
  expect_error(ranef_like(fit), 'exactly one grouping factor')
})

test_that('scoring and M-MOM together recover the OLS no-predictor-error limit', {
  skip_if_not_installed('lme4')
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  scores <- ranef_like(fit)
  # The balanced design gives equal outcome precision. With an observed
  # predictor, REPAIR must reduce to OLS of the likelihood-only slope scores.
  x <- seq_along(scores$levels)
  result <- fit_mmom(list(
    outcome = scores$estimates[, 'Days'], predictors = cbind(observed = x),
    err_x = matrix(0, 1, 1), err_y = scores$err_cov['Days', 'Days', ]
  ))
  ordinary <- lm(scores$estimates[, 'Days'] ~ x)
  expect_true(result@converged)
  expect_equal(result@coefficients$estimate, unname(coef(ordinary)), tolerance = 1e-9)
  expect_equal(unname(result@vcov), unname(vcov(ordinary)), tolerance = 1e-9)
})
