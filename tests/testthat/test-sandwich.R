if_example <- function(REML = TRUE) {
  d <- repair_example
  m <- repair_lmer_refitter(score ~ day + (day | person), REML = REML)(d$repeated)
  list(model = m, repeated = d$repeated, people = d$people,
    fit = repair(re("day") ~ practice, m, d$people, "person", backend = "r"))
}

test_that("ML and REML likelihood scores and frequency derivatives match independent differences", {
  for (reml in c(FALSE, TRUE)) {
    ex <- if_example(reml); ctx <- repair_if_context(ex$fit); psi <- ctx$psi
    val <- repair_if_likelihood(ctx, psi, case_scores = TRUE)
    expect_equal(val$loglik, as.numeric(logLik(ex$model)), tolerance = 1e-8)
    grad <- vapply(seq_along(psi), function(j) {
      h <- 1e-5 * max(1, abs(psi[j])); a <- b <- psi; a[j] <- a[j] + h; b[j] <- b[j] - h
      (repair_if_likelihood(ctx, a)$loglik - repair_if_likelihood(ctx, b)$loglik) / (2 * h)
    }, numeric(1))
    expect_equal(unname(val$score), grad, tolerance = 2e-5)
    for (i in c(1L, 12L)) {
      h <- 1e-4; f <- g <- rep(1, nrow(ex$people)); f[i] <- f[i] + h; g[i] <- g[i] - h
      numeric <- (repair_if_likelihood(ctx, psi, f)$score - repair_if_likelihood(ctx, psi, g)$score) / (2 * h)
      expect_equal(unname(val$case_scores[i, ]), numeric, tolerance = 1e-6)
    }
  }
})

test_that("fractional-frequency M-MOM reproduces explicit whole-person duplication", {
  set.seed(724)
  n <- 40; x <- matrix(rnorm(n * 2), n, 2); y <- 1 + x %*% c(.4, -.2) + rnorm(n)
  omega <- array(0, c(2, 2, n)); for (i in 1:n) omega[, , i] <- diag(c(.15, .12)) * (1 + i / n)
  data <- list(outcome = drop(y), predictors = x, err_x = omega, err_y = rep(.1, n), err_xy = matrix(0, n, 2))
  for (correction in c("minimum", "none", "custom")) {
    control <- list(correction = correction, alpha = if (correction == "custom") 7 else NULL)
    inputs <- mmom_prepare(data, control, control, "r")
    draw <- c(2:40, 2L); f <- tabulate(draw, nbins = n)
    fractional <- repair_if_structural(inputs, f)
    manual <- fit_mmom(list(outcome = data$outcome[draw], predictors = x[draw, ], err_x = omega[, , draw],
      err_y = data$err_y[draw], err_xy = data$err_xy[draw, ]), control, control, backend = "r")
    expect_true(manual@converged)
    expect_equal(unname(fractional$value), unname(manual@coefficients$estimate), tolerance = 1e-8)
    expect_equal(fractional$steps[[2]]@resid_var, manual@diagnostics$resid_var, tolerance = 1e-9)
    expect_equal(fractional$steps[[3]]@correction_final, manual@diagnostics$correction_final, tolerance = 1e-9)
  }
})

test_that("full covariance includes signed cross-stage terms and preserves fitted results", {
  ex <- if_example(); before <- list(coef = coef(ex$fit), vcov = vcov(ex$fit), scores = repair_scores(ex$fit))
  set.seed(421); seed <- .Random.seed
  out <- repair_sandwich(ex$fit)
  expect_true(out$diagnostics$full_available); expect_true(out$diagnostics$stage2_available)
  expect_identical(.Random.seed, seed)
  expect_equal(coef(out), coef(ex$fit)); expect_equal(nobs(out), nrow(ex$people))
  expect_equal(vcov(out), vcov(out, component = "stage2") + vcov(out, component = "measurement") + vcov(out, component = "cross"), tolerance = 1e-12)
  expect_lt(vcov(out, component = "cross")[1, 1], 0)
  expect_equal(vcov(out), crossprod(out$influence$full) / nobs(out)^2, tolerance = 1e-12)
  expect_lt(max(abs(colMeans(out$influence$full))), 1e-7)
  expect_identical(rownames(out$influence$full), ex$people$person)
  table <- repair_sandwich_coefficients(out)
  expect_equal(table$std.error, sqrt(diag(vcov(out))), ignore_attr = TRUE)
  expect_equal(unname(confint(out)), unname(cbind(table$conf.low, table$conf.high)))
  expect_output(print(out), "influence-function")
  path <- tempfile(); on.exit(unlink(path)); saveRDS(out, path)
  expect_equal(vcov(readRDS(path)), vcov(out))
  expect_identical(before, list(coef = coef(ex$fit), vcov = vcov(ex$fit), scores = repair_scores(ex$fit)))
})

test_that("measurement influence agrees with independently reoptimized frequency likelihoods", {
  ex <- if_example(); ctx <- repair_if_context(ex$fit); psi <- ctx$psi; n <- nrow(ex$people)
  measurement <- repair_if_measurement(ctx, 1e-4, .03)
  # Optimization of the frequency likelihood is independent of the implicit
  # derivative solve. Central frequency shifts preserve the N-person population.
  for (person in c(3L, 18L)) {
    epsilon <- .02
    solve_at <- function(delta) {
      f <- rep(1 - delta / n, n); f[person] <- f[person] + delta
      optim(psi, function(p) -repair_if_likelihood(ctx, p, f)$loglik,
        gr = function(p) -repair_if_likelihood(ctx, p, f)$score, method = "BFGS",
        control = list(reltol = 1e-13, maxit = 500))$par
    }
    derivative <- (solve_at(epsilon) - solve_at(-epsilon)) / (2 * epsilon)
    expect_equal(unname(measurement$influence[person, ] / n), unname(derivative), tolerance = 3e-3)
  }
})

test_that("old OLS-outcome sandwich is an exact special-case anchor", {
  set.seed(104)
  n <- 60; times <- seq(-1, 1, length.out = 8)
  x <- rnorm(n); b0 <- rnorm(n); b1 <- .5 * x + rnorm(n)
  d <- do.call(rbind, lapply(1:n, function(i) data.frame(person = sprintf("p%03d", i), time = times,
    y = .4 * times + b0[i] + b1[i] * times + rnorm(8, sd = .5))))
  people <- data.frame(person = sprintf("p%03d", 1:n), x = x)
  m <- repair_lmer_refitter(y ~ time + (time | person))(d)
  f <- repair(re("time") ~ x, m, people, "person", backend = "r")
  out <- repair_sandwich(f)
  expect_true(out$diagnostics$full_available)
  y <- f@analysis_data@outcome; X <- cbind(1, x); b <- solve(crossprod(X), crossprod(X, y))
  r <- drop(y - X %*% b); inverse <- solve(crossprod(X))
  manual <- inverse %*% crossprod(r * X) %*% inverse
  expect_equal(vcov(out, component = "stage2")[2, 2], manual[2, 2], tolerance = 1e-7)
  expect_equal(vcov(out)[2, 2], manual[2, 2], tolerance = 1e-7)
  expect_lt(max(abs(out$influence$measurement[, 2])), 1e-5)
})

test_that("weighted offset nested designs reconstruct scores and the native likelihood", {
  ex <- nested_example(); d <- ex$case$data[match(rownames(model.frame(ex$model)), rownames(ex$case$data)), ]
  m <- repair_lmer_refitter(y ~ x + offset(off) + (x | person) + (1 | person:day), weights = "prior", REML = FALSE)(d)
  f <- repair(ex$formula, list(growth = m), ex$case$people, "person", backend = "r")
  ctx <- repair_if_context(f)
  expect_equal(repair_if_likelihood(ctx, ctx$psi)$loglik, as.numeric(logLik(m)), tolerance = 1e-7)
  data <- repair_if_inputs(f, ctx, ctx$psi)
  for (field in c("outcome", "predictors", "err_x", "err_y", "err_xy"))
    expect_equal(data[[field]], S7::prop(f@analysis_data, field), tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("unavailable inference is explicit and never falls back silently", {
  ex <- if_example(); ctx <- repair_if_context(ex$fit); ctx$singular <- TRUE
  out <- repair_sandwich_compute(ex$fit, 1e-4, .03, context = ctx)
  expect_false(out$diagnostics$full_available); expect_true(out$diagnostics$stage2_available)
  expect_error(vcov(out), "singular")
  expect_true(all(is.finite(vcov(out, component = "stage2")))); expect_null(summary(out)$coefficients)
  expect_error(repair_sandwich(ex$fit, step = 0), "positive finite")
  expect_error(confint(out, component = "stage2", parm = "missing"), "existing coefficients")
  expect_error(repair_if_difference(function(j, d) list(value = abs(d), signature = if (d > 0) "right" else "left"),
    list(value = 0, signature = "left"), 1L, .001, .03), "unstable|branch")
})

test_that("full noisy-predictor influences match complete frequency perturbations", {
  ex <- if_example()
  f <- repair(practice ~ re("(Intercept)") + re("day"), ex$model, ex$people, "person", backend = "r")
  ctx <- repair_if_context(f); out <- repair_sandwich(f); n <- nrow(ex$people)
  expect_true(out$diagnostics$full_available)
  expect_lt(max(abs(out$measurement$structural_sensitivity[-1, 1:2])), 1e-5)
  epsilon <- .02; person <- 11L
  solve_all <- function(delta) {
    frequency <- rep(1 - delta / n, n); frequency[person] <- frequency[person] + delta
    opt <- optim(ctx$psi, function(p) -repair_if_likelihood(ctx, p, frequency)$loglik,
      gr = function(p) -repair_if_likelihood(ctx, p, frequency)$score, method = "BFGS",
      control = list(reltol = 1e-14, maxit = 1000))
    expect_equal(opt$convergence, 0L)
    repair_if_structural(repair_if_inputs(f, ctx, opt$par), frequency)$value
  }
  numerical <- (solve_all(epsilon) - solve_all(-epsilon)) / (2 * epsilon)
  expect_equal(unname(out$influence$full[person, ] / n), unname(numerical), tolerance = 5e-4)
})

test_that("sample-dependent measurement recipes are not silently held fixed", {
  d <- repair_example
  m <- lme4::lmer(score ~ scale(day) + (scale(day) | person), d$repeated)
  f <- repair(re("scale(day)") ~ practice, m, d$people, "person")
  expect_error(repair_sandwich(f), "pointwise")
  m <- lme4::lmer(score ~ day + (day | person), d$repeated, weights = (score * 0 + 1) / sd(score)^2)
  f <- repair(re("day") ~ practice, m, d$people, "person")
  expect_error(repair_sandwich(f), "pointwise")
})
