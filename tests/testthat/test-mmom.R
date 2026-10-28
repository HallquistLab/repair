test_that("M-MOM fits the zero-error OLS limit with both backends", {
  set.seed(1)
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
  y <- as.vector(1 + x %*% c(0.5, -0.25) + rnorm(n))
  data_list <- lapply(seq_len(n), function(row) list(
    outcome = y[[row]], predictors = x[row, ], err_x = matrix(0, 2, 2), err_y = 0, err_xy = matrix(0, 2, 1)
  ))

  fit_r <- fit_mmom(data_list, backend = "r")
  fit_cpp <- fit_mmom(data_list, backend = "rcpp")
  expect_true(fit_r@converged)
  expect_true(fit_cpp@converged)
  expect_equal(
    unname(fit_r@coefficients$estimate),
    unname(coef(lm(y ~ x))),
    tolerance = 1e-7
  )
  expect_equal(
    fit_r@coefficients$estimate,
    fit_cpp@coefficients$estimate,
    tolerance = 1e-7
  )
  expect_equal(
    fit_r@vcov,
    fit_cpp@vcov,
    tolerance = 1e-10
  )
  expect_equal(fit_r@steps$step1@step, "prelim_unweighted_coef")
  expect_equal(fit_r@steps$step2@step, "resid_var_and_weights")
  expect_equal(fit_r@steps$step3@step, "final_weighted_coef")
})

test_that("M-MOM fits an errors-in-variables model with both backends", {
  set.seed(4)
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
  y <- as.vector(1 + x %*% c(0.5, -0.25) + rnorm(n))
  err_cov <- lapply(seq_len(n), function(row) {
    A <- matrix(rnorm(9, sd = 0.1), 3, 3)
    crossprod(A)
  })
  err_x <- lapply(err_cov, function(covariance) {
    covariance[1:2, 1:2, drop = FALSE]
  })
  err_y <- lapply(err_cov, function(covariance) {
    covariance[3, 3]
  })
  err_xy <- lapply(err_cov, function(covariance) {
    covariance[1:2, 3, drop = FALSE]
  })
  data_list <- lapply(seq_len(n), function(row) {
    list(
      outcome = y[[row]],
      predictors = x[row, ],
      err_x = err_x[[row]],
      err_y = err_y[[row]],
      err_xy = err_xy[[row]]
    )
  })

  fit_r <- fit_mmom(data_list, backend = "r")
  fit_cpp <- fit_mmom(data_list, backend = "rcpp")
  expect_true(fit_r@converged)
  expect_true(fit_cpp@converged)
  expect_equal(
    fit_r@coefficients$estimate,
    fit_cpp@coefficients$estimate,
    tolerance = 1e-7
  )
  expect_equal(
    fit_r@vcov,
    fit_cpp@vcov,
    tolerance = 1e-10
  )
  expect_equal(fit_r@steps$step1@step, "prelim_unweighted_coef")
  expect_equal(fit_r@steps$step2@step, "resid_var_and_weights")
  expect_equal(fit_r@steps$step3@step, "final_weighted_coef")
})

test_that("rank-deficient moments return specific solver diagnostics", {
  set.seed(3)
  predictor <- rnorm(20)
  data_list <- lapply(seq_along(predictor), function(row) list(
    y = predictor[[row]] + rnorm(1),
    predictors = c(predictor[[row]], predictor[[row]]),
    err_x = matrix(0, 2, 2),
    err_y = 0,
    err_xy = matrix(0, 2, 1)
  ))

  fit <- fit_mmom(data_list, backend = "r")
  expect_false(fit@converged)
  expect_match(fit@message, "prelim_coef_solve")
  expect_equal(fit@diagnostics$step, "prelim_coef_solve")
  expect_equal(fit@diagnostics$solver, "qr")
  expect_lt(fit@diagnostics$rank, fit@diagnostics$dimension)
})

test_that("class_repair_mod is an S7 object", {
  expect_true(inherits(class_repair_mod(list(), list()), "S7_object"))
})

test_that("ranef_like uses the closed-form within-cluster score", {
  skip_if_not_installed("lme4")
  set.seed(2)
  data <- do.call(rbind, lapply(seq_len(10), function(cluster) {
    z <- seq(-1, 1, length.out = 5)
    data.frame(id = factor(cluster), z = z, y = 0.4 + 0.7 * z +
      rnorm(1, sd = 0.5) + rnorm(length(z), sd = 0.2))
  }))
  fit <- lme4::lmer(y ~ z + (1 + z | id), data = data, REML = FALSE)
  scored <- ranef_like(fit)
  fixed_residual <- data$y - as.numeric(
    lme4::getME(fit, "X") %*% lme4::fixef(fit)
  )
  expected <- lapply(split(fixed_residual, data$id), function(residual) {
    z <- seq(-1, 1, length.out = 5)
    unname(coef(lm(residual ~ z)))
  })
  expect_equal(
    unname(scored$estimates[, 1]),
    unname(vapply(expected, `[[`, numeric(1), 1)),
    tolerance = 1e-7
  )
  expect_equal(
    unname(scored$estimates[, 2]),
    unname(vapply(expected, `[[`, numeric(1), 2)),
    tolerance = 1e-7
  )
  expect_equal(dim(scored$err_cov), c(2, 2, 10))
  expect_equal(scored$diagnostics$solver, "chol")
  expect_true(all(vapply(scored$diagnostics$clusters, function(item) {
    identical(item$solver, "chol")
  }, logical(1))))
})