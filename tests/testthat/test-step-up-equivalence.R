test_that("optimized searches match the frozen estimator exactly on independent cases", {
  reference <- readRDS(test_path("fixtures", "step-up-reference.rds"))
  for (case in reference$cases) for (backend in c("r", "rcpp")) {
    control <- list(correction = "step_up", step_up = case$settings)
    expected <- step_up_legacy$fit_mmom(case$inputs, control, control, backend)
    actual <- fit_mmom(case$inputs, control, control, backend)
    expect_step_up_unchanged(actual, expected)
  }
})

test_that("cache invalidation preserves mixed policies, grids and refinement boundaries", {
  input <- readRDS(test_path("fixtures", "step-up-reference.rds"))$cases$strict_correlated$inputs
  policies <- list(
    list(correction = "step_up", step_up = list(min_relative_eigenvalue = .65)),
    list(correction = "step_up", step_up = list(min_relative_eigenvalue = .65,
      candidate_alphas = c(80, 4, 20, 20), max_refinements = 0L)),
    list(correction = "step_up", step_up = list(min_relative_eigenvalue = .65,
      coarse_grid_size = 3L, search_tolerance = .0001, max_refinements = 12L)),
    list(correction = "minimum"), list(correction = "none"), list(correction = "custom", alpha = 8))
  for (backend in c("r", "rcpp")) for (preliminary in policies) for (final in policies) {
    if (preliminary$correction != "step_up" && final$correction != "step_up") next
    expect_step_up_unchanged(fit_mmom(input, preliminary, final, backend),
      step_up_legacy$fit_mmom(input, preliminary, final, backend))
  }
  base <- step_up_legacy$fit_mmom(input, policies[[1L]], policies[[1L]])
  boundary <- base@diagnostics$step_up$preliminary$selected$relative_min_eigenvalue
  for (floor in boundary + c(-1, 0, 1) * 10 * .Machine$double.eps) {
    control <- list(correction = "step_up", step_up = list(min_relative_eigenvalue = floor))
    expect_step_up_unchanged(fit_mmom(input, control, control),
      step_up_legacy$fit_mmom(input, control, control))
  }
})

test_that("zero errors, determinant branches and numerical failures retain their behavior", {
  input <- readRDS(test_path("fixtures", "step-up-reference.rds"))$cases$weak_predictor$inputs
  control <- list(correction = "step_up", step_up = list(min_relative_eigenvalue = .6))
  lambda <- step_up_legacy$fit_mmom(input)@diagnostics$lambda_prelim
  n <- length(input$outcome)
  for (target in c(.99, 1, 1 + 1/n - 1e-12, 1 + 1/n + 1e-12, 1.05)) {
    candidate <- input; candidate$err_x <- candidate$err_x * lambda / target
    expect_step_up_unchanged(fit_mmom(candidate, control, control),
      step_up_legacy$fit_mmom(candidate, control, control))
  }
  input$err_x[] <- 0
  expect_step_up_unchanged(fit_mmom(input, control, control), step_up_legacy$fit_mmom(input, control, control))
  input$outcome[] <- 0
  expect_step_up_unchanged(fit_mmom(input, control, control), step_up_legacy$fit_mmom(input, control, control))
  input$predictors <- cbind(input$predictors, input$predictors)
  colnames(input$predictors) <- c("x1", "x2")
  input$err_x <- array(0, c(2, 2, n)); input$err_xy <- matrix(0, n, 2)
  expect_step_up_unchanged(fit_mmom(input, control, control), step_up_legacy$fit_mmom(input, control, control))
})

test_that("additional predictor dimensions and heterogeneous errors preserve complete fits", {
  set.seed(43822)
  for (q in c(1L, 3L, 4L)) {
    n <- 48L
    latent <- matrix(rnorm(n * q), n, q)
    covariance <- crossprod(matrix(rnorm((q + 1L)^2, sd = .1), q + 1L))
    h <- rep(c(.3, 1, 2), length.out = n)
    error <- (matrix(rnorm(n * (q + 1L)), n) %*% chol(covariance)) * sqrt(h)
    input <- list(outcome = .4 + drop(latent %*% seq(.2, .6, length.out = q)) +
      rnorm(n, sd = .5) + error[, 1], predictors = latent + error[, -1, drop = FALSE],
      err_x = array(0, c(q, q, n)), err_y = h * covariance[1, 1],
      err_xy = outer(h, covariance[1, -1]))
    for (i in seq_len(n)) input$err_x[, , i] <- h[i] * covariance[-1, -1, drop = FALSE]
    for (backend in c("r", "rcpp")) for (floor in c(.05, .65)) {
      control <- list(correction = "step_up", step_up = list(min_relative_eigenvalue = floor))
      expect_step_up_unchanged(fit_mmom(input, control, control, backend),
        step_up_legacy$fit_mmom(input, control, control, backend))
    }
  }
})

test_that("complete candidate checks are retained while invariant work is cached", {
  reference <- readRDS(test_path("fixtures", "step-up-reference.rds"))
  prepared <- mmom_fit_prepared; root <- mmom_smallest_det_root
  first <- mmom_step1; second <- mmom_step2; final_moments <- mmom_final_moments
  counts <- integer(5L)
  local_mocked_bindings(
    mmom_fit_prepared = function(...) { counts[1] <<- counts[1] + 1L; prepared(...) },
    mmom_smallest_det_root = function(...) { counts[2] <<- counts[2] + 1L; root(...) },
    mmom_step1 = function(...) { counts[3] <<- counts[3] + 1L; first(...) },
    mmom_step2 = function(...) { counts[4] <<- counts[4] + 1L; second(...) },
    mmom_final_moments = function(...) { counts[5] <<- counts[5] + 1L; final_moments(...) })
  control <- list(correction = "step_up")
  fit <- fit_mmom(reference$cases$both_sides_correlated$inputs, control, control)
  expect_true(fit@converged)
  expect_identical(counts, c(1L, 2L, 1L, 1L, 1L))
  expect_identical(vapply(fit@diagnostics$step_up, function(x) nrow(x$candidates), integer(1)),
    c(preliminary = 1L, final = 1L))
  counts[] <- 0L
  final <- list(correction = "step_up", step_up = reference$cases$strict_correlated$settings)
  fit <- fit_mmom(reference$cases$strict_correlated$inputs, list(correction = "minimum"), final)
  expect_gt(nrow(fit@diagnostics$step_up$final$candidates), 1L)
  expect_identical(counts[2:5], c(2L, 1L, 1L, 1L))
  # All evaluated candidates still go through the complete fit, including SEs.
  expect_identical(counts[1], as.integer(nrow(fit@diagnostics$step_up$final$candidates)))
})
