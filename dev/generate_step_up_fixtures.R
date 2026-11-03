# Maintenance only: freeze numerical results from the separate research
# implementation. No research source is copied into the MIT package, and normal
# tests neither load that repository nor regenerate these expectations.
# Rscript dev/generate_step_up_fixtures.R /path/to/blup_unshrinking
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else "../blup_unshrinking")
reference <- new.env(parent = globalenv())
source_path <- file.path(root, "R/stage2_estimators.R")
sys.source(source_path, reference)
baseline <- readRDS("tests/testthat/fixtures/repair-reference.rds")
cases <- lapply(baseline$cases, function(x) list(inputs = x$inputs, settings = list()))
# A stricter threshold exercises bisection in both stages; an impossible guard
# exercises the documented return-to-minimum fallback.
cases$strict_correlated <- list(inputs = baseline$cases$both_sides_correlated$inputs,
  settings = list(min_relative_eigenvalue = 0.65))
cases$impossible_guard <- list(inputs = baseline$cases$both_sides_correlated$inputs,
  settings = list(min_relative_eigenvalue = 2))
cases <- lapply(cases, function(case) {
  input <- case$inputs; q <- ncol(input$predictors)
  dat <- data.frame(y = input$outcome, x1 = input$predictors[, q],
    s22 = input$err_x[q, q, ], syy = input$err_y, c1 = input$err_xy[, q])
  call <- list(stage2_df = dat, outcome = "y", predictor_u0 = NULL,
    predictor_u1 = "x1", meas11 = NULL, meas12 = NULL, meas22 = "s22",
    outcome_meas_var = "syy", predictor_outcome_meas_cov_u1 = "c1",
    skip_internal_scaling = TRUE)
  if (q == 2L) {
    call$stage2_df$x0 <- input$predictors[, 1L]
    call$stage2_df$s11 <- input$err_x[1L, 1L, ]
    call$stage2_df$s12 <- input$err_x[1L, 2L, ]
    call$stage2_df$c0 <- input$err_xy[, 1L]
    call$predictor_u0 <- "x0"; call$meas11 <- "s11"; call$meas12 <- "s12"
    call$predictor_outcome_meas_cov_u0 <- "c0"
  }
  if (length(case$settings)) {
    call$min_sx1_star_relative_eigen <- case$settings$min_relative_eigenvalue
    call$min_sx_star_relative_eigen <- case$settings$min_relative_eigenvalue
  }
  result <- do.call(reference$fit_fuller_dual_alpha_stepdown, call)
  stopifnot(result$status_code == 0L)
  # Independently confirm that the selected constants yield the same fixed-alpha
  # core result. This does not use repair or its search implementation.
  call$min_sx1_star_relative_eigen <- call$min_sx_star_relative_eigen <- NULL
  call$alpha_step1 <- result$fuller_alpha_step1_used
  call$alpha_step3 <- result$fuller_alpha_step3_used
  fixed <- do.call(reference$fit_fuller_dual_core, call)
  stopifnot(isTRUE(all.equal(result$estimate, fixed$estimate, tolerance = 1e-12)),
    isTRUE(all.equal(result$se, fixed$se, tolerance = 1e-12)))
  case$expected <- as.data.frame(result[, c("estimate", "se", "ci_low", "ci_high",
    "fuller_alpha_step1_used", "fuller_alpha_step3_used", "fuller_lambda1", "fuller_lambda2",
    "fuller_sigma2", "fuller_correction1", "fuller_correction_c",
    "fuller_sx1_star_relative_min_eigen", "fuller_sx_star_relative_min_eigen")])
  case
})
saveRDS(list(provenance = list(
  source = "R/stage2_estimators.R: fit_fuller_dual_alpha_stepdown; skip_internal_scaling=TRUE",
  source_md5 = unname(tools::md5sum(source_path)),
  research_commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
  R = R.version.string,
  note = "Frozen per-dataset algorithm agreement, not coverage validation."), cases = cases),
  "tests/testthat/fixtures/step-up-reference.rds", version = 2)
cat("Saved", length(cases), "independent step-up cases.\n")
