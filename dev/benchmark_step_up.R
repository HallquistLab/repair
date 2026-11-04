# Run against an installed optimized candidate from the repository root:
# R_LIBS=/path/to/checked/library Rscript dev/benchmark_step_up.R [repetitions]
# Compares the exact frozen 60c46ec estimator with the optimized estimator in
# the same process, dependency versions, backend and inputs. This is a structural
# microbenchmark, not a claim about end-to-end bootstrap speed or coverage.
library(repair)
ns <- asNamespace("repair")
stopifnot(exists("mmom_preliminary_moments", ns, inherits = FALSE))
legacy <- new.env(parent = ns)
sys.source("tests/testthat/fixtures/step-up-60c46ec.R", legacy)
reference <- readRDS("tests/testthat/fixtures/step-up-reference.rds")
args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args)) as.integer(args[1L]) else 5L
stopifnot(is.finite(repetitions), repetitions >= 1L)

profile_calls <- function(input, control, backend, use_legacy) {
  fits <- roots <- 0L
  root <- get("mmom_smallest_det_root", ns)
  prepared <- get("mmom_fit_prepared", ns)
  previous <- legacy$mmom_fit_prepared
  on.exit(legacy$mmom_fit_prepared <- previous, add = TRUE)
  legacy$mmom_fit_prepared <- function(...) { fits <<- fits + 1L; previous(...) }
  testthat::local_mocked_bindings(
    mmom_fit_prepared = function(...) { fits <<- fits + 1L; prepared(...) },
    mmom_smallest_det_root = function(...) { roots <<- roots + 1L; root(...) }, .package = "repair")
  fun <- if (use_legacy) legacy$fit_mmom else repair::fit_mmom
  fun(input, control, control, backend)
  c(fits = fits, roots = roots)
}

cases <- list(minimum_passes = reference$cases$both_sides_correlated,
  both_stages_refine = reference$cases$strict_correlated,
  no_passing_candidate = reference$cases$impossible_guard)
# Repeated synthetic inputs increase the computational workload, not the amount
# of independent information. No statistical interpretation is attached to them.
large <- reference$cases$strict_correlated
index <- rep(seq_along(large$inputs$outcome), 10L)
large$inputs <- with(large$inputs, list(outcome = outcome[index],
  predictors = predictors[index, , drop = FALSE], err_x = err_x[, , index, drop = FALSE],
  err_y = err_y[index], err_xy = err_xy[index, , drop = FALSE]))
cases$repeated_inputs <- large
results <- list()
for (name in names(cases)) for (backend in c("r", "rcpp")) {
  case <- cases[[name]]; input <- case$inputs
  control <- list(correction = "step_up", step_up = case$settings)
  old <- legacy$fit_mmom(input, control, control, backend)
  new <- repair::fit_mmom(input, control, control, backend)
  stopifnot(identical(old@coefficients, new@coefficients), identical(old@vcov, new@vcov),
    identical(old@diagnostics$alpha_prelim, new@diagnostics$alpha_prelim),
    identical(old@diagnostics$alpha_final, new@diagnostics$alpha_final))
  before <- profile_calls(input, control, backend, TRUE)
  after <- profile_calls(input, control, backend, FALSE)
  times <- matrix(NA_real_, repetitions, 2L)
  for (i in seq_len(repetitions)) for (which in if (i %% 2) c(1L, 2L) else c(2L, 1L)) {
    fun <- if (which == 1L) legacy$fit_mmom else repair::fit_mmom
    times[i, which] <- system.time(fun(input, control, control, backend), gcFirst = FALSE)[["elapsed"]]
  }
  medians <- apply(times, 2L, median)
  results[[length(results) + 1L]] <- data.frame(case = name, backend = backend,
    n = length(input$outcome), repetitions = repetitions,
    fits_before = unname(before["fits"]), fits_after = unname(after["fits"]),
    roots_before = unname(before["roots"]), roots_after = unname(after["roots"]),
    median_seconds_before = medians[1], median_seconds_after = medians[2],
    speed_ratio = medians[1] / medians[2], exact_estimates_SEs_alphas = TRUE)
}
result <- do.call(rbind, results)
write.csv(result, "dev/step-up-efficiency.csv", row.names = FALSE)
writeLines(trimws(c("Baseline: 60c46ece17071f7401717a2eefe1a169e5a598d6",
  "Same-process structural microbenchmark; timings are environment-specific.",
  capture.output(sessionInfo())), which = "right"), "dev/step-up-efficiency-session.txt")
print(result, row.names = FALSE)
