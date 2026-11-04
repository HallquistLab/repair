#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: summarize.R checkpoint_directory [output_directory] [--allow-partial]")
root <- normalizePath(args[1L])
dest <- if (length(args) >= 2L) args[2L] else "dev/covariance-screen/results"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
manifest <- readRDS(file.path(root, "manifest.rds"))
design <- manifest$design
paths <- list.files(root, "^cell-[0-9]+\\.rds$", full.names = TRUE)
parts <- lapply(paths, readRDS)
stopifnot(length(parts) > 0L)
cells <- vapply(parts, function(x) x$config$cell, integer(1))
stopifnot(!anyDuplicated(cells))
if (!"--allow-partial" %in% args && !setequal(cells, design$cell)) stop("Missing design cells")
for (x in parts) {
  stopifnot(identical(x$settings$reps, manifest$settings$reps),
    identical(x$settings$seed, manifest$settings$seed),
    identical(x$settings$audit_reps, manifest$settings$audit_reps),
    nrow(x$results) == manifest$settings$reps * 12L,
    nrow(x$errors) == manifest$settings$reps,
    length(unique(x$results$rep)) == manifest$settings$reps)
}
results <- do.call(rbind, lapply(parts, `[[`, "results"))
errors <- do.call(rbind, lapply(parts, `[[`, "errors"))
audits <- do.call(rbind, lapply(parts, `[[`, "audits"))
results <- merge(results, design, by = "cell", sort = FALSE)
results$scenario_role <- ifelse(results$assumed_rho == results$rho, "correct",
  ifelse(results$assumed_rho == 0, "independence", "misspecified_nonzero"))
results$covered <- results$ok & results$lower <= results$beta & results$upper >= results$beta
results$rejected <- results$ok & (results$lower > 0 | results$upper < 0)

summarize_estimates <- function(x) {
  y <- x[x$ok, ]; n <- nrow(x); b <- nrow(y); truth <- unique(x$beta)
  stopifnot(length(truth) == 1L)
  error <- y$estimate - truth
  sd_emp <- if (b > 1L) sd(y$estimate) else NA_real_
  rmse <- if (b) sqrt(mean(error^2)) else NA_real_
  rate <- if (b) mean(y$covered) else NA_real_
  rejection <- if (b) mean(y$rejected) else NA_real_
  mn <- function(x) if (length(x) && any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
  wilson <- function(p, n) {
    if (!n) return(c(NA_real_, NA_real_))
    z <- qnorm(.975); center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(center - half, center + half)
  }
  coverage_bounds <- wilson(rate, b); rejection_bounds <- wilson(rejection, b)
  data.frame(n_attempted = n, n_usable = b, n_failed = n - b,
    failure_rate = (n - b) / n, bias = mn(error), bias_mcse = sd_emp / sqrt(b),
    rmse = rmse, rmse_mcse = if (b > 1L && rmse > 0) sd(error^2) / (2 * rmse * sqrt(b)) else NA_real_,
    empirical_sd = sd_emp, mean_se = mn(y$se), rms_se = sqrt(mn(y$se^2)),
    mean_se_over_sd = mn(y$se) / sd_emp, rms_se_over_sd = sqrt(mn(y$se^2)) / sd_emp,
    coverage_usable = rate, coverage_mcse = sqrt(rate * (1 - rate) / b),
    coverage_mc_lower = coverage_bounds[1L], coverage_mc_upper = coverage_bounds[2L],
    coverage_all_attempts = mean(x$covered), reject_usable = rejection,
    reject_mcse = sqrt(rejection * (1 - rejection) / b),
    reject_mc_lower = rejection_bounds[1L], reject_mc_upper = rejection_bounds[2L],
    type1_error = if (truth == 0) rejection else NA_real_,
    power = if (truth != 0) rejection else NA_real_,
    median_abs_error = if (b) median(abs(error)) else NA_real_,
    q95_abs_error = if (b) unname(quantile(abs(error), .95)) else NA_real_,
    max_abs_error = if (b) max(abs(error)) else NA_real_,
    mean_alpha_preliminary = mn(x$alpha_preliminary), mean_alpha_final = mn(x$alpha_final),
    max_alpha_preliminary = if (any(is.finite(x$alpha_preliminary))) max(x$alpha_preliminary, na.rm = TRUE) else NA_real_,
    max_alpha_final = if (any(is.finite(x$alpha_final))) max(x$alpha_final, na.rm = TRUE) else NA_real_,
    alpha_changed_preliminary = mn(x$alpha_preliminary > 3 + 1e-12),
    alpha_changed_final = mn(x$alpha_final > 3 + 1e-12),
    fallback_preliminary = mn(x$fallback_preliminary), fallback_final = mn(x$fallback_final),
    warning_rate = mean(nzchar(trimws(x$warning))))
}

group_summary <- function(data, keys, fun) {
  groups <- split(seq_len(nrow(data)), interaction(data[keys], drop = TRUE, lex.order = TRUE))
  out <- lapply(groups, function(i) cbind(data[i[1L], keys, drop = FALSE], fun(data[i, , drop = FALSE])))
  value <- do.call(rbind, out); rownames(value) <- NULL; value
}
summary <- group_summary(results, c("cell", "n", "m", "layout", "rho", "beta", "mode", "assumed_rho", "scenario_role", "policy"), summarize_estimates)
pooled <- group_summary(results, c("mode", "policy", "scenario_role", "rho", "assumed_rho", "beta"), summarize_estimates)
write.csv(summary, file.path(dest, "performance.csv"), row.names = FALSE)
write.csv(pooled, file.path(dest, "pooled-performance.csv"), row.names = FALSE)

covariance <- do.call(rbind, lapply(cells, function(cell) {
  e <- errors[errors$cell == cell, ]
  do.call(rbind, lapply(c("known", "fitted"), function(mode) {
    do.call(rbind, lapply(c("x", "y", "xy"), function(component) {
      prediction <- e[[paste0(mode, "_", component)]]
      empirical <- e[[paste0("empirical_", component)]]
      keep <- is.finite(prediction) & is.finite(empirical)
      diff <- empirical[keep] - prediction[keep]
      mcse <- sd(diff) / sqrt(length(diff))
      cbind(design[design$cell == cell, ], data.frame(mode = mode, component = component,
        n_reps = sum(keep), empirical_mean = mean(empirical[keep]), propagated_mean = mean(prediction[keep]),
        difference = mean(diff), difference_mcse = mcse, difference_z = if (mcse > 0) mean(diff) / mcse else NA_real_))
    }))
  }))
}))
write.csv(covariance, file.path(dest, "prediction-covariance.csv"), row.names = FALSE)

measurement <- group_summary(errors, "cell", function(x) data.frame(
  reps = nrow(x), measurement_failures = sum(x$measurement_failure),
  singular = sum(x$singular, na.rm = TRUE), warnings = sum(x$measurement_warning),
  oracle_bias = mean(x$oracle_estimate - design$beta[design$cell == x$cell[1L]]),
  oracle_empirical_sd = sd(x$oracle_estimate), oracle_mean_se = mean(x$oracle_se),
  oracle_coverage = mean(x$oracle_covered), oracle_coverage_mcse = sqrt(mean(x$oracle_covered) * (1 - mean(x$oracle_covered)) / nrow(x)),
  oracle_reject = mean(x$oracle_reject)))
measurement <- merge(design, measurement, by = "cell")
write.csv(measurement, file.path(dest, "measurement-and-oracle.csv"), row.names = FALSE)
write.csv(audits, file.path(dest, "interface-audits.csv"), row.names = FALSE)
write.csv(design[design$cell %in% cells, ], file.path(dest, "design.csv"), row.names = FALSE)
problems <- results[!results$ok | nzchar(trimws(results$warning)), ]
if (nrow(problems)) {
  counts <- group_summary(problems, c("cell", "mode", "assumed_rho", "policy", "ok", "message", "warning"),
    function(x) data.frame(count = nrow(x)))
} else counts <- data.frame(cell = integer(), mode = character(), assumed_rho = numeric(),
  policy = character(), ok = logical(), message = character(), warning = character(), count = integer())
write.csv(counts, file.path(dest, "failures-and-warnings.csv"), row.names = FALSE)
writeLines(manifest$session, file.path(dest, "session-info.txt"))
manifest_lines <- c("REPAIR cross-model covariance screen",
  paste("Started UTC:", manifest$started),
  paste("Package source commit:", manifest$package_commit),
  paste("Simulation script MD5:", unname(manifest$source_md5)),
  paste("Seed:", manifest$settings$seed),
  "RNG: L'Ecuyer-CMRG; one stream per cell and substream per replicate",
  paste("Replicates per condition:", manifest$settings$reps),
  paste("Workers:", manifest$settings$workers),
  paste("Public audits per cell:", manifest$settings$audit_reps),
  paste("R:", manifest$R), paste(names(manifest$packages), manifest$packages, sep = ": "),
  paste("Dataset count:", nrow(errors)), paste("Structural-fit count:", nrow(results)),
  paste("Measurement failures:", sum(errors$measurement_failure)),
  paste("Structural failures:", sum(!results$ok)),
  paste("Public scenario audits:", sum(audits$public_scenarios_checked)),
  paste("R/C++ comparisons:", sum(audits$backend_fits_checked)))
writeLines(manifest_lines, file.path(dest, "run-manifest.txt"))
hashes <- data.frame(file = basename(paths), md5 = unname(tools::md5sum(paths)))
write.csv(hashes, file.path(dest, "checkpoint-md5.csv"), row.names = FALSE)
cat("Datasets:", nrow(errors), "structural fits:", nrow(results), "\n")
cat("Measurement failures:", sum(errors$measurement_failure), "structural failures:", sum(!results$ok), "\n")
cat("Public fits checked:", sum(audits$public_scenarios_checked), "R/C++ fits checked:", sum(audits$backend_fits_checked), "\n")
cat("Correct-scenario coverage range:", range(summary$coverage_usable[summary$scenario_role == "correct"], na.rm = TRUE), "\n")
print(pooled[, c("mode", "policy", "scenario_role", "rho", "assumed_rho", "beta", "bias", "empirical_sd", "mean_se", "coverage_usable", "type1_error", "n_failed")], row.names = FALSE)
