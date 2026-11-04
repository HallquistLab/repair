#!/usr/bin/env Rscript
# Repeated-sampling screen for explicit cross-model covariance assumptions.
# Run against an installed PR #24 package; does not modify the estimator.
options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
settings <- list(reps = 500L, workers = 4L, seed = 20261104L,
  output = "check/covariance-screen", cells = "all", audit_reps = 5L)
for (arg in args) {
  pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
  if (length(pair) != 2L || !pair[1L] %in% names(settings)) stop("Unknown argument: ", arg)
  settings[[pair[1L]]] <- pair[2L]
}
for (name in c("reps", "workers", "seed", "audit_reps")) {
  settings[[name]] <- as.integer(settings[[name]])
  stopifnot(length(settings[[name]]) == 1L, is.finite(settings[[name]]), settings[[name]] >= 0L)
}
stopifnot(settings$reps >= 2L, settings$workers >= 1L)
suppressPackageStartupMessages(library(repair))
stopifnot(requireNamespace("lme4", quietly = TRUE))
dir.create(settings$output, recursive = TRUE, showWarnings = FALSE)
settings$output <- normalizePath(settings$output)
design <- expand.grid(n = c(50L, 200L), m = c(6L, 12L),
  layout = c("balanced", "unbalanced"), rho = c(-.5, 0, .5), beta = c(0, .5))
design$cell <- seq_len(nrow(design))
cells <- if (settings$cells == "all") design$cell else as.integer(strsplit(settings$cells, ",", fixed = TRUE)[[1L]])
stopifnot(length(cells) > 0L, !anyNA(cells), !anyDuplicated(cells), all(cells %in% design$cell))
write.csv(design, file.path(settings$output, "design.csv"), row.names = FALSE)
RNGkind("L'Ecuyer-CMRG"); set.seed(settings$seed)
streams <- vector("list", nrow(design)); streams[[1L]] <- .Random.seed
for (i in 2:nrow(design)) streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])

screen_capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = value, warnings = paste(unique(warnings), collapse = "\n"))
}

screen_fit <- function(input, policy) {
  control <- list(correction = policy)
  fit_mmom(input, control, control, backend = "rcpp")
}

screen_record <- function(fit, cell, rep, mode, assumed, policy, warning = "", mlm_error = "") {
  failed <- inherits(fit, "error") || is.null(fit)
  ok <- !failed && isTRUE(fit@converged) &&
    all(is.finite(fit@coefficients$estimate)) && all(is.finite(fit@coefficients$se)) &&
    all(fit@coefficients$se > 0)
  d <- if (failed) list() else fit@diagnostics
  searches <- d$step_up
  get_alpha <- function(stage) {
    value <- d[[if (stage == "preliminary") "alpha_prelim" else "alpha_final"]]
    if (is.null(value)) value <- d[[paste0(stage, "_moment_control")]]$alpha
    if (is.null(value)) NA_real_ else value
  }
  data.frame(cell = cell, rep = rep, mode = mode, assumed_rho = assumed, policy = policy,
    ok = ok, estimate = if (ok) tail(fit@coefficients$estimate, 1L) else NA_real_,
    se = if (ok) tail(fit@coefficients$se, 1L) else NA_real_,
    lower = if (ok) tail(fit@coefficients$ci_low, 1L) else NA_real_,
    upper = if (ok) tail(fit@coefficients$ci_high, 1L) else NA_real_,
    alpha_preliminary = get_alpha("preliminary"), alpha_final = get_alpha("final"),
    fallback_preliminary = if (is.null(searches$preliminary)) NA else searches$preliminary$fallback,
    fallback_final = if (is.null(searches$final)) NA else searches$final$fallback,
    message = if (nzchar(mlm_error)) mlm_error else if (failed) {
      if (inherits(fit, "error")) conditionMessage(fit) else "Unavailable measurement fit"
    } else fit@message,
    warning = warning)
}

screen_equal <- function(x, y, label) {
  if (!isTRUE(all.equal(x, y, tolerance = 1e-9, check.attributes = FALSE))) {
    stop("Public-interface audit failed: ", label)
  }
}

screen_cell <- function(cell) {
  cfg <- design[cell, ]; n <- cfg$n
  counts <- if (cfg$layout == "balanced") rep(cfg$m, n) else
    as.integer(cfg$m * rep(c(.5, 1, 1.5), length.out = n))
  times <- lapply(counts, function(m) seq(-1, 1, length.out = m))
  h <- vapply(times, function(t) 1 / sum(t^2), numeric(1))
  ids <- sprintf("p%03d", seq_len(n))
  index <- rep(seq_len(n), counts); time <- unlist(times, use.names = FALSE)
  rows <- split(seq_along(index), index)
  dat <- data.frame(id = factor(ids[index], levels = ids), time = time,
    x = 0, y = 0, occasion = unlist(lapply(counts, seq_len)))
  rownames(dat) <- sprintf("obs%05d", seq_len(nrow(dat)))
  keys <- setNames(as.character(dat$occasion), rownames(dat))
  keys <- list(x = keys, y = keys)
  equation <- re(time, "y") ~ re(time, "x")
  people <- data.frame(id = ids)
  sigma_x <- 1.5; sigma_y <- 1.5; sigma_eta <- .8
  stream <- streams[[cell]]
  results <- vector("list", settings$reps)
  errors <- vector("list", settings$reps)
  audits <- vector("list", settings$reps)
  started <- proc.time()[[3L]]
  for (r in seq_len(settings$reps)) {
    assign(".Random.seed", stream, envir = .GlobalEnv)
    stream <- parallel::nextRNGSubStream(stream)
    bx <- rnorm(n)
    by <- cfg$beta * bx + rnorm(n, sd = sigma_eta)
    b0x <- .2 * bx + sqrt(.96) * rnorm(n)
    b0y <- .2 * by / sqrt(cfg$beta^2 + sigma_eta^2) + sqrt(.96) * rnorm(n)
    ex <- rnorm(nrow(dat)); ey <- cfg$rho * ex + sqrt(1 - cfg$rho^2) * rnorm(nrow(dat))
    dat$x <- 2 + .4 * time + b0x[index] + bx[index] * time + sigma_x * ex
    dat$y <- 1 - .2 * time + b0y[index] + by[index] * time + sigma_y * ey
    # Independent closed-form OLS extraction, with known fixed slopes removed.
    score_x <- vapply(rows, function(i) sum(time[i] * dat$x[i]) / sum(time[i]^2) - .4, numeric(1))
    score_y <- vapply(rows, function(i) sum(time[i] * dat$y[i]) / sum(time[i]^2) + .2, numeric(1))
    known <- list(outcome = score_y, predictors = matrix(score_x, ncol = 1L,
      dimnames = list(ids, 're("time", model = "x")')),
      err_x = array(sigma_x^2 * h, c(1L, 1L, n)), err_y = sigma_y^2 * h,
      err_xy = matrix(0, n, 1L))
    # Marginal covariance parameters are estimated; cross-process rho never is.
    fitted <- screen_capture({
      control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
      models <- list(x = lme4::lmer(x ~ time + (time | id), dat, control = control),
        y = lme4::lmer(y ~ time + (time | id), dat, control = control))
      # For this centered iid design, GLS is OLS and h=1/sum(t^2).
      # Evaluate the closed form efficiently; audit against the public scorer,
      # propagation helper and complete sensitivity calls in every design cell.
      sx <- sigma(models$x); sy <- sigma(models$y)
      input <- list(outcome = score_y - .2 - lme4::fixef(models$y)["time"],
        predictors = matrix(score_x + .4 - lme4::fixef(models$x)["time"], ncol = 1L,
          dimnames = list(ids, 're("time", model = "x")')),
        err_x = array(sx^2 * h, c(1L, 1L, n)), err_y = sy^2 * h,
        err_xy = matrix(0, n, 1L))
      unit_cross <- sx * sy * h
      if (r <= settings$audit_reps) {
        analysis <- repair_data(equation, models, people, "id", cross_model_error = "independent")
        block <- repair_cross_covariance(equation, models, keys, rho = .5, ids = ids)
        screen_equal(input$outcome, analysis@outcome, "closed-form outcome")
        screen_equal(input$predictors, analysis@predictors, "closed-form predictors")
        screen_equal(input$err_x, analysis@err_x, "closed-form predictor covariance")
        screen_equal(input$err_y, analysis@err_y, "closed-form outcome covariance")
        screen_equal(unit_cross, as.numeric(block[[1L]]$covariance[1L, 1L, ids]) / .5,
          "closed-form cross-covariance")
      }
      list(models = models, input = input, unit_cross = unit_cross,
        singular = any(vapply(models, lme4::isSingular, logical(1))))
    })
    modes <- list(known = known)
    good <- !inherits(fitted$value, "error")
    if (!good && grepl("Public-interface audit failed", conditionMessage(fitted$value), fixed = TRUE)) {
      stop(fitted$value)
    }
    if (good) {
      modes$fitted <- fitted$value$input
    }
    rep_results <- list(); k <- 0L
    for (mode in c("known", "fitted")) {
      for (assumed in c(-.5, 0, .5)) {
        for (policy in c("minimum", "step_up")) {
          k <- k + 1L
          if (mode == "fitted" && !good) {
            fit <- NULL; warning <- fitted$warnings
          } else {
            input <- modes[[mode]]
            input$err_xy <- matrix(assumed * if (mode == "known") sigma_x * sigma_y * h else
              fitted$value$unit_cross, n, 1L)
            captured <- screen_capture(screen_fit(input, policy))
            fit <- captured$value; warning <- captured$warnings
            if (mode == "fitted" && nzchar(fitted$warnings)) warning <- paste(fitted$warnings, warning, sep = "\n")
          }
          rep_results[[k]] <- screen_record(fit, cell, r, mode, assumed, policy, warning,
            if (mode == "fitted" && !good) conditionMessage(fitted$value) else "")
          if (r <= settings$audit_reps && !inherits(fit, "error") && !is.null(fit)) {
            # Both computation engines agree on the exact same score inputs.
            control <- list(correction = policy)
            ref <- fit_mmom(input, control, control, backend = "r")
            screen_equal(fit@coefficients, ref@coefficients, "R/C++ coefficient tables")
            screen_equal(fit@vcov, ref@vcov, "R/C++ covariance")
          }
        }
      }
    }
    results[[r]] <- do.call(rbind, rep_results)
    err_x <- score_x - bx; err_y <- score_y - by
    # Sample covariance removes common fixed-slope shifts; its expectation is
    # the average within-person observation-noise covariance, even when m varies.
    oracle <- lm(by ~ bx)
    oracle_estimate <- unname(coef(oracle)[2L]); oracle_se <- sqrt(vcov(oracle)[2L, 2L])
    errors[[r]] <- data.frame(cell = cell, rep = r,
      empirical_x = var(err_x), empirical_y = var(err_y), empirical_xy = cov(err_x, err_y),
      known_x = mean(known$err_x), known_y = mean(known$err_y),
      known_xy = mean(cfg$rho * sigma_x * sigma_y * h),
      fitted_x = if (good) mean(modes$fitted$err_x) else NA_real_,
      fitted_y = if (good) mean(modes$fitted$err_y) else NA_real_,
      fitted_xy = if (good) mean(cfg$rho * fitted$value$unit_cross) else NA_real_,
      singular = if (good) fitted$value$singular else NA,
      measurement_failure = !good, measurement_warning = nzchar(fitted$warnings),
      oracle_estimate = oracle_estimate, oracle_se = oracle_se,
      oracle_covered = abs(oracle_estimate - cfg$beta) <= qnorm(.975) * oracle_se,
      oracle_reject = abs(oracle_estimate) > qnorm(.975) * oracle_se)
    if (good && r <= settings$audit_reps) {
      scenarios <- setNames(lapply(c(-.5, 0, .5), function(rho) {
        repair_cross_covariance(equation, fitted$value$models, keys, rho = rho, ids = ids)
      }), c("negative", "zero", "positive"))
      for (policy in c("minimum", "step_up")) {
        public <- repair_sensitivity(equation, fitted$value$models, people, "id", scenarios,
          regularization = policy, backend = "rcpp")
        for (j in seq_along(scenarios)) {
          manual <- results[[r]][results[[r]]$mode == "fitted" &
            results[[r]]$assumed_rho == c(-.5, 0, .5)[j] & results[[r]]$policy == policy, ]
          actual <- public$fits[[j]]@struc_fit
          stopifnot(manual$ok == isTRUE(actual@converged))
          if (manual$ok) {
            screen_equal(manual$estimate, tail(actual@coefficients$estimate, 1L), "public estimate")
            screen_equal(manual$se, tail(actual@coefficients$se, 1L), "public SE")
            screen_equal(manual$lower, tail(actual@coefficients$ci_low, 1L), "public interval")
            screen_equal(manual$alpha_preliminary, public$scenarios$alpha_preliminary[j], "preliminary alpha")
            screen_equal(manual$alpha_final, public$scenarios$alpha_final[j], "final alpha")
          }
        }
      }
      audits[[r]] <- data.frame(cell = cell, rep = r, public_scenarios_checked = 6L,
        backend_fits_checked = 12L)
    }
  }
  value <- list(config = cfg, settings = settings, results = do.call(rbind, results),
    errors = do.call(rbind, errors), audits = do.call(rbind, audits),
    elapsed = proc.time()[[3L]] - started)
  path <- file.path(settings$output, sprintf("cell-%02d.rds", cell))
  saveRDS(value, paste0(path, ".tmp"), version = 2)
  if (!file.rename(paste0(path, ".tmp"), path)) stop("Could not save checkpoint")
  data.frame(cell = cell, reps = settings$reps, elapsed = value$elapsed)
}

manifest <- list(settings = settings, design = design, streams = streams,
  R = R.version.string, packages = vapply(c("repair", "lme4", "S7", "Rcpp"),
    function(x) as.character(packageVersion(x)), character(1)),
  session = capture.output(sessionInfo()), started = format(Sys.time(), tz = "UTC"),
  source_md5 = tools::md5sum("dev/covariance-screen/simulate.R"),
  package_commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_))
saveRDS(manifest, file.path(settings$output, "manifest.rds"), version = 2)
writeLines(manifest$session, file.path(settings$output, "session-info.txt"))
cat("Starting", length(cells), "cells x", settings$reps, "replicates with", settings$workers, "workers\n")
if (settings$workers == 1L) {
  timing <- do.call(rbind, lapply(cells, screen_cell))
} else {
  cluster <- parallel::makePSOCKcluster(settings$workers,
    outfile = file.path(settings$output, "workers.log"))
  parallel::clusterExport(cluster, c("settings", "design", "streams", "screen_capture", "screen_fit",
    "screen_record", "screen_equal", "screen_cell"), envir = globalenv())
  parallel::clusterEvalQ(cluster, { library(repair); NULL })
  timing <- tryCatch(do.call(rbind, parallel::parLapplyLB(cluster, cells, screen_cell)),
    finally = parallel::stopCluster(cluster))
}
write.csv(timing, file.path(settings$output, "timing.csv"), row.names = FALSE)
cat("Completed", nrow(timing), "cells; worker seconds", sum(timing$elapsed), "\n")
