# From the repair package root:
# Rscript dev/generate_nlme_fixtures.R /path/to/blup_unshrinking
# Never load repair here: frozen expectations come from the paper's production
# prediction extraction/estimation functions, cross-checked against independent algebra.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else "../blup_unshrinking")
ref <- new.env(parent = globalenv())
paths <- c("R/blup_helpers.R", "R/stage2_estimators.R", "R/fuller_mlm.R",
  "vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R")
for (path in paths) sys.source(file.path(root, path), ref)
close <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-8, check.attributes = FALSE)))

make_case <- function(seed, kind) {
  set.seed(seed)
  n <- 60L
  b <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, .2, .2, .64), 2))
  people <- data.frame(id = sprintf("p%02d", seq_len(n)), observed = rnorm(n),
    distal = .4 - .3 * b[, 1] + .8 * b[, 2] + rnorm(n, sd = .7))
  # Outcome case has an observed predictor related to the latent random slope.
  if (kind == "outcome_unbalanced") b[, 2] <- .5 * people$observed + b[, 2]
  data <- do.call(rbind, lapply(seq_len(n), function(i) {
    size <- if (kind == "predictor_balanced") 8L else 5L + i %% 7L
    time <- seq_len(size) - 1L
    z <- time / 4
    phase <- factor(time %% 2)
    r <- if (kind == "both_ar1_ident") .35^abs(outer(time, time, "-")) else diag(size)
    sd <- if (kind == "both_ar1_ident") ifelse(phase == "1", .9, .5) else rep(.6, size)
    error <- drop(t(chol(r * outer(sd, sd))) %*% rnorm(size))
    data.frame(id = people$id[i], time = time, z = z, phase = phase,
      y = 2 + .7 * z + b[i, 1] + b[i, 2] * z + error)
  }))
  data$id <- factor(data$id)
  model <- nlme::lme(y ~ z, random = ~ z | id, data = data,
    correlation = if (kind == "both_ar1_ident") nlme::corAR1(form = ~ time | id) else NULL,
    weights = if (kind == "both_ar1_ident") nlme::varIdent(form = ~ 1 | phase) else NULL,
    control = nlme::lmeControl(msMaxIter = 300))
  # The production prediction extraction helper receives the complete, sorted observations
  # in these fixtures. Missingness/permutation are separate extraction tests.
  predictions <- ref$get_closed_form_corrected_scores(model, data, "id", "y", "z")
  predictions <- as.data.frame(predictions)[match(people$id, predictions$id), ]
  if (kind == "predictor_balanced") {
    formula <- 'distal ~ re("(Intercept)", "growth") + re("z", "growth")'
    x <- cbind(predictions$corrected_intercept_full, predictions$corrected_slope_full)
    y <- people$distal; ey <- rep(0, n); exy <- matrix(0, n, 2)
    ex <- array(0, c(2, 2, n))
    ex[1, 1, ] <- predictions$gls_var11; ex[2, 2, ] <- predictions$gls_var22
    ex[1, 2, ] <- ex[2, 1, ] <- predictions$gls_var12
    reporting_scale <- .8 # generating population slope SD, as in VH Study 2
  } else if (kind == "outcome_unbalanced") {
    formula <- 're("z", "growth") ~ observed'
    x <- cbind(people$observed); y <- predictions$corrected_slope_full
    ex <- array(0, c(1, 1, n)); ey <- predictions$gls_var22; exy <- matrix(0, n, 1)
    reporting_scale <- 1
  } else {
    formula <- 're("z", "growth") ~ re("(Intercept)", "growth")'
    x <- cbind(predictions$corrected_intercept_full); y <- predictions$corrected_slope_full
    ex <- array(predictions$gls_var11, c(1, 1, n)); ey <- predictions$gls_var22
    exy <- cbind(predictions$gls_var12); reporting_scale <- 1
  }
  q <- ncol(x)
  colnames(x) <- paste0("x", seq_len(q))
  full <- ref$fit_fuller_mlm(y, x, ex, ey, exy, backend = "r")
  stage2 <- data.frame(y = y, x1 = x[, q], s22 = ex[q, q, ], syy = ey, c1 = exy[, q])
  arguments <- list(stage2_df = stage2, outcome = "y", predictor_u1 = "x1",
    meas22 = "s22", outcome_meas_var = "syy", predictor_outcome_meas_cov_u1 = "c1",
    skip_internal_scaling = TRUE)
  if (q == 2L) {
    arguments$stage2_df$x0 <- x[, 1]; arguments$stage2_df$s11 <- ex[1, 1, ]
    arguments$stage2_df$s12 <- ex[1, 2, ]; arguments$stage2_df$c0 <- exy[, 1]
    arguments$predictor_u0 <- "x0"; arguments$meas11 <- "s11"; arguments$meas12 <- "s12"
    arguments$predictor_outcome_meas_cov_u0 <- "c0"
  }
  production <- do.call(ref$fit_fuller_dual_core, arguments)
  omega <- array(0, c(n, q + 1L, q + 1L)); omega[, -1, -1] <- aperm(ex, c(3, 1, 2))
  audit <- ref$fit_fuller_audit(y, cbind(1, x), omega, ey, cbind(0, exy),
    preliminary = "current_cross", se_bread = "current_final")
  fields <- c("estimate", "se", "ci_low", "ci_high")
  stopifnot(full$converged, production$status_code == 0L, audit$status == "ok")
  close(unlist(tail(full$coefficients, 1)[fields]), unlist(production[fields]))
  close(unlist(production[fields]), unlist(as.data.frame(audit)[fields]))
  list(data = data, people = people, kind = kind, formula = formula,
    fitted_beta = nlme::fixef(model), predictions = predictions,
    input = list(outcome = y, predictors = x, err_x = ex, err_y = ey, err_xy = exy),
    coefficients = full$coefficients, vcov = full$vcov,
    focal = as.data.frame(production[fields]), reporting_scale = reporting_scale,
    reported_focal = as.data.frame(production[fields]) * reporting_scale)
}
fixture <- list(provenance = list(
  source_commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
  source_md5 = stats::setNames(unname(tools::md5sum(file.path(root, paths))), paths),
  nlme_version = as.character(packageVersion("nlme")), R_version = R.version.string,
  preliminary_moment = "modified", variance_bread = "modified", skip_internal_scaling = TRUE,
  note = "Deterministic synthetic end-to-end reproduction cases; not published Monte Carlo table values."),
  cases = list(predictor_balanced = make_case(202610291, "predictor_balanced"),
    outcome_unbalanced = make_case(202610292, "outcome_unbalanced"),
    both_ar1_ident = make_case(202610293, "both_ar1_ident")))
saveRDS(fixture, "inst/extdata/nlme-workflow-reference.rds", version = 2)
cat("Saved", length(fixture$cases), "independently cross-checked end-to-end fixtures.\n")
