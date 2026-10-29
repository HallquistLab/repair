# Rscript dev/generate_multimodel_fixtures.R /path/to/blup_unshrinking
# Do not load repair: expectations use the paper's production functions and
# independent algebra. This script generates synthetic reproduction datasets.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else "../blup_unshrinking")
ref <- new.env(parent = globalenv())
paths <- c("R/blup_helpers.R", "R/stage2_estimators.R", "R/fuller_mlm.R",
  "vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R")
for (path in paths) sys.source(file.path(root, path), ref)
close <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-8, check.attributes = FALSE)))
make_case <- function(seed, rho) {
  set.seed(seed)
  n <- 70L
  b <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, .2, .2, .64), 2))
  q0 <- .4 * b[, 1] + rnorm(n, sd = .8)
  q1 <- -.25 * b[, 1] + .65 * b[, 2] + rnorm(n, sd = .6)
  people <- data.frame(person = sprintf("p%02d", seq_len(n)))
  observations <- lapply(seq_len(n), function(i) {
    t <- seq_len(6L + i %% 5L) - 1L
    z <- t / 4; w <- .5 + t / 5
    e1 <- rnorm(length(t)); e2 <- rho * e1 + sqrt(1 - rho^2) * rnorm(length(t))
    list(a = data.frame(id = people$person[i], z = z,
      y = 2 + .7 * z + b[i, 1] + b[i, 2] * z + .6 * e1),
      b = data.frame(subject = people$person[i], w = w,
      q = 1 - .3 * w + q0[i] + q1[i] * w + .7 * e2))
  })
  a <- do.call(rbind, lapply(observations, `[[`, "a"))
  b <- do.call(rbind, lapply(observations, `[[`, "b"))
  a$id <- factor(a$id); b$subject <- factor(b$subject)
  fit_a <- nlme::lme(y ~ z, random = ~ z | id, data = a)
  fit_b <- nlme::lme(q ~ w, random = ~ w | subject, data = b)
  sa <- as.data.frame(ref$get_closed_form_corrected_scores(fit_a, a, "id", "y", "z"))
  sb <- as.data.frame(ref$get_closed_form_corrected_scores(fit_b, b, "subject", "q", "w"))
  sa <- sa[match(people$person, sa$id), ]; sb <- sb[match(people$person, sb$id), ]
  y <- sb$corrected_slope_full
  x <- cbind(x0 = sa$corrected_intercept_full, x1 = sa$corrected_slope_full)
  ex <- array(0, c(2, 2, n)); ex[1, 1, ] <- sa$gls_var11; ex[2, 2, ] <- sa$gls_var22
  ex[1, 2, ] <- ex[2, 1, ] <- sa$gls_var12
  ey <- sb$gls_var22; exy <- matrix(0, n, 2)
  # Known illustrative residual correlation, conditional on fitted marginal
  # residual scales. L_a R_ab L_b' propagates observation errors into scores.
  # Cross-model covariance is NOT the correlation between latent slopes.
  for (i in seq_len(n)) {
    za <- cbind(1, a$z[a$id == people$person[i]])
    zb <- cbind(1, b$w[b$subject == people$person[i]])
    la <- solve(crossprod(za), t(za)); lb <- solve(crossprod(zb), t(zb))
    cross <- rho * sigma(fit_a) * sigma(fit_b) * tcrossprod(la, lb)
    exy[i, ] <- cross[, 2]
  }
  cross <- array(t(exy), c(2, 1, n),
    dimnames = list(c("(Intercept)", "z"), "w", people$person))
  full <- ref$fit_fuller_mlm(y, x, ex, ey, exy, backend = "r")
  stage2 <- data.frame(y = y, x0 = x[, 1], x1 = x[, 2], s11 = ex[1, 1, ],
    s12 = ex[1, 2, ], s22 = ex[2, 2, ], syy = ey, c0 = exy[, 1], c1 = exy[, 2])
  production <- ref$fit_fuller_dual_core(stage2, outcome = "y", predictor_u0 = "x0",
    predictor_u1 = "x1", meas11 = "s11", meas12 = "s12", meas22 = "s22",
    outcome_meas_var = "syy", predictor_outcome_meas_cov_u0 = "c0",
    predictor_outcome_meas_cov_u1 = "c1", skip_internal_scaling = TRUE)
  omega <- array(0, c(n, 3, 3)); omega[, -1, -1] <- aperm(ex, c(3, 1, 2))
  audit <- ref$fit_fuller_audit(y, cbind(1, x), omega, ey, cbind(0, exy),
    preliminary = "current_cross", se_bread = "current_final")
  fields <- c("estimate", "se", "ci_low", "ci_high")
  stopifnot(full$converged, production$status_code == 0L, audit$status == "ok")
  close(unlist(tail(full$coefficients, 1)[fields]), unlist(production[fields]))
  close(unlist(production[fields]), unlist(as.data.frame(audit)[fields]))
  # Dual-process reporting convention: population predictor SD / outcome SD.
  reporting_scale <- .8 / sqrt((-.25)^2 + .65^2 * .64 + 2 * (-.25) * .65 * .2 + .6^2)
  list(a = a, b = b, people = people, rho = rho,
    formula = 're("w", "response") ~ re("(Intercept)", "growth") + re("z", "growth")',
    cross_model_error = if (rho == 0) "independent" else
      list(list(models = c("growth", "response"), covariance = cross)),
    fitted_beta = list(a = nlme::fixef(fit_a), b = nlme::fixef(fit_b)),
    scores = list(a = sa, b = sb),
    input = list(outcome = y, predictors = x, err_x = ex, err_y = ey, err_xy = exy),
    coefficients = full$coefficients, vcov = full$vcov, reporting_scale = reporting_scale,
    reported_focal = as.data.frame(production[fields]) * reporting_scale)
}
fixture <- list(provenance = list(
  source_commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
  source_md5 = stats::setNames(unname(tools::md5sum(file.path(root, paths))), paths),
  nlme_version = as.character(packageVersion("nlme")), R_version = R.version.string,
  preliminary_moment = "modified", variance_bread = "modified", skip_internal_scaling = TRUE,
  note = paste("Synthetic dual-process cases; not published Monte Carlo table values.",
    "Supplied cross-error case uses known residual correlation with fitted marginal scales.")),
  cases = list(independent = make_case(202610294, 0), supplied = make_case(202610295, .25)))
saveRDS(fixture, "inst/extdata/multimodel-workflow-reference.rds", version = 2)
cat("Saved", length(fixture$cases), "independently cross-checked dual-process fixtures.\n")
