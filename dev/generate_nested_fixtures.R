# From the package root:
# Rscript dev/generate_nested_fixtures.R /path/to/blup_unshrinking
# Do not load repair: build GLS/posterior references directly from fitted lme4
# quantities, then cross-check structural inference against the research code.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else "../blup_unshrinking")
ref <- new.env(parent = globalenv())
paths <- c("R/stage2_estimators.R", "R/fuller_mlm.R",
  "vig_hallquist_2026/analysis/fuller_algebra_audit_simulation.R")
for (path in paths) sys.source(file.path(root, path), ref)
close <- function(x, y) stopifnot(isTRUE(all.equal(x, y,
  tolerance = 1e-8, check.attributes = FALSE)))

make_case <- function(seed, kind) {
  set.seed(seed)
  n <- 48L
  people <- data.frame(person = sprintf("P%03d", seq_len(n)), practice = rnorm(n))
  b0 <- rnorm(n, sd = 1.1)
  b1 <- .25 * b0 + .4 * people$practice + rnorm(n, sd = .45)
  data <- do.call(rbind, lapply(seq_len(n), function(i) {
    ndays <- 3L + i %% 3L
    day_effect <- matrix(rnorm(ndays * 2), ndays) %*% chol(matrix(c(.49, .06, .06, .16), 2))
    block_effect <- rnorm(2, sd = .45)
    do.call(rbind, lapply(seq_len(ndays), function(day) {
      size <- 4L + (i + day) %% 3L
      within <- seq(-1, 1, length.out = size)
      x <- (day - 2.5) / 2 + within
      prior <- runif(size, .6, 1.8)
      off <- .1 * sin(x)
      block <- 1L + day %% 2L
      nuisance <- day_effect[day, 1] + if (kind == "nested_slope")
        day_effect[day, 2] * within + block_effect[block] else 0
      data.frame(person = people$person[i], day = day, block = block,
        within = within, x = x, prior = prior, off = off,
        y = 2 + .8 * x + off + b0[i] + b1[i] * x + nuisance + rnorm(size, sd = .6 / sqrt(prior)))
    }))
  }))
  data$person <- factor(data$person, levels = c(rev(people$person), "unused"))
  data$y[7L] <- NA_real_
  data <- data[sample(nrow(data)), ]
  formula <- if (kind == "nested_slope")
    "y ~ x + offset(off) + (1 + x | person) + (1 + within | person:day) + (1 | person:block)" else
    "y ~ x + offset(off) + (1 + x | person) + (1 | person:day)"
  model <- lme4::lmer(as.formula(formula), data, weights = prior, na.action = na.exclude,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  stopifnot(!lme4::isSingular(model), is.null(model@optinfo$conv$lme4$messages))
  mf <- model.frame(model)
  residual <- model.response(mf) - model.offset(mf) -
    drop(lme4::getME(model, "X") %*% lme4::fixef(model))
  vc <- lme4::VarCorr(model)
  g <- as.matrix(vc$person)
  scores <- matrix(0, n, 2, dimnames = list(people$person, c("(Intercept)", "x")))
  error <- posterior <- array(0, c(2, 2, n), dimnames = list(colnames(scores), colnames(scores), people$person))
  covariance <- vector("list", n); names(covariance) <- people$person
  counts <- integer(n)
  for (i in seq_len(n)) {
    rows <- which(as.character(mf$person) == people$person[i]); counts[i] <- length(rows)
    z <- cbind(1, mf$x[rows])
    r <- diag(sigma(model)^2 / model.weights(mf)[rows])
    same_day <- outer(mf$day[rows], mf$day[rows], "==")
    if (kind == "nested_slope") {
      h <- cbind(1, mf$within[rows])
      r <- r + same_day * (h %*% as.matrix(vc[["person:day"]]) %*% t(h)) +
        as.numeric(vc[["person:block"]]) * outer(mf$block[rows], mf$block[rows], "==")
    } else r <- r + as.numeric(vc[["person:day"]]) * same_day
    covariance[[i]] <- r
    error[, , i] <- solve(crossprod(z, solve(r, z)))
    scores[i, ] <- error[, , i] %*% crossprod(z, solve(r, residual[rows]))
    posterior[, , i] <- g - g %*% t(z) %*% solve(r + z %*% g %*% t(z), z %*% g)
  }
  # Joint random outcome/predictor errors plus an observed error-free predictor.
  y <- scores[, "x"]; x <- cbind(scores[, "(Intercept)"], people$practice)
  ex <- array(0, c(2, 2, n)); ex[1, 1, ] <- error[1, 1, ]
  ey <- error[2, 2, ]; exy <- cbind(error[2, 1, ], 0)
  full <- ref$fit_fuller_mlm(y, x, ex, ey, exy, backend = "r")
  production <- ref$fit_fuller_dual_core(
    stage2_df = data.frame(y = y, x0 = x[, 1], x1 = x[, 2], s11 = ex[1, 1, ],
      s22 = 0, s12 = 0, syy = ey, c0 = exy[, 1], c1 = 0),
    outcome = "y", predictor_u0 = "x0", predictor_u1 = "x1", meas11 = "s11",
    meas22 = "s22", meas12 = "s12", outcome_meas_var = "syy",
    predictor_outcome_meas_cov_u0 = "c0", predictor_outcome_meas_cov_u1 = "c1",
    skip_internal_scaling = TRUE)
  omega <- array(0, c(n, 3, 3)); omega[, -1, -1] <- aperm(ex, c(3, 1, 2))
  audit <- ref$fit_fuller_audit(y, cbind(1, x), omega, ey, cbind(0, exy),
    preliminary = "current_cross", se_bread = "current_final")
  fields <- c("estimate", "se", "ci_low", "ci_high")
  stopifnot(full$converged, production$status_code == 0L, audit$status == "ok")
  close(unlist(tail(full$coefficients, 1)[fields]), unlist(production[fields]))
  close(unlist(production[fields]), unlist(as.data.frame(audit)[fields]))
  list(kind = kind, data = data, people = people, measurement_formula = formula,
    formula = 're("x", "growth", group = "person") ~ re("(Intercept)", "growth", group = "person") + practice',
    fitted_beta = lme4::fixef(model), target_g = g, n_observations = counts,
    covariance = covariance, scores = scores, error = error, posterior = posterior,
    input = list(outcome = y, predictors = x, err_x = ex, err_y = ey, err_xy = exy),
    coefficients = full$coefficients, vcov = full$vcov)
}
fixture <- list(provenance = list(
  source_commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
  source_md5 = stats::setNames(unname(tools::md5sum(file.path(root, paths))), paths),
  generator_md5 = unname(tools::md5sum("dev/generate_nested_fixtures.R")),
  lme4_version = as.character(packageVersion("lme4")), R_version = R.version.string,
  scoring = "Independent within-person GLS with all fitted nuisance contributions; target G excluded",
  posterior = "Independent G - GZ'(R_target+ZGZ')^-1 ZG",
  inference = "Conditional modified M-MOM; production Fuller and independent algebra audit agree",
  note = "Synthetic per-dataset implementation references, not a coverage study"),
  cases = list(nested_intercept = make_case(202610301, "nested_intercept"),
    nested_slope = make_case(202610302, "nested_slope")))
saveRDS(fixture, "inst/extdata/nested-workflow-reference.rds", version = 2, compress = "xz")
cat("Saved", length(fixture$cases), "independent nested-nuisance references.\n")
