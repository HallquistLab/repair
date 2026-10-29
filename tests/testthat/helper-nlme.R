nlme_example_data <- function(seed = 203L) {
  set.seed(seed)
  dat <- do.call(rbind, lapply(seq_len(35), function(i) {
    time <- 0:(5 + i %% 5)
    z <- time / 4
    b <- rnorm(2, sd = c(1, .6))
    phase <- factor(time %% 2)
    error <- drop(t(chol(.4^abs(outer(time, time, "-")))) %*% rnorm(length(time))) *
      ifelse(phase == "1", 1.4, .6)
    data.frame(id = sprintf("p%02d", i), time = time, z = z, phase = phase,
      scale_var = 1 + z, y = 10 + .7 * z + b[1] + b[2] * z + error)
  }))
  rownames(dat) <- paste0("observation", seq_len(nrow(dat)))
  dat$id <- factor(dat$id, levels = c(rev(unique(dat$id)), "unused"))
  dat <- dat[sample(nrow(dat)), ]
  dat$y[which(dat$time > 0)[1L]] <- NA_real_
  dat
}

# Test oracle: construct residual covariance directly from the fitted parameters
# and observed covariates; do not use repair's extractor or nlme::getVarCov(R).
independent_nlme_gls <- function(model, data, structure = "iid") {
  ids <- rownames(nlme::ranef(model))
  data <- data[match(rownames(model$groups), rownames(data)), , drop = FALSE]
  residual <- data$y - drop(model.matrix(~ z + phase, data,
    contrasts.arg = model$contrasts) %*% nlme::fixef(model))
  estimates <- matrix(NA_real_, length(ids), 2, dimnames = list(ids, c("(Intercept)", "z")))
  errors <- array(NA_real_, c(2, 2, length(ids)), dimnames = list(colnames(estimates), colnames(estimates), ids))
  residual_blocks <- list()
  for (id in ids) {
    rows <- which(as.character(data$id) == id)
    d <- data[rows, ]
    sd <- rep(model$sigma, nrow(d))
    if (grepl("ident", structure)) {
      multipliers <- coef(model$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)
      sd <- sd * multipliers[as.character(d$phase)]
    }
    if (grepl("power", structure)) {
      power <- as.numeric(coef(model$modelStruct$varStruct, unconstrained = FALSE))
      if (structure == "fitted_power") {
        # The fitted varPower weights, even their stored covariate/power, can
        # differ at the last convergence update. nlme::getVarCov uses these
        # weights; this oracle does not use the residual-SD attribute.
        fitted_weights <- nlme::varWeights(model$modelStruct$varStruct)
        sd <- model$sigma / fitted_weights[order(order(model$groups[[1L]]))][rows]
      } else {
        sd <- sd * abs(d$scale_var)^power
      }
    }
    corr <- diag(nrow(d))
    if (grepl("ar1|car1", structure)) {
      rho <- as.numeric(coef(model$modelStruct$corStruct, unconstrained = FALSE))
      times <- if (grepl("car1", structure)) d$z else d$time
      corr <- rho^abs(outer(times, times, "-"))
    }
    r <- outer(sd, sd) * corr
    z <- cbind(1, d$z)
    inverse <- solve(crossprod(z, solve(r, z)))
    estimates[id, ] <- drop(inverse %*% crossprod(z, solve(r, residual[rows])))
    errors[, , id] <- inverse
    residual_blocks[[id]] <- r
  }
  list(estimates = estimates, errors = errors, residual_blocks = residual_blocks)
}
