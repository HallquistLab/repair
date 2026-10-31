# Independent-person marginal Gaussian likelihood, including contained nuisance
# blocks. Native lme4 theta is log-transformed on its diagonal only.
repair_if_context <- function(object) {
  a <- object@analysis_data
  if (length(object@mlm_fit) != 1L || !inherits(object@mlm_fit[[1L]], "lmerMod"))
    repair_error("Influence inference currently supports one Gaussian lmer model; nlme and multiple-model derivatives are not implemented.")
  model <- object@mlm_fit[[1L]]; target <- unique(a@effects$group)
  if (length(target) != 1L) repair_error("Select one target grouping factor.")
  # Stored X/Z are sufficient only for pointwise transformations. A callback
  # that recomputes scale(), sd(), sample knots, etc. changes those designs under
  # case resampling and needs additional derivatives that are not implemented.
  tryCatch(repair_replay_formula(stats::formula(model), all.vars(stats::formula(model))),
    error = function(e) repair_error("Influence inference requires supported pointwise measurement formulas; precompute transformations or implement their sampling derivatives."))
  for (field in c("weights", "offset")) {
    expression <- stats::getCall(model)[[field]]
    if (is.null(expression)) next
    variables <- all.vars(expression)
    dummy <- as.data.frame(stats::setNames(rep(list(rep(1, 2)), length(variables)), variables))
    value <- tryCatch(repair_replay_eval(expression, dummy), error = function(e) NULL)
    if (is.null(value) || !is.numeric(value) ||
        (!length(variables) && length(value) > 1L && length(unique(value)) > 1L))
      repair_error("Influence inference requires pointwise {.val {field}} columns or expressions; sample-dependent or position-specific vector recipes are not supported.")
  }
  repair_validate_target(model, target, names(object@mlm_fit))
  e <- repair_get_mlm_elements(model)
  if (!setequal(e@group_levels[[target]], a@ids))
    repair_error("Influence inference requires coinciding measurement and analyzed cohorts; fit an explicit restricted reference first.")
  theta <- lme4::getME(model, "theta"); diagonal <- lme4::getME(model, "lower") == 0
  beta <- lme4::fixef(model); sigma <- stats::sigma(model)
  blocks <- lme4::getME(model, "cnms"); ptr <- lme4::getME(model, "Gp")
  Z <- as.matrix(lme4::getME(model, "Z"))
  pairs <- lapply(blocks, function(x) which(lower.tri(matrix(0, length(x), length(x)), diag = TRUE), arr.ind = TRUE))
  sizes <- vapply(pairs, nrow, integer(1)); offsets <- c(0L, cumsum(sizes))
  stopifnot(sum(sizes) == length(theta))
  people <- lapply(a@ids, function(id) {
    rows <- e@group_rows[[target]][[id]]
    basis <- lapply(seq_along(blocks), function(j) {
      k <- length(blocks[[j]]); z <- Z[rows, seq.int(ptr[j] + 1L, ptr[j + 1L]), drop = FALSE]
      lapply(seq_len(nrow(pairs[[j]])), function(t) {
        pair <- pairs[[j]][t, ]; za <- z[, seq.int(pair[1], ncol(z), by = k), drop = FALSE]
        zb <- z[, seq.int(pair[2], ncol(z), by = k), drop = FALSE]
        value <- tcrossprod(za, zb)
        if (pair[1] != pair[2]) value <- value + t(value)
        value
      })
    })
    list(y = e@response[rows] - e@offset[rows], X = e@fixed_design[rows, , drop = FALSE],
      R = diag(1 / e@observation_weights[rows], length(rows)), basis = basis,
      target = get_ranef_design(e, target, rows, id))
  })
  psi <- c(unname(beta), ifelse(diagonal, log(pmax(theta, .Machine$double.xmin)), theta), log(sigma))
  names(psi) <- c(paste0("fixed:", names(beta)), paste0("theta:", seq_along(theta)), "log_sigma")
  list(people = people, psi = psi, beta_size = length(beta), diagonal = diagonal,
    block_sizes = lengths(blocks), pairs = pairs, offsets = offsets,
    nuisance = names(blocks) != target, reml = lme4::isREML(model), ids = a@ids,
    terms = e@random_terms[[target]], model = model,
    singular = lme4::isSingular(model), converged = repair_measurement_converged(model) &&
      (is.null(attr(model, "repair_lmer_optimization")) || isTRUE(attr(model, "repair_lmer_optimization")$usable)))
}

repair_if_parameters <- function(context, psi) {
  p <- context$beta_size; q <- length(context$diagonal)
  raw <- psi[p + seq_len(q)]; theta <- raw; theta[context$diagonal] <- exp(raw[context$diagonal])
  covariance <- derivative <- vector("list", length(context$block_sizes))
  for (j in seq_along(covariance)) {
    k <- context$block_sizes[j]; indices <- context$offsets[j] + seq_len(nrow(context$pairs[[j]]))
    L <- matrix(0, k, k); L[context$pairs[[j]]] <- theta[indices]
    covariance[[j]] <- tcrossprod(L)[context$pairs[[j]]]
    derivative[[j]] <- lapply(seq_along(indices), function(t) {
      E <- matrix(0, k, k); E[context$pairs[[j]][t, , drop = FALSE]] <-
        if (context$diagonal[indices[t]]) theta[indices[t]] else 1
      (tcrossprod(E, L) + tcrossprod(L, E))[context$pairs[[j]]]
    })
  }
  list(beta = psi[seq_len(p)], sigma2 = exp(2 * psi[p + q + 1L]),
    covariance = covariance, derivative = derivative)
}

repair_if_covariance <- function(person, pars, context, nuisance_only = FALSE) {
  V <- person$R
  for (j in seq_along(person$basis)) if (!nuisance_only || context$nuisance[j])
    for (k in seq_along(person$basis[[j]])) V <- V + pars$covariance[[j]][k] * person$basis[[j]][[k]]
  pars$sigma2 * V
}

# Returns exact Gaussian score and case-frequency derivatives of that score.
# REML includes the derivative of log|sum_i f_i X_i' V_i^-1 X_i|; a plain
# cluster ML score at a REML fit would not reproduce the fitted procedure.
repair_if_likelihood <- function(context, psi, frequency = rep(1, length(context$people)), case_scores = FALSE, deletions = FALSE) {
  if (deletions && any(frequency != 1)) repair_error("Deletion scores require unit reference frequencies.")
  pars <- repair_if_parameters(context, psi); p <- context$beta_size
  q <- length(context$diagonal); n <- length(context$people)
  scores <- matrix(0, n, length(psi)); H <- matrix(0, p, p)
  Hs <- vector("list", n); Ds <- vector("list", n); D <- replicate(q + 1L, matrix(0, p, p), simplify = FALSE)
  loglik <- 0
  for (i in seq_len(n)) {
    person <- context$people[[i]]; V <- repair_if_covariance(person, pars, context)
    U <- chol(V); inv <- chol2inv(U); residual <- person$y - drop(person$X %*% pars$beta)
    ir <- drop(inv %*% residual); ix <- inv %*% person$X
    Hi <- crossprod(person$X, ix); H <- H + frequency[i] * Hi
    Hs[[i]] <- Hi
    loglik <- loglik - frequency[i] * (length(residual) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(residual * ir)) / 2
    scores[i, seq_len(p)] <- drop(crossprod(person$X, ir))
    dVs <- vector("list", q + 1L)
    for (j in seq_along(person$basis)) for (t in seq_along(pars$derivative[[j]])) {
      v <- matrix(0, nrow(V), ncol(V))
      for (k in seq_along(person$basis[[j]])) v <- v + pars$derivative[[j]][[t]][k] * person$basis[[j]][[k]]
      dVs[[context$offsets[j] + t]] <- pars$sigma2 * v
    }
    dVs[[q + 1L]] <- 2 * V
    Ds[[i]] <- vector("list", q + 1L)
    for (j in seq_along(dVs)) {
      dV <- dVs[[j]]
      scores[i, p + j] <- (sum(ir * drop(dV %*% ir)) - sum(inv * dV)) / 2
      Dij <- crossprod(ix, dV %*% ix); Ds[[i]][[j]] <- Dij
      D[[j]] <- D[[j]] + frequency[i] * Dij
    }
  }
  score <- colSums(scores * frequency)
  deleted <- if (deletions) matrix(score, n, length(score), byrow = TRUE) - scores else NULL
  if (context$reml) {
    Hinv <- chol2inv(chol(H))
    loglik <- loglik + p * log(2 * pi) / 2 - sum(log(diag(chol(H))))
    for (j in seq_len(q + 1L)) {
      score[p + j] <- score[p + j] + sum(Hinv * D[[j]]) / 2
      if (case_scores) for (i in seq_len(n)) scores[i, p + j] <- scores[i, p + j] +
        (sum(Hinv * Ds[[i]][[j]]) - sum((Hinv %*% Hs[[i]] %*% Hinv) * D[[j]])) / 2
    }
    if (deletions) for (i in seq_len(n)) {
      inverse <- chol2inv(chol(H - Hs[[i]]))
      for (j in seq_len(q + 1L)) deleted[i, p + j] <- deleted[i, p + j] +
        sum(inverse * (D[[j]] - Ds[[i]][[j]])) / 2
    }
  }
  list(loglik = loglik, score = score, case_scores = scores, deletion_scores = deleted)
}

repair_if_inputs <- function(object, context, psi) {
  a <- object@analysis_data; pars <- repair_if_parameters(context, psi)
  h <- matrix(NA_real_, length(context$people), length(context$terms))
  omega <- array(NA_real_, c(ncol(h), ncol(h), nrow(h)))
  for (i in seq_len(nrow(h))) {
    person <- context$people[[i]]; R <- repair_if_covariance(person, pars, context, nuisance_only = TRUE)
    irz <- solve(R, person$target); W <- solve(crossprod(person$target, irz))
    h[i, ] <- drop(W %*% crossprod(irz, person$y - drop(person$X %*% pars$beta)))
    omega[, , i] <- W
  }
  outcome <- a@outcome; predictors <- a@predictors
  err_x <- a@err_x; err_y <- a@err_y; err_xy <- a@err_xy
  e <- a@effects; terms <- match(e$term, context$terms)
  columns <- match(e$label, colnames(predictors))
  for (j in seq_len(nrow(e))) {
    if (e$role[j] == "outcome") { outcome <- h[, terms[j]]; err_y <- omega[terms[j], terms[j], ] }
    else { predictors[, columns[j]] <- h[, terms[j]] }
  }
  for (j in seq_len(nrow(e))) for (k in seq_len(nrow(e))) {
    if (e$role[j] == "predictor" && e$role[k] == "predictor") err_x[columns[j], columns[k], ] <- omega[terms[j], terms[k], ]
    if (e$role[j] == "outcome" && e$role[k] == "predictor") err_xy[, columns[k]] <- omega[terms[j], terms[k], ]
  }
  controls <- object@struc_fit@diagnostics
  mmom_prepare(list(outcome = outcome, predictors = predictors, err_x = err_x, err_y = err_y, err_xy = err_xy),
    controls$preliminary_moment_control, controls$final_moment_control, controls$backend)
}
