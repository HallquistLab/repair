# Augmented estimating system. Raw and weighted augmented moments are explicit
# nuisance parameters; the final coefficient equation remains person additive.
# This representation is deliberately specified: leverage corrections are not
# invariant to arbitrary re-expression of nonlinear estimating equations.
repair_lev_pack <- function(x) x[upper.tri(x, diag = TRUE)]
repair_lev_unpack <- function(x, k) {
  out <- matrix(0, k, k); out[upper.tri(out, diag = TRUE)] <- x
  out + t(out) - diag(diag(out), k)
}

repair_lev_rows <- function(inputs) {
  n <- inputs$nobs; k <- inputs$ncoef + 1L; p <- inputs$npredictors
  z <- cbind(inputs$outcome, inputs$design)
  M <- E <- matrix(0, n, k * (k + 1L) / 2L)
  for (i in seq_len(n)) {
    M[i, ] <- repair_lev_pack(tcrossprod(z[i, ]))
    E[i, ] <- repair_lev_pack(mmom_error_moment(matrix(inputs$err_x[, , i], p, p), inputs$err_y[i], c(0, inputs$err_xy[i, ])))
  }
  diagonals <- which(upper.tri(diag(k), diag = TRUE) & diag(k) == 1)
  noisy <- which(colSums(abs(E[, match(diagonals, which(upper.tri(diag(k), diag = TRUE))), drop = FALSE])) > 0)
  allowed <- outer(seq_len(k) %in% noisy, seq_len(k) %in% noisy, "&")
  positions <- which(allowed[upper.tri(allowed, diag = TRUE)])
  list(M = M, E = E[, positions, drop = FALSE], error_positions = positions)
}

repair_lev_state <- function(inputs) {
  base <- repair_if_structural(inputs); rows <- repair_lev_rows(inputs)
  w <- base$steps[[2]]@case_weight
  list(value = c(colMeans(rows$M), colMeans(rows$E), colMeans(rows$M * w), colMeans(rows$E * w), log(mean(w)), base$value),
    inputs = inputs, rows = rows, signature = base$signature)
}

repair_lev_equations <- function(inputs, state, rows = repair_lev_rows(inputs)) {
  n <- inputs$nobs; p <- inputs$ncoef; k <- p + 1L; m <- ncol(rows$M); e <- ncol(rows$E)
  error_matrix <- function(x) {
    full <- numeric(k * (k + 1L) / 2L); full[rows$error_positions] <- x
    repair_lev_unpack(full, k)
  }
  M <- repair_lev_unpack(state[seq_len(m)], k)
  E <- error_matrix(state[m + seq_len(e)])
  W <- repair_lev_unpack(state[m + e + seq_len(m)], k)
  F <- error_matrix(state[2L * m + e + seq_len(e)])
  average_weight <- exp(state[2L * (m + e) + 1L])
  beta <- utils::tail(state, p); index <- seq.int(2L, k)
  correction <- function(M, E, control) {
    lambda <- mmom_smallest_det_root(M, E)
    value <- if (control$correction == "none" || max(abs(E)) <= .Machine$double.eps) 1 else
      mmom_correction(lambda, n, control$alpha)
    branch <- if (control$correction == "none") "none" else if (!is.finite(lambda)) "fallback" else
      if (lambda <= 1 + 1 / n) "root" else "alpha"
    list(value = value, lambda = lambda, branch = branch)
  }
  first <- correction(M, E, inputs$preliminary_moment_control)
  b0 <- drop(solve(M[index, index, drop = FALSE] - first$value * E[index, index, drop = FALSE],
    M[index, 1L] - first$value * E[index, 1L]))
  contrast <- c(1, -b0)
  resid <- if (is.finite(first$lambda) && first$lambda < 1) 0 else
    max(0, n / max(1, n - p) * drop(crossprod(contrast, M %*% contrast)) -
      drop(crossprod(contrast, E %*% contrast)))
  variance <- mmom_composite_err_y(inputs, b0[-1L])
  if (any(variance < -1e-8)) repair_error("Negative composite error variance in augmented system.")
  w <- 1 / (resid + pmax(variance, 0))
  if (any(!is.finite(w)) || any(1 / w <= sqrt(.Machine$double.eps)))
    repair_error("Invalid augmented-system case weights.")
  last <- correction(W, F, inputs$final_moment_control)
  coefficients <- matrix(0, n, p)
  for (i in seq_len(n)) {
    corrected <- repair_lev_unpack(rows$M[i, ], k) - last$value * error_matrix(rows$E[i, ])
    coefficients[i, ] <- w[i] / average_weight * (corrected[index, 1L] - drop(corrected[index, index, drop = FALSE] %*% beta))
  }
  equations <- cbind(sweep(rows$M, 2L, repair_lev_pack(M)), sweep(rows$E, 2L, repair_lev_pack(E)[rows$error_positions]),
    sweep(rows$M * w, 2L, repair_lev_pack(W)), sweep(rows$E * w, 2L, repair_lev_pack(F)[rows$error_positions]),
    w / average_weight - 1, coefficients)
  signature <- paste(first$branch, if (is.finite(first$lambda) && first$lambda < 1) "forced_zero" else
    if (resid == 0) "clipped_zero" else "positive", last$branch, sep = ":")
  list(value = equations, signature = signature)
}

# Scale both rows and columns; a raw condition number depends arbitrarily on
# units of moment equations. No ridge, truncation or pseudoinverse is used.
repair_lev_solve <- function(J, rhs, scales, tolerance = 1e-9) {
  scaled <- sweep(J, 2L, scales, "*")
  rows <- sqrt(rowSums(scaled^2))
  if (any(!is.finite(rows)) || any(rows == 0)) repair_error("Zero/invalid sensitivity row.")
  scaled <- scaled / rows; condition <- rcond(scaled)
  if (!is.finite(condition) || condition < tolerance) repair_error("Weak or singular deletion sensitivity.")
  value <- solve(scaled, rhs / rows) * scales
  list(value = value, rcond = condition)
}

repair_leverage_compute <- function(object, step, tolerance, context, baseline) {
  result <- baseline; result$method <- "Full-system HC3-style leverage correction (experimental)"
  result$inference <- "Experimental full-system HC3-style covariance; normal Wald inference"
  class(result) <- c("repair_leverage", "repair_sandwich")
  result$covariance <- lapply(result$covariance, function(x) { x[] <- NA_real_; x })
  result$influence <- list()
  result$diagnostics$full_available <- result$diagnostics$stage2_available <- FALSE
  result$diagnostics$full_reason <- result$diagnostics$stage2_reason <- "Leverage computation incomplete."
  n <- length(context$ids); q <- length(context$psi)
  inputs <- repair_if_inputs(object, context, context$psi)
  state <- repair_lev_state(inputs); s <- state$value; r <- length(s); p <- inputs$ncoef
  selected <- r - p + seq_len(p); parameter <- c(context$psi, s)
  scales <- pmax(1, abs(parameter)); d <- length(parameter)
  records <- data.frame(id = context$ids, usable = FALSE, reason = "", rcond = NA_real_)
  result$deletions <- matrix(NA_real_, n, p, dimnames = list(context$ids, names(coef(object))))
  result$diagnostics$correction <- "HC3-style; uncentered deletion-shift outer products"
  result$diagnostics$finite_N_policy <- "Original N and regularization policy held fixed in analytic correction."
  computed <- tryCatch({
    # The same finite-difference validation is applied to every block. Sparse
    # moment/state directions require no measurement likelihood evaluation.
    equations <- repair_lev_equations(inputs, s, state$rows)
    if (!identical(equations$signature, state$signature) || max(abs(colMeans(equations$value))) > 1e-7 * (1 + max(abs(s))))
      repair_error("Augmented system does not reproduce the canonical estimator.")
    native <- repair_if_likelihood(context, context$psi, case_scores = TRUE, deletions = TRUE)
    score <- rbind(native$score, native$deletion_scores)
    derivative <- repair_if_difference(function(j, delta) {
      psi <- context$psi; psi[j] <- psi[j] + delta
      value <- repair_if_likelihood(context, psi, deletions = TRUE)
      data <- repair_if_inputs(object, context, psi)
      eq <- repair_lev_equations(data, s)
      list(value = c(c(t(rbind(value$score, value$deletion_scores))), c(eq$value)), signature = eq$signature)
    }, list(value = c(c(t(score)), c(equations$value)), signature = state$signature), q,
      step * scales[seq_len(q)], tolerance)
    state_derivative <- repair_if_difference(function(j, delta) {
      changed <- s; changed[j] <- changed[j] + delta
      eq <- repair_lev_equations(inputs, changed, state$rows)
      list(value = c(eq$value), signature = eq$signature)
    }, list(value = c(equations$value), signature = state$signature), r,
      step * scales[q + seq_len(r)], tolerance)
    nscore <- (n + 1L) * q
    # Per-person augmented equations: rows people, columns equations, layers parameters.
    structural <- array(0, c(n, r, d))
    for (j in seq_len(q)) structural[, , j] <- matrix(derivative$jacobian[nscore + seq_len(n * r), j], n, r)
    for (j in seq_len(r)) structural[, , q + j] <- matrix(state_derivative$jacobian[, j], n, r)
    J <- matrix(0, d, d)
    J[seq_len(q), seq_len(q)] <- derivative$jacobian[seq_len(q), ] / n
    J[q + seq_len(r), ] <- apply(structural, c(2, 3), mean)
    cases <- cbind(native$case_scores, equations$value)
    cases <- scale(cases, scale = FALSE)
    plain <- -repair_lev_solve(J, t(cases) / n, scales)$value
    plain <- t(plain[q + selected, , drop = FALSE]) * n
    if (baseline$diagnostics$full_available) {
      gap <- max(abs(plain - baseline$influence$full)) / (1 + max(abs(baseline$influence$full)))
      if (gap > .001) repair_error("Augmented HC0 does not reproduce the canonical full influence (relative gap {signif(gap, 3)}).")
    } else gap <- NA_real_
    direct <- propagated <- result$deletions
    for (i in seq_len(n)) {
      value <- tryCatch({
        # Exact deletion of the REML determinant at the original parameters.
        rows <- i * q + seq_len(q)
        Ji <- J
        Ji[seq_len(q), seq_len(q)] <- derivative$jacobian[rows, ] / n
        Ji[q + seq_len(r), ] <- (apply(structural, c(2, 3), sum) - structural[i, , ]) / n
        residual <- c((score[i + 1L, ] - score[1L, ]) / n, -equations$value[i, ] / n)
        one <- two <- residual
        one[seq_len(q)] <- 0; two[q + seq_len(r)] <- 0
        a <- -repair_lev_solve(Ji, one, scales)$value
        b <- -repair_lev_solve(Ji, two, scales)$value
        list(direct = a[q + selected], propagated = b[q + selected],
          rcond = repair_lev_solve(Ji, residual, scales)$rcond)
      }, error = function(e) e)
      if (inherits(value, "error")) records$reason[i] <- conditionMessage(value) else {
        direct[i, ] <- value$direct; propagated[i, ] <- value$propagated
        records$usable[i] <- TRUE; records$rcond[i] <- value$rcond
      }
    }
    list(direct = direct, propagated = propagated, records = records, gap = gap,
      derivative_error = max(derivative$error, state_derivative$error), dimension = d)
  }, error = function(e) e)
  if (inherits(computed, "error")) {
    result$diagnostics$full_reason <- result$diagnostics$stage2_reason <- conditionMessage(computed)
    records$reason <- conditionMessage(computed)
  } else {
    records <- computed$records
    result$diagnostics$augmented_dimension <- computed$dimension
    result$diagnostics$hc0_relative_gap <- computed$gap
    result$diagnostics$leverage_derivative_error <- computed$derivative_error
    if (all(records$usable) && baseline$diagnostics$full_available) {
      total <- computed$direct + computed$propagated
      result$deletions <- total
      result$covariance$full <- crossprod(total)
      result$covariance$stage2 <- crossprod(computed$direct)
      result$covariance$measurement <- crossprod(computed$propagated)
      cross <- crossprod(computed$direct, computed$propagated)
      result$covariance$cross <- cross + t(cross)
      result$influence <- list(full = -n * total, stage2 = -n * computed$direct, measurement = -n * computed$propagated)
      result$diagnostics$full_available <- result$diagnostics$stage2_available <- TRUE
      result$diagnostics$full_reason <- result$diagnostics$stage2_reason <- ""
    } else {
      reason <- if (!baseline$diagnostics$full_available) baseline$diagnostics$full_reason else
        paste(sum(!records$usable), "person deletion sensitivities unavailable; no omissions permitted.")
      result$diagnostics$full_reason <- result$diagnostics$stage2_reason <- reason
    }
  }
  result$diagnostics$deletions <- records
  result
}

#' Full-system HC3-style leverage correction
#' @param object A fitted [repair()] analysis supported by [repair_sandwich()].
#' @param step,tolerance Validated numerical sensitivity controls as in [repair_sandwich()].
#' @details Experimental deletion-sensitivity correction for an explicitly
#'   augmented system: Gaussian measurement equations, unweighted observed/error
#'   moments, weighted observed/error moments, and final corrected coefficient
#'   equations. An additional log mean-weight equation normalizes the final
#'   coefficient equation, cancelling irrelevant common weights in OLS.
#'   Structurally zero error-moment entries remain fixed. All nuisance weights
#'   and correction multipliers are differentiated.
#'   For each person, solve the deleted system's Jacobian at the original fit.
#'   REML uses the exact deleted fixed-design determinant. Original N and user
#'   policies are held fixed in the analytic structural equations. Sum the
#'   uncentered coefficient deletion-shift outer products; this reduces to HC3
#'   in ordinary OLS. This is not a derived CR2/Satterthwaite correction and is
#'   representation dependent in nonlinear systems. No fitted models are refitted.
#'   All deletion sensitivities and the original full influence must be usable;
#'   otherwise inference is withheld, with no omitted cases or fallback.
#' @return A `repair_leverage` result inheriting from `repair_sandwich`, supported
#'   by [repair_sandwich_coefficients()], `summary()`, `vcov()` and `confint()`.
#'   `deletions` and diagnostics retain person shifts, conditioning and failures.
#'   Original point estimates are unchanged; intervals use normal reference.
#' @export
repair_leverage <- function(object, step = 1e-4, tolerance = .03) {
  baseline <- repair_sandwich(object, step, tolerance)
  context <- repair_if_context(object)
  if (!baseline$diagnostics$full_available) {
    baseline$method <- "Full-system HC3-style leverage correction (experimental)"
    baseline$inference <- "Experimental full-system HC3-style covariance; normal Wald inference"
    baseline$diagnostics$stage2_available <- FALSE
    baseline$diagnostics$stage2_reason <- baseline$diagnostics$full_reason
    baseline$covariance <- lapply(baseline$covariance, function(x) { x[] <- NA_real_; x })
    baseline$influence <- list()
    class(baseline) <- c("repair_leverage", "repair_sandwich")
    baseline$call <- match.call()
    return(baseline)
  }
  result <- repair_leverage_compute(object, step, tolerance, context, baseline)
  result$call <- match.call(); result
}
