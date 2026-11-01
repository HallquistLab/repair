repair_if_structural <- function(inputs, frequency = rep(1, inputs$nobs)) {
  if (length(frequency) != inputs$nobs || any(!is.finite(frequency)) || any(frequency < 0) ||
      abs(sum(frequency) - inputs$nobs) > 1e-8) repair_error("Case frequencies must be nonnegative and sum to the original person count.")
  inputs$case_frequency <- frequency
  s1 <- mmom_step1(inputs)
  if (mmom_is_failure(s1)) repair_error("Influence preliminary moment failed: {s1@message}")
  s2 <- mmom_step2(inputs, s1@coef0, s1@predictor_idx, s1@lambda_prelim)
  if (mmom_is_failure(s2)) repair_error("Influence weighting failed: {s2@message}")
  s3 <- mmom_step3(inputs, s2@case_weight)
  if (mmom_is_failure(s3)) repair_error("Influence final moment failed: {s3@message}")
  branch <- function(lambda, control) {
    if (control$correction == "none") return("none")
    if (!is.finite(lambda)) return("fallback")
    if (lambda <= 1 + 1 / inputs$nobs) "root" else "alpha"
  }
  signature <- paste(branch(s1@lambda_prelim, inputs$preliminary_moment_control),
    if (is.finite(s1@lambda_prelim) && s1@lambda_prelim < 1) "forced_zero" else if (s2@resid_var == 0) "clipped_zero" else "positive",
    branch(s3@lambda_final, inputs$final_moment_control), sep = ":")
  list(value = s3@coef, signature = signature, steps = list(s1, s2, s3))
}

# Central sensitivities with step halving, one-sided consistency and unchanged
# active branches. This differentiates the finite-sample estimator; it does not
# resample subjects or freeze its estimated structural weights/corrections.
repair_if_difference <- function(fn, baseline, dimension, steps, tolerance, max_refine = 3L) {
  jacobian <- matrix(NA_real_, length(baseline$value), dimension)
  errors <- rep(NA_real_, dimension); used <- steps
  for (j in seq_len(dimension)) {
    step <- steps[j]; accepted <- FALSE
    for (refine in 0:max_refine) {
      eval <- function(delta) fn(j, delta)
      plus <- eval(step); minus <- eval(-step); ph <- eval(step / 2); mh <- eval(-step / 2)
      if (!all(vapply(list(plus, minus, ph, mh), function(x) identical(x$signature, baseline$signature), logical(1)))) {
        step <- step / 2; next
      }
      coarse <- (plus$value - minus$value) / (2 * step)
      fine <- (ph$value - mh$value) / step
      scale <- max(1e-5, max(abs(fine)))
      change <- max(abs(fine - coarse)) / scale
      one_sided <- max(abs((ph$value - baseline$value) * 2 / step - (baseline$value - mh$value) * 2 / step)) / scale
      if (all(is.finite(fine)) && max(change, one_sided) <= tolerance) {
        jacobian[, j] <- fine; errors[j] <- max(change, one_sided); used[j] <- step / 2; accepted <- TRUE; break
      }
      step <- step / 2
    }
    if (!accepted) repair_error("Local derivative unstable or crossed an estimation branch at coordinate {j}; no influence interval is available.")
  }
  list(jacobian = jacobian, error = errors, step = used)
}

repair_if_measurement <- function(context, step, tolerance) {
  if (!context$converged) repair_error("Measurement fit failed its recorded convergence policy.")
  if (context$singular) repair_error("Full influence inference is unavailable at a singular random-effect covariance boundary.")
  psi <- context$psi; n <- length(context$ids)
  base <- repair_if_likelihood(context, psi, case_scores = TRUE)
  scores <- scale(base$case_scores, scale = FALSE)
  stationarity <- max(abs(base$score) / (1 + sqrt(colSums(scores^2))))
  if (!is.finite(stationarity) || stationarity > .001)
    repair_error("Measurement estimating equations are not sufficiently centered at this fit (scaled score {signif(stationarity, 3)}).")
  difference <- repair_if_difference(function(j, delta) {
    changed <- psi; changed[j] <- changed[j] + delta
    list(value = repair_if_likelihood(context, changed)$score / n, signature = "measurement")
  }, list(value = base$score / n, signature = "measurement"), length(psi), step * pmax(1, abs(psi)), tolerance)
  bread <- -difference$jacobian; asymmetry <- max(abs(bread - t(bread))) / (1 + max(abs(bread)))
  if (asymmetry > tolerance) repair_error("Measurement sensitivity matrix is not numerically symmetric.")
  bread <- (bread + t(bread)) / 2
  if (any(diag(bread) <= 0)) repair_error("Measurement sensitivity has nonpositive diagonal curvature.")
  scaling <- sqrt(diag(bread)); standardized <- bread / outer(scaling, scaling)
  condition <- rcond(standardized)
  if (condition < 1e-9 || min(eigen(standardized, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8)
    repair_error("Measurement sensitivity is singular, indefinite or weakly identified.")
  influence <- t(solve(bread, t(scores)))
  colnames(influence) <- colnames(scores) <- names(psi); rownames(influence) <- rownames(scores) <- context$ids
  list(influence = influence, scores = scores, bread = bread, psi = psi,
    stationarity = stationarity, rcond = condition, derivative_error = max(difference$error),
    loglik = base$loglik, criterion = if (context$reml) "REML" else "ML")
}

repair_sandwich_compute <- function(object, step, tolerance, context = NULL, measurement = NULL) {
  a <- object@analysis_data; n <- length(a@ids); fit <- object@struc_fit
  controls <- fit@diagnostics
  inputs <- mmom_prepare(list(outcome = a@outcome, predictors = a@predictors, err_x = a@err_x,
    err_y = a@err_y, err_xy = a@err_xy), controls$preliminary_moment_control, controls$final_moment_control, controls$backend)
  estimate <- stats::coef(object); baseline <- repair_if_structural(inputs)
  if (max(abs(baseline$value - estimate)) > 1e-7 * (1 + max(abs(estimate))))
    repair_error("Stored structural estimates do not reproduce the current estimation pipeline.")
  blank <- matrix(NA_real_, length(estimate), length(estimate), dimnames = list(names(estimate), names(estimate)))
  result <- structure(list(estimate = estimate, ids = a@ids,
    covariance = list(full = blank, stage2 = blank, measurement = blank, cross = blank),
    influence = list(), diagnostics = list(full_available = FALSE, stage2_available = FALSE,
      full_reason = "", stage2_reason = "", branch = baseline$signature, step = step, tolerance = tolerance,
      criterion = if (lme4::isREML(object@mlm_fit[[1L]])) "REML" else "ML",
      measurement_singular = lme4::isSingular(object@mlm_fit[[1L]])),
    measurement = NULL, versions = list(R = R.version.string,
      repair = as.character(utils::packageVersion("repair")), lme4 = as.character(utils::packageVersion("lme4"))),
    call = match.call()), class = "repair_sandwich")
  stage2 <- tryCatch({
    difference <- repair_if_difference(function(j, delta) {
      f <- rep(1 - delta / n, n); f[j] <- f[j] + delta
      repair_if_structural(inputs, f)
    }, baseline, n, rep(step * 10, n), tolerance)
    influence <- t(difference$jacobian) * n
    influence <- scale(influence, scale = FALSE)
    colnames(influence) <- names(estimate); rownames(influence) <- a@ids
    list(influence = influence, error = max(difference$error))
  }, error = function(e) e)
  if (inherits(stage2, "error")) {
    result$diagnostics$stage2_reason <- result$diagnostics$full_reason <- conditionMessage(stage2)
    return(result)
  }
  result$diagnostics$stage2_available <- TRUE
  result$diagnostics$frequency_derivative_error <- stage2$error
  result$influence$stage2 <- stage2$influence
  result$covariance$stage2 <- crossprod(stage2$influence) / n^2
  full <- tryCatch({
    if (is.null(context)) context <- repair_if_context(object)
    if (!identical(context$ids, a@ids)) repair_error("Influence context has a different person order.")
    if (is.null(measurement)) measurement <- repair_if_measurement(context, step, tolerance)
    if (inherits(measurement, "error")) stop(measurement)
    rebuilt <- repair_if_inputs(object, context, context$psi)
    for (field in c("outcome", "predictors", "err_x", "err_y", "err_xy"))
      if (!isTRUE(all.equal(rebuilt[[field]], inputs[[field]], tolerance = 1e-7, check.attributes = FALSE)))
        repair_error("Likelihood-only scoring reconstruction does not reproduce {.val {field}}.")
    difference <- repair_if_difference(function(j, delta) {
      psi <- context$psi; psi[j] <- psi[j] + delta
      repair_if_structural(repair_if_inputs(object, context, psi))
    }, baseline, length(context$psi), step * pmax(1, abs(context$psi)), tolerance)
    nuisance <- measurement$influence %*% t(difference$jacobian)
    colnames(nuisance) <- names(estimate); rownames(nuisance) <- a@ids
    influence <- stage2$influence + nuisance
    Vn <- crossprod(nuisance) / n^2
    Vcross <- crossprod(stage2$influence, nuisance) / n^2
    list(influence = influence, nuisance = nuisance, measurement = measurement,
      covariance = crossprod(influence) / n^2, nuisance_covariance = Vn,
      cross = Vcross + t(Vcross), derivative = difference$jacobian, error = max(difference$error))
  }, error = function(e) e)
  if (inherits(full, "error")) result$diagnostics$full_reason <- conditionMessage(full) else {
    result$diagnostics$full_available <- TRUE
    result$diagnostics$measurement_derivative_error <- full$error
    result$covariance$full <- full$covariance
    result$covariance$measurement <- full$nuisance_covariance
    result$covariance$cross <- full$cross
    result$influence$measurement <- full$nuisance
    result$influence$full <- full$influence
    result$measurement <- full$measurement
    result$measurement$structural_sensitivity <- full$derivative
  }
  result
}

#' Full person-level influence-function covariance for REPAIR
#'
#' Experimental empirical influence covariance for the complete fitted M-MOM
#' pipeline. The original coefficient estimates and ordinary conditional SEs are
#' unchanged. A separate result records full and frozen-measurement (stage2)
#' covariance, person influences, cross-stage contributions and availability.
#' @param object A fitted [repair()] analysis with one Gaussian lmer model.
#' @param step Positive numerical sensitivity step; default 1e-4. Person-frequency
#'   steps start at ten times this value. Central differences, step halving,
#'   one-sided agreement and active-branch checks validate local derivatives.
#' @param tolerance Positive derivative-agreement tolerance; default .03.
#' @details Supports a single Gaussian lmer with independent target people and
#'   contained nuisance groups, prior weights, offsets and stored factor coding.
#'   Measurement transformations and weight/offset recipes must be supported
#'   pointwise expressions or explicit columns. Recomputed sample-wide scaling,
#'   data-adaptive bases and position-specific literal vectors are rejected;
#'   their additional sampling derivatives are not implemented.
#'   Analyzed and measurement cohorts must coincide. nlme and multiple-model
#'   derivatives are not implemented. Unsupported scope errors explicitly.
#'
#'   The structural case influence differentiates all preliminary/final moments,
#'   score-error corrections, residual variance, case weights and data-dependent
#'   correction multipliers. Measurement influence uses marginal Gaussian scores
#'   and a numerically differentiated score Jacobian. REML includes the global
#'   restricted-likelihood determinant and its case-frequency derivative, rather
#'   than substituting ML scores. Total influence is the structural influence plus
#'   the measurement-parameter influence propagated through that same pipeline.
#'   Its outer product includes both cross-stage covariance terms automatically.
#'
#'   This is an asymptotic covariance computed with validated numerical
#'   sensitivities, not a resampling SE or a closed-form symbolic implementation.
#'   The full result is withheld at singular/weakly identified measurement fits
#'   or unstable/branch-crossing derivatives. Stage2 inference can remain available
#'   separately; there is no silent conditional fallback. A stable local active
#'   regularization/floor branch does not establish uniform coverage near a boundary.
#'   No finite-sample HC correction, bias correction, or coverage guarantee is made.
#' @return A separate S3 `repair_sandwich` result. `estimate` retains original
#'   coefficients; `covariance` has full, stage2, measurement and cross blocks;
#'   `influence` has ID-keyed matrices; `measurement` retains the score/bread,
#'   parameter influences and structural sensitivity; `diagnostics` retains
#'   availability/reasons, active branches and derivative agreement. Failed local
#'   linearizations return an inspectable result with inference unavailable.
#'   Use [repair_sandwich_coefficients()], `vcov()`, `confint()` or `summary()`.
#' @examples
#' data(repair_example)
#' model <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, model, repair_example$people, "person")
#' inf <- repair_sandwich(fit)
#' summary(inf)
#' @export
repair_sandwich <- function(object, step = 1e-4, tolerance = .03) {
  repair_require_model(object); repair_require_fit(object)
  for (value in list(step, tolerance))
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      repair_error("`step` and `tolerance` must be positive finite scalars.")
  # Scope checks are programming/input errors; numerical unavailability is data.
  context <- repair_if_context(object)
  result <- repair_sandwich_compute(object, step, tolerance, context)
  result$call <- match.call()
  result
}
