#' Internal fitting function for the M-MOM structural model.
#'
#' @param data_list A list. Each element should correspond to an observation and itself
#'  be a list with fields `outcome`, `predictors`, `err_x`, `err_y`, and `err_xy`.
#' @param preliminary_moment_control Correction controls with `correction` and `alpha`.
#' @param final_moment_control Correction controls with `correction` and `alpha`.
#' @param backend Computation backend: `"auto"`, `"r"`, or `"rcpp"`.
#' @return An S7 `class_mmom_fit` object. Numerical failures use the same class.
#' @export
fit_mmom <- function(data_list,
                     preliminary_moment_control = list(
                       correction = "minimum", alpha = NULL
                     ),
                     final_moment_control = list(
                       correction = "minimum", alpha = NULL
                     ),
                     backend = c("auto", "r", "rcpp")) {
  set_caller_context()

  # Normalize inputs once so all three steps share the same matrices.
  backend <- mmom_backend(backend)
  inputs <- mmom_prepare(data_list, preliminary_moment_control, final_moment_control, backend)

  # Stop at the first numerical failure and return it as an mmom_fit.
  step1 <- mmom_step1(inputs)
  if (mmom_is_failure(step1)) {
    return(step1)
  }
  step2 <- mmom_step2(
    inputs,
    coef0 = step1@coef0,
    predictor_idx = step1@predictor_idx,
    lambda_prelim = step1@lambda_prelim
  )
  if (mmom_is_failure(step2)) {
    return(step2)
  }
  step3 <- mmom_step3(inputs, case_weight = step2@case_weight)
  if (mmom_is_failure(step3)) {
    return(step3)
  }
  step_se <- mmom_step_se(inputs, step1, step2, step3)
  if (mmom_is_failure(step_se)) {
    return(step_se)
  }
  # return mmom fit object
  construct_mmom_fit(inputs, step1, step2, step3, step_se)
}

# Step 1: estimate preliminary coefficients, residual variance, and weights.
mmom_step1 <- function(inputs) {
  # Build the unweighted measurement-error moment.
  err_x_sum <- mmom_weighted_err_cov_sum(
    inputs$err_x, rep(1, inputs$nobs), inputs$backend
  )
  err_moment <- mmom_error_moment(
    err_x_sum, sum(inputs$err_y), colSums(cbind(0, inputs$err_xy))
  )
  observed_moment <- crossprod(cbind(inputs$outcome, inputs$design))
  lambda_prelim <- mmom_smallest_det_root(observed_moment, err_moment)
  correction_prelim <- if (
    identical(inputs$preliminary_moment_control$correction, "none") ||
      max(abs(err_moment)) <= .Machine$double.eps
  ) {
    1
  } else {
    mmom_correction(lambda_prelim, inputs$nobs, inputs$preliminary_moment_control$alpha)
  }

  # Solve the preliminary corrected coefficient equation.
  err_cov_design <- matrix(0, inputs$ncoef, inputs$ncoef)
  err_cov_design[inputs$predictor_index, inputs$predictor_index] <- err_x_sum
  coef_moment <- crossprod(inputs$design) - correction_prelim * err_cov_design
  outcome_moment <- as.vector(
    crossprod(inputs$design, inputs$outcome) -
      correction_prelim * colSums(cbind(0, inputs$err_xy))
  )
  coef0_solve <- repair_solve(
    coef_moment,
    outcome_moment,
    method = "qr"
  )
  if (!coef0_solve@success) {
    fail <- mmom_failure(
      "prelim_coef_solve",
      coef0_solve@message,
      diagnostics = coef0_solve@diagnostics
    )
    return(fail)
  }
  coef0_diagnostics <- coef0_solve@diagnostics
  coef0 <- as.vector(coef0_solve@solution)
  colnames(coef_moment) <- rownames(coef_moment) <- colnames(inputs$design)
  names(coef0) <- colnames(inputs$design)

  class_mmom_step1(
    step = "prelim_unweighted_coef",
    converged = TRUE,
    status_code = 0,
    message = "ok",
    err_x_sum = err_x_sum,
    err_cov_design = err_cov_design,
    lambda_prelim = lambda_prelim,
    correction_prelim = correction_prelim,
    coef0 = coef0,
    predictor_idx = inputs$predictor_index,
    diagnostics = list(
      preliminary_solver = coef0_diagnostics,
      lambda_prelim = lambda_prelim,
      correction_prelim = correction_prelim
    )
  )
}

# Step 2: Estimate structural residual variance and feasible case weights.
mmom_step2 <- function(inputs, coef0, predictor_idx, lambda_prelim) {
  err_resid_var <- mmom_composite_err_y(inputs, coef0[predictor_idx])
  if (any(err_resid_var < -1e-8)) {
    return(mmom_failure(
      "resid_var_meas_err_subtract",
      "The measurement model produced a negative composite error variance."
    ))
  }
  err_resid_var <- pmax(err_resid_var, 0)
  residual <- inputs$outcome - as.vector(inputs$design %*% coef0)
  resid_var <- if (is.finite(lambda_prelim) && lambda_prelim < 1) {
    0
  } else {
    max(
      0,
      sum(residual^2) / max(1, inputs$nobs - inputs$ncoef) -
        mean(err_resid_var)
    )
  }
  if (!is.finite(resid_var)) {
    return(mmom_failure(
      "resid_var",
      "The corrected structural residual variance was non-finite.",
      diagnostics = list(computation_backend = inputs$backend)
    ))
  }
  if (resid_var < 0) {
    return(mmom_failure(
      "resid_var",
      "The corrected structural residual variance was non-positive"
    ))
  }
  case_var <- resid_var + err_resid_var
  if (any(!is.finite(case_var)) || any(case_var <= sqrt(.Machine$double.eps))) {
    return(mmom_failure(
      "case_weights",
      "Case weights were nonpositive or non-finite."
    ))
  }
  class_mmom_step2(
    step = "resid_var_and_weights",
    converged = TRUE,
    status_code = 0,
    message = "ok",
    err_resid_var = err_resid_var,
    resid_var = resid_var,
    case_var = case_var,
    case_weight = 1 / case_var,
    diagnostics = list(
      case_variance_min = min(case_var),
      case_variance_max = max(case_var)
    )
  )
}

# Step 3: apply feasible case weights and solve the final corrected moment.
mmom_step3 <- function(inputs, case_weight) {
  # Recompute moments using the inverse case-variance weights.
  weighted_design <- inputs$design * sqrt(case_weight)
  weighted_outcome <- inputs$outcome * sqrt(case_weight)
  weighted_moment <- crossprod(cbind(weighted_outcome, weighted_design))
  weighted_x_sum <- mmom_weighted_err_cov_sum(
    inputs$err_x, case_weight, inputs$backend
  )
  weighted_cov_design <- matrix(0, inputs$ncoef, inputs$ncoef)
  weighted_cov_design[
    inputs$predictor_index, inputs$predictor_index
  ] <- weighted_x_sum
  weighted_xy_sum <- colSums(
    cbind(0, inputs$err_xy) * case_weight
  )
  weighted_err_moment <- mmom_error_moment(
    weighted_x_sum,
    sum(case_weight * inputs$err_y),
    weighted_xy_sum
  )
  lambda_final <- mmom_smallest_det_root(weighted_moment, weighted_err_moment)
  correction_final <- if (
    identical(inputs$final_moment_control$correction, "none") ||
      max(abs(weighted_err_moment)) <= .Machine$double.eps
  ) {
    1
  } else {
    mmom_correction(
      lambda_final,
      inputs$nobs,
      inputs$final_moment_control$alpha
    )
  }

  # Solve the final corrected coefficient equation.
  corrected_moment <- weighted_moment - correction_final * weighted_err_moment
  coef_moment <- corrected_moment[
    inputs$design_index, inputs$design_index,
    drop = FALSE
  ]
  outcome_moment <- corrected_moment[
    inputs$design_index, 1,
    drop = FALSE
  ]
  coef_solve <- repair_solve(
    coef_moment,
    outcome_moment,
    method = "qr"
  )
  if (!coef_solve@success) {
    fail <- mmom_failure(
      "final_coef_solve",
      coef_solve@message,
      diagnostics = coef_solve@diagnostics
    )
    return(fail)
  }
  coef_diagnostics <- coef_solve@diagnostics
  coef <- as.vector(coef_solve@solution)
  names(coef) <- colnames(inputs$design)
  class_mmom_step3(
    step = "final_weighted_coef",
    converged = TRUE,
    status_code = 0,
    message = "ok",
    weighted_design = weighted_design,
    weighted_moment = weighted_moment,
    weighted_x_sum = weighted_x_sum,
    weighted_cov_design = weighted_cov_design,
    weighted_xy_sum = weighted_xy_sum,
    weighted_err_moment = weighted_err_moment,
    lambda_final = lambda_final,
    correction_final = correction_final,
    corrected_moment = corrected_moment,
    coef_moment = coef_moment,
    coef = coef,
    diagnostics = list(
      final_solver = coef_diagnostics,
      lambda_final = lambda_final,
      correction_final = correction_final
    )
  )
}

# Compute the conditional model-based coefficient covariance.
mmom_step_se <- function(inputs, step1, step2, step3) {
  bread <- step3@coef_moment
  bread_inverse_solve <- repair_solve(
    bread / inputs$nobs, diag(inputs$ncoef),
    method = "qr"
  )
  if (!bread_inverse_solve@success) {
    fail <- mmom_failure(
      "std_err_bread_inverse",
      bread_inverse_solve@message,
      diagnostics = bread_inverse_solve@diagnostics
    )
    return(fail)
  }
  bread_inverse <- as.matrix(bread_inverse_solve@solution)

  tilde_crossprod <- mmom_tilde_crossprod(
    inputs, step1@coef0[step1@predictor_idx], step2@case_weight
  )
  meat <- crossprod(step3@weighted_design) + tilde_crossprod
  bread_diagnostics <- bread_inverse_solve@diagnostics
  coef_vcov <- mmom_symmetrize(
    bread_inverse %*% meat %*% bread_inverse / inputs$nobs^2
  )
  dimnames(coef_vcov) <- list(colnames(inputs$design), colnames(inputs$design))
  coef_var <- diag(coef_vcov)
  if (any(!is.finite(coef_var)) || any(coef_var < 0)) {
    return(mmom_failure(
      "final_coef_covariance",
      "The M-MOM coefficient variance was negative or non-finite."
    ))
  }
  standard_error <- sqrt(coef_var)
  class_mmom_step_se(
    step = "standard_errors",
    converged = TRUE,
    status_code = 0L,
    message = "ok",
    standard_error = standard_error,
    coef_vcov = coef_vcov,
    diagnostics = list(
      bread_solver = bread_diagnostics$solver,
      covariance_dimension = nrow(coef_vcov)
    )
  )
}

# Package estimates, covariance, and diagnostics as one S7 object.
construct_mmom_fit <- function(inputs, step1, step2, step3, step_se, alpha = 0.05) {
  if (alpha <= 0 || alpha >= 1) {
    repair_error("The alpha level must be between 0 and 1.")
  }
  coef_table <- data.frame(
    term = colnames(inputs$design),
    estimate = unname(step3@coef),
    se = unname(step_se@standard_error),
    zstat = unname(step3@coef / step_se@standard_error),
    pval = 2 * (1 - pnorm(abs(unname(step3@coef / step_se@standard_error)))),
    ci_low = unname(step3@coef - stats::qnorm(1 - alpha / 2) * step_se@standard_error),
    ci_high = unname(step3@coef + stats::qnorm(1 - alpha / 2) * step_se@standard_error),
    row.names = NULL,
    check.names = FALSE
  )
  class_mmom_fit(
    converged = TRUE,
    status_code = 0L,
    message = "ok",
    coefficients = coef_table,
    vcov = step_se@coef_vcov,
    preliminary_coef = step1@coef0,
    corrected_design_moment = step3@coef_moment,
    analysis_weights = step2@case_weight,
    steps = list(
      step1 = step1,
      step2 = step2,
      step3 = step3,
      step_se = step_se
    ),
    diagnostics = list(
      nobs = inputs$nobs,
      npredictors = inputs$npredictors,
      backend = inputs$backend,
      lambda_prelim = step1@lambda_prelim,
      lambda_final = step3@lambda_final,
      correction_prelim = step1@correction_prelim,
      correction_final = step3@correction_final,
      resid_var = step2@resid_var,
      preliminary_moment_control = inputs$preliminary_moment_control,
      final_moment_control = inputs$final_moment_control
    )
  )
}
