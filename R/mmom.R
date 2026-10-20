#' Fit the stabilized M-MOM structural model.
#'
#' @param data_list Level-2 observations or a normalized list containing
#'   `y`, `predictors`, and `err_cov`.
#' @param preliminary_moment Correction controls with `regularize` and `alpha`.
#' @param variance_bread Correction controls for the variance bread.
#' @param backend Computation backend: `"auto"`, `"r"`, or `"rcpp"`.
#' @return An S7 `mmom_fit` object. Numerical failures use the same class.
#' @export
fit_mmom <- function(data_list,
                     preliminary_moment = list(
                       regularize = FALSE, alpha = NULL
                     ),
                     variance_bread = list(
                       regularize = FALSE, alpha = NULL
                     ),
                     backend = c("auto", "r", "rcpp")) {
  # Normalize inputs once so all three steps share the same matrices.
  backend <- mmom_backend(backend)
  inputs <- mmom_prepare(data_list, preliminary_moment, variance_bread, backend)

  # Stop at the first numerical failure and return it as an mmom_fit.
  step1 <- mmom_step1(inputs)
  if (mmom_is_failure(step1)) {
    return(step1)
  }
  step2 <- mmom_step2(inputs, step1)
  if (mmom_is_failure(step2)) {
    return(step2)
  }
  mmom_step3(inputs, step1, step2)
}

# Step 1: estimate preliminary coefficients, residual variance, and weights.
mmom_step1 <- function(inputs) {
  # Build the unweighted measurement-error moment.
  err_cov_sum <- mmom_weighted_err_cov_sum(
    inputs$err_cov, rep(1, inputs$nobs), inputs$backend
  )
  err_moment <- mmom_error_moment(
    err_cov_sum, sum(inputs$err_var), colSums(cbind(0, inputs$err_cross))
  )
  observed_moment <- crossprod(cbind(inputs$outcome, inputs$design))
  lambda1 <- mmom_smallest_det_root(observed_moment, err_moment)
  correction1 <- mmom_correction(
    lambda1, inputs$nobs, inputs$preliminary$alpha
  )
  moment_scale <- if (inputs$preliminary$regularize) {
    1
  } else {
    correction1
  }

  # Solve the preliminary corrected coefficient equation.
  err_cov_design <- matrix(0, inputs$ncoef, inputs$ncoef)
  err_cov_design[inputs$predictor_index, inputs$predictor_index] <- err_cov_sum
  coef_moment <- crossprod(inputs$design) - moment_scale * err_cov_design
  outcome_moment <- as.vector(
    crossprod(inputs$design, inputs$outcome) -
      moment_scale * colSums(cbind(0, inputs$err_cross))
  )
  coef0 <- mmom_solve(
    coef_moment,
    outcome_moment,
    "step1_preliminary",
    "The preliminary corrected moment was singular."
  )
  if (mmom_is_failure(coef0)) {
    return(coef0)
  }
  coef0_diagnostics <- attr(coef0, "mmom_solver_diagnostics")
  coef0 <- as.vector(coef0)
  colnames(coef_moment) <- rownames(coef_moment) <- colnames(inputs$design)
  names(coef0) <- colnames(inputs$design)
  coef0_predictors <- coef0[inputs$predictor_index]

  # Estimate structural residual variance and feasible case weights.
  err_resid_var <- mmom_composite_err_var(inputs, coef0_predictors)
  if (any(err_resid_var < -1e-8)) {
    return(mmom_failure(
      "step1_measurement_error",
      "The measurement model produced a negative composite error variance."
    ))
  }
  err_resid_var <- pmax(err_resid_var, 0)
  residual <- inputs$outcome - as.vector(inputs$design %*% coef0)
  resid_var <- if (is.finite(lambda1) && lambda1 < 1) {
    0
  } else {
    max(
      0,
      sum(residual^2) / max(1, inputs$nobs - inputs$ncoef) -
        mean(err_resid_var)
    )
  }
  case_var <- resid_var + err_resid_var
  if (any(!is.finite(case_var)) || any(case_var <= sqrt(.Machine$double.eps))) {
    return(mmom_failure(
      "step1_weights",
      "Feasible case weights were nonpositive or non-finite."
    ))
  }
  list(
    err_cov_sum = err_cov_sum,
    err_cov_design = err_cov_design,
    lambda1 = lambda1,
    correction1 = correction1,
    coef0 = coef0,
    coef0_predictors = coef0_predictors,
    err_resid_var = err_resid_var,
    resid_var = resid_var,
    case_var = case_var,
    case_weight = 1 / case_var,
    diagnostics = list(
      step = "step1",
      preliminary_solver = coef0_diagnostics,
      lambda1 = lambda1,
      correction1 = correction1,
      case_variance_min = min(case_var),
      case_variance_max = max(case_var)
    )
  )
}

# Step 2: apply feasible case weights and solve the final corrected moment.
mmom_step2 <- function(inputs, step1) {
  # Recompute moments using the inverse case-variance weights.
  weighted_design <- inputs$design * sqrt(step1$case_weight)
  weighted_outcome <- inputs$outcome * sqrt(step1$case_weight)
  weighted_moment <- crossprod(cbind(weighted_outcome, weighted_design))
  weighted_cov_sum <- mmom_weighted_err_cov_sum(
    inputs$err_cov, step1$case_weight, inputs$backend
  )
  weighted_cov_design <- matrix(0, inputs$ncoef, inputs$ncoef)
  weighted_cov_design[
    inputs$predictor_index, inputs$predictor_index
  ] <- weighted_cov_sum
  weighted_cross_sum <- colSums(
    cbind(0, inputs$err_cross) * step1$case_weight
  )
  weighted_err_moment <- mmom_error_moment(
    weighted_cov_sum,
    sum(step1$case_weight * inputs$err_var),
    weighted_cross_sum
  )
  lambda2 <- mmom_smallest_det_root(weighted_moment, weighted_err_moment)
  correction3 <- mmom_correction(
    lambda2, inputs$nobs, inputs$variance$alpha
  )

  # Solve the final corrected coefficient equation.
  corrected_moment <- weighted_moment - correction3 * weighted_err_moment
  coef_moment <- corrected_moment[
    inputs$design_index, inputs$design_index,
    drop = FALSE
  ]
  outcome_moment <- corrected_moment[
    inputs$design_index, 1,
    drop = FALSE
  ]
  coef <- mmom_solve(
    coef_moment,
    outcome_moment,
    "step2_final",
    "The final corrected M-MOM moment was singular."
  )
  if (mmom_is_failure(coef)) {
    return(coef)
  }
  coef_diagnostics <- attr(coef, "mmom_solver_diagnostics")
  coef <- as.vector(coef)
  names(coef) <- colnames(inputs$design)
  list(
    weighted_design = weighted_design,
    weighted_moment = weighted_moment,
    weighted_cov_sum = weighted_cov_sum,
    weighted_cov_design = weighted_cov_design,
    weighted_cross_sum = weighted_cross_sum,
    weighted_err_moment = weighted_err_moment,
    lambda2 = lambda2,
    correction3 = correction3,
    corrected_moment = corrected_moment,
    coef_moment = coef_moment,
    coef = coef,
    diagnostics = list(
      step = "step2",
      final_solver = coef_diagnostics,
      lambda2 = lambda2,
      correction3 = correction3
    )
  )
}

# Step 3: compute the conditional model-based coefficient covariance.
mmom_step3 <- function(inputs, step1, step2) {
  # Assemble the variance bread and estimating-equation meat.
  bread <- if (inputs$variance$regularize) {
    step2$weighted_moment[
      inputs$design_index, inputs$design_index,
      drop = FALSE
    ] - step2$weighted_cov_design
  } else {
    step2$coef_moment
  }
  bread_inverse <- mmom_solve(
    bread / inputs$nobs, diag(inputs$ncoef),
    "step3_bread",
    "The M-MOM variance bread was singular."
  )
  if (mmom_is_failure(bread_inverse)) {
    return(bread_inverse)
  }

  tilde_crossprod <- mmom_tilde_crossprod(
    inputs, step1$coef0_predictors, step1$case_weight
  )
  meat <- crossprod(step2$weighted_design) + tilde_crossprod
  bread_diagnostics <- attr(bread_inverse, "mmom_solver_diagnostics")
  coef_vcov <- mmom_symmetrize(
    bread_inverse %*% meat %*% bread_inverse / inputs$nobs^2
  )
  coef_var <- diag(coef_vcov)
  if (any(!is.finite(coef_var)) || any(coef_var < 0)) {
    return(mmom_failure(
      "step3_covariance",
      "The M-MOM coefficient variance was negative or non-finite."
    ))
  }
  standard_error <- sqrt(coef_var)

  # Package estimates, covariance, and diagnostics as one S7 object.
  coef_table <- data.frame(
    term = colnames(inputs$design),
    estimate = unname(step2$coef),
    se = unname(standard_error),
    ci_low = unname(step2$coef - stats::qnorm(.975) * standard_error),
    ci_high = unname(step2$coef + stats::qnorm(.975) * standard_error),
    check.names = FALSE
  )
  mmom_fit(
    converged = TRUE,
    status_code = 0L,
    message = "ok",
    coefficients = coef_table,
    vcov = coef_vcov,
    preliminary_coefficients = step1$coef0,
    corrected_design_moment = step2$coef_moment,
    analysis_weights = step1$case_weight,
    diagnostics = list(
      nobs = inputs$nobs,
      npredictors = inputs$npredictors,
      backend = inputs$backend,
      lambda1 = step1$lambda1,
      lambda2 = step2$lambda2,
      correction1 = step1$correction1,
      correction3 = step2$correction3,
      resid_var = step1$resid_var,
      preliminary_moment = inputs$preliminary,
      variance_bread = inputs$variance,
      steps = list(
        step1 = step1$diagnostics,
        step2 = step2$diagnostics,
        step3 = list(
          step = "step3",
          bread_solver = bread_diagnostics,
          covariance_dimension = nrow(coef_vcov)
        )
      )
    )
  )
}

#' Extract likelihood-only random-effect scores.
#'
#' The direct method computes the closed-form GLS score for each cluster.
#'
#' @param mlm An `lme4` mixed-model fit.
#' @param method Scoring method. `"direct"` is currently implemented.
#' @return A list containing cluster-level scores and covariance blocks.
#' @export
ranef_like <- function(mlm, method = c("direct", "unshrink")) {
  method <- match.arg(method)
  if (identical(method, "unshrink")) {
    repair_error("`method = \"unshrink\"` is not implemented yet.")
  }
  if (!inherits(mlm, "merMod")) {
    repair_error("`mlm` must be an lme4 mixed-model fit.")
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    repair_error("The `lme4` package is required for direct scoring.")
  }

  # Extract fixed effects, random-effect designs, and cluster membership.
  random_effects <- lme4::ranef(mlm)
  group <- names(random_effects)[[1L]]
  grouping <- as.character(lme4::getME(mlm, "flist")[[group]])
  cluster_ids <- rownames(random_effects[[group]])
  cluster_rows <- split(seq_along(grouping), grouping, drop = TRUE)[cluster_ids]
  z_terms <- lapply(lme4::getME(mlm, "Ztlist"), as.matrix)
  fixed_design <- lme4::getME(mlm, "X")
  fixed_coefs <- lme4::fixef(mlm)
  response <- stats::model.response(stats::model.frame(mlm))
  resid_var <- stats::sigma(mlm)^2
  if (!is.finite(resid_var) || resid_var <= 0) {
    repair_error(
      "Direct scoring requires a finite, positive residual variance."
    )
  }
  random_names <- colnames(random_effects[[group]])
  nrandom <- length(random_names)
  scores <- matrix(
    NA_real_,
    length(cluster_ids),
    nrandom,
    dimnames = list(cluster_ids, random_names)
  )
  err_cov <- array(
    NA_real_,
    c(nrandom, nrandom, length(cluster_ids)),
    dimnames = list(random_names, random_names, cluster_ids)
  )
  score_diagnostics <- vector("list", length(cluster_ids))

  # Apply the closed-form GLS score independently within each cluster.
  for (cluster_pos in seq_along(cluster_ids)) {
    rows <- cluster_rows[[cluster_ids[[cluster_pos]]]]
    random_design <- do.call(
      cbind,
      lapply(z_terms, function(term) {
        as.numeric(term[cluster_pos, rows, drop = TRUE])
      })
    )
    residual <- response[rows] - as.vector(
      fixed_design[rows, , drop = FALSE] %*% fixed_coefs
    )

    score <- tryCatch(
      {
        information <- crossprod(random_design) / resid_var
        score_rhs <- crossprod(random_design, residual) / resid_var
        chol_result <- mmom_chol_solve(information, score_rhs)
        list(
          score = as.vector(chol_result$solution),
          covariance = chol_result$inverse,
          diagnostics = chol_result$diagnostics
        )
      }, error = function(error) {
        list(
          score = rep(NA_real_, nrandom),
          covariance = matrix(NA_real_, nrandom, nrandom),
          diagnostics = list(
            solver = "chol",
            rank = NA_integer_,
            dimension = nrandom,
            condition_number = Inf,
            error = conditionMessage(error)
          )
        )
      }
    )
    scores[cluster_pos, ] <- score$score
    err_cov[, , cluster_pos] <- score$covariance
    score_diagnostics[[cluster_pos]] <- score$diagnostics
  }
  list(
    group = group,
    levels = cluster_ids,
    estimates = scores,
    err_cov = err_cov,
    score_cov = err_cov,
    method = method,
    diagnostics = list(
      step = "ranef_like",
      solver = "chol",
      clusters = stats::setNames(score_diagnostics, cluster_ids)
    )
  )
}
