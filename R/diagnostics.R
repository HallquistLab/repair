#' Inspect REPAIR sample, prediction extraction, and numerical diagnostics
#'
#' Assemble readable tables without refitting or changing the estimator.
#' @param object A fitted [repair()] result, including a numerical failure.
#' @param table Optional exact table name: `"status"`, `"scope"`, `"excluded"`,
#'   `"mlm"`, `"prediction_failures"`, `"prediction_error_variances"`, `"regularization"`,
#'   `"conditioning"`, `"weights"`, `"weight_summary"`, or `"reliability"`.
#'   Default NULL returns the existing list. A table name returns a plain data
#'   frame with automatic row names, including typed zero-row tables on failure.
#'   Reliability is computed only for NULL or `"reliability"`. The cross-model
#'   covariance policy is metadata in the list, not a selectable table.
#' @details Error variances describe likelihood-only predictions on the structural
#'   input scale; they are not posterior variances. Observed variables have zero
#'   error under the model's assumptions. Weight fractions are normalized inverse
#'   case-variance weights. `weight_equivalent_n = 1 / sum(fraction^2)` describes
#'   concentration only: it is not an inferential sample size or degrees of freedom.
#'
#'   Conditioning uses eigenvalues of the preliminary and final corrected predictor
#'   moments, including the error-free structural intercept. The reported spectral
#'   condition number depends on variable units; no automatic cutoff is applied.
#'   A nonpositive minimum eigenvalue yields an infinite condition number.
#'   The subtraction multiplier equals one for full measurement-error subtraction;
#'   it is not a measure of the amount of measurement error. `regularized` indicates
#'   that the multiplier differs from one beyond numerical tolerance (1e-8).
#'
#'   Prediction extraction failures include failures outside the structural scope and identify
#'   their MLM and inclusion status. Exclusions retain all reasons.
#'   On numerical structural failure, MLM diagnostics remain available,
#'   while unavailable structural tables have zero rows and `status` records why.
#' @return With `table = NULL`, a list with data frames `status`, `scope`, `excluded`, `mlm`,
#'   `prediction_failures`, `prediction_error_variances`, `regularization`, `conditioning`,
#'   `weights`, `weight_summary`, and `reliability`. Reliability uses the analyzed
#'   scope and the definitions in [repair_reliability()]. `cross_model_error`
#'   records the covariance policy and original supplied blocks.
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' repair_diagnostics(fit, table = "scope")
#' repair_diagnostics(fit, table = "weight_summary")
#' repair_diagnostics(fit, table = "excluded")
#' @export
repair_diagnostics <- function(object, table = NULL) {
  repair_require_model(object)
  tables <- c("status", "scope", "excluded", "mlm", "prediction_failures",
    "prediction_error_variances", "regularization", "conditioning", "weights", "weight_summary", "reliability")
  if (!is.null(table) && (!is.character(table) || length(table) != 1L || is.na(table) || !table %in% tables)) {
    repair_error("`table` must be NULL or one exact diagnostic table name: {.val {tables}}.")
  }
  if (identical(table, "reliability")) {
    result <- as.data.frame(repair_reliability(object))
    rownames(result) <- NULL
    return(result)
  }
  a <- object@analysis_data
  d <- a@diagnostics
  predictions <- if (length(d$model) == 1L) stats::setNames(list(a@predictions), d$model) else a@predictions
  failures <- data.frame(model = character(), group = character(), id = character(),
    analyzed = logical(), reason = character())
  for (name in names(predictions)) {
    group <- d$mlm_models[[name]]$group
    details <- predictions[[name]]@diagnostics[[group]]
    for (id in details$failed_levels) {
      failures <- rbind(failures, data.frame(model = name, group = group, id = id,
        analyzed = id %in% a@ids, reason = details$per_level[[id]]$error))
    }
  }
  mlm <- do.call(rbind, lapply(names(d$mlm_models), function(name) {
    m <- d$mlm_models[[name]]
    data.frame(model = name, group = m$group, backend = m$mlm_backend,
      singular = m$mlm_singular, n_levels = length(predictions[[name]]@levels[[m$group]]),
      n_outside_requested = length(m$unused_prediction_ids),
      messages = paste(m$mlm_messages, collapse = "; "))
  }))
  outcome_label <- paste(deparse(a@formula[[2L]]), collapse = "")
  errors <- data.frame(id = a@ids, term = outcome_label, role = "outcome", variance = a@err_y,
                       row.names = NULL)
  for (j in seq_len(ncol(a@predictors))) {
    errors <- rbind(errors, data.frame(id = a@ids, term = colnames(a@predictors)[j],
      role = "predictor", variance = a@err_x[j, j, ], row.names = NULL))
  }
  result <- c(list(status = data.frame(converged = isTRUE(object@struc_fit@converged),
       step = object@struc_fit@diagnostics$step %||% "complete", message = object@struc_fit@message),
    scope = data.frame(n_requested = d$n_requested, n_analyzed = length(a@ids),
                        n_excluded = nrow(d$excluded)),
    excluded = d$excluded, mlm = mlm, prediction_failures = failures,
    prediction_error_variances = errors), repair_structural_diagnostics(object),
    list(cross_model_error = d$cross_model_error))
  if (is.null(table)) {
    # Preserve the existing list's names, order, and reliability class.
    return(c(result[setdiff(names(result), "cross_model_error")],
      list(reliability = repair_reliability(object), cross_model_error = d$cross_model_error)))
  }
  result <- as.data.frame(result[[table]])
  rownames(result) <- NULL
  result
}

repair_structural_diagnostics <- function(object) {
  fit <- object@struc_fit; a <- object@analysis_data
  regularization <- data.frame(step = character(), policy = character(), alpha = numeric(),
    subtraction_multiplier = numeric(), regularized = logical())
  conditioning <- data.frame(step = character(), dimension = integer(), min_eigenvalue = numeric(),
    max_eigenvalue = numeric(), condition_number = numeric())
  weights <- data.frame(id = character(), weight = numeric(), fraction = numeric())
  weight_summary <- data.frame(n_units = integer(), max_fraction = numeric(),
                               weight_equivalent_n = numeric())
  if (isTRUE(fit@converged)) {
    design <- cbind(1, a@predictors)
    preliminary <- crossprod(design) - fit@steps$step1@correction_prelim * fit@steps$step1@err_cov_design
    matrices <- list(preliminary = preliminary, final = fit@corrected_design_moment)
    for (step in names(matrices)) {
      control <- fit@diagnostics[[paste0(step, "_moment_control")]]
      multiplier <- if (step == "preliminary") fit@diagnostics$correction_prelim else fit@diagnostics$correction_final
      regularization <- rbind(regularization, data.frame(step = step, policy = control$correction,
        alpha = control$alpha %||% NA_real_, subtraction_multiplier = multiplier,
        regularized = abs(multiplier - 1) > 1e-8))
      m <- matrices[[step]]
      values <- eigen((m + t(m)) / 2, symmetric = TRUE, only.values = TRUE)$values
      conditioning <- rbind(conditioning, data.frame(step = step, dimension = nrow(m),
        min_eigenvalue = min(values), max_eigenvalue = max(values),
        condition_number = if (min(values) > 0) max(values) / min(values) else Inf))
    }
    w <- as.numeric(fit@analysis_weights)
    fraction <- w / max(w); fraction <- fraction / sum(fraction)
    weights <- data.frame(id = a@ids, weight = w, fraction = fraction)
    weight_summary <- data.frame(n_units = length(w), max_fraction = max(fraction),
                                 weight_equivalent_n = 1 / sum(fraction^2))
  }
  list(regularization = regularization, conditioning = conditioning,
       weights = weights, weight_summary = weight_summary)
}
