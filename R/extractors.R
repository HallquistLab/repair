# Resolve public fit/assembly inputs in one place; extractors never rescore or refit.
repair_analysis_object <- function(object) {
  if (S7::S7_inherits(object, class_repair_mod)) return(object@analysis_data)
  if (S7::S7_inherits(object, class_repair_data)) return(object)
  repair_error("Supply a `repair()` result or an assembled `repair_data()` object.")
}

#' Extract a structural coefficient table
#'
#' Obtain an ordinary data frame for reporting, joins, and export without using
#' internal S7 properties. No fitting or reliability computation is performed.
#' @param object A fitted [repair()] result.
#' @param conf.level Confidence level strictly between zero and one; default 0.95.
#' @details Estimates retain their input scale. Standard errors, two-sided p-values,
#'   and normal confidence intervals condition on fitted MLM quantities
#'   and any supplied cross-model error covariance. Changing `conf.level` changes
#'   interval endpoints only. Failed structural fits stop with an informative error;
#'   inspect [repair_diagnostics()] for their status and MLM diagnostics.
#' @return A plain data frame with columns, in order: `term` (character),
#'   `estimate`, `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`, and
#'   `conf.level` (numeric). Rows follow `coef(object)`: the structural intercept,
#'   then model-matrix predictor columns. Factor contrasts retain their column
#'   names; random-coefficient predictors retain model-qualified labels. Row names
#'   are automatic, not identifiers. `vcov(object)` remains the public accessor
#'   for structural coefficient covariance.
#' @examples
#' growth <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
#' people <- data.frame(id = levels(lme4::sleepstudy$Subject), observed = 1:18)
#' fit <- repair(re("Days") ~ observed, growth, people, "id")
#' repair_coefficients(fit, conf.level = 0.90)
#' @export
repair_coefficients <- function(object, conf.level = 0.95) {
  repair_require_model(object)
  fit <- repair_require_fit(object)
  intervals <- stats::confint(object, level = conf.level)
  b <- fit@coefficients
  data.frame(term = as.character(b$term), estimate = as.numeric(b$estimate),
    std.error = as.numeric(b$se), statistic = as.numeric(b$zstat), p.value = as.numeric(b$pval),
    conf.low = as.numeric(intervals[, 1]), conf.high = as.numeric(intervals[, 2]),
    conf.level = rep(as.numeric(conf.level), nrow(b)), row.names = NULL)
}

#' Extract selected likelihood-only predictions as an ID-keyed table
#'
#' Read saved predictions and their sampling-error variances without refitting or
#' navigating the different single- and multiple-model prediction containers.
#' @param object A [repair()] result or an assembled [repair_data()] object.
#'   A numerical structural failure is accepted because MLM predictions exist.
#' @param scope `"structural"` (default) returns retained structural IDs in their
#'   analysis order. `"all"` returns each model's effective grouping levels
#'   in saved scorer order, including units outside or excluded from the analysis.
#' @details Only effects selected in the structural formula are returned. Predictions
#'   are likelihood-only random deviations, not BLUPs or fixed-plus-random totals.
#'   `std.error` is the square root of the prediction-error variance conditional on
#'   fitted MLM parameters; it is not a posterior SD or a structural SE.
#'   Marginal variances alone do not describe all prediction uncertainty: use
#'   [repair_prediction_covariance()] for retained within- and cross-model covariances.
#'
#'   A unit's structural exclusion does not itself indicate a failed prediction.
#'   Failed predictions remain as rows, with missing estimates/uncertainty where saved,
#'   `status = "failed"`, and a reason. Absent MLM IDs have no prediction row;
#'   use `repair_diagnostics(fit, table = "excluded")` for their exclusion reasons.
#' @return A plain data frame with columns in order: `id`, `model`, `group`, `term`,
#'   `label`, `role` (character), `prediction`, `std.error`, `prediction.error.variance` (numeric),
#'   `scope` (character), `analyzed` (logical), `status`, `reason` (character).
#'   Rows follow selected-effect order (outcome first when random, then predictor
#'   references), then IDs within each effect. The key is `(model, group, term, id)`;
#'   IDs are character and never stored only in row names. `label` is the qualified
#'   formula label, and `role` is `"outcome"` or `"predictor"`. Successful rows have
#'   status `"ok"` and an empty reason. The schema is identical for one or many models.
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' head(repair_predictions(fit))
#' head(repair_predictions(fit, scope = "all"))
#' @export
repair_predictions <- function(object, scope = c("structural", "all")) {
  a <- repair_analysis_object(object)
  scope <- match.arg(scope)
  predictions <- if (length(a@diagnostics$model) == 1L)
    stats::setNames(list(a@predictions), a@diagnostics$model) else a@predictions
  rows <- lapply(seq_len(nrow(a@effects)), function(j) {
    e <- a@effects[j, ]; prediction <- predictions[[e$model]]
    ids <- if (scope == "structural") a@ids else prediction@levels[[e$group]]
    estimate <- as.numeric(prediction@predictions[[e$group]][ids, e$term])
    variance <- as.numeric(prediction@prediction_error_covariance[[e$group]][e$term, e$term, ids])
    ok <- is.finite(estimate) & is.finite(variance) & variance >= 0
    reason <- vapply(seq_along(ids), function(i) {
      if (ok[i]) "" else prediction@diagnostics[[e$group]]$per_level[[ids[i]]]$error %||%
        "Nonfinite prediction or invalid sampling-error variance."
    }, character(1))
    data.frame(id = as.character(ids), model = e$model, group = e$group, term = e$term,
      label = e$label, role = e$role, prediction = estimate,
      std.error = sqrt(ifelse(is.finite(variance) & variance >= 0, variance, NA_real_)),
      prediction.error.variance = variance, scope = scope, analyzed = ids %in% a@ids,
      status = ifelse(ok, "ok", "failed"), reason = reason, row.names = NULL)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Extract prediction-error covariance for the structural scope
#'
#' Return a long table of the covariance actually assembled for selected random
#' coefficients, including cross-model blocks. This is distinct from `vcov(fit)`,
#' which describes uncertainty in structural regression coefficients.
#' @inheritParams repair_predictions
#' @details Only structural units and formula-selected random coefficients are
#'   included. MLM-scope cross-model covariance may not be available,
#'   so this accessor has no MLM-scope option. Observed structural columns
#'   have zero error by assumption and are not included. Within-model entries come
#'   from prediction extraction; cross-model entries reflect the explicit policy or supplied blocks.
#'   Values are conditional on fitted MLM quantities and are not posterior
#'   covariances. Extraction does not change inputs or assume new zero covariances.
#' @return A plain data frame with columns in order: `id`, `model1`, `group1`,
#'   `term1`, `label1`, `role1`, `model2`, `group2`, `term2`, `label2`, `role2`
#'   (character), `covariance` (numeric), and `source` (character). There is one row
#'   per unordered pair (including the diagonal) per analyzed ID. Pairs follow the
#'   lower triangle in selected-effect order: (1,1), (2,1), (2,2), etc.; IDs within
#'   each pair follow analysis order. Mirror off-diagonal rows to reconstruct a
#'   symmetric matrix. Source is `"within_model"`, `"assumed_independent"`, or
#'   `"supplied"`. The key is ID plus the two model/group/term identities.
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' head(repair_prediction_covariance(fit))
#' # vcov(fit) instead describes structural coefficient uncertainty.
#' @export
repair_prediction_covariance <- function(object) {
  a <- repair_analysis_object(object)
  e <- a@effects
  indices <- match(e$label, colnames(a@predictors))
  rows <- list()
  for (i in seq_len(nrow(e))) {
    for (j in seq_len(i)) {
      y1 <- e$role[i] == "outcome"; y2 <- e$role[j] == "outcome"
      covariance <- if (y1 && y2) a@err_y else if (y1) a@err_xy[, indices[j]] else
        if (y2) a@err_xy[, indices[i]] else a@err_x[indices[i], indices[j], ]
      source <- if (e$model[i] == e$model[j]) "within_model" else
        if (a@diagnostics$cross_model_error$policy == "independent") "assumed_independent" else "supplied"
      rows[[length(rows) + 1L]] <- data.frame(id = a@ids,
        model1 = e$model[i], group1 = e$group[i], term1 = e$term[i], label1 = e$label[i], role1 = e$role[i],
        model2 = e$model[j], group2 = e$group[j], term2 = e$term[j], label2 = e$label[j], role2 = e$role[j],
        covariance = as.numeric(covariance), source = source, row.names = NULL)
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
