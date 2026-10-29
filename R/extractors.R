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
#'   and normal confidence intervals condition on fitted measurement-model quantities
#'   and any supplied cross-model error covariance. Changing `conf.level` changes
#'   interval endpoints only. Failed structural fits stop with an informative error;
#'   inspect [repair_diagnostics()] for their status and measurement diagnostics.
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

#' Extract selected likelihood-only scores as an ID-keyed table
#'
#' Read saved scores and their sampling-error variances without refitting or
#' navigating the different single- and multiple-model score containers.
#' @param object A [repair()] result or an assembled [repair_data()] object.
#'   A numerical structural failure is accepted because measurement scores exist.
#' @param cohort `"analyzed"` (default) returns retained structural IDs in their
#'   analysis order. `"measurement"` returns each model's effective grouping levels
#'   in saved scorer order, including people outside or excluded from the analysis.
#' @details Only effects selected in the structural formula are returned. Scores
#'   are likelihood-only random deviations, not BLUPs or fixed-plus-random totals.
#'   `std.error` is the square root of the scoring-error variance conditional on
#'   fitted measurement parameters; it is not a posterior SD or a structural SE.
#'   Marginal variances alone do not describe all score uncertainty: use
#'   [repair_score_covariance()] for retained within- and cross-model covariances.
#'
#'   A person's structural exclusion does not itself indicate a failed score.
#'   Failed scores remain as rows, with missing estimates/uncertainty where saved,
#'   `status = "failed"`, and a reason. Absent measurement IDs have no score row;
#'   use `repair_diagnostics(fit, table = "excluded")` for their exclusion reasons.
#' @return A plain data frame with columns in order: `id`, `model`, `group`, `term`,
#'   `label`, `role` (character), `score`, `std.error`, `error.variance` (numeric),
#'   `cohort` (character), `analyzed` (logical), `status`, `reason` (character).
#'   Rows follow selected-effect order (outcome first when random, then predictor
#'   references), then IDs within each effect. The key is `(model, group, term, id)`;
#'   IDs are character and never stored only in row names. `label` is the qualified
#'   formula label, and `role` is `"outcome"` or `"predictor"`. Successful rows have
#'   status `"ok"` and an empty reason. The schema is identical for one or many models.
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' head(repair_scores(fit))
#' head(repair_scores(fit, cohort = "measurement"))
#' @export
repair_scores <- function(object, cohort = c("analyzed", "measurement")) {
  a <- repair_analysis_object(object)
  cohort <- match.arg(cohort)
  scores <- if (length(a@diagnostics$model) == 1L)
    stats::setNames(list(a@scores), a@diagnostics$model) else a@scores
  rows <- lapply(seq_len(nrow(a@effects)), function(j) {
    e <- a@effects[j, ]; score <- scores[[e$model]]
    ids <- if (cohort == "analyzed") a@ids else score@levels[[e$group]]
    estimate <- as.numeric(score@predictions[[e$group]][ids, e$term])
    variance <- as.numeric(score@error_covariance[[e$group]][e$term, e$term, ids])
    ok <- is.finite(estimate) & is.finite(variance) & variance >= 0
    reason <- vapply(seq_along(ids), function(i) {
      if (ok[i]) "" else score@diagnostics[[e$group]]$per_level[[ids[i]]]$error %||%
        "Nonfinite score or invalid sampling-error variance."
    }, character(1))
    data.frame(id = as.character(ids), model = e$model, group = e$group, term = e$term,
      label = e$label, role = e$role, score = estimate,
      std.error = sqrt(ifelse(is.finite(variance) & variance >= 0, variance, NA_real_)),
      error.variance = variance, cohort = cohort, analyzed = ids %in% a@ids,
      status = ifelse(ok, "ok", "failed"), reason = reason, row.names = NULL)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Extract scoring-error covariance for the analyzed cohort
#'
#' Return a long table of the covariance actually assembled for selected random
#' coefficients, including cross-model blocks. This is distinct from `vcov(fit)`,
#' which describes uncertainty in structural regression coefficients.
#' @inheritParams repair_scores
#' @details Only analyzed people and formula-selected random coefficients are
#'   included. Measurement-cohort cross-model covariance may not be available,
#'   so this accessor has no measurement-cohort option. Observed structural columns
#'   have zero error by assumption and are not included. Within-model entries come
#'   from scoring; cross-model entries reflect the explicit policy or supplied blocks.
#'   Values are conditional on fitted measurement quantities and are not posterior
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
#' head(repair_score_covariance(fit))
#' # vcov(fit) instead describes structural coefficient uncertainty.
#' @export
repair_score_covariance <- function(object) {
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
