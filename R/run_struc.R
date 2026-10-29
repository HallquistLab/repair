#' Assemble and inspect structural analysis inputs
#'
#' Match likelihood-only scores and their full within-person error covariance
#' to a person-level data frame by ID. The data frame defines the analysis cohort
#' and its order. No aggregation, row-order matching, or silent omission occurs.
#' @param formula An additive structural formula with at least one [re()] marker.
#'   Observed predictors may be numeric, logical, or factors. The outcome must
#'   be numeric or a random coefficient. An intercept is always included.
#' @param models A fitted Gaussian `lme4::lmer()` model or named list. All effects
#'   in an equation must come from one model with one grouping factor.
#' @param data A data frame with one row per analysis ID, containing the observed
#'   structural variables. Its ID values must match the model's grouping levels.
#' @param id Name of the ID column in `data`; need not equal the grouping name.
#' @param na_action `"fail"` (default) or `"omit"`. Omission covers missing or
#'   nonfinite observed values, unmatched IDs, and failed score blocks; it warns
#'   and records every excluded ID and reason.
#' @return An S7 analysis-data object. Use `@ids`, `@outcome`, `@predictors`,
#'   `@err_x`, `@err_y`, `@err_xy`, `@effects`, and `@diagnostics` to inspect it.
#'   Error arrays are ordered predictor by predictor by ID. Observed variables
#'   have zero sampling error. Scores are deviations from fitted fixed effects.
#' @export
repair_data <- function(formula, models, data, id, na_action = c("fail", "omit")) {
  na_action <- match.arg(na_action)
  models <- normalize_repair_models(models)
  parsed <- parse_struc_formula(formula, repair_effects(models))
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      anyNA(names(data)) || any(!nzchar(names(data)))) {
    repair_error("`data` must be a nonempty data frame with unique column names.")
  }
  if (!is.character(id) || length(id) != 1L || is.na(id) || !id %in% names(data)) {
    repair_error("`id` must name one column in `data`.")
  }
  if (any(startsWith(names(data), ".repair_"))) repair_error("Names beginning with `.repair_` are reserved.")
  if (!(is.character(data[[id]]) || is.factor(data[[id]]) || is.numeric(data[[id]]))) {
    repair_error("IDs must be character, factor, or numeric values.")
  }
  ids <- as.character(data[[id]])
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids) ||
      (is.numeric(data[[id]]) && any(!is.finite(data[[id]])))) {
    repair_error("IDs must be nonmissing, nonempty, finite, and unique: supply one row per person.")
  }
  missing_vars <- setdiff(parsed$observed, names(data))
  if (length(missing_vars)) repair_error("Observed variables missing from `data`: {.val {missing_vars}}.")
  if (id %in% parsed$observed) repair_error("The ID column cannot be a structural variable.")
  model <- models[[parsed$model]]
  if (length(lme4::getME(model, "flist")) != 1L) {
    repair_error(paste0("Structural fitting currently requires one grouping factor. ",
      "Multi-group scores can have cross-person errors that `fit_mmom()` cannot represent."))
  }
  group <- parsed$effects[[1L]]$group
  score <- repair_ranef_gls(model, target_group = group)
  positions <- match(ids, score@levels[[group]])
  reason <- rep("", length(ids))
  add_reason <- function(bad, label) {
    reason[bad] <<- ifelse(nzchar(reason[bad]), paste(reason[bad], label, sep = "; "), label)
  }
  add_reason(is.na(positions), "ID absent from measurement model")
  work <- data
  for (variable in parsed$observed) {
    value <- data[[variable]]
    if (!(is.numeric(value) || is.logical(value) || is.factor(value)) || !is.null(dim(value))) {
      repair_error("Observed structural variables must be numeric, logical, or factors: {.val {variable}}.")
    }
    bad <- if (is.numeric(value)) !is.finite(value) else is.na(value)
    add_reason(bad, paste("missing/nonfinite observed variable", variable))
  }
  for (effect in parsed$effects) {
    work[[effect$column]] <- score@predictions[[group]][positions, effect$term]
  }
  selected_terms <- vapply(parsed$effects, `[[`, character(1), "term")
  covariances <- score@error_covariance[[group]][selected_terms, selected_terms, positions, drop = FALSE]
  bad_scores <- !vapply(seq_along(ids), function(i) {
    all(is.finite(covariances[, , i])) &&
      all(vapply(parsed$effects, function(e) is.finite(work[[e$column]][i]), logical(1)))
  }, logical(1))
  add_reason(bad_scores & !is.na(positions), "likelihood-only scoring failed")
  keep <- !nzchar(reason)
  excluded <- data.frame(id = ids[!keep], reason = reason[!keep], row.names = NULL)
  if (any(!keep) && na_action == "fail") {
    examples <- paste(utils::head(paste(excluded$id, excluded$reason, sep = ": "), 5L), collapse = "; ")
    repair_error("{nrow(excluded)} analysis row(s) are incomplete ({examples}). Use `na_action = \"omit\"` to explicitly exclude and record them.")
  }
  if (!any(keep)) repair_error("No complete analysis rows remain.")
  if (any(!keep)) repair_warn("Omitting {sum(!keep)} analysis row(s); inspect `@diagnostics$excluded`.")
  work <- work[keep, , drop = FALSE]
  mf <- stats::model.frame(parsed$formula, data = work, na.action = stats::na.fail,
                          drop.unused.levels = TRUE)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || !is.null(dim(y))) repair_error("The structural outcome must be numeric.")
  design <- stats::model.matrix(parsed$formula, data = mf)
  x <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  n <- nrow(x)
  q <- ncol(x)
  if (n < max(8L, q + 2L)) repair_error("M-MOM requires at least max(8, npredictors + 2) complete analysis rows.")
  if (qr(design)$rank < ncol(design)) repair_error("The structural design is rank-deficient; remove redundant or constant predictors.")
  err_x <- array(0, c(q, q, n))
  err_y <- numeric(n)
  err_xy <- matrix(0, n, q)
  covariances <- covariances[, , keep, drop = FALSE]
  columns <- vapply(parsed$effects, `[[`, character(1), "column")
  predictor_effects <- match(colnames(x), columns)
  outcome_effect <- match(parsed$outcome, columns)
  for (i in seq_len(n)) {
    for (a in which(!is.na(predictor_effects))) {
      for (b in which(!is.na(predictor_effects))) {
        err_x[a, b, i] <- covariances[predictor_effects[a], predictor_effects[b], i]
      }
      if (!is.na(outcome_effect)) err_xy[i, a] <- covariances[outcome_effect, predictor_effects[a], i]
    }
    if (!is.na(outcome_effect)) err_y[i] <- covariances[outcome_effect, outcome_effect, i]
  }
  labels <- vapply(parsed$effects, `[[`, character(1), "label")
  replace <- !is.na(predictor_effects)
  colnames(x)[replace] <- labels[predictor_effects[replace]]
  rownames(x) <- ids[keep]
  dimnames(err_x) <- list(colnames(x), colnames(x), ids[keep])
  dimnames(err_xy) <- list(ids[keep], colnames(x))
  names(y) <- names(err_y) <- ids[keep]
  effects <- do.call(rbind, lapply(parsed$effects, function(e) {
    data.frame(model = e$model, group = e$group, term = e$term, label = e$label,
               role = if (e$column == parsed$outcome) "outcome" else "predictor")
  }))
  class_repair_data(
    formula = formula, ids = ids[keep], outcome = y, predictors = x,
    err_x = err_x, err_y = err_y, err_xy = err_xy, effects = effects,
    scores = score,
    diagnostics = list(n_requested = length(ids), n_analyzed = n,
      excluded = excluded, unused_score_ids = setdiff(score@levels[[group]], ids),
      id_column = id, group = group, model = parsed$model, na_action = na_action,
      score_semantics = "deviations", conditional_on_stage1 = TRUE,
      measurement_singular = lme4::isSingular(model),
      measurement_messages = model@optinfo$conv$lme4$messages %||% character())
  )
}

run_struc <- function(struc_data, regularization, alpha, backend) {
  control <- list(correction = regularization, alpha = alpha)
  fit_mmom(list(outcome = struc_data@outcome, predictors = struc_data@predictors,
    err_x = struc_data@err_x, err_y = struc_data@err_y, err_xy = struc_data@err_xy),
    preliminary_moment_control = control, final_moment_control = control, backend = backend)
}
