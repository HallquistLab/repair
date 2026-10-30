#' Assemble and inspect structural analysis inputs
#'
#' Match likelihood-only scores and their full within-person error covariance
#' to a person-level data frame by ID. The data frame defines the analysis cohort
#' and its order. No aggregation, row-order matching, or silent omission occurs.
#' @param formula An additive structural formula with at least one [re()] marker.
#'   Observed predictors may be numeric, logical, or factors. The outcome must
#'   be numeric or a random coefficient. An intercept is always included.
#' @param models A fitted Gaussian `lme4::lmer()` or linear `nlme::lme()`
#'   model, or named list. Select one target grouping factor per model. Gaussian
#'   lmer models may also contain nuisance groups whose every fitted level belongs
#'   to exactly one target level; nlme models must have one grouping factor.
#'   Grouping column names may differ, but their values must identify the same people.
#' @param data A data frame with one row per analysis ID, containing the observed
#'   structural variables. Its ID values must match the model's grouping levels.
#' @param id Name of the ID column in `data`; need not equal the grouping name.
#' @param na_action `"fail"` (default) or `"omit"`. Omission covers missing or
#'   nonfinite observed values, unmatched IDs, and failed score blocks; it warns
#'   and records every excluded ID and reason.
#' @param cross_model_error Cross-model score-error covariance policy. Default
#'   `NULL` is valid only when the formula references one model. For multiple
#'   models, use `"independent"` to explicitly assume zero cross-model scoring
#'   error covariance, or a list of blocks, one for every unordered model pair.
#'   Each block is `list(models = c("a", "b"), covariance = array(...))`.
#'   The numeric array must have three named dimensions: coefficient names from
#'   model a, coefficient names from model b, and person IDs. Term names must
#'   exactly cover the effects selected from those models. Every retained ID
#'   must be present; extra IDs are allowed. Order does not matter. Values are
#'   covariances of likelihood-only scoring errors, not latent coefficients or
#'   raw observations. The assembled covariance must be positive semidefinite.
#'   Missing/invalid supplied blocks are errors even with `na_action = "omit"`.
#' @details Use `re("term", model = "name", group = "person")` to resolve
#'   coefficients shared by several groups. All references to one model must
#'   select the same target group. Nesting is checked on all effective fitted
#'   observations before structural omission; shared nuisance levels are rejected
#'   even when only one of their target people is requested. Group names alone
#'   do not establish nesting. Nested random intercepts/slopes and multiple
#'   nuisance groups are allowed, provided every nuisance group passes the check.
#'   Target-specific scoring-error covariance includes every nuisance contribution.
#'   See `vignette("nested-nuisance", package = "repair")` for a complete example.
#'
#'   Within-model covariance always comes from the scorer and cannot be
#'   overridden. Supplied blocks are conditional on the fitted measurement models;
#'   uncertainty in their estimation is not propagated. Independent errors do not
#'   imply independent latent processes. The package does not estimate cross-model
#'   covariance from separate fits. See `vignette("multiple-models", package = "repair")`.
#' @return An S7 analysis-data object. Use `@ids`, `@outcome`, `@predictors`,
#'   `@err_x`, `@err_y`, `@err_xy`, `@effects`, and `@diagnostics` to inspect it.
#'   Error arrays are ordered predictor by predictor by ID. Observed variables
#'   have zero sampling error. Scores are deviations from fitted fixed effects.
#'   `@scores` holds a score object for one referenced model, or a named list for
#'   multiple models. `@diagnostics$measurement_models` always provides diagnostics
#'   by model, including the verified `nuisance_groups`;
#'   `@diagnostics$cross_model_error` records the applied policy and blocks.
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' analysis <- repair_data(re("day") ~ practice, growth, repair_example$people, "person")
#' analysis
#' head(repair_scores(analysis))
#' head(repair_score_covariance(analysis))
#' @export
repair_data <- function(formula, models, data, id, na_action = c("fail", "omit"),
                        cross_model_error = NULL) {
  repair_assemble(formula, models, data, id, na_action, cross_model_error)
}

repair_assemble <- function(formula, models, data, id, na_action = c("fail", "omit"),
                            cross_model_error = NULL, structural_contrasts = NULL) {
  na_action <- match.arg(na_action)
  models <- normalize_repair_models(models)
  parsed <- parse_struc_formula(formula, repair_effects(models))
  cross_error <- repair_cross_error_policy(cross_model_error, parsed)
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
  reason <- rep("", length(ids))
  add_reason <- function(bad, label) {
    reason[bad] <<- ifelse(nzchar(reason[bad]), paste(reason[bad], label, sep = "; "), label)
  }
  work <- data
  for (variable in parsed$observed) {
    value <- data[[variable]]
    if (!(is.numeric(value) || is.logical(value) || is.factor(value)) || !is.null(dim(value))) {
      repair_error("Observed structural variables must be numeric, logical, or factors: {.val {variable}}.")
    }
    bad <- if (is.numeric(value)) !is.finite(value) else is.na(value)
    add_reason(bad, paste("missing/nonfinite observed variable", variable))
  }
  scores <- model_diagnostics <- list()
  covariances <- array(0, c(length(parsed$effects), length(parsed$effects), length(ids)))
  effect_models <- vapply(parsed$effects, `[[`, character(1), "model")
  for (name in parsed$model) {
    model <- models[[name]]
    selected <- which(effect_models == name)
    effects <- parsed$effects[selected]
    group <- effects[[1L]]$group
    nuisance <- repair_validate_target(model, group, name)
    score <- repair_ranef_gls(model, target_group = group)
    scores[[name]] <- score
    positions <- match(ids, score@levels[[group]])
    add_reason(is.na(positions), paste("ID absent from measurement model", name))
    for (effect in effects) {
      work[[effect$column]] <- score@predictions[[group]][positions, effect$term]
    }
    terms <- vapply(effects, `[[`, character(1), "term")
    block <- score@error_covariance[[group]][terms, terms, positions, drop = FALSE]
    covariances[selected, selected, ] <- block
    bad <- !vapply(seq_along(ids), function(i) {
      all(is.finite(block[, , i])) &&
        all(vapply(effects, function(e) is.finite(work[[e$column]][i]), logical(1)))
    }, logical(1))
    add_reason(bad & !is.na(positions), paste("likelihood-only scoring failed in model", name))
    model_diagnostics[[name]] <- c(list(group = group, nuisance_groups = nuisance,
      unused_score_ids = setdiff(score@levels[[group]], ids)), repair_measurement_diagnostics(model))
  }
  keep <- !nzchar(reason)
  excluded <- data.frame(id = ids[!keep], reason = reason[!keep], row.names = NULL)
  if (any(!keep) && na_action == "fail") {
    examples <- paste(utils::head(paste(excluded$id, excluded$reason, sep = ": "), 5L), collapse = "; ")
    repair_error("{nrow(excluded)} analysis row(s) are incomplete ({examples}). Use `na_action = \"omit\"` to explicitly exclude and record them.")
  }
  if (!any(keep)) repair_error("No complete analysis rows remain.")
  covariances <- repair_apply_cross_error(covariances[, , keep, drop = FALSE],
                                         cross_error, parsed, ids[keep])
  if (any(!keep)) repair_warn("Omitting {sum(!keep)} analysis row(s); inspect `@diagnostics$excluded`.")
  work <- work[keep, , drop = FALSE]
  mf <- stats::model.frame(parsed$formula, data = work, na.action = stats::na.fail,
                          drop.unused.levels = TRUE)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || !is.null(dim(y))) repair_error("The structural outcome must be numeric.")
  design <- stats::model.matrix(parsed$formula, data = mf, contrasts.arg = structural_contrasts)
  x <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  n <- nrow(x)
  q <- ncol(x)
  if (n < max(8L, q + 2L)) repair_error("M-MOM requires at least max(8, npredictors + 2) complete analysis rows.")
  if (qr(design)$rank < ncol(design)) repair_error("The structural design is rank-deficient; remove redundant or constant predictors.")
  err_x <- array(0, c(q, q, n))
  err_y <- numeric(n)
  err_xy <- matrix(0, n, q)
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
    scores = if (length(scores) == 1L) scores[[1L]] else scores,
    diagnostics = c(list(n_requested = length(ids), n_analyzed = n,
      excluded = excluded,
      unused_score_ids = if (length(scores) == 1L) model_diagnostics[[1L]]$unused_score_ids else
        lapply(model_diagnostics, `[[`, "unused_score_ids"),
      structural_factors = repair_refit_factors(mf[intersect(parsed$observed, names(mf))], attr(design, "contrasts"), strict = FALSE),
      id_column = id, group = if (length(scores) == 1L) model_diagnostics[[1L]]$group else
        vapply(model_diagnostics, `[[`, character(1), "group"),
      model = parsed$model, na_action = na_action, measurement_models = model_diagnostics,
      cross_model_error = cross_error,
      score_semantics = "deviations", conditional_on_stage1 = TRUE),
      if (length(scores) == 1L) repair_measurement_diagnostics(models[[parsed$model]]) else list(
        measurement_backend = vapply(model_diagnostics, `[[`, character(1), "measurement_backend"),
        measurement_singular = vapply(model_diagnostics, `[[`, logical(1), "measurement_singular"),
        measurement_messages = unlist(lapply(model_diagnostics, `[[`, "measurement_messages"))))
  )
}

run_struc <- function(struc_data, regularization, alpha, backend) {
  control <- list(correction = regularization, alpha = alpha)
  fit_mmom(list(outcome = struc_data@outcome, predictors = struc_data@predictors,
    err_x = struc_data@err_x, err_y = struc_data@err_y, err_xy = struc_data@err_xy),
    preliminary_moment_control = control, final_moment_control = control, backend = backend)
}
