# Fitted measurement models are supplied explicitly, preserving their controls,
# contrasts, weights, offsets, and effective observations.
normalize_repair_models <- function(models) {
  if (inherits(models, "lmerMod") || inherits(models, "lme")) models <- list(measurement = models)
  if (!is.list(models) || !length(models) || is.null(names(models)) ||
      anyNA(names(models)) || any(!nzchar(names(models))) || anyDuplicated(names(models))) {
    repair_error("`models` must be a fitted lmer/lme model or a uniquely named list of fitted lmer/lme models.")
  }
  if (!all(vapply(models, function(x) inherits(x, "lmerMod") ||
      (inherits(x, "lme") && !inherits(x, "nlme")), logical(1)))) {
    repair_error("The structural workflow currently supports Gaussian `lme4::lmer()` or linear `nlme::lme()` models.")
  }
  models
}

#' Inspect available random coefficients
#' @param models A fitted Gaussian `lme4::lmer()` or linear `nlme::lme()`
#'   model, or a named list of models.
#' @return A data frame with `model`, `group`, `term`, and `n_levels` columns.
#'   A bare model is named `measurement`. Terms are scalar coefficient names.
#'   The catalog includes nuisance groups; discovery does not establish eligibility
#'   for structural fitting. Select the target with the `group` argument of [re()].
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' repair_effects(list(growth = growth))
#' @export
repair_effects <- function(models) {
  models <- normalize_repair_models(models)
  do.call(rbind, unname(lapply(names(models), function(name) {
    model <- models[[name]]
    if (inherits(model, "lmerMod")) {
      groups <- lme4::ranef(model, condVar = FALSE)
    } else {
      repair_check_pkg("nlme")
      groups <- stats::setNames(lapply(seq_len(ncol(model$groups)), function(level) {
        nlme::ranef(model, level = level)
      }), names(model$groups))
    }
    do.call(rbind, lapply(names(groups), function(group) {
      data.frame(model = name, group = group, term = colnames(groups[[group]]),
                 n_levels = nrow(groups[[group]]), row.names = NULL)
    }))
  })))
}

repair_measurement_diagnostics <- function(model) {
  if (inherits(model, "lmerMod")) {
    return(list(measurement_backend = "lme4", measurement_singular = lme4::isSingular(model),
      measurement_messages = model@optinfo$conv$lme4$messages %||% character()))
  }
  list(measurement_backend = "nlme", measurement_singular = NA,
    measurement_messages = if (is.character(model$apVar))
      paste("nlme variance-parameter uncertainty:", model$apVar) else character())
}

# Qualification uses all effective fitted rows, before structural cohort selection.
# A name such as person:day is neither required nor sufficient: actual level
# membership determines whether any nuisance effect is shared between targets.
repair_validate_target <- function(model, group, name) {
  groups <- if (inherits(model, "lmerMod")) lme4::getME(model, "flist") else model$groups
  if (inherits(model, "lme") && length(groups) != 1L) {
    repair_error("nlme structural fitting requires exactly one grouping factor.")
  }
  nuisance <- setdiff(names(groups), group)
  target <- as.character(groups[[group]])
  for (other in nuisance) {
    owners <- split(target, as.character(groups[[other]]), drop = TRUE)
    shared <- names(owners)[vapply(owners, function(x) length(unique(x)) > 1L, logical(1))]
    if (length(shared)) {
      examples <- paste(utils::head(shared, 3L), collapse = ", ")
      repair_error(paste0("Model {.val {name}} has nuisance group {.val {other}} with levels ",
        "shared across target group {.val {group}} (e.g. {examples}). ",
        "Cross-person scoring-error covariance is not supported, even with `na_action = \"omit\"`."))
    }
  }
  nuisance
}
