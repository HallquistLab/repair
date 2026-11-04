# Fitted MLMs are supplied explicitly, preserving their controls,
# contrasts, weights, offsets, and effective observations.
normalize_repair_models <- function(models) {
  if (inherits(models, "brmsfit") || inherits(models, "stanreg"))
    repair_error("Bayesian fits require explicit opt-in with `repair_bayes_plugin()`; see its supported scope and inference limits.")
  if (inherits(models, "lmerMod") || inherits(models, "lme") || repair_is_bayes_plugin(models)) models <- list(mlm = models)
  if (!is.list(models) || !length(models) || is.null(names(models)) ||
      anyNA(names(models)) || any(!nzchar(names(models))) || anyDuplicated(names(models))) {
    repair_error("`models` must be a fitted lmer/lme model, a `repair_bayes_plugin()` snapshot, or a uniquely named list of these.")
  }
  if (!all(vapply(models, function(x) inherits(x, "lmerMod") ||
      (inherits(x, "lme") && !inherits(x, "nlme")) || repair_is_bayes_plugin(x), logical(1)))) {
    repair_error("Supply Gaussian `lme4::lmer()` or linear `nlme::lme()` models, or explicitly opt in with `repair_bayes_plugin()` for supported Bayesian fits.")
  }
  models
}

#' Inspect available random coefficients
#' @param models A fitted Gaussian `lme4::lmer()` or linear `nlme::lme()`
#'   model, an experimental [repair_bayes_plugin()] snapshot, or a named list.
#' @return A data frame with `model`, `group`, `term`, and `n_levels` columns.
#'   A bare model is named `mlm`. Terms are scalar coefficient names.
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
    } else if (repair_is_bayes_plugin(model)) {
      groups <- model$elements@ranef
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

repair_mlm_diagnostics <- function(model) {
  if (repair_is_bayes_plugin(model)) return(model$diagnostics)
  if (inherits(model, "lmerMod")) {
    return(list(mlm_backend = "lme4", mlm_singular = lme4::isSingular(model),
      mlm_messages = model@optinfo$conv$lme4$messages %||% character()))
  }
  list(mlm_backend = "nlme", mlm_singular = NA,
    mlm_messages = if (is.character(model$apVar))
      paste("nlme variance-parameter uncertainty:", model$apVar) else character())
}

# Qualification uses all effective fitted rows, before structural scope selection.
# A name such as person:day is neither required nor sufficient: actual level
# membership determines whether any nuisance effect is shared between targets.
repair_validate_target <- function(model, group, name) {
  groups <- if (repair_is_bayes_plugin(model)) model$elements@grouping else
    if (inherits(model, "lmerMod")) lme4::getME(model, "flist") else model$groups
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
        "Cross-unit prediction-error covariance is not supported, even with `na_action = \"omit\"`."))
    }
  }
  nuisance
}
