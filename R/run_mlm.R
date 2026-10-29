# Fitted measurement models are supplied explicitly, preserving their controls,
# contrasts, weights, offsets, and effective observations.
normalize_repair_models <- function(models) {
  if (inherits(models, "lmerMod")) models <- list(measurement = models)
  if (!is.list(models) || !length(models) || is.null(names(models)) ||
      anyNA(names(models)) || any(!nzchar(names(models))) || anyDuplicated(names(models))) {
    repair_error("`models` must be a fitted lmer model or a uniquely named list of fitted lmer models.")
  }
  if (!all(vapply(models, inherits, logical(1), "lmerMod"))) {
    repair_error("The structural workflow currently supports Gaussian `lme4::lmer()` models only.")
  }
  models
}

#' Inspect available random coefficients
#' @param models A fitted Gaussian `lme4::lmer()` model or named list of models.
#' @return A data frame with `model`, `group`, `term`, and `n_levels` columns.
#'   A bare model is named `measurement`. Terms are scalar coefficient names.
#' @export
repair_effects <- function(models) {
  models <- normalize_repair_models(models)
  do.call(rbind, unname(lapply(names(models), function(name) {
    groups <- lme4::ranef(models[[name]], condVar = FALSE)
    do.call(rbind, lapply(names(groups), function(group) {
      data.frame(model = name, group = group, term = colnames(groups[[group]]),
                 n_levels = nrow(groups[[group]]), row.names = NULL)
    }))
  })))
}
