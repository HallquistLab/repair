#' Fit a REPAIR structural regression from a measurement model
#'
#' Select random coefficients with [re()], match them to person-level data by
#' ID, and fit the modified method-of-moments (M-MOM) regression. Measurement
#' models are fitted beforehand so that their fitting choices remain explicit.
#' @inheritParams repair_data
#' @param regularization Moment regularization: `"minimum"` (default), `"none"`,
#'   or `"custom"`. `"none"` still subtracts measurement error; it disables only
#'   stabilization. The same setting is used for preliminary and final moments.
#' @param alpha Custom regularization constant, at least the number of structural
#'   predictor columns plus two. Required for `regularization = "custom"`.
#' @param backend M-MOM computation backend: `"auto"`, `"r"`, or `"rcpp"`.
#' @details The supported structural grammar is additive with an intercept.
#'   Effects on either side must come from one Gaussian lmer model with one
#'   grouping factor. Full within-person score-error covariance is retained.
#'   Observed variables are treated as measured without error. Scores are random
#'   deviations, not fixed-plus-random individual coefficients.
#'
#'   Standard errors condition on fitted measurement-model quantities and use
#'   the modified final predictor moment (the manuscript M-MOM calculation).
#'   They assume the model-based normal-error covariance; measurement-model
#'   parameter uncertainty is not propagated. Intervals and tests use normal
#'   reference distributions. Numerical fitting failures return a result with
#'   `@struc_fit@converged = FALSE` and issue a warning.
#' @return An S7 `class_repair_mod` with `@mlm_fit`, `@struc_fit`, and
#'   inspectable `@analysis_data`. Use `summary()`, `coef()`, `vcov()`,
#'   `confint()`, and `nobs()` for standard output.
#' @examples
#' growth <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
#' people <- data.frame(Subject = levels(lme4::sleepstudy$Subject),
#'                      observed = seq_len(18))
#' repair_effects(list(growth = growth))
#' fit <- repair(re("Days", model = "growth") ~ observed,
#'               models = list(growth = growth), data = people, id = "Subject")
#' summary(fit)
#' confint(fit)
#' @export
repair <- function(formula, models, data, id, na_action = c("fail", "omit"),
                   regularization = c("minimum", "none", "custom"),
                   alpha = NULL, backend = c("auto", "r", "rcpp")) {
  regularization <- match.arg(regularization)
  if (regularization == "custom" && is.null(alpha)) repair_error("Supply `alpha` for custom regularization.")
  if (regularization != "custom" && !is.null(alpha)) repair_error("`alpha` is used only with custom regularization.")
  models <- normalize_repair_models(models)
  analysis_data <- repair_data(formula, models, data, id, na_action)
  struc_fit <- run_struc(analysis_data, regularization, alpha, backend)
  if (!isTRUE(struc_fit@converged)) repair_warn("Structural fit failed: {struc_fit@message}")
  class_repair_mod(mlm_fit = models[analysis_data@diagnostics$model],
    struc_fit = struc_fit, analysis_data = analysis_data, call = match.call())
}
