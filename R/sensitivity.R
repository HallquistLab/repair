#' Compare explicit cross-model error assumptions
#'
#' Run the same structural analysis under named covariance scenarios. Measurement
#' models are fitted at most once per call. Each scenario uses [repair()] and its
#' existing estimator, regularization, and conditional inference unchanged.
#' @inheritParams repair
#' @param scenarios Nonempty uniquely named list. Each entry is a
#'   `cross_model_error` value: `"independent"` or complete prediction-error
#'   covariance blocks, for example from [repair_cross_covariance()].
#' @param conf.level Confidence level for successful coefficient tables.
#' @details Scenarios are assumptions supplied by the analyst, not estimates,
#'   model selection candidates, confidence limits for a correlation, or a
#'   probability distribution. Choose a scientifically defensible range before
#'   inspecting results and report all scenarios. Intervals are conditional on
#'   each scenario and fitted measurement parameters; they do not include
#'   uncertainty in cross-model covariance or replace bootstrap inference.
#'
#'   Common formula, model, scope, and fitting-argument errors stop the call.
#'   Scenario-specific validation errors (including inadmissible covariance)
#'   are recorded as `"error"`; numerical structural failures as `"failed"`.
#'   Neither yields coefficient rows. Successful minimum-alpha fallbacks remain
#'   `"ok"` with explicit fallback flags. Scenario warnings are collected in the
#'   summary rather than signaled repeatedly. No scenario is silently discarded.
#'   Covariance construction errors occur before this function is called: retain
#'   those separately when generating a grid, as shown in the multiple-model guide.
#' @return A plain list with `scenarios` (one summary row per assumption),
#'   `coefficients` (successful coefficient tables prefixed by `scenario`),
#'   `fits` (named results, NULL for errors), `assumptions` (original scenario
#'   inputs), and `conf.level`. Summary columns include `status`, `message`,
#'   `warnings`, `nobs`, selected preliminary/final alphas, and step-up fallback
#'   flags (NA when not applicable). Full diagnostic histories are in `fits`.
#' @seealso [repair_cross_covariance()], `vignette("multiple-models")`
#' @examples
#' if (requireNamespace("nlme", quietly = TRUE)) {
#'   case <- readRDS(system.file("extdata", "multimodel-workflow-reference.rds",
#'                              package = "repair"))$cases$supplied
#'   models <- list(growth = nlme::lme(y ~ z, random = ~ z | id, data = case$a),
#'     response = nlme::lme(q ~ w, random = ~ w | subject, data = case$b))
#'   result <- repair_sensitivity(re(w, "response") ~ re(z, "growth"),
#'     models, case$people, "person", scenarios = list(
#'       independent = "independent",
#'       assumed_nonzero = list(list(models = c("growth", "response"),
#'         covariance = case$cross_model_error[[1]]$covariance["z", , , drop = FALSE]))))
#'   result$scenarios
#'   result$coefficients
#' }
#' @export
repair_sensitivity <- function(formula, models, data, id, scenarios,
                                na_action = c("fail", "omit"),
                                regularization = c("minimum", "none", "step_up", "custom"),
                                alpha = NULL, backend = c("auto", "r", "rcpp"),
                                step_up_control = list(), conf.level = .95) {
  if (!is.list(scenarios) || !repair_cov_names(names(scenarios))) {
    repair_error("`scenarios` must be a nonempty uniquely named list.")
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    repair_error("`conf.level` must be strictly between zero and one.")
  }
  na_action <- match.arg(na_action); regularization <- match.arg(regularization)
  backend <- match.arg(backend)
  models <- normalize_repair_models(models)
  # Validate shared settings outside the scenario error handler. Use the same
  # moment-control validator as the estimator, without fitting a new estimator.
  if (regularization == "custom" && is.null(alpha)) repair_error("Supply `alpha` for custom regularization.")
  if (regularization != "custom" && !is.null(alpha)) repair_error("`alpha` is used only with custom regularization.")
  if (regularization != "step_up" && length(step_up_control)) repair_error("`step_up_control` is used only with step-up regularization.")
  baseline <- repair_data(formula, models, data, id, na_action,
    cross_model_error = "independent")
  control <- list(correction = regularization, alpha = alpha)
  if (regularization == "step_up") control$step_up <- step_up_control
  mmom_control(control, ncol(baseline@predictors), "moment_control")
  fits <- stats::setNames(vector("list", length(scenarios)), names(scenarios))
  rows <- tables <- vector("list", length(scenarios))
  for (j in seq_along(scenarios)) {
    warnings <- character(); error <- ""
    fit <- tryCatch(withCallingHandlers(
      repair(formula, models, data, id, na_action, regularization, alpha, backend,
        cross_model_error = scenarios[[j]], step_up_control = step_up_control),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) { error <<- conditionMessage(e); NULL })
    fits[j] <- list(fit)
    ok <- !is.null(fit) && isTRUE(fit@struc_fit@converged)
    d <- if (is.null(fit)) list() else fit@struc_fit@diagnostics
    searches <- d$step_up %||% list()
    rows[[j]] <- data.frame(scenario = names(scenarios)[j],
      status = if (is.null(fit)) "error" else if (ok) "ok" else "failed",
      message = if (is.null(fit)) error else fit@struc_fit@message,
      warnings = paste(unique(warnings), collapse = "\n"),
      nobs = if (is.null(fit)) NA_integer_ else length(fit@analysis_data@ids),
      alpha_preliminary = d$alpha_prelim %||% d$preliminary_moment_control$alpha %||% NA_real_,
      alpha_final = d$alpha_final %||% d$final_moment_control$alpha %||% NA_real_,
      fallback_preliminary = searches$preliminary$fallback %||% NA,
      fallback_final = searches$final$fallback %||% NA,
      row.names = NULL)
    if (ok) tables[[j]] <- cbind(data.frame(scenario = names(scenarios)[j]),
      repair_coefficients(fit, conf.level))
  }
  coefficients <- do.call(rbind, tables)
  if (is.null(coefficients)) coefficients <- data.frame(scenario = character(),
    term = character(), estimate = numeric(), std.error = numeric(), statistic = numeric(),
    p.value = numeric(), conf.low = numeric(), conf.high = numeric(), conf.level = numeric())
  rownames(coefficients) <- NULL
  list(scenarios = do.call(rbind, rows), coefficients = coefficients,
    fits = fits, assumptions = scenarios, conf.level = conf.level)
}
