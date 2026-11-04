#' Specify a first-stage multilevel model for REPAIR
#'
#' Describe a Gaussian MLM to fit when [repair()], [repair_data()], or
#' [repair_effects()] consumes it. Construction does not fit a model.
#' @param formula A two-sided formula. For `engine = "lmer"`, include random
#'   effects using lme4 syntax. For `engine = "lme"`, this is the fixed-effects
#'   formula; supply `random` separately in `...`.
#' @param data A nonempty repeated-observation data frame. This is separate from
#'   the one-row-per-unit structural `data` supplied to [repair()].
#' @param engine Native fitter: `"lmer"` (default) or `"lme"`.
#' @param ... Fully named arguments to the selected native fitter, such as
#'   `REML`, `method`, `random`, `control`, `contrasts`, `weights`, `offset`,
#'   `correlation`, `subset`, or `na.action`, as supported by that fitter.
#'   For lmer, `weights`, `offset`, and `subset` use native data expressions;
#'   for lme, `subset` does. Other arguments are evaluated at construction.
#' @details Pass a specification directly as `models`, or combine specifications
#'   and fitted models in a uniquely named list. A bare specification is named
#'   `mlm`. Each consuming call fits its specifications once. Reuse the fitted
#'   models in a result's `@mlm_fit` to avoid fitting again during discovery or
#'   another structural analysis.
#'
#'   Fitting uses [lme4::lmer()] or [nlme::lme()] with native defaults, controls,
#'   warnings and errors; no optimizer retries or model selection are added.
#'   Formula and data are retained in the specification. For lme, `random` is
#'   required and `keep.data = FALSE` is rejected because prediction extraction
#'   needs the data. The optional nlme package is required for that engine.
#'   Existing Gaussian-model, grouping, residual-structure, structural-scope and
#'   cross-model covariance restrictions still apply. This helper changes neither
#'   the structural estimator nor its inference.
#'
#'   Results contain ordinary fitted models with explicit fitting arguments and
#'   a retained data binding. [repair_refit()] and [repair_bootstrap()] still
#'   require the original repeated and structural data explicitly, and retain
#'   their single-model and reference-replay checks. Data-only expressions and
#'   supported formulas can use native replay. External vectors, custom formula
#'   dependencies and other unsupported recipes still require an explicit
#'   `refitter`. This helper does not turn an arbitrary R expression into a
#'   portable replay recipe.
#' @return An unfitted `repair_mlm_spec` specification, accepted wherever REPAIR
#'   accepts measurement models. Its print method displays the formula, engine,
#'   observation count and supplied argument names without fitting.
#' @examples
#' data(repair_example)
#' growth <- repair_mlm(score ~ day + (day | person),
#'   repair_example$repeated, REML = TRUE)
#' fit <- repair(re("day", "growth") ~ practice, list(growth = growth),
#'   repair_example$people, "person")
#' repair_coefficients(fit)
#' if (requireNamespace("nlme", quietly = TRUE)) {
#'   growth_lme <- repair_mlm(score ~ day, repair_example$repeated,
#'     engine = "lme", random = ~ day | person)
#'   repair_effects(list(growth = growth_lme))
#' }
#' @export
repair_mlm <- function(formula, data, engine = c("lmer", "lme"), ...) {
  engine <- match.arg(engine)
  if (!inherits(formula, "formula") || length(formula) != 3L)
    repair_error("`formula` must be a two-sided MLM formula.")
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      anyNA(names(data)) || any(!nzchar(names(data))))
    repair_error("MLM `data` must be a nonempty data frame with unique column names.")
  expressions <- as.list(substitute(list(...)))[-1L]
  if (length(expressions) && (is.null(names(expressions)) ||
      anyNA(names(expressions)) || any(!nzchar(names(expressions))) ||
      anyDuplicated(names(expressions))))
    repair_error("All native fitting arguments in `...` must be uniquely and fully named.")
  if (engine == "lme") repair_check_pkg("nlme")
  fitter <- if (engine == "lmer") lme4::lmer else nlme::lme
  allowed <- setdiff(names(formals(fitter)), c("formula", "fixed", "data", "..."))
  unknown <- setdiff(names(expressions), allowed)
  if (length(unknown)) repair_error("Unsupported {.val {engine}} fitting argument(s): {.val {unknown}}. Use full native argument names.")
  nse <- intersect(names(expressions), if (engine == "lmer") c("weights", "offset", "subset") else "subset")
  caller <- parent.frame()
  arguments <- lapply(expressions[setdiff(names(expressions), nse)], eval, envir = caller)
  if (engine == "lmer" && "devFunOnly" %in% names(arguments) &&
      !identical(arguments$devFunOnly, FALSE))
    repair_error("The lmer engine requires `devFunOnly = FALSE` to return a fitted MLM.")
  if (engine == "lme") {
    if (is.null(arguments$random)) repair_error("Supply an explicit `random` argument for `engine = \"lme\"`.")
    if ("keep.data" %in% names(arguments) && !isTRUE(arguments$keep.data))
      repair_error("The lme engine requires `keep.data = TRUE`.")
    arguments$keep.data <- TRUE
  }
  structure(list(formula = formula, data = data, engine = engine,
    arguments = arguments, expressions = expressions[nse]), class = "repair_mlm_spec")
}

#' @export
#' @noRd
print.repair_mlm_spec <- function(x, ...) {
  cat("REPAIR first-stage MLM specification (", x$engine, ")\n", sep = "")
  print(x$formula)
  cat("Observations:", nrow(x$data), "\n")
  arguments <- c(names(x$arguments), names(x$expressions))
  if (length(arguments)) cat("Native arguments:", paste(arguments, collapse = ", "), "\n")
  invisible(x)
}

repair_fit_mlm_spec <- function(spec, name) {
  if (spec$engine == "lme") repair_check_pkg("nlme")
  env <- new.env(parent = environment(spec$formula))
  env$.repair_mlm_data <- spec$data
  formula <- spec$formula; environment(formula) <- env
  env$.repair_mlm_formula <- formula
  arguments <- list(.repair_mlm_formula = as.name(".repair_mlm_formula"))
  names(arguments) <- if (spec$engine == "lmer") "formula" else "fixed"
  arguments$data <- as.name(".repair_mlm_data")
  for (argument in names(spec$arguments)) {
    key <- paste0(".repair_mlm_arg_", argument)
    assign(key, spec$arguments[[argument]], envir = env)
    arguments[[argument]] <- as.name(key)
  }
  arguments <- c(arguments, spec$expressions)
  call <- as.call(c(list(if (spec$engine == "lmer") quote(lme4::lmer) else quote(nlme::lme)), arguments))
  model <- tryCatch(eval(call, env), error = function(e)
    repair_error("Fitting MLM {.val {name}} with {.val {spec$engine}} failed: {conditionMessage(e)}"))
  # Resolve helper-local symbols in the stored native call. Data expressions
  # remain expressions so replay can evaluate them on replacement observations.
  recorded <- stats::getCall(model)
  recorded[[1L]] <- if (spec$engine == "lmer") quote(lme4::lmer) else quote(nlme::lme)
  recorded[[if (spec$engine == "lmer") "formula" else "fixed"]] <- formula
  for (argument in names(spec$arguments)) recorded[argument] <- spec$arguments[argument]
  if (spec$engine == "lmer") model@call <- recorded else model$call <- recorded
  model
}
