# Keep the convergence policy shared by robust fitting and ordinary replay.
repair_measurement_converged <- function(model, warnings = character()) {
  codes <- if (inherits(model, "lmerMod"))
    c(unlist(model@optinfo$conv$opt), unlist(model@optinfo$conv$lme4$code)) else numeric()
  !anyNA(codes) && !any(codes != 0) &&
    !any(grepl("convergence error|iteration limit|maximum number of iterations|nlminb problem",
      warnings, ignore.case = TRUE))
}

repair_lmer_attempt <- function(args, label) {
  warnings <- messages <- character()
  model <- tryCatch(withCallingHandlers(do.call(lme4::lmer, args),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }),
    error = function(e) e)
  error <- inherits(model, "error")
  objective <- if (error) NA_real_ else -2 * as.numeric(stats::logLik(model))
  usable <- !error && is.finite(objective) && is.finite(stats::sigma(model)) && stats::sigma(model) > 0 &&
    all(is.finite(lme4::fixef(model))) && repair_measurement_converged(model, warnings)
  record <- data.frame(attempt = label, usable = usable, objective = objective,
    singular = if (error) NA else lme4::isSingular(model),
    reason = if (error) conditionMessage(model) else paste(c(model@optinfo$conv$lme4$messages, warnings), collapse = "; "),
    warnings = I(list(warnings)), messages = I(list(messages)), stringsAsFactors = FALSE)
  list(model = if (error) NULL else model, record = record)
}

repair_lmer_fit <- function(data, recipe) {
  if (!is.data.frame(data) || !nrow(data)) repair_error("Supply a nonempty measurement data frame.")
  data <- repair_refit_factor_data(data, recipe$factors)
  formula <- repair_replay_formula(recipe$formula, names(data))
  env <- new.env(parent = asNamespace("stats")); env$.repair_refit_data <- data
  environment(formula) <- env
  args <- list(formula = formula, data = data, REML = recipe$REML, na.action = stats::na.fail)
  for (name in c("weights", "offset")) if (!is.null(recipe[[name]])) {
    column <- recipe[[name]]
    if (!column %in% names(data)) repair_error("Measurement data lack {.val {column}}.")
    args[[name]] <- data[[column]]
  }
  if (length(recipe$contrasts)) args$contrasts <- recipe$contrasts
  args$control <- lme4::lmerControl(optimizer = "nloptwrap",
    optCtrl = list(maxeval = recipe$maxfun, xtol_abs = 1e-8, ftol_abs = 1e-8))
  attempts <- list(repair_lmer_attempt(args, "nlopt_start"))
  start <- attempts[[1L]]$model
  theta <- if (!is.null(start)) lme4::getME(start, "theta") else NULL
  control <- lme4::lmerControl(optimizer = "bobyqa",
    optCtrl = list(maxfun = recipe$maxfun, rhobeg = 2e-3, rhoend = 1e-8))
  args$control <- control
  args$start <- if (!is.null(theta) && all(is.finite(theta))) theta else NULL
  attempts[[2L]] <- repair_lmer_attempt(args, "bobyqa_warm")
  acceptable <- function() {
    records <- do.call(rbind, lapply(attempts, `[[`, "record"))
    best <- suppressWarnings(min(records$objective, na.rm = TRUE))
    which(seq_len(nrow(records)) > 1L & records$usable &
      is.finite(records$objective) & records$objective <= best + recipe$objective_tolerance)
  }
  if (!length(acceptable())) {
    args$start <- NULL
    attempts[[3L]] <- repair_lmer_attempt(args, "bobyqa_cold")
  }
  if (!length(acceptable()) && !is.null(theta) && all(is.finite(theta))) {
    lower <- lme4::getME(start, "lower")
    # Deterministic interior perturbation; no random starts or RNG consumption.
    args$start <- pmax(theta + .01 * pmax(1, abs(theta)), lower + .01)
    attempts[[length(attempts) + 1L]] <- repair_lmer_attempt(args, "bobyqa_perturbed")
  }
  records <- do.call(rbind, lapply(attempts, `[[`, "record")); rownames(records) <- NULL
  eligible <- acceptable(); usable <- length(eligible) > 0L
  if (!usable) eligible <- which(seq_len(nrow(records)) > 1L & is.finite(records$objective))
  if (!length(eligible)) {
    e <- simpleError("Robust lmer fitting produced no finite BOBYQA fit.")
    e$repair_optimization <- records; stop(e)
  }
  chosen <- eligible[which.min(records$objective[eligible])]
  model <- attempts[[chosen]]$model
  # Literal controls/columns and an explicit data binding allow safe saved replay.
  model@call[[1L]] <- quote(lme4::lmer); model@call$formula <- formula
  model@call$data <- as.name(".repair_refit_data"); model@call$control <- control
  model@call$start <- NULL
  for (name in c("weights", "offset")) model@call[name] <- list(
    if (!is.null(recipe[[name]])) as.name(recipe[[name]]) else NULL)
  coding <- attr(lme4::getME(model, "X"), "contrasts")
  frame <- stats::model.frame(model)
  recipe$factors <- repair_refit_factors(frame[intersect(names(coding), names(frame))], coding)
  attr(model, "repair_lmer_recipe") <- recipe
  attr(model, "repair_lmer_optimization") <- list(attempts = records, selected = chosen, usable = usable,
    objective_tolerance = recipe$objective_tolerance)
  if (!usable) warning("Robust lmer fitting remains unresolved; inspect attr(model, 'repair_lmer_optimization'). Refitting/bootstrap inference will reject this fit.", call. = FALSE)
  model
}

#' Create an audited, deterministic Gaussian lmer fitting recipe
#'
#' Use the returned function for the original measurement fit. Models retain the
#' recipe, so [repair_refit()] and [repair_bootstrap()] apply it automatically to
#' reference and replacement data, including after saving/reloading a model.
#'
#' The recipe fits nloptwrap to obtain a starting point, then BOBYQA with tighter
#' tolerances and a larger evaluation budget. If needed it retries BOBYQA from
#' its default start and a deterministic interior perturbation. A selected fit
#' must pass the ordinary convergence checks and have REML/ML deviance within
#' `objective_tolerance` of the best finite attempted value. Singular-but-usable
#' fits are retained. No convergence checks are disabled, no random starts are
#' drawn, and the bootstrap failure gate is unchanged.
#'
#' This changes the numerical fitting recipe. Refit the original measurement
#' model and rebuild `repair()` before bootstrapping; do not attach this callback
#' to an analysis fitted with a different recipe. Nonlinear formulas and implicit
#' omissions unsupported by native replay remain unsupported. For categorical
#' predictors supply explicit contrasts to insulate a callback from global options.
#'
#' @param formula A supported data-only lmer formula.
#' @param REML Logical; retain REML (default) or use ML consistently.
#' @param weights,offset Optional names of explicit numeric measurement columns.
#'   Formula offsets are also supported. Expressions/global vectors are not read.
#' @param contrasts Optional explicit contrast list passed to lmer.
#' @param maxfun Positive integer evaluation budget for each optimizer attempt.
#' @param objective_tolerance Nonnegative finite absolute deviance tolerance.
#' @return A function of one measurement data frame. It returns an lmer model
#'   with `repair_lmer_recipe` and `repair_lmer_optimization` attributes. The latter
#'   records every attempted fit, warnings/messages, objective, singularity and
#'   selection. If finite fits remain unresolved, the best BOBYQA model is returned
#'   with a warning and `usable = FALSE`; refitting/bootstrap inference rejects it.
#'   If no finite BOBYQA fit exists, an error retains the attempt history.
#' @examples
#' data(repair_example)
#' fitter <- repair_lmer_refitter(score ~ day + (day | person))
#' model <- fitter(repair_example$repeated)
#' fit <- repair(re("day") ~ practice, model, repair_example$people, "person")
#' replay <- repair_refit(fit, repair_example$repeated, repair_example$people)
#' attr(model, "repair_lmer_optimization")
#' @export
repair_lmer_refitter <- function(formula, REML = TRUE, weights = NULL, offset = NULL,
                                 contrasts = NULL, maxfun = 200000L, objective_tolerance = 1e-6) {
  if (!inherits(formula, "formula")) repair_error("`formula` must be a formula.")
  formula <- repair_replay_formula(formula, all.vars(formula))
  if (!is.logical(REML) || length(REML) != 1L || is.na(REML)) repair_error("`REML` must be TRUE or FALSE.")
  for (name in c("weights", "offset")) {
    value <- get(name)
    if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)))
      repair_error("`{name}` must name one explicit data column.")
  }
  if (!is.null(contrasts) && (!is.list(contrasts) || is.null(names(contrasts)))) repair_error("`contrasts` must be a named list.")
  if (!is.numeric(maxfun) || length(maxfun) != 1L || !is.finite(maxfun) || maxfun < 1 || maxfun %% 1 || maxfun > .Machine$integer.max)
    repair_error("`maxfun` must be a positive integer.")
  if (!is.numeric(objective_tolerance) || length(objective_tolerance) != 1L || !is.finite(objective_tolerance) || objective_tolerance < 0)
    repair_error("`objective_tolerance` must be nonnegative and finite.")
  recipe <- list(formula = formula, REML = REML, weights = weights, offset = offset,
    contrasts = contrasts, maxfun = as.integer(maxfun), objective_tolerance = objective_tolerance, factors = list())
  force(recipe)
  function(data) repair_lmer_fit(data, recipe)
}
