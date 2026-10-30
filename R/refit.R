# Replay expressions are interpreted against explicit columns and a small set of
# constructors. The original call's data/subset/environment are never evaluated.
repair_replay_formula <- function(expr, columns, dynamic = FALSE) {
  allowed <- c("~", "+", "-", "*", "/", "^", ":", "|", "||", "(",
    "I", "offset", "log", "log1p", "sqrt", "exp", "abs", "interaction")
  if (dynamic) allowed <- c(allowed, "fitted", "resid", "residuals")
  walk <- function(x) {
    if (is.call(x)) {
      if (!is.symbol(x[[1L]]) || !as.character(x[[1L]]) %in% allowed)
        repair_error("Replay requires data-only formulas with supported operators; supply an explicit `refitter` for other transformations.")
      lapply(as.list(x)[-1L], walk)
    }
    invisible(NULL)
  }
  walk(expr)
  if (length(setdiff(all.vars(expr), c(columns, if (dynamic) "."))))
    repair_error("A replay formula depends on columns absent from the explicit data. Supply those columns or an explicit `refitter`.")
  stats::as.formula(expr, env = asNamespace("stats"))
}

repair_replay_eval <- function(expr, data, literal = FALSE) {
  if (is.null(expr) || is.atomic(expr)) return(expr)
  if (is.symbol(expr)) {
    name <- as.character(expr)
    if (!literal && name %in% names(data)) return(data[[name]])
    repair_error("Replay argument {.val {name}} is not an explicit data column. Supply an explicit `refitter`; global bindings are not replayed.")
  }
  if (!is.call(expr)) repair_error("Replay needs a literal fitting recipe or an explicit `refitter`.")
  head <- expr[[1L]]
  if (identical(head, as.name("~"))) return(repair_replay_formula(expr, names(data), dynamic = TRUE))
  if (is.call(head) && identical(head[[1L]], as.name("::"))) {
    pkg <- as.character(head[[2L]]); name <- as.character(head[[3L]])
    if (!pkg %in% c("base", "stats", "lme4", "nlme")) repair_error("Unsupported replay constructor; supply an explicit `refitter`.")
  } else if (is.symbol(head)) name <- as.character(head) else
    repair_error("Unsupported replay expression; supply an explicit `refitter`.")
  base <- c("list", "c", "+", "-", "*", "/", "^", "(", "I", "log", "log1p", "sqrt", "exp", "abs")
  nlme <- c("lmeControl", "corAR1", "corCAR1", "corARMA", "varIdent", "varPower",
    "pdDiag", "pdSymm", "pdLogChol", "pdIdent", "pdCompSymm")
  fun <- if (name %in% base) get(name, baseenv()) else if (name == "lmerControl")
    lme4::lmerControl else if (name %in% nlme) getExportedValue("nlme", name) else
      repair_error("Unsupported replay constructor {.val {name}}; supply an explicit `refitter`.")
  args <- lapply(as.list(expr)[-1L], repair_replay_eval, data = data,
    literal = literal || name %in% c(nlme, "lmerControl"))
  do.call(fun, args)
}

repair_refit_close <- function(x, y, tolerance, what) {
  if (!isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
    repair_error("Reference replay does not reproduce {.val {what}}. Check reference data, coding, and the fitting recipe; no replacement fit was attempted.")
}

repair_refit_factors <- function(frame, coding = NULL, strict = TRUE) {
  fields <- names(frame)[vapply(frame, function(x) is.factor(x) || is.logical(x) || is.character(x), logical(1))]
  stats::setNames(lapply(fields, function(n) {
    logical_column <- is.logical(frame[[n]])
    value <- if (logical_column) factor(frame[[n]], levels = c(FALSE, TRUE)) else
      if (is.character(frame[[n]])) factor(frame[[n]]) else frame[[n]]
    if (nlevels(value) < 2L) return(NULL)
    code <- coding[[n]]
    if (is.character(code)) {
      if (code %in% c("contr.treatment", "contr.sum", "contr.helmert", "contr.poly", "contr.SAS")) {
        matrix <- getExportedValue("stats", code)(levels(value))
      } else if (!strict) matrix <- stats::contrasts(value) else
        repair_error("An unrecorded custom measurement contrast function requires an explicit contrast matrix in the data.")
    } else matrix <- if (is.matrix(code)) code else stats::contrasts(value)
    list(levels = levels(value), contrasts = matrix, ordered = is.ordered(value), logical = logical_column)
  }), fields)
}

repair_refit_factor_data <- function(data, spec, require_levels = TRUE) {
  for (name in names(spec)) {
    s <- spec[[name]]
    if (is.null(s)) next
    if (!name %in% names(data)) repair_error("Refit data lack factor column {.val {name}}.")
    values <- as.character(data[[name]])
    if (any(!is.na(values) & !values %in% s$levels) ||
        (require_levels && !all(s$levels %in% values)))
      repair_error("Refit factor {.val {name}} has new or absent levels; the original coefficient schema must remain estimable.")
    if (isTRUE(s$logical)) data[[name]] <- ifelse(is.na(values), NA, values == "TRUE") else {
      data[[name]] <- factor(values, levels = s$levels, ordered = s$ordered)
      stats::contrasts(data[[name]]) <- s$contrasts
    }
  }
  data
}

repair_refit_rows <- function(model) {
  if (inherits(model, "lmerMod")) rownames(stats::model.frame(model)) else rownames(model$groups)
}

repair_refit_frame <- function(data, name) {
  if (is.list(data) && !is.data.frame(data)) {
    if (!identical(names(data), name) || length(data) != 1L)
      repair_error("Supply a data frame or a one-element measurement-data list named {.val {name}}.")
    data <- data[[1L]]
  }
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      anyNA(names(data)) || any(!nzchar(names(data))))
    repair_error("Measurement data must be a nonempty data frame with unique column names.")
  data
}

repair_refit_spec <- function(model, data) {
  lmer <- inherits(model, "lmerMod")
  call <- stats::getCall(model)
  formula <- repair_replay_formula(stats::formula(model), names(data))
  frame <- if (lmer) stats::model.frame(model) else {
    stats::model.frame(model$terms, data, na.action = stats::na.fail)
  }
  coding <- if (lmer) attr(lme4::getME(model, "X"), "contrasts") else model$contrasts
  fixed <- all.vars(stats::delete.response(stats::terms(if (lmer) stats::formula(stats::terms(model, fixed.only = TRUE)) else formula)))
  if (length(setdiff(names(coding), names(data))))
    repair_error("Transformed factor coding requires an explicit `refitter` or precomputed factor columns.")
  factors <- repair_refit_factors(frame[intersect(names(coding), names(frame))], coding)
  arguments <- list()
  if (lmer) {
    arguments$REML <- lme4::isREML(model)
    arguments$control <- if (is.null(call$control)) lme4::lmerControl(
      optimizer = model@optinfo$optimizer, optCtrl = model@optinfo$control) else
        repair_replay_eval(call$control, data)
    if (!inherits(arguments$control, "merControl")) repair_error("The lmer replay control must be a lmerControl object.")
    if (!is.null(call$start)) arguments$start <- repair_replay_eval(call$start, data, literal = TRUE)
    slopes <- character()
    find_slopes <- function(x) {
      if (is.call(x)) {
        if (is.symbol(x[[1L]]) && as.character(x[[1L]]) %in% c("|", "||")) slopes <<- union(slopes, all.vars(x[[2L]]))
        else lapply(as.list(x)[-1L], find_slopes)
      }
      invisible(NULL)
    }
    find_slopes(formula)
    random_factors <- setdiff(intersect(slopes, names(frame)[vapply(frame, function(x) is.factor(x) || is.logical(x) || is.character(x), logical(1))]), fixed)
    if (length(random_factors)) repair_error("Random-only factor coding requires an explicit `refitter`.")
  } else {
    random_vars <- unique(unlist(lapply(stats::formula(model$modelStruct$reStruct, asList = TRUE), all.vars)))
    categorical <- names(data)[vapply(data, function(x) is.factor(x) || is.logical(x) || is.character(x), logical(1))]
    if (length(setdiff(intersect(random_vars, categorical), fixed)))
      repair_error("Random-only factor coding requires an explicit `refitter`.")
    arguments$method <- model$method
    if (is.null(call$random)) repair_error("Implicit nlme random formulas require an explicit `refitter`.")
    arguments$random <- repair_replay_eval(call$random, data)
    arguments$control <- if (is.null(call$control)) nlme::lmeControl() else repair_replay_eval(call$control, data)
    arguments$correlation <- repair_replay_eval(call$correlation, data)
    arguments$weights <- repair_replay_eval(call$weights, data)
    arguments$keep.data <- TRUE
  }
  list(lmer = lmer, formula = formula, arguments = arguments, factors = factors,
    weights = if (lmer) call$weights else NULL, offset = if (lmer) call$offset else NULL,
    expressions = if (!lmer) list(random = call$random, weights = call$weights, correlation = call$correlation) else NULL)
}

repair_refit_native <- function(spec, data) {
  data <- repair_refit_factor_data(data, spec$factors)
  args <- spec$arguments
  env <- new.env(parent = asNamespace("stats")); env$.repair_refit_data <- data
  formula <- repair_replay_formula(spec$formula, names(data)); environment(formula) <- env
  args[[if (spec$lmer) "formula" else "fixed"]] <- formula
  args$data <- as.name(".repair_refit_data"); args$na.action <- stats::na.fail
  if (!spec$lmer) for (name in names(spec$expressions))
    args[name] <- list(repair_replay_eval(spec$expressions[[name]], data, literal = TRUE))
  if (length(spec$factors)) args$contrasts <- lapply(spec$factors, `[[`, "contrasts")
  if (spec$lmer) {
    if (!is.null(spec$weights)) args$weights <- repair_replay_eval(spec$weights, data)
    if (!is.null(spec$offset)) args$offset <- repair_replay_eval(spec$offset, data)
  }
  model <- do.call(if (spec$lmer) lme4::lmer else nlme::lme, args, envir = env)
  # Keep a readable, self-contained data binding for ordinary model inspection.
  if (spec$lmer) {
    f <- stats::formula(model); environment(f) <- env
    model@call[[1L]] <- quote(lme4::lmer); model@call$formula <- f
    model@call$data <- as.name(".repair_refit_data")
    model@call$weights <- spec$weights; model@call$offset <- spec$offset
  } else {
    f <- model$call$fixed; environment(f) <- env
    model$call[[1L]] <- quote(nlme::lme); model$call$fixed <- f
    model$call$data <- as.name(".repair_refit_data")
  }
  model
}

repair_refit_control_declaration <- function(model) {
  raw <- stats::getCall(model)$control
  if (is.null(raw)) return(if (inherits(model, "lme")) nlme::lmeControl() else NULL)
  value <- if (is.list(raw)) raw else tryCatch(
    repair_replay_eval(raw, data.frame(), literal = TRUE), error = function(e) NULL)
  if (inherits(model, "lme") && !is.null(value)) {
    defaults <- nlme::lmeControl(); defaults[names(value)] <- value
    value <- defaults
  }
  value
}

repair_refit_configuration <- function(model) {
  e <- repair_get_mlm_elements(model)
  fixed_formula <- if (inherits(model, "lmerMod")) stats::formula(stats::terms(model, fixed.only = TRUE)) else stats::formula(model)
  frame <- if (inherits(model, "lmerMod")) stats::model.frame(model) else
    stats::model.frame(model$terms, droplevels(model$data[match(repair_refit_rows(model), rownames(model$data)), , drop = FALSE]), na.action = stats::na.fail)
  coding <- if (inherits(model, "lmerMod")) attr(lme4::getME(model, "X"), "contrasts") else model$contrasts
  group_vars <- unique(unlist(lapply(names(e@grouping), function(n) {
    if (n %in% names(frame)) n else tryCatch(all.vars(str2lang(n)), error = function(e) n)
  })))
  factor_fields <- names(frame)[vapply(frame, function(x) is.factor(x) || is.logical(x) || is.character(x), logical(1))]
  fields <- union(intersect(names(coding), names(frame)), setdiff(factor_fields, group_vars))
  out <- list(backend = e@backend, formula = paste(deparse(stats::formula(model)), collapse = " "),
    factors = repair_refit_factors(frame[fields], coding), fixed = colnames(e@fixed_design), random = e@random_terms,
    method = if (inherits(model, "lmerMod")) if (lme4::isREML(model)) "REML" else "ML" else model$method)
  if (inherits(model, "lmerMod")) {
    out$optimizer <- model@optinfo$optimizer
    out$optimizer_control <- model@optinfo$control
  } else {
    cs <- model$modelStruct$corStruct; vs <- model$modelStruct$varStruct
    out$random_class <- lapply(model$modelStruct$reStruct, class)
    out$correlation <- if (is.null(cs)) "independent" else if (inherits(cs, "corAR1") ||
      (inherits(cs, "corARMA") && identical(as.integer(attr(cs, "p")), 1L) && identical(as.integer(attr(cs, "q")), 0L))) "AR1" else class(cs)[1L]
    out$correlation_formula <- if (!is.null(cs)) deparse(stats::formula(cs))
    out$correlation_fixed <- if (!is.null(cs)) attr(cs, "fixed")
    out$variance <- if (is.null(vs)) "constant" else class(vs)[1L]
    out$variance_formula <- if (!is.null(vs)) deparse(stats::formula(vs))
    out$variance_fixed <- if (!is.null(vs)) attr(vs, "fixed")
    if (inherits(vs, "varIdent")) {
      groups <- names(stats::coef(vs, unconstrained = FALSE, allCoef = TRUE))
      out$variance_groups <- sort(groups)
      if (length(attr(vs, "fixed"))) out$variance_reference <- groups[1L]
    }
  }
  out
}

repair_refit_stage <- function(fitter, data, original, stage, tolerance) {
  before <- repair_refit_configuration(original)
  declared_control <- repair_refit_control_declaration(original)
  warnings <- messages <- character()
  model <- tryCatch(withCallingHandlers(fitter(data),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }),
    error = function(e) repair_error("{stage} measurement refit failed: {conditionMessage(e)}"))
  if (!(inherits(model, "lmerMod") || (inherits(model, "lme") && !inherits(model, "nlme"))))
    repair_error("{stage} refitter must return a supported Gaussian lmer or linear lme model.")
  if (!identical(repair_refit_rows(model), rownames(data)))
    repair_error("{stage} refit changed effective observations or their order. Supply exactly the intended effective rows; implicit subset/omission is not allowed.")
  if (inherits(model, "lmerMod") &&
      (any(unlist(model@optinfo$conv$opt) != 0) || any(unlist(model@optinfo$conv$lme4$code) != 0)))
    repair_error("{stage} measurement optimizer did not converge.")
  if (any(grepl("convergence error|iteration limit|maximum number of iterations|nlminb problem", warnings, ignore.case = TRUE)))
    repair_error("{stage} measurement optimizer did not converge: {paste(warnings, collapse = '; ')}")
  after <- repair_refit_configuration(model)
  if (!is.null(declared_control) && !isTRUE(all.equal(declared_control,
      repair_refit_control_declaration(model), tolerance = 0)))
    repair_error("{stage} refit changed the recorded optimizer control recipe.")
  if (!isTRUE(all.equal(before, after, tolerance = 0)))
    repair_error("{stage} refit changed the measurement specification, optimizer controls, or coefficient schema.")
  list(model = model, warnings = warnings, messages = messages,
    controls = if (!is.null(declared_control)) "recorded_recipe" else if (inherits(original, "lmerMod"))
      "recorded_optimizer" else "callback_declared")
}

repair_refit_identity <- function(original, replay, tolerance) {
  before <- repair_get_mlm_elements(original@mlm_fit[[1L]])
  after <- repair_get_mlm_elements(replay@mlm_fit[[1L]])
  for (field in c("response", "fixed_design", "offset", "observation_weights", "random_design", "grouping", "random_covariance"))
    repair_refit_close(S7::prop(before, field), S7::prop(after, field), tolerance, paste("measurement", field))
  if (!identical(before@group_levels, after@group_levels) ||
      !identical(lapply(before@grouping, as.character), lapply(after@grouping, as.character)))
    repair_error("Reference replay changed measurement IDs or grouping membership.")
  group <- original@analysis_data@effects$group[1L]
  for (id in before@group_levels[[group]]) {
    rows <- before@group_rows[[group]][[id]]
    repair_refit_close(get_ranef_target_covariance(before, group, rows),
      get_ranef_target_covariance(after, group, rows), tolerance, "measurement residual/nuisance covariance")
  }
  a <- original@analysis_data; b <- replay@analysis_data
  if (!identical(a@ids, b@ids) || !identical(a@effects, b@effects) ||
      !identical(a@diagnostics$excluded, b@diagnostics$excluded) || a@diagnostics$n_requested != b@diagnostics$n_requested)
    repair_error("Reference replay changed effect identities or the requested/analyzed cohort.")
  for (field in c("outcome", "predictors", "err_x", "err_y", "err_xy"))
    repair_refit_close(S7::prop(a, field), S7::prop(b, field), tolerance, paste("structural", field))
  repair_refit_close(repair_scores(original, "measurement"), repair_scores(replay, "measurement"), tolerance, "scores")
  repair_refit_close(stats::coef(original), stats::coef(replay), tolerance, "coefficients")
  repair_refit_close(stats::vcov(original), stats::vcov(replay), tolerance, "coefficient covariance")
}

# No random draws are made here. Re-key each requested copy of a complete person.
repair_refit_draw <- function(object, measurement, people, id, draw) {
  a <- object@analysis_data; group <- a@effects$group[1L]
  if (!is.character(draw) || !length(draw) || anyNA(draw) || any(!draw %in% a@ids))
    repair_error("`draw` must be a nonempty character vector of analyzed source IDs; repeats are allowed.")
  e <- repair_get_mlm_elements(object@mlm_fit[[1L]])
  if (!setequal(e@group_levels[[group]], a@ids))
    repair_error("Deterministic draws require the measurement and analyzed cohorts to coincide; first create an explicit restricted reference fit.")
  if (!group %in% names(measurement)) repair_error("Draws require a bare target grouping column; use explicit replacement data for computed targets.")
  model <- object@mlm_fit[[1L]]
  forbidden <- if (inherits(model, "lmerMod")) all.vars(stats::formula(stats::terms(model, fixed.only = TRUE))) else
    all.vars(stats::formula(model))
  if (inherits(model, "lmerMod")) {
    visit <- function(x) {
      if (is.call(x)) {
        if (is.symbol(x[[1L]]) && as.character(x[[1L]]) %in% c("|", "||")) forbidden <<- union(forbidden, all.vars(x[[2L]]))
        else lapply(as.list(x)[-1L], visit)
      }
      invisible(NULL)
    }
    visit(stats::formula(model))
    forbidden <- union(forbidden, c(all.vars(stats::getCall(model)$weights), all.vars(stats::getCall(model)$offset)))
  } else {
    forbidden <- union(forbidden, unlist(lapply(stats::formula(model$modelStruct$reStruct, asList = TRUE), all.vars)))
    for (part in c("corStruct", "varStruct")) {
      component <- model$modelStruct[[part]]
      if (!is.null(component)) forbidden <- union(forbidden, all.vars(nlme::getCovariateFormula(component)))
    }
  }
  if (group %in% forbidden) repair_error("A grouping ID used as a covariate cannot be re-keyed by `draw`.")
  nuisance <- setdiff(names(e@grouping), group)
  rekey <- character()
  for (n in nuisance) {
    vars <- if (n %in% names(measurement)) n else tryCatch(all.vars(str2lang(n)),
      error = function(e) repair_error("Computed nuisance groups require explicit replacement data."))
    if (group %in% vars) next
    if (!n %in% names(measurement) || n %in% forbidden)
      repair_error("Draws require nuisance grouping expressions containing the target ID or bare nuisance ID columns unused as covariates; otherwise supply explicit replacement data.")
    rekey <- c(rekey, n)
  }
  mapping <- data.frame(draw = seq_along(draw), source_id = draw,
    refit_id = sprintf("r%07d", seq_along(draw)), stringsAsFactors = FALSE)
  blocks <- lapply(seq_along(draw), function(i) {
    block <- measurement[as.character(measurement[[group]]) == draw[i], , drop = FALSE]
    block[[group]] <- mapping$refit_id[i]
    for (n in rekey) block[[n]] <- paste(mapping$refit_id[i], as.character(block[[n]]), sep = "::")
    block
  })
  measurement <- do.call(rbind, blocks); rownames(measurement) <- sprintf("observation%09d", seq_len(nrow(measurement)))
  people <- people[match(draw, as.character(people[[id]])), , drop = FALSE]
  people[[id]] <- mapping$refit_id; rownames(people) <- NULL
  list(measurement = measurement, people = people, mapping = mapping)
}

#' Refit a REPAIR analysis after verifying its reference data
#'
#' Replay one measurement model and the complete structural analysis using explicit
#' reference inputs. Optional replacement inputs or a deterministic whole-person
#' draw are fitted only after identity replay succeeds. This is a refitting tool,
#' not a bootstrap inference engine.
#' @param object A converged [repair()] fit using one measurement model.
#' @param measurement_data Reference repeated-observation data, as a data frame or
#'   one-element list named like the fitted model. Retain original fitted row names.
#'   Rows excluded from the original measurement fit are filtered by those names.
#' @param data Original person-level data, including requested/excluded rows.
#' @param new_measurement_data Optional replacement effective repeated observations.
#'   Every supplied row must be used; implicit subset or omission is rejected.
#' @param new_data Optional replacement person-level data. Every requested person
#'   must be complete; the original coefficient schema must remain estimable.
#' @param id Structural ID column; defaults to the original ID column and must match it.
#' @param refitter Optional deterministic function of one data frame returning a
#'   fitted lmer/lme model. Required for nonliteral control objects or unsupported
#'   formula/argument dependencies. Its reference replay is subject to the same checks.
#' @param draw Optional character vector of analyzed source IDs, with repeats
#'   allowed. Copy whole people and assign fresh target/nuisance identities. Cannot
#'   be combined with replacement data; no randomness or intervals are generated.
#' @param tolerance Positive finite numerical tolerance for identity replay;
#'   default 1e-6. Identity/schema/cohort checks are exact.
#' @details Native adapters rebuild supported Gaussian lmer and single-group lme
#'   fits using data-only formulas, weights/offsets, literal fitting constructors,
#'   recorded ML/REML, and fixed factor coding. They do not evaluate original data
#'   or subset expressions in the caller's environment. Reference fitting uses
#'   the original effective rows; replacements require explicit complete rows.
#'   Symbolic controls and other external dependencies require `refitter`.
#'   Recorded control recipes are compared exactly; lmer also retains optimizer
#'   settings. Unrecoverable nlme options are the callback author's explicit
#'   declaration, identified as `callback_declared` in the audit.
#'
#'   Original-data replay checks effective observations, designs, weights/offsets,
#'   covariance, scores, structural inputs, coefficients and their covariance.
#'   Measurement convergence, specification, and coefficient identities are checked
#'   again for replacements. Singular but usable measurement fits are retained.
#'   Warnings/messages are recorded in `summary(result)$refit`; fitting errors stop.
#'   Original structural formula, regularization and computation backend are retained.
#'
#'   `draw` requires matching measurement/analyzed cohorts and a bare target ID.
#'   Nuisance IDs must be bare columns or grouping expressions containing that ID.
#'   Shared nuisance groups are rejected by the ordinary structural validator.
#'   See `vignette("refitting", package = "repair")` for scope and examples.
#' @return An ordinary fitted REPAIR object, with a reference/replacement audit
#'   in `summary(result)$refit`. Original conditional SEs remain conditional;
#'   no bootstrap uncertainty or coverage claim is added.
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' replay <- repair_refit(fit, repair_example$repeated, repair_example$people)
#' repair_coefficients(replay)
#' summary(replay)$refit$checks
#' @export
repair_refit <- function(object, measurement_data, data, new_measurement_data = NULL,
                         new_data = NULL, id = NULL, refitter = NULL, draw = NULL,
                         tolerance = 1e-6) {
  repair_require_model(object); repair_require_fit(object)
  if (length(object@mlm_fit) != 1L) repair_error("Initial refitting supports one measurement model; multiple-model covariance replay is not yet supported.")
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0)
    repair_error("`tolerance` must be one positive finite number.")
  if (!is.null(refitter) && !is.function(refitter)) repair_error("`refitter` must be a function of one data frame.")
  if (!is.data.frame(data) || !nrow(data)) repair_error("`data` must contain the original person-level reference rows.")
  if (!is.null(draw) && (!is.null(new_measurement_data) || !is.null(new_data)))
    repair_error("Use either `draw` or explicit replacement data, not both.")
  a <- object@analysis_data; name <- names(object@mlm_fit); original <- object@mlm_fit[[1L]]
  if (is.null(id)) id <- a@diagnostics$id_column
  if (!identical(id, a@diagnostics$id_column)) repair_error("Use the original structural ID column.")
  supplied <- repair_refit_frame(measurement_data, name)
  rows <- repair_refit_rows(original)
  positions <- match(rows, rownames(supplied))
  if (anyNA(positions)) repair_error("Reference measurement data lack fitted observation row names; restore the original row identities.")
  reference <- supplied[positions, , drop = FALSE]
  spec <- if (is.null(refitter)) a@diagnostics$refit$spec %||% repair_refit_spec(original, reference) else NULL
  fitter <- if (is.null(refitter)) function(d) repair_refit_native(spec, d) else refitter
  controls <- object@struc_fit@diagnostics
  preliminary <- controls$preliminary_moment_control; final <- controls$final_moment_control
  if (!identical(preliminary, final)) repair_error("Replay requires matching preliminary/final regularization controls from `repair()`.")
  factors <- a@diagnostics$structural_factors %||% list()
  structural_warnings <- list()
  structural <- function(model, people, action, stage) {
    models <- stats::setNames(list(model), name)
    assembled <- withCallingHandlers(repair_assemble(a@formula, models, people, id, na_action = action,
      structural_contrasts = lapply(factors, `[[`, "contrasts")), warning = function(w) {
        structural_warnings[[stage]] <<- c(structural_warnings[[stage]], conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    fitted <- run_struc(assembled, preliminary$correction, preliminary$alpha, controls$backend)
    if (!isTRUE(fitted@converged)) repair_error("Structural refit failed: {fitted@message}")
    class_repair_mod(mlm_fit = models, struc_fit = fitted, analysis_data = assembled, call = match.call())
  }
  baseline <- repair_refit_stage(fitter, reference, original, "Reference", tolerance)
  replay <- structural(baseline$model, data, a@diagnostics$na_action, "reference")
  repair_refit_identity(object, replay, tolerance)
  if (!length(factors)) factors <- replay@analysis_data@diagnostics$structural_factors
  replacement <- !is.null(draw) || !is.null(new_measurement_data) || !is.null(new_data)
  mapping <- NULL; result <- replay
  checks <- data.frame(stage = "reference", n_measurement = nrow(reference),
    n_requested = a@diagnostics$n_requested, n_analyzed = length(a@ids), controls = baseline$controls, verified = TRUE)
  warning_log <- list(reference = c(baseline$warnings, structural_warnings$reference)); message_log <- list(reference = baseline$messages)
  if (replacement) {
    d <- if (is.null(new_measurement_data)) reference else repair_refit_frame(new_measurement_data, name)
    people <- if (is.null(new_data)) data else new_data
    if (!is.null(draw)) {
      copied <- repair_refit_draw(object, reference, data, id, draw)
      d <- copied$measurement; people <- copied$people; mapping <- copied$mapping
    }
    people <- repair_refit_factor_data(people, factors)
    changed <- repair_refit_stage(fitter, d, original, "Replacement", tolerance)
    result <- structural(changed$model, people, "fail", "replacement")
    if (!identical(names(stats::coef(object)), names(stats::coef(result))) ||
        !identical(a@effects, result@analysis_data@effects))
      repair_error("Replacement refit changed the structural coefficient/effect schema.")
    checks <- rbind(checks, data.frame(stage = "replacement", n_measurement = nrow(d),
      n_requested = nrow(people), n_analyzed = stats::nobs(result), controls = changed$controls, verified = TRUE))
    warning_log$replacement <- c(changed$warnings, structural_warnings$replacement); message_log$replacement <- changed$messages
  }
  audit <- list(mode = if (replacement) "replacement" else "identity", adapter = if (is.null(refitter)) "native" else "callback",
    tolerance = tolerance, checks = checks, warnings = warning_log, messages = message_log,
    draw = mapping, reference_excluded_rows = setdiff(rownames(supplied), rows),
    versions = list(R = R.version.string, repair = as.character(utils::packageVersion("repair")),
      measurement = as.character(utils::packageVersion(if (inherits(original, "lmerMod")) "lme4" else "nlme"))),
    spec = spec)
  result@analysis_data@diagnostics$refit <- audit
  result
}
