#' Refer to a measurement-model random coefficient
#'
#' Use only inside the structural formula passed to [repair()] or
#' [repair_data()]. Each reference selects one scalar coefficient, not a
#' multi-column formula term. Inspect available names with [repair_effects()].
#' @param term Random coefficient name; use `"(Intercept)"` for an intercept.
#' @param model Name in the `models` list. May be omitted if unambiguous.
#' @param group Target grouping-factor name from [repair_effects()]. May be
#'   omitted if the coefficient is unambiguous. Select one target group per model;
#'   other lmer groups must have every fitted level contained within that target.
#' @return A formula marker; calling it outside a structural formula is an error.
#' @examples
#' equation <- re("day", model = "growth") ~ practice
#' equation
#' # Use this formula with models = list(growth = your_fitted_model).
#' @export
re <- function(term, model = NULL, group = NULL) {
  repair_error("`re()` is a formula marker; use it inside `repair()` or `repair_data()`.")
}

# Deliberately parse a small grammar rather than evaluating arbitrary re() calls.
parse_struc_formula <- function(formula, catalog) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    repair_error("Supply one two-sided structural formula.")
  }
  effects <- list()
  observed <- character()
  parse_atom <- function(expr) {
    if (is.symbol(expr) && as.character(expr) != ".") {
      name <- as.character(expr)
      if (startsWith(name, ".repair_")) {
        repair_error("Names beginning with `.repair_` are reserved.")
      }
      observed <<- union(observed, name)
      return(expr)
    }
    if (!is.call(expr) || !identical(expr[[1L]], as.name("re"))) {
      repair_error(paste0(
        "Use additive bare variable names and `re(\"term\", model = \"name\", group = \"target\")`. ",
        "Interactions, transformations, offsets, and `.` are not supported; ",
        "precompute transformations of observed variables in `data`."
      ))
    }
    matched <- tryCatch(match.call(definition = re, call = expr),
      error = function(e) repair_error("Use `re(\"term\", model = \"name\", group = \"target\")` with literal strings."))
    args <- list(term = matched$term)
    if (!is.null(matched$model)) args$model <- matched$model
    if (!is.null(matched$group)) args$group <- matched$group
    if (!all(vapply(args, function(x) is.character(x) && length(x) == 1L &&
          !is.na(x) && nzchar(x), logical(1)))) {
      repair_error("Use `re(\"term\", model = \"name\", group = \"target\")` with literal strings.")
    }
    term <- args$term
    if (term %in% c("1", "Intercept")) term <- "(Intercept)"
    candidates <- catalog[catalog$term == term, , drop = FALSE]
    if (!is.null(args$model)) candidates <- candidates[candidates$model == args$model, , drop = FALSE]
    if (!is.null(args$group)) candidates <- candidates[candidates$group == args$group, , drop = FALSE]
    if (nrow(candidates) != 1L) {
      repair_error("Random-effect reference is unknown or ambiguous. Use `repair_effects(models)` and specify `model` and, where needed, `group`.")
    }
    entry <- candidates[1L, , drop = FALSE]
    key <- paste(entry$model, entry$group, entry$term, sep = "\r")
    existing <- which(vapply(effects, function(x) identical(x$key, key), logical(1)))
    if (length(existing)) return(as.name(effects[[existing]]$column))
    column <- paste0(".repair_re", length(effects) + 1L)
    label_call <- call("re", entry$term, model = entry$model)
    if (length(unique(catalog$group[catalog$model == entry$model])) > 1L) {
      label_call$group <- entry$group
    }
    label <- paste(deparse(label_call), collapse = "")
    effects[[length(effects) + 1L]] <<- list(
      key = key, column = column, label = label,
      model = entry$model, group = entry$group, term = entry$term
    )
    as.name(column)
  }
  parse_rhs <- function(expr) {
    if (is.call(expr) && identical(expr[[1L]], as.name("+")) && length(expr) == 3L) {
      return(call("+", parse_rhs(expr[[2L]]), parse_rhs(expr[[3L]])))
    }
    if (is.numeric(expr) && identical(as.numeric(expr), 1)) return(1)
    parse_atom(expr)
  }
  lhs <- parse_atom(formula[[2L]])
  rhs <- parse_rhs(formula[[3L]])
  parsed <- stats::as.formula(call("~", lhs, rhs), env = baseenv())
  tt <- stats::terms(parsed)
  if (!length(attr(tt, "term.labels"))) repair_error("At least one structural predictor is required.")
  if (as.character(lhs) %in% all.vars(rhs)) repair_error("The outcome cannot also be a predictor.")
  if (!length(effects)) repair_error("The structural formula must contain at least one `re()` reference.")
  used_models <- unique(vapply(effects, `[[`, character(1), "model"))
  for (name in used_models) {
    selected <- effects[vapply(effects, function(e) e$model == name, logical(1))]
    if (length(unique(vapply(selected, `[[`, character(1), "group"))) != 1L) {
      repair_error("All selected random coefficients from model {.val {name}} must use one target grouping factor.")
    }
  }
  list(formula = parsed, original = formula, effects = effects, observed = observed,
       model = used_models, outcome = as.character(lhs))
}
