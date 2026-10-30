#' Extract and summarize a REPAIR fit
#'
#' Standard extraction methods operate on [repair()] results. `summary()` returns
#' a list containing a coefficient table, effect sources, sample accounting,
#' regularization, measurement-model messages, posterior reliability, weight
#' concentration, and the inference convention. See [repair_reliability()] and
#' [repair_diagnostics()] for definitions and individual diagnostic tables.
#' Use [repair_coefficients()] for a plain coefficient table with a chosen
#' confidence level, and [repair_scores()] for ID-keyed measurement scores.
#' Results from [repair_refit()] also expose the replay audit in `summary(fit)$refit`.
#' `nobs()` counts analyzed people, not repeated observations. `confint()` uses
#' normal quantiles and supports named or integer coefficient selection.
#' Failed numerical fits remain inspectable with `summary()`; coefficient and
#' covariance accessors stop rather than presenting failed estimates.
#' @name repair-methods
#' @importFrom stats coef vcov nobs confint
#' @param object,x A fitted REPAIR result.
#' @param parm Coefficient names or integer positions; missing selects all.
#' @param level Confidence level strictly between zero and one.
#' @param ... Reserved for generic compatibility.
#' @return `coef()` returns a named numeric vector; `vcov()` a matrix;
#'   `confint()` a two-column matrix; `nobs()` an integer; `summary()` a list.
NULL

repair_require_fit <- function(object) {
  if (!isTRUE(object@struc_fit@converged)) {
    message <- object@struc_fit@message
    repair_error("Structural fit did not converge: {message}")
  }
  object@struc_fit
}

S7::method(coef, class_repair_mod) <- function(object, ...) {
  fit <- repair_require_fit(object)
  stats::setNames(fit@coefficients$estimate, fit@coefficients$term)
}

S7::method(vcov, class_repair_mod) <- function(object, ...) {
  repair_require_fit(object)@vcov
}

S7::method(nobs, class_repair_mod) <- function(object, ...) {
  length(object@analysis_data@ids)
}

S7::method(confint, class_repair_mod) <- function(object, parm, level = 0.95, ...) {
  beta <- stats::coef(object)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    repair_error("`level` must be a finite number strictly between zero and one.")
  }
  if (missing(parm)) parm <- names(beta)
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(!is.finite(parm)) || any(parm %% 1 != 0) ||
        any(parm < 1 | parm > length(beta))) repair_error("Invalid coefficient positions in `parm`.")
    parm <- names(beta)[parm]
  }
  if (!is.character(parm) || anyNA(parm) || any(!parm %in% names(beta))) {
    repair_error("Unknown coefficient names in `parm`.")
  }
  width <- stats::qnorm((1 + level) / 2) * sqrt(diag(stats::vcov(object)))[parm]
  result <- cbind(beta[parm] - width, beta[parm] + width)
  colnames(result) <- paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), " %")
  result
}

S7::method(summary, class_repair_mod) <- function(object, ...) {
  fit <- object@struc_fit
  data <- object@analysis_data
  structure(list(call = object@call, formula = data@formula,
    converged = isTRUE(fit@converged), message = fit@message,
    coefficients = fit@coefficients, effects = data@effects,
    n_requested = data@diagnostics$n_requested, n_analyzed = length(data@ids),
    excluded = data@diagnostics$excluded,
    unused_score_ids = data@diagnostics$unused_score_ids,
    measurement_backend = data@diagnostics$measurement_backend,
    measurement_singular = data@diagnostics$measurement_singular,
    measurement_messages = data@diagnostics$measurement_messages,
    measurement_models = data@diagnostics$measurement_models,
    cross_model_error = data@diagnostics$cross_model_error,
    regularization = fit@diagnostics[c("preliminary_moment_control", "final_moment_control",
      "correction_prelim", "correction_final")],
    diagnostics = fit@diagnostics,
    reliability = repair_reliability(object),
    structural_diagnostics = repair_structural_diagnostics(object),
    refit = data@diagnostics$refit[setdiff(names(data@diagnostics$refit), "spec")],
    inference = "Conditional on fitted measurement-model parameters; model-based normal SEs and intervals."),
    class = "summary_repair")
}

#' @exportS3Method base::print
print.summary_repair <- function(x, ...) {
  cat("REPAIR structural regression\n")
  cat(paste(deparse(x$formula), collapse = " "), "\n")
  cat("People:", x$n_analyzed, "analyzed /", x$n_requested, "requested;",
      nrow(x$excluded), "excluded\n")
  for (name in names(x$measurement_models)) {
    d <- x$measurement_models[[name]]
    cat("Measurement model: ", name, " (", d$measurement_backend, "; group ", d$group,
        "; ", length(d$unused_score_ids), " IDs outside cohort)\n", sep = "")
    if (length(d$nuisance_groups)) {
      cat("Nested nuisance groups: ", paste(d$nuisance_groups, collapse = ", "),
        " (each fitted level contained within ", d$group, ")\n", sep = "")
    }
    if (identical(d$measurement_backend, "nlme")) cat("nlme singularity is not assessed; inspect the measurement fit and its fitting warnings.\n")
    if (isTRUE(d$measurement_singular)) cat("Measurement model: singular random-effect covariance.\n")
    if (length(d$measurement_messages)) cat("Measurement-model messages:", paste(d$measurement_messages, collapse = "; "), "\n")
  }
  cat("Cross-model scoring errors: ", switch(x$cross_model_error$policy,
      independent = "assumed independent (zero covariance)",
      supplied = "user-supplied covariance blocks (conditional inference)",
      not_applicable = "not applicable (one model)"), "\n", sep = "")
  if (!is.null(x$refit)) {
    cat("Refit: original-data replay verified; ", x$refit$mode,
      " (", x$refit$adapter, ")\n", sep = "")
    if (length(unlist(x$refit$warnings))) cat("Refit warnings recorded: inspect summary(fit)$refit$warnings.\n")
  }
  cat("Scores: likelihood-only random deviations\n")
  print(x$effects[, c("model", "group", "term", "role"), drop = FALSE], row.names = FALSE)
  cat("\nPosterior reliability (analyzed cohort; equal person weights):\n")
  print(as.data.frame(x$reliability)[, c("model", "term", "role", "measure", "reliability",
    "n_people", "n_defined"), drop = FALSE], row.names = FALSE, digits = 3)
  if (any(x$reliability$n_undefined > 0)) cat("Undefined reliability: inspect repair_reliability(fit) for reasons.\n")
  if (x$converged) {
    cat("\n")
    table <- x$coefficients
    rownames(table) <- table$term
    table$term <- NULL
    print(table, digits = 4)
    cat("\nRegularization (preliminary / final):",
      x$regularization$preliminary_moment_control$correction, "/",
      x$regularization$final_moment_control$correction, "\n")
    cat("Error-moment subtraction multipliers:",
      format(x$regularization$correction_prelim, digits = 4), "/",
      format(x$regularization$correction_final, digits = 4), "\n")
    weight <- x$structural_diagnostics$weight_summary
    cat("Largest analysis-weight fraction:", format(weight$max_fraction, digits = 3),
        "; weight-equivalent count:", format(weight$weight_equivalent_n, digits = 3), "\n")
    cat("Moment conditioning and individual weights: repair_diagnostics(fit).\n")
  } else {
    cat("Structural fit FAILED:", x$message, "\n")
  }
  if (nrow(x$excluded)) cat("Excluded IDs and reasons: inspect summary(fit)$excluded.\n")
  cat(x$inference, "\n")
  invisible(x)
}

S7::method(print, class_repair_mod) <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

S7::method(print, class_repair_data) <- function(x, ...) {
  cat("REPAIR analysis data:", length(x@ids), "people,", ncol(x@predictors), "predictor columns\n")
  cat("Excluded:", nrow(x@diagnostics$excluded), "of", x@diagnostics$n_requested, "requested IDs\n")
  print(x@effects, row.names = FALSE)
  invisible(x)
}
