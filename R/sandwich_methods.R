repair_sandwich_covariance <- function(object, component) {
  if (!inherits(object, "repair_sandwich")) repair_error("Supply a `repair_sandwich()` result.")
  key <- if (component == "stage2") "stage2" else "full"
  if (!isTRUE(object$diagnostics[[paste0(key, "_available")]]))
    repair_error("Influence inference unavailable: {object$diagnostics[[paste0(key, '_reason')]]}")
  object$covariance[[component]]
}

#' Extract influence-function standard errors and normal intervals
#' @param object A [repair_sandwich()] result.
#' @param conf.level Confidence level, strictly between zero and one.
#' @param component `"full"` (default) or frozen-measurement `"stage2"`.
#' @return A plain data frame with term, estimate, std.error, conf.low, conf.high,
#'   conf.level, component and inference columns. Original point estimates are
#'   unchanged; intervals are experimental normal Wald intervals, not bootstrap
#'   percentiles. Unavailable inference stops explicitly.
#' @export
repair_sandwich_coefficients <- function(object, conf.level = .95, component = c("full", "stage2")) {
  component <- match.arg(component)
  if (!is.numeric(conf.level) || length(conf.level) != 1L || !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    repair_error("`conf.level` must be strictly between zero and one.")
  V <- repair_sandwich_covariance(object, component); se <- sqrt(diag(V)); z <- stats::qnorm((1 + conf.level) / 2)
  data.frame(term = names(object$estimate), estimate = unname(object$estimate), std.error = unname(se),
    conf.low = unname(object$estimate - z * se), conf.high = unname(object$estimate + z * se),
    conf.level = conf.level, component = component,
    inference = object$inference %||% "Experimental person influence-function sandwich; local asymptotic normal inference", row.names = NULL)
}

#' Methods for influence-function inference
#' @name repair-sandwich-methods
#' @param object,x A [repair_sandwich()] result.
#' @param component Covariance component: full, stage2, measurement or cross for
#'   `vcov()`; full or stage2 for intervals. The cross component is a signed
#'   contribution and is not a standalone positive-semidefinite covariance.
#' @param parm Coefficient names or indices; missing selects all coefficients.
#' @param level Confidence level strictly between zero and one.
#' @param ... Reserved for generic compatibility.
#' @return `coef()` returns the original coefficient vector; `nobs()` the number
#'   of people; `vcov()` a named covariance/contribution matrix; `confint()` normal
#'   intervals. `summary()` remains inspectable when inference is unavailable,
#'   returning full/stage2 coefficient tables or NULL and diagnostic reasons.
#'   Printing returns the object invisibly.
NULL

#' @rdname repair-sandwich-methods
#' @exportS3Method stats::coef
coef.repair_sandwich <- function(object, ...) object$estimate

#' @rdname repair-sandwich-methods
#' @exportS3Method stats::nobs
nobs.repair_sandwich <- function(object, ...) length(object$ids)

#' @rdname repair-sandwich-methods
#' @exportS3Method stats::vcov
vcov.repair_sandwich <- function(object, ..., component = c("full", "stage2", "measurement", "cross")) {
  repair_sandwich_covariance(object, match.arg(component))
}

#' @rdname repair-sandwich-methods
#' @exportS3Method stats::confint
confint.repair_sandwich <- function(object, parm, level = .95, ..., component = c("full", "stage2")) {
  table <- repair_sandwich_coefficients(object, level, match.arg(component))
  interval <- cbind(table$conf.low, table$conf.high); rownames(interval) <- table$term
  colnames(interval) <- paste0(100 * c((1 - level) / 2, (1 + level) / 2), "%")
  if (!missing(parm)) {
    if (is.character(parm)) parm <- match(parm, table$term)
    if (!is.numeric(parm) || anyNA(parm) || any(parm < 1 | parm > nrow(interval) | parm %% 1))
      repair_error("`parm` must identify existing coefficients.")
    interval <- interval[parm, , drop = FALSE]
  }
  interval
}

#' @rdname repair-sandwich-methods
#' @exportS3Method base::summary
summary.repair_sandwich <- function(object, ...) {
  table <- function(component) if (object$diagnostics[[paste0(component, "_available")]])
    repair_sandwich_coefficients(object, component = component) else NULL
  structure(list(coefficients = table("full"), stage2 = table("stage2"), diagnostics = object$diagnostics,
    criterion = object$diagnostics$criterion, n_people = length(object$ids), method = object$method), class = "summary_repair_sandwich")
}

#' @rdname repair-sandwich-methods
#' @exportS3Method base::print
print.summary_repair_sandwich <- function(x, ...) {
  cat(x$method %||% "REPAIR person influence-function sandwich (experimental)", "\n", sep = "")
  cat("People:", x$n_people, "\n")
  if (!is.null(x$coefficients)) print(x$coefficients) else cat("Full unavailable:", x$diagnostics$full_reason, "\n")
  if (is.null(x$stage2)) cat("Stage2 unavailable:", x$diagnostics$stage2_reason, "\n")
  cat("Local linearization; unchanged original estimates; coverage is not guaranteed.\n")
  invisible(x)
}

#' @rdname repair-sandwich-methods
#' @exportS3Method base::print
print.repair_sandwich <- function(x, ...) {
  print(summary(x)); invisible(x)
}
