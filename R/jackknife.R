#' Refit REPAIR after deleting each complete target person
#' @inheritParams repair_refit
#' @details Experimental serial delete-one person jackknife through
#'   [repair_refit()]. Reference replay is verified first. Each deletion refits
#'   the measurement model and complete structural estimator on N-1 people,
#'   updating finite-N corrections and all nuisance parameters. Original
#'   coefficients remain the point estimates; the estimated jackknife bias is
#'   reported separately and is not subtracted from confidence intervals.
#'   All N complete coefficient vectors must be usable. Failed deletions are
#'   retained and invalidate covariance/interval inference; no successes-only
#'   option, omitted cases or replacement attempts are provided. Singular but
#'   converged measurement fits follow repair_refit's policy and are recorded.
#'   Scope matches the single-model deterministic draws supported by refitting,
#'   including contained nuisance effects. Callbacks must be deterministic.
#'   Standard caller RNG kind/state are preserved on success/error; hidden-state
#'   Box-Muller and user-supplied RNGs are rejected, as for person bootstrapping.
#' @return Separate S3 `repair_jackknife` result containing original estimates,
#'   ID-keyed N-by-p delete-one coefficients, attempt diagnostics/optimizer logs,
#'   reference audit, bias estimate and versions. `vcov()` uses
#'   (N-1)/N times the centered deletion cross-product; `confint()` and
#'   [repair_jackknife_coefficients()] use normal Wald intervals by default.
#'   Neither branch switching nor biased point estimation is automatically fixed.
#' @export
repair_jackknife <- function(object, measurement_data, data, id = NULL, refitter = NULL, tolerance = 1e-6) {
  repair_require_model(object); repair_require_fit(object)
  kind <- RNGkind()
  if (kind[1L] == "user-supplied" || kind[2L] %in% c("Box-Muller", "user-supplied"))
    repair_error("Jackknife RNG isolation cannot restore hidden RNG state; select a standard RNG with Inversion normals.")
  present <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  seed <- if (present) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (present) assign(".Random.seed", seed, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  replay <- repair_refit(object, measurement_data, data, id = id, refitter = refitter, tolerance = tolerance)
  cohort <- object@analysis_data@ids; n <- length(cohort)
  if (n < 3L) repair_error("Delete-one inference requires at least three people.")
  if (is.null(id)) id <- object@analysis_data@diagnostics$id_column
  measurement <- repair_refit_frame(measurement_data, names(object@mlm_fit))
  measurement <- measurement[match(repair_refit_rows(object@mlm_fit[[1L]]), rownames(measurement)), , drop = FALSE]
  repair_refit_draw(object, measurement, data, id, cohort)
  estimate <- stats::coef(object)
  values <- matrix(NA_real_, n, length(estimate), dimnames = list(cohort, names(estimate)))
  records <- vector("list", n)
  for (i in seq_len(n)) {
    warnings <- messages <- character(); started <- proc.time()[[3L]]
    value <- tryCatch(withCallingHandlers({
      fitted <- repair_refit(object, measurement_data, data, id = id, refitter = refitter,
        draw = cohort[-i], tolerance = tolerance)
      beta <- stats::coef(fitted)
      if (!identical(names(beta), names(estimate)) || any(!is.finite(beta)))
        repair_error("Deletion returned incompatible or nonfinite coefficients.")
      if (stats::nobs(fitted) != n - 1L) repair_error("Deletion did not retain exactly N-1 analyzed people.")
      list(fit = fitted, beta = beta)
    }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
       message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }), error = function(e) e)
    ok <- !inherits(value, "error"); singular <- NA; optimization <- NULL; residual <- preliminary <- final <- NA_real_
    if (ok) {
      values[i, ] <- value$beta; fit <- value$fit
      audit <- fit@analysis_data@diagnostics$refit
      warnings <- c(warnings, unlist(audit$warnings, use.names = FALSE))
      messages <- c(messages, unlist(audit$messages, use.names = FALSE))
      optimization <- audit$optimization$replacement
      if (inherits(fit@mlm_fit[[1L]], "lmerMod")) singular <- lme4::isSingular(fit@mlm_fit[[1L]])
      dg <- fit@struc_fit@diagnostics
      residual <- dg$resid_var; preliminary <- dg$correction_prelim; final <- dg$correction_final
    } else {
      warnings <- c(warnings, value$repair_warnings)
      messages <- c(messages, value$repair_messages)
      optimization <- value$repair_optimization
    }
    row <- data.frame(deleted_id = cohort[i], n_retained = n - 1L, usable = ok,
      stage = if (ok) "complete" else value$repair_stage %||% "refit",
      reason = if (ok) "" else conditionMessage(value), measurement_singular = singular,
      residual_variance = residual, correction_preliminary = preliminary, correction_final = final,
      elapsed = proc.time()[[3L]] - started)
    row$warnings <- I(list(unique(warnings))); row$messages <- I(list(unique(messages)))
    row$optimization <- I(list(optimization)); records[[i]] <- row
  }
  diagnostics <- do.call(rbind, records)
  bias <- if (all(diagnostics$usable)) (n - 1) * (colMeans(values) - estimate) else
    stats::setNames(rep(NA_real_, length(estimate)), names(estimate))
  structure(list(estimate = estimate, ids = cohort, replicates = values, diagnostics = diagnostics,
    bias = bias, reference = replay@analysis_data@diagnostics$refit,
    versions = replay@analysis_data@diagnostics$refit$versions,
    call = match.call()), class = "repair_jackknife")
}

repair_jackknife_covariance <- function(object) {
  if (!inherits(object, "repair_jackknife")) repair_error("Supply a repair_jackknife result.")
  if (!all(object$diagnostics$usable)) repair_error("Jackknife inference withheld: {sum(!object$diagnostics$usable)} deletion attempts failed; all N vectors are required.")
  n <- length(object$ids); centered <- scale(object$replicates, scale = FALSE)
  (n - 1) / n * crossprod(centered)
}

#' Extract person-jackknife coefficients and intervals
#' @param object A [repair_jackknife()] result.
#' @param conf.level Confidence level strictly between zero and one.
#' @param reference `"normal"` (default) or experimental `"t"` with N-1 degrees
#'   of freedom. The latter is a heuristic, not derived Satterthwaite inference.
#' @return Data frame retaining original estimates, jackknife standard errors,
#'   bias estimates, interval endpoints, confidence level and reference.
#' @export
repair_jackknife_coefficients <- function(object, conf.level = .95, reference = c("normal", "t")) {
  reference <- match.arg(reference); repair_bootstrap_level(conf.level)
  se <- sqrt(diag(repair_jackknife_covariance(object))); n <- length(object$ids)
  critical <- if (reference == "normal") stats::qnorm((1 + conf.level) / 2) else stats::qt((1 + conf.level) / 2, n - 1)
  data.frame(term = names(object$estimate), estimate = unname(object$estimate), std.error = unname(se),
    bias = unname(object$bias), conf.low = unname(object$estimate - critical * se),
    conf.high = unname(object$estimate + critical * se), conf.level = conf.level, reference = reference,
    df = if (reference == "t") n - 1 else Inf, row.names = NULL)
}

#' Methods for person delete-one inference
#' @name repair-jackknife-methods
#' @param object,x A [repair_jackknife()] result.
#' @param parm Coefficient names or indices; missing selects all.
#' @param level Confidence level strictly between zero and one.
#' @param reference Normal or experimental t reference as in [repair_jackknife_coefficients()].
#' @param ... Reserved for generic compatibility.
#' @return Original coefficients, person count, centered jackknife covariance,
#'   Wald intervals or inspectable summary with all attempt diagnostics. Printing
#'   returns its input invisibly. No fallback inference is used on failures.
NULL

#' @rdname repair-jackknife-methods
#' @exportS3Method stats::coef
coef.repair_jackknife <- function(object, ...) object$estimate
#' @rdname repair-jackknife-methods
#' @exportS3Method stats::nobs
nobs.repair_jackknife <- function(object, ...) length(object$ids)
#' @rdname repair-jackknife-methods
#' @exportS3Method stats::vcov
vcov.repair_jackknife <- function(object, ...) repair_jackknife_covariance(object)
#' @rdname repair-jackknife-methods
#' @exportS3Method stats::confint
confint.repair_jackknife <- function(object, parm, level = .95, ..., reference = c("normal", "t")) {
  table <- repair_jackknife_coefficients(object, level, match.arg(reference))
  out <- cbind(table$conf.low, table$conf.high)
  dimnames(out) <- list(table$term, paste0(100 * c((1 - level) / 2, (1 + level) / 2), "%"))
  if (!missing(parm)) {
    if (is.character(parm)) parm <- match(parm, table$term)
    if (!is.numeric(parm) || anyNA(parm) || any(parm < 1 | parm > nrow(out) | parm %% 1))
      repair_error("parm must identify existing coefficients.")
    out <- out[parm, , drop = FALSE]
  }
  out
}
#' @rdname repair-jackknife-methods
#' @exportS3Method base::summary
summary.repair_jackknife <- function(object, ...) {
  structure(list(method = "Experimental whole-person delete-one jackknife", n_people = length(object$ids),
    coefficients = if (all(object$diagnostics$usable)) repair_jackknife_coefficients(object) else NULL,
    bias = object$bias, diagnostics = object$diagnostics), class = "summary_repair_jackknife")
}
#' @rdname repair-jackknife-methods
#' @exportS3Method base::print
print.summary_repair_jackknife <- function(x, ...) {
  cat(x$method, "\n")
  cat(sum(x$diagnostics$usable), "/", x$n_people, "complete deletion vectors\n")
  if (!is.null(x$coefficients)) print(x$coefficients) else
    cat("Inference withheld: every deletion vector is required.\n")
  invisible(x)
}
#' @rdname repair-jackknife-methods
#' @exportS3Method base::print
print.repair_jackknife <- function(x, ...) {
  print(summary(x))
  invisible(x)
}
