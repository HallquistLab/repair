repair_bootstrap_require <- function(object) {
  if (!inherits(object, "repair_bootstrap")) repair_error("Supply a `repair_bootstrap()` result.")
  invisible(object)
}

repair_bootstrap_eligible <- function(object, allow_failed) {
  repair_bootstrap_require(object)
  if (!is.logical(allow_failed) || length(allow_failed) != 1L || is.na(allow_failed))
    repair_error("`allow_failed` must be TRUE or FALSE.")
  usable <- object$diagnostics$usable
  if (any(!usable) && !allow_failed)
    repair_error("Bootstrap inference withheld: {sum(!usable)} of {length(usable)} attempts failed. Inspect `repair_bootstrap_diagnostics()`; use `allow_failed = TRUE` only for inference conditional on successful replicates.")
  if (sum(usable) < 2L)
    repair_error("Bootstrap inference needs at least two complete usable coefficient vectors; {sum(usable)} are available.")
  object$replicates[usable, , drop = FALSE]
}

repair_bootstrap_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    repair_error("`level` must be a finite number strictly between zero and one.")
}

repair_bootstrap_metadata <- function(value, object, method) {
  attr(value, "method") <- method
  attr(value, "attempted") <- nrow(object$replicates)
  attr(value, "usable") <- sum(object$diagnostics$usable)
  attr(value, "failed") <- sum(!object$diagnostics$usable)
  attr(value, "conditional_on_success") <- any(!object$diagnostics$usable)
  attr(value, "inference") <- "Experimental person bootstrap"
  value
}

#' Bootstrap whole people through the complete REPAIR analysis
#'
#' Experimental serial cases bootstrap: sample people with replacement, copy
#' their complete effective histories with fresh target/nuisance IDs, and rerun
#' [repair_refit()]. Estimate uncertainty of the fitted procedure, including
#' measurement-model estimation. Coverage has not yet been established by simulation.
#' @inheritParams repair_refit
#' @param B Positive integer number of attempted resamples; default 1999.
#'   Failed attempts are retained and never replaced. One attempt is permitted
#'   for diagnostics, but inference requires at least two usable vectors.
#' @param seed Integer between zero and .Machine$integer.max; default 1.
#'   Uses saved L'Ecuyer-CMRG streams, Inversion normals, and Rejection sampling.
#'   Caller RNG kind and .Random.seed are restored on success or error. Box-Muller
#'   and user-supplied caller RNGs are rejected before changing state because
#'   their hidden state cannot be restored from .Random.seed.
#' @details Supports the single-model scope of [repair_refit()], including
#'   contained lmer nuisance effects and supported linear lme covariance structures.
#'   Measurement and analyzed cohorts must coincide. Reference replay and draw
#'   eligibility are checked before sampling. Each attempt samples N source IDs
#'   uniformly with replacement and refits measurement and structural models;
#'   no observations are resampled within a person. The reference cohort and its
#'   selection/missingness pattern define the sampling population.
#'
#'   Each attempt has a separate recorded stream, so fitting RNG use/failure in
#'   one attempt cannot change subsequent draws. Preflight has its own stream.
#'   Seeded runs share their first B attempts when extended to larger B in the
#'   same numerical environment. Cross-version/optimizer bitwise equivalence is
#'   not promised. Callbacks must be deterministic recipes as in [repair_refit()].
#'
#'   Errors invalidate the entire coefficient vector. Informational warnings and
#'   singular-but-usable lmer fits are recorded and retained. No retries, trimming,
#'   silent omissions, or replacement attempts are made. Inference is withheld
#'   when any attempt fails unless `allow_failed = TRUE` is explicitly passed to
#'   an inference accessor; that changes the interpretation to conditional on
#'   successful replicates. It does not make failure selection ignorable.
#'
#'   `vcov()` uses the sample covariance of the same complete usable vectors for
#'   every term. `confint()` uses marginal percentile intervals with R quantile
#'   type 7. Bootstrap variance is not added to the original conditional variance.
#'   No bootstrap p-values, BCa/studentized/simultaneous intervals, parallelism,
#'   checkpoint/resume, stratification, or multiple-model sampling are provided.
#'   Small B is for exercising code; increase B and assess Monte Carlo stability
#'   before interpreting tails. Software correctness is not a coverage guarantee.
#' @return A separate S3 `repair_bootstrap` result containing the original estimate
#'   and conditional coefficient table, a B-by-p coefficient matrix (failed rows
#'   are all NA), replicate diagnostics, source-ID draw indices/cohort, reference
#'   audit, static `specification` (formula text, effects, ID, measurement settings,
#'   structural backend and regularization controls), seed streams, and software
#'   versions. Replicate fitted models and copied
#'   raw datasets are not retained. The original fit is unchanged. Use
#'   [repair_bootstrap_coefficients()] and [repair_bootstrap_diagnostics()] for
#'   stable tables; `coef()`, `vcov()`, `confint()`, `nobs()`, `summary()` and
#'   `print()` are supported. See `vignette("person-bootstrap", package = "repair")`.
#' @examples
#' data(repair_example)
#' model <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, model, repair_example$people, "person")
#' # Tiny B demonstrates the API, not adequate tail precision.
#' boot <- repair_bootstrap(fit, repair_example$repeated,
#'   repair_example$people, B = 3, seed = 41)
#' summary(boot)
#' confint(boot, level = .90)
#' repair_bootstrap_coefficients(boot, conf.level = .90)
#' repair_bootstrap_diagnostics(boot)
#' @export
repair_bootstrap <- function(object, measurement_data, data, id = NULL,
                             refitter = NULL, B = 1999L, seed = 1L, tolerance = 1e-6) {
  repair_require_model(object); repair_require_fit(object)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1 ||
      B %% 1 != 0 || B >= .Machine$integer.max)
    repair_error("`B` must be one positive integer smaller than .Machine$integer.max.")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed %% 1 != 0 || seed > .Machine$integer.max)
    repair_error("`seed` must be one integer between zero and .Machine$integer.max.")
  B <- as.integer(B); seed <- as.integer(seed)
  old_kind <- RNGkind()
  if (old_kind[1L] == "user-supplied" || old_kind[2L] %in% c("Box-Muller", "user-supplied"))
    repair_error("Bootstrap RNG isolation cannot restore hidden Box-Muller or user-supplied RNG state. Select a standard RNG with Inversion normals first.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  streams <- vector("list", B + 1L)
  streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
  for (i in seq_len(B)) streams[[i + 1L]] <- parallel::nextRNGStream(streams[[i]])
  # The public refitter validates the original estimator before any sampling.
  replay <- repair_refit(object, measurement_data, data, id = id,
    refitter = refitter, tolerance = tolerance)
  a <- object@analysis_data; cohort <- a@ids
  if (is.null(id)) id <- a@diagnostics$id_column
  measurement <- repair_refit_frame(measurement_data, names(object@mlm_fit))
  measurement <- measurement[match(repair_refit_rows(object@mlm_fit[[1L]]), rownames(measurement)), , drop = FALSE]
  repair_refit_draw(object, measurement, data, id, cohort)
  counts <- vapply(cohort, function(x) sum(as.character(measurement[[a@effects$group[1L]]]) == x), integer(1))
  estimate <- stats::coef(object)
  draws <- matrix(NA_integer_, B, length(cohort))
  coefficients <- matrix(NA_real_, B, length(estimate),
    dimnames = list(as.character(seq_len(B)), names(estimate)))
  records <- vector("list", B)
  for (i in seq_len(B)) {
    assign(".Random.seed", streams[[i + 1L]], envir = .GlobalEnv)
    draws[i, ] <- sample.int(length(cohort), length(cohort), replace = TRUE)
    warnings <- messages <- character()
    started <- proc.time()[[3L]]
    result <- tryCatch(withCallingHandlers({
      fitted <- repair_refit(object, measurement_data, data, id = id,
        refitter = refitter, draw = cohort[draws[i, ]], tolerance = tolerance)
      beta <- stats::coef(fitted)
      if (!identical(names(beta), names(estimate)) || any(!is.finite(beta)))
        repair_error("Bootstrap refit produced an incompatible or nonfinite coefficient vector.")
      list(fit = fitted, beta = beta)
    }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
       message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }),
      error = function(e) e)
    ok <- !inherits(result, "error")
    singular <- NA; controls <- ""; n_analyzed <- NA_integer_
    preliminary <- final <- NA_real_
    if (ok) {
      coefficients[i, ] <- result$beta
      audit <- result$fit@analysis_data@diagnostics$refit
      warnings <- c(warnings, unlist(audit$warnings, use.names = FALSE))
      messages <- c(messages, unlist(audit$messages, use.names = FALSE))
      model <- result$fit@mlm_fit[[1L]]
      if (inherits(model, "lmerMod")) singular <- lme4::isSingular(model)
      controls <- audit$checks$controls[audit$checks$stage == "replacement"]
      n_analyzed <- as.integer(stats::nobs(result$fit))
      preliminary <- result$fit@struc_fit@diagnostics$correction_prelim
      final <- result$fit@struc_fit@diagnostics$correction_final
    } else {
      warnings <- c(warnings, result$repair_warnings)
      messages <- c(messages, result$repair_messages)
    }
    records[[i]] <- data.frame(replicate = i, usable = ok,
      stage = if (ok) "complete" else result$repair_stage %||% "refit",
      reason = if (ok) "" else conditionMessage(result),
      n_people = length(cohort), n_measurement = sum(counts[draws[i, ]]),
      n_analyzed = n_analyzed, measurement_singular = singular, controls = controls,
      preliminary_multiplier = preliminary, final_multiplier = final,
      elapsed = proc.time()[[3L]] - started,
      warnings = I(list(as.character(warnings))), messages = I(list(as.character(messages))),
      stringsAsFactors = FALSE)
  }
  diagnostics <- do.call(rbind, records); rownames(diagnostics) <- NULL
  reference <- replay@analysis_data@diagnostics$refit
  reference$spec <- NULL
  answer <- structure(list(estimate = estimate, conditional = repair_coefficients(object),
    replicates = coefficients, diagnostics = diagnostics, draws = draws, cohort = cohort,
    specification = list(structural_formula = paste(deparse(a@formula), collapse = " "),
      measurement_model = names(object@mlm_fit),
      measurement_formula = paste(deparse(stats::formula(object@mlm_fit[[1L]])), collapse = " "),
      id_column = id, effects = a@effects,
      measurement = repair_refit_configuration(object@mlm_fit[[1L]]),
      backend = object@struc_fit@diagnostics$backend,
      regularization = object@struc_fit@diagnostics[c("preliminary_moment_control", "final_moment_control")]),
    reference = reference, rng = list(seed = seed, kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
      preflight = streams[[1L]], replicates = streams[-1L]), versions = reference$versions,
    call = match.call()), class = "repair_bootstrap")
  if (any(!diagnostics$usable)) cli::cli_warn(
    "{sum(!diagnostics$usable)} of {B} bootstrap attempts failed; all attempts retained. Inference is withheld by default; inspect `repair_bootstrap_diagnostics()`.")
  answer
}

#' Extract experimental bootstrap inference as a table
#'
#' This separate table leaves the ordinary [repair_coefficients()] contract intact.
#' No conditional Wald statistics or p-values are assigned bootstrap meanings.
#' @param object A [repair_bootstrap()] result.
#' @param conf.level Confidence level strictly between zero and one.
#' @param allow_failed FALSE (default) withholds inference after any failed attempt.
#'   TRUE explicitly requests inference conditional on successful replicates.
#' @return A plain data frame with columns `term`, `estimate` (original fit),
#'   `std.error` (bootstrap SD), `conf.low`, `conf.high`, `conf.level`, `method`
#'   (percentile), `attempted`, `usable`, `failed`, `conditional_on_success`, and
#'   `inference` (experimental label). Rows follow original coefficient order.
#'   All terms use the same complete usable vectors. At least two are required;
#'   this numerical minimum does not establish adequate precision or validity.
#' @inherit repair_bootstrap examples
#' @export
repair_bootstrap_coefficients <- function(object, conf.level = .95, allow_failed = FALSE) {
  vectors <- repair_bootstrap_eligible(object, allow_failed)
  intervals <- stats::confint(object, level = conf.level, allow_failed = allow_failed)
  failed <- sum(!object$diagnostics$usable)
  data.frame(term = names(object$estimate), estimate = unname(object$estimate),
    std.error = sqrt(diag(stats::cov(vectors))), conf.low = intervals[, 1L],
    conf.high = intervals[, 2L], conf.level = conf.level, method = "percentile",
    attempted = nrow(object$replicates), usable = nrow(vectors), failed = failed,
    conditional_on_success = failed > 0L,
    inference = if (failed) "Experimental person bootstrap, conditional on successful replicates" else
      "Experimental person bootstrap", row.names = NULL)
}

#' Extract bootstrap attempts and source-ID maps
#'
#' Read saved results without fitting or sampling.
#' @param object A [repair_bootstrap()] result.
#' @param table Exact table name: `"replicates"` (default) or `"draws"`.
#' @return A plain data frame. `replicates` has one row per attempt with columns
#'   `replicate` (integer), `usable` (logical), `stage`, `reason` (character),
#'   `n_people`, `n_measurement`, `n_analyzed` (integer), `measurement_singular`
#'   (logical; NA for lme or unavailable), `controls` (character),
#'   `preliminary_multiplier`, `final_multiplier`, `elapsed` (numeric seconds),
#'   `warnings`, `messages` (list columns of character vectors).
#'   Failed fits have unavailable diagnostics marked NA; stages identify reference
#'   or replacement measurement/structural/data validation when available.
#'   `draws` has `replicate`, `draw` (integer), `source_id`, `refit_id` (character).
#'   Its key is (replicate, draw); generated IDs repeat across independently fitted
#'   replicates and must be interpreted with that key. Source IDs are never parsed.
#' @inherit repair_bootstrap examples
#' @export
repair_bootstrap_diagnostics <- function(object, table = "replicates") {
  repair_bootstrap_require(object)
  if (!is.character(table) || length(table) != 1L || is.na(table) ||
      !table %in% c("replicates", "draws")) repair_error("`table` must be 'replicates' or 'draws'.")
  if (table == "replicates") return(object$diagnostics)
  B <- nrow(object$draws); n <- length(object$cohort)
  data.frame(replicate = rep(seq_len(B), each = n), draw = rep(seq_len(n), B),
    source_id = object$cohort[as.vector(t(object$draws))],
    refit_id = rep(sprintf("r%07d", seq_len(n)), B), row.names = NULL)
}

#' Methods for experimental person-bootstrap results
#'
#' `coef()` returns original estimates; `vcov()` returns empirical bootstrap
#' covariance; `confint()` returns marginal percentile intervals (quantile type 7).
#' `summary()` stays inspectable when inference is withheld. `nobs()` counts
#' source people, not attempted fits or repeated observations.
#' @name repair-bootstrap-methods
#' @param object,x A [repair_bootstrap()] result.
#' @param parm Coefficient names or integer positions; missing selects all.
#' @param level Confidence level strictly between zero and one.
#' @param type Only `"percentile"` is supported.
#' @param allow_failed Explicitly permit inference conditional on successful
#'   replicates. FALSE by default; see [repair_bootstrap_coefficients()].
#' @param ... Reserved for generic compatibility.
#' @return `coef()` returns a named numeric vector; `vcov()` a p-by-p matrix;
#'   `confint()` a two-column matrix; `nobs()` an integer; `summary()` an
#'   inspectable list with counts, availability/reason, and an inference table
#'   (NULL when unavailable). `print()` returns its input invisibly.
#'   Covariance/interval matrices carry method, attempted/usable/failed counts,
#'   conditional-on-success status, and an experimental inference label as attributes.
NULL

#' @rdname repair-bootstrap-methods
#' @exportS3Method stats::coef
coef.repair_bootstrap <- function(object, ...) object$estimate

#' @rdname repair-bootstrap-methods
#' @exportS3Method stats::nobs
nobs.repair_bootstrap <- function(object, ...) length(object$cohort)

#' @rdname repair-bootstrap-methods
#' @exportS3Method stats::vcov
vcov.repair_bootstrap <- function(object, ..., allow_failed = FALSE) {
  repair_bootstrap_metadata(stats::cov(repair_bootstrap_eligible(object, allow_failed)),
    object, "empirical covariance")
}

#' @rdname repair-bootstrap-methods
#' @exportS3Method stats::confint
confint.repair_bootstrap <- function(object, parm, level = .95, type = "percentile",
                                     ..., allow_failed = FALSE) {
  repair_bootstrap_level(level)
  if (!identical(type, "percentile")) repair_error("Only `type = 'percentile'` is supported.")
  vectors <- repair_bootstrap_eligible(object, allow_failed)
  terms <- names(object$estimate)
  if (missing(parm)) parm <- terms
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(!is.finite(parm)) || any(parm %% 1 != 0) || any(parm < 1 | parm > length(terms)))
      repair_error("Invalid coefficient positions in `parm`.")
    parm <- terms[parm]
  }
  if (!is.character(parm) || anyNA(parm) || any(!parm %in% terms)) repair_error("Unknown coefficient names in `parm`.")
  probs <- c((1 - level) / 2, (1 + level) / 2)
  result <- matrix(NA_real_, length(parm), 2L,
    dimnames = list(parm, paste0(format(100 * probs, trim = TRUE), " %")))
  for (i in seq_along(parm)) result[i, ] <- stats::quantile(vectors[, parm[i]], probs, type = 7, names = FALSE)
  repair_bootstrap_metadata(result, object, "percentile")
}

#' @rdname repair-bootstrap-methods
#' @exportS3Method base::summary
summary.repair_bootstrap <- function(object, ..., level = .95, allow_failed = FALSE) {
  repair_bootstrap_level(level)
  # A valid summary remains useful even when the inference gate refuses access.
  table <- tryCatch(repair_bootstrap_coefficients(object, level, allow_failed), error = function(e) e)
  if (!is.logical(allow_failed) || length(allow_failed) != 1L || is.na(allow_failed))
    repair_error("`allow_failed` must be TRUE or FALSE.")
  available <- !inherits(table, "error")
  structure(list(n_people = length(object$cohort), attempted = nrow(object$replicates),
    usable = sum(object$diagnostics$usable), failed = sum(!object$diagnostics$usable),
    inference_available = available, reason = if (available) "" else conditionMessage(table),
    coefficients = if (available) table else NULL, level = level,
    conditional_on_success = any(!object$diagnostics$usable) && allow_failed,
    reference = object$reference, specification = object$specification,
    inference = "Experimental person-level cases bootstrap; marginal percentile intervals; coverage not yet validated."),
    class = "summary_repair_bootstrap")
}

#' @rdname repair-bootstrap-methods
#' @exportS3Method base::print
print.summary_repair_bootstrap <- function(x, ...) {
  cat("REPAIR person bootstrap (experimental)\n")
  cat(x$specification$structural_formula, "\n")
  cat("People:", x$n_people, "; attempts:", x$attempted, "; usable:", x$usable, "; failed:", x$failed, "\n")
  if (x$inference_available) print(x$coefficients, row.names = FALSE) else cat(x$reason, "\n")
  cat(x$inference, "\n")
  cat("Assess Monte Carlo stability before interpreting interval tails.\n")
  invisible(x)
}

#' @rdname repair-bootstrap-methods
#' @exportS3Method base::print
print.repair_bootstrap <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}
