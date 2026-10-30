# Posterior covariance uses the marginal observation covariance, avoiding G^-1.
# This permits singular G and rank-deficient within-unit designs.
repair_posterior_covariance <- function(g, z, r) {
  zg <- z %*% g
  ch <- chol(r + tcrossprod(zg, z))
  solved <- backsolve(ch, forwardsolve(t(ch), zg))
  v <- g - crossprod(zg, solved)
  (v + t(v)) / 2
}

repair_require_model <- function(object) {
  if (!S7::S7_inherits(object, class_repair_mod)) repair_error("Supply a fitted `repair()` result.")
  invisible(object)
}

#' Inspect posterior reliability of selected random coefficients
#'
#' Report plug-in posterior reliability using the manuscript's marginal and
#' intercept-residualized definitions. These diagnostics do not change the fit.
#' @param object A fitted [repair()] result, including a numerical structural failure.
#' @param detail `"summary"` (default) gives equally weighted scope summaries;
#'   `"person"` gives one row per unit, selected effect, and reliability measure.
#' @param scope `"structural"` (default) uses the structural analysis IDs for every
#'   model. `"all"` uses all effective grouping levels separately for each
#'   referenced model, including units excluded from structural analysis.
#' @param type `"both"` (default), `"marginal"`, or `"residualized"`. Residualized
#'   reliability is available for a selected slope only in a two-coefficient
#'   random-intercept/random-slope model. Other selected effects are returned
#'   with an explicit not-applicable status when this measure is requested alone.
#' @details For fitted random-effect covariance G and posterior covariance V_i,
#'   marginal reliability of coefficient k is `1 - V_i[k,k] / G[k,k]`.
#'   Posterior covariance is computed as
#'   `G - G Z_i' (R_i + Z_i G Z_i')^-1 Z_i G` using fitted effective observations
#'   and target-specific residual covariance. For nested lmer models, R_i includes
#'   every nuisance-group covariance contribution, as in likelihood prediction extraction.
#'   This supports singular G without inverting G.
#'
#'   For a slope residualized on the intercept, the contrast is
#'   `a = c(-G[1,2] / G[1,1], 1)` in intercept/slope order and reliability is
#'   `1 - (a' V_i a) / (a' G a)`. This uses the full covariance, not just diagonals.
#'   It does not residualize on observed structural predictors or other models.
#'
#'   Cohort reliability is the equally weighted mean across units, matching the
#'   manuscript. If any requested unit has undefined reliability, the overall
#'   `reliability` is NA; `mean_defined` summarizes only defined values and the
#'   counts make that denominator explicit. Zero or numerically degenerate contrast
#'   variance is undefined, not zero reliability. Singular G does not itself make
#'   all marginal reliabilities undefined. Tiny boundary deviations within 1e-8
#'   are rounded to the interval from zero to one; larger invalid values receive an explanatory status.
#'
#'   Posterior reliability describes uncertainty reduction under the fitted
#'   MLM. It is distinct from the likelihood-only prediction-error
#'   covariance used by REPAIR. A rank-deficient cluster can have defined posterior
#'   reliability while failing likelihood-only prediction extraction. Reliability is calculated
#'   within each MLM, without borrowing information across models
#'   or applying the cross-model prediction-error policy. It conditions on fitted
#'   parameters and supplies no confidence interval or automatic adequacy cutoff.
#' @return A data frame of class `repair_reliability`, with model, group, selected
#'   term, structural role, measure, scope, reliability, and denominator/status
#'   columns. Person rows also include ID, effective observation count, prior and
#'   posterior variance, and a reason for undefined results. Summary rows include
#'   `n_units`, `n_defined`, `n_undefined`, `mean_defined`, and distribution summaries.
#' @seealso [repair_diagnostics()], [plot.repair_reliability()]
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' repair_reliability(fit)
#' head(repair_reliability(fit, detail = "unit"))
#' @export
repair_reliability <- function(object, detail = c("summary", "unit"),
                               scope = c("structural", "all"),
                               type = c("both", "marginal", "residualized")) {
  repair_require_model(object)
  detail <- match.arg(detail); scope <- match.arg(scope); type <- match.arg(type)
  effects <- object@analysis_data@effects
  rows <- list()
  for (name in unique(effects$model)) {
    model <- object@mlm_fit[[name]]
    selected <- effects[effects$model == name, , drop = FALSE]
    elements <- repair_get_mlm_elements(model)
    group <- selected$group[1L]
    terms <- elements@random_terms[[group]]
    k <- length(terms)
    g <- as.matrix(elements@random_covariance[[group]][seq_len(k), seq_len(k), drop = FALSE])
    ids <- if (scope == "structural") object@analysis_data@ids else elements@group_levels[[group]]
    posterior <- lapply(ids, function(id) {
      obs <- elements@group_rows[[group]][[id]]
      tryCatch(list(v = repair_posterior_covariance(g,
        get_ranef_design(elements, group, obs, id), get_ranef_target_covariance(elements, group, obs)),
        n = length(obs), reason = ""),
        error = function(e) list(v = NULL, n = length(obs), reason = conditionMessage(e)))
    })
    for (j in seq_len(nrow(selected))) {
      term <- selected$term[j]; index <- match(term, terms)
      applicable <- k == 2L && "(Intercept)" %in% terms && term != "(Intercept)"
      measures <- if (type == "both") c("marginal", if (applicable) "residualized") else type
      for (measure in measures) {
        a <- numeric(k); a[index] <- 1
        status <- "ok"; reason <- ""
        if (measure == "residualized") {
          intercept <- match("(Intercept)", terms)
          if (!applicable) {
            status <- "not_applicable"; reason <- "Requires a slope in a two-coefficient intercept/slope model."
          } else if (g[intercept, intercept] <= 0) {
            status <- "undefined_variance"; reason <- "Intercept variance is zero; residualization is undefined."
          } else a[intercept] <- -g[intercept, index] / g[intercept, intercept]
        }
        prior <- if (status == "ok") drop(crossprod(a, g %*% a)) else NA_real_
        # Cancellation tolerance has the same units as the contrast variance.
        tolerance <- 100 * .Machine$double.eps * sum(abs(outer(a, a) * g))
        if (status == "ok" && (!is.finite(prior) || prior <= tolerance)) {
          status <- "undefined_variance"; reason <- "Prior contrast variance is zero or numerically degenerate."
        }
        for (i in seq_along(ids)) {
          s <- status; why <- reason; post <- value <- NA_real_
          if (s == "ok") {
            if (is.null(posterior[[i]]$v)) {
              s <- "posterior_failure"; why <- posterior[[i]]$reason
            } else {
              post <- drop(crossprod(a, posterior[[i]]$v %*% a))
              value <- 1 - post / prior
              if (!is.finite(value) || value < -1e-8 || value > 1 + 1e-8) {
                value <- NA_real_; s <- "invalid_reliability"
                why <- "Posterior/prior variance ratio is outside its numerical bounds."
              } else value <- min(1, max(0, value))
            }
          }
          rows[[length(rows) + 1L]] <- data.frame(model = name, group = group, term = term,
            role = selected$role[j], measure = measure, scope = scope, id = ids[i],
            n_observations = posterior[[i]]$n, prior_variance = prior,
            posterior_variance = post, reliability = value, status = s, reason = why)
        }
      }
    }
  }
  units <- do.call(rbind, rows)
  result <- if (detail == "unit") units else repair_summarize_reliability(units)
  structure(result, class = c("repair_reliability", "data.frame"), detail = detail, scope = scope)
}

repair_summarize_reliability <- function(units) {
  keys <- unique(units[c("model", "group", "term", "role", "measure", "scope")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    selected <- Reduce(`&`, lapply(names(key), function(n) units[[n]] == key[[n]]))
    p <- units[selected, , drop = FALSE]; ok <- is.finite(p$reliability)
    values <- p$reliability[ok]
    data.frame(key, n_units = nrow(p), n_defined = sum(ok), n_undefined = sum(!ok),
      reliability = if (all(ok)) mean(values) else NA_real_,
      mean_defined = if (any(ok)) mean(values) else NA_real_,
      min = if (any(ok)) min(values) else NA_real_,
      median = if (any(ok)) stats::median(values) else NA_real_,
      max = if (any(ok)) max(values) else NA_real_,
      status = if (all(ok)) "ok" else "undefined_values",
      reason = paste(unique(p$reason[nzchar(p$reason)]), collapse = "; "), row.names = NULL)
  }))
}

#' Plot unit-specific posterior reliability
#' @param x A unit-level result from [repair_reliability()].
#' @param type `"distribution"` for histograms or `"observations"` for reliability
#'   against effective observation count. Each model/effect/measure gets a panel.
#' @param ... Graphical parameters passed to `hist()` or `plot()`.
#' @details Undefined rows are omitted from drawing and counted in panel subtitles.
#'   No adequacy thresholds are drawn. Observation count alone does not determine
#'   reliability: design, residual covariance, and random-effect covariance matter.
#' @return The plotted data frame invisibly. Graphics parameters are restored.
#' @exportS3Method graphics::plot
plot.repair_reliability <- function(x, type = c("distribution", "observations"), ...) {
  type <- match.arg(type)
  if (!identical(attr(x, "detail"), "unit")) repair_error("Plot unit-level reliability: use `repair_reliability(fit, detail = \"unit\")`.")
  if (!nrow(x)) repair_error("No reliability rows to plot.")
  keys <- unique(x[c("model", "term", "measure")])
  old <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old)
    # Restoring mfrow resets cex/mex after par() has already visited them.
    graphics::par(cex = old$cex, mex = old$mex)
  }, add = TRUE)
  graphics::par(mfrow = c(ceiling(nrow(keys) / 2), min(2, nrow(keys))))
  for (i in seq_len(nrow(keys))) {
    p <- x[x$model == keys$model[i] & x$term == keys$term[i] & x$measure == keys$measure[i], ]
    keep <- is.finite(p$reliability)
    title <- paste0(keys$model[i], ": ", keys$term[i], "\n", keys$measure[i])
    subtitle <- paste(sum(keep), "defined /", nrow(p), "units;", attr(x, "scope"), "scope")
    if (!any(keep)) {
      graphics::plot.new(); graphics::title(main = title, sub = subtitle)
      graphics::text(.5, .5, "Reliability undefined; inspect status/reason")
    } else if (type == "distribution") {
      args <- utils::modifyList(list(x = p$reliability[keep], breaks = seq(0, 1, .1),
        xlim = c(0, 1), main = title, sub = subtitle, xlab = "Posterior reliability",
        cex.main = .85, cex.sub = .7), list(...))
      do.call(graphics::hist, args)
    } else {
      args <- utils::modifyList(list(x = p$n_observations[keep], y = p$reliability[keep],
        ylim = c(0, 1), main = title, sub = subtitle, xlab = "Effective observations",
        ylab = "Posterior reliability", cex.main = .85, cex.sub = .7), list(...))
      do.call(graphics::plot, args)
    }
  }
  invisible(x)
}
