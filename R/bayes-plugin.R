# Deliberately wrap supported fits rather than implicitly accepting every model
# with a brmsfit/stanreg class. The snapshot is a plug-in Gaussian likelihood,
# never a collection of marginal posterior random-effect moments.
repair_is_bayes_plugin <- function(model) inherits(model, "repair_bayes_plugin")

repair_findbars <- function(formula) {
  if (requireNamespace("reformulas", quietly = TRUE)) reformulas::findbars(formula) else lme4::findbars(formula)
}

repair_nobars <- function(formula) {
  if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(formula) else lme4::nobars(formula)
}

repair_bayes_formula <- function(formula, data) {
  fail <- function() repair_error(paste(
    "Experimental Bayesian plug-in support requires a numeric response, additive",
    "numeric columns, and exactly one `(terms | group)` block with a bare grouping column.",
    "Precompute transformations; factors, interactions, offsets and special terms are unsupported."
  ))
  if (!inherits(formula, "formula") || length(formula) != 3L || !is.symbol(formula[[2L]])) fail()
  walk <- function(x) {
    if (is.symbol(x)) return(invisible(NULL))
    if (is.numeric(x) && length(x) == 1L && x %in% c(0, 1)) return(invisible(NULL))
    if (!is.call(x) || !is.symbol(x[[1L]]) ||
        !as.character(x[[1L]]) %in% c("+", "-", "(", "|")) fail()
    lapply(as.list(x)[-1L], walk)
    invisible(NULL)
  }
  walk(formula[[3L]])
  bars <- repair_findbars(formula)
  if (length(bars) != 1L || !is.symbol(bars[[1L]][[3L]])) fail()
  group <- as.character(bars[[1L]][[3L]])
  fields <- all.vars(formula)
  if (!is.data.frame(data) || anyDuplicated(names(data)) || !all(fields %in% names(data))) {
    repair_error("Bayesian plug-in extraction requires the fitted data columns to be retained in the model.")
  }
  numeric_fields <- setdiff(fields, group)
  if (group %in% all.vars(repair_nobars(formula)) ||
      group %in% all.vars(bars[[1L]][[2L]]) ||
      !all(vapply(data[numeric_fields], function(x) is.numeric(x) && is.null(dim(x)), logical(1)))) fail()
  if (!(is.factor(data[[group]]) || is.character(data[[group]]) || is.numeric(data[[group]]))) fail()
  list(formula = formula, group = group, fields = fields)
}

repair_bayes_name <- function(x) replace(x, x == "Intercept", "(Intercept)")

repair_bayes_depth_limit <- function(fit) {
  if (methods::is(fit, "stanfit") && length(fit@stan_args)) {
    limits <- vapply(fit@stan_args, function(x) x$control$max_treedepth %||% NA_real_, numeric(1))
    if (all(is.finite(limits)) && length(unique(limits)) == 1L) return(limits[[1L]])
  }
  NA_real_
}

repair_bayes_brms <- function(model) {
  repair_check_pkg("brms")
  f <- model$formula
  if (!inherits(f, "brmsformula") || inherits(f, "mvbrmsformula") ||
      !identical(model$family$family, "gaussian") || !identical(model$family$link, "identity") ||
      !identical(model$algorithm, "sampling") || isTRUE(attr(f$formula, "nl")) ||
      length(f$pforms) || length(f$pfix) || length(model$autocor) ||
      length(model$stanvars) || length(model$data2) || any(nzchar(model$stan_funs))) {
    repair_error(paste("Experimental brms support requires a sampled, univariate Gaussian identity-link MLM",
      "with constant independent residual variance; distributional/nonlinear models, autocorrelation,",
      "custom Stan code and additional data are unsupported."))
  }
  spec <- repair_bayes_formula(f$formula, model$data)
  re <- model$ranef
  if (!is.data.frame(re) || !nrow(re) || !all(c("id", "group", "coef", "by", "cov", "dist") %in% names(re)) ||
      !identical(unique(re$group), spec$group) || length(unique(re$id)) != 1L ||
      any(re$by != "") || any(re$cov != "") || any(re$dist != "gaussian")) {
    repair_error("Experimental brms support requires exactly one ordinary Gaussian random-effects block.")
  }
  s <- brms::standata(model)
  block <- unique(re$id)
  z_names <- paste0("Z_", block, "_", seq_len(nrow(re)))
  allowed <- c("N", "Y", "K", "Kc", "X", "prior_only", z_names,
    paste0(c("J_", "N_", "M_", "NC_"), block))
  if (length(setdiff(names(s), allowed)) || !identical(as.integer(s$prior_only), 0L) ||
      !all(c("N", "Y", "X", z_names, paste0("J_", block)) %in% names(s))) {
    repair_error(paste("Unsupported brms likelihood features or Stan data layout:",
      "weights, censoring, truncation, measurement error and residual extensions are unsupported."))
  }
  fixed <- brms::fixef(model, summary = FALSE)
  if (!is.matrix(fixed) || !ncol(fixed)) repair_error("Bayesian plug-in extraction requires at least one fixed coefficient.")
  x <- s$X
  colnames(x) <- repair_bayes_name(colnames(x))
  colnames(fixed) <- repair_bayes_name(colnames(fixed))
  if (!identical(colnames(x), colnames(fixed))) {
    repair_error("brms fitted fixed-effect names do not match its native design matrix.")
  }
  levels <- attr(re, "levels")[[spec$group]]
  index <- s[[paste0("J_", block)]]
  if (is.null(levels) || anyNA(index) || any(index < 1L | index > length(levels))) {
    repair_error("Cannot align brms native grouping indices with fitted grouping IDs.")
  }
  terms <- repair_bayes_name(re$coef)
  z <- do.call(cbind, s[z_names]); colnames(z) <- terms
  variance <- brms::VarCorr(model, summary = FALSE)
  g <- variance[[spec$group]]$cov
  if (length(dim(g)) != 3L ||
      !identical(repair_bayes_name(dimnames(g)[[2L]]), terms) ||
      !identical(repair_bayes_name(dimnames(g)[[3L]]), terms)) {
    repair_error("brms covariance draws do not match its fitted random coefficients.")
  }
  sigma <- variance$residual__$sd
  if (!is.matrix(sigma) || ncol(sigma) != 1L) repair_error("brms must retain one residual-SD draw per posterior draw.")
  list(spec = spec, x = x, beta = colMeans(fixed), y = s$Y, z = z,
    grouping = levels[index], levels = levels, terms = terms,
    g = apply(g, c(2L, 3L), mean), sigma2 = mean(sigma^2),
    priors = brms::prior_summary(model), draws = brms::as_draws_array(model),
    sampler = brms::nuts_params(model), max_treedepth = repair_bayes_depth_limit(model$fit))
}

repair_bayes_rstanarm <- function(model) {
  repair_check_pkg("rstanarm")
  if (inherits(model, "stanmvreg") ||
      !isTRUE(model$stan_function %in% c("stan_lmer", "stan_glmer")) ||
      !identical(model$family$family, "gaussian") || !identical(model$family$link, "identity") ||
      !identical(model$algorithm, "sampling") || is.null(model$glmod)) {
    repair_error("Experimental rstanarm support requires a sampled Gaussian identity-link `stan_lmer()` or `stan_glmer()` fit.")
  }
  if (length(model$weights) || length(model$prior.weights) || length(model$offset) ||
      !is.null(model$call$weights) || !is.null(model$call$offset)) {
    repair_error("Weights and offsets are unsupported for experimental Bayesian plug-in fits.")
  }
  if (!is.null(model$call$prior_PD) && !identical(model$call$prior_PD, FALSE)) {
    repair_error("Prior-predictive rstanarm fits (or nonliteral `prior_PD` settings) are unsupported.")
  }
  frame <- model$glmod$fr
  spec <- repair_bayes_formula(stats::formula(model), frame)
  re <- model$glmod$reTrms
  if (length(re$cnms) != 1L || length(re$flist) != 1L ||
      !identical(names(re$flist), spec$group)) {
    repair_error("Experimental rstanarm support requires exactly one random-effects block and grouping factor.")
  }
  x <- model$glmod$X
  if (!is.matrix(x) || !ncol(x)) repair_error("Bayesian plug-in extraction requires at least one fixed coefficient.")
  terms <- re$cnms[[1L]]
  grouping <- re$flist[[1L]]
  levels <- levels(grouping)
  # Collapse the stored sparse, level-specific design without rebuilding formulas
  # or applying the caller's current factor/contrast settings.
  zt <- as.matrix(re$Zt)
  if (nrow(zt) != length(terms) * length(levels)) repair_error("Invalid rstanarm random-effects design dimensions.")
  z <- vapply(seq_along(terms), function(k) colSums(zt[seq.int(k, nrow(zt), by = length(terms)), , drop = FALSE]), numeric(ncol(zt)))
  colnames(z) <- terms
  draws <- as.matrix(model)
  needed <- c(colnames(x), "sigma")
  if (!all(needed %in% colnames(draws))) repair_error("rstanarm must retain all shared fixed-effect and residual draws.")
  g <- matrix(0, length(terms), length(terms), dimnames = list(terms, terms))
  for (j in seq_along(terms)) for (k in seq_len(j)) {
    name <- paste0("Sigma[", spec$group, ":", terms[j], ",", terms[k], "]")
    if (!name %in% colnames(draws)) repair_error("Missing rstanarm covariance draws: {.val {name}}.")
    g[j, k] <- g[k, j] <- mean(draws[, name])
  }
  sp <- rstan::get_sampler_params(model$stanfit, inc_warmup = FALSE)
  sampler <- do.call(rbind, lapply(seq_along(sp), function(chain) {
    data.frame(Chain = chain, Parameter = rep(colnames(sp[[chain]]), each = nrow(sp[[chain]])),
      Value = as.vector(sp[[chain]]))
  }))
  list(spec = spec, x = x, beta = colMeans(draws[, colnames(x), drop = FALSE]),
    y = model$y, z = z, grouping = as.character(grouping), levels = levels, terms = terms,
    g = g, sigma2 = mean(draws[, "sigma"]^2), priors = model$prior.info,
    draws = rstanarm::as_draws_array(model), sampler = sampler,
    max_treedepth = repair_bayes_depth_limit(model$stanfit))
}

repair_bayes_elements <- function(parts, model, backend) {
  x <- as.matrix(parts$x); z <- as.matrix(parts$z)
  n <- length(parts$y); group <- parts$spec$group; terms <- parts$terms
  levels <- as.character(parts$levels); g <- as.matrix(parts$g)
  if (!n || nrow(x) != n || nrow(z) != n || length(parts$grouping) != n ||
      !identical(colnames(x), names(parts$beta)) || anyDuplicated(terms) ||
      anyDuplicated(levels) || anyNA(parts$grouping) || any(!nzchar(trimws(levels))) ||
      !setequal(as.character(parts$grouping), levels) || ncol(z) != length(terms) ||
      !identical(dim(g), rep(length(terms), 2L)) ||
      any(!is.finite(c(x, z, parts$beta, parts$y, g, parts$sigma2))) ||
      length(parts$sigma2) != 1L || parts$sigma2 <= 0) {
    repair_error("Extracted Bayesian plug-in components have invalid values, identities or dimensions.")
  }
  if (qr(x)$rank != ncol(x) || any(eigen(g, symmetric = TRUE, only.values = TRUE)$values < -1e-8)) {
    repair_error("Bayesian plug-in extraction requires a full-rank fixed design and positive-semidefinite random covariance.")
  }
  dimnames(g) <- list(terms, terms)
  sigma2 <- parts$sigma2
  operator <- list(type = "independent", group = NULL,
    matrix = function(group = NULL, level = NULL, rows) diag(sigma2, length(rows)))
  key <- function(value) stats::setNames(list(value), group)
  # Only the catalogue shape is needed; these are not posterior predictions.
  catalogue <- matrix(NA_real_, length(levels), length(terms), dimnames = list(levels, terms))
  class_mlm_elements(backend = backend, model = model, ranef = key(catalogue),
    fixed_design = x, fixed_coefs = parts$beta, response = as.numeric(parts$y),
    offset = numeric(), observation_weights = numeric(), random_design = key(z),
    random_terms = key(terms), grouping = key(parts$grouping), group_levels = key(levels),
    group_rows = key(repair_group_rows(parts$grouping, levels)), group_ptr = NULL,
    random_covariance = key(kronecker(diag(length(levels)), g)), residual_covariance = operator,
    diagnostics = list(nobs = n, groups = group, residual_variance = sigma2,
      conditional_on_stage1 = TRUE, experimental = TRUE,
      summary = "posterior means of fixed coefficients and complete covariance matrices; E(sigma^2 | y)"))
}

#' Opt in to experimental Bayesian measurement-model plug-in support
#'
#' Snapshot shared parameter summaries from a narrowly supported Bayesian Gaussian
#' MLM for use by [repair()], [repair_data()], and [repair_ranef_gls()].
#' @param model An already sampled `brms::brm()` Gaussian identity-link fit or
#'   `rstanarm::stan_lmer()` / Gaussian identity-link `rstanarm::stan_glmer()` fit.
#' @details This is an explicit experimental opt-in, not Bayesian REPAIR.
#'   The initial scope is a univariate numeric response, additive numeric columns,
#'   exactly one ordinary `(terms | group)` random-effects block, a bare grouping
#'   column, and constant independent Gaussian residual variance. Factors in the
#'   measurement equation, transformations, interactions, offsets, weights,
#'   independent multiple blocks (`||`), nuisance groups, nonidentity links,
#'   distributional/nonlinear models, residual extensions and custom Stan code
#'   are unsupported. Precompute numeric transformations before fitting.
#'
#'   The snapshot uses posterior means of fixed coefficients and complete
#'   random-effect covariance matrices, and `mean(sigma^2)`, not `mean(sigma)^2`.
#'   Likelihood-only GLS predictions are reconstructed from native designs;
#'   marginal posterior random-effect means/SDs never enter the correction.
#'   REPAIR scores remain deviations from fitted fixed effects. Each unit must
#'   have an identifiable likelihood-only random-coefficient design even when
#'   its Bayesian posterior is proper.
#'
#'   Structural SEs, tests and confidence intervals remain conditional and
#'   model-based. Step-1 parameter uncertainty, its induced cross-unit dependence,
#'   and prior sensitivity are not propagated. Full refitted-workflow coverage
#'   has not been established. Sampler diagnostics are recorded and warnings
#'   retained; they cannot establish structural coverage or model adequacy.
#'   Priors on shared parameters still affect the inputs. [repair_refit()] and
#'   [repair_bootstrap()] do not yet support Bayesian plug-in snapshots.
#'   See `vignette("bayesian-plugin", package = "repair")`.
#' @return A `repair_bayes_plugin` snapshot containing the original fit,
#'   Gaussian likelihood components, backend/version, priors, sampler diagnostics,
#'   summary convention and explicit inference limitations. Recreate the snapshot
#'   after changing a fit. The brms/rstanarm packages are optional dependencies.
#' @examples
#' \dontrun{
#' library(rstanarm)
#' growth <- rstanarm::stan_lmer(Reaction ~ Days + (Days | Subject),
#'   data = lme4::sleepstudy, chains = 4, seed = 42)
#' measurement <- repair_bayes_plugin(growth)
#' people <- data.frame(Subject = levels(lme4::sleepstudy$Subject), x = 1:18)
#' fit <- repair(re("Days") ~ x, measurement, people, "Subject")
#' summary(fit)
#' }
#' @export
repair_bayes_plugin <- function(model) {
  repair_check_pkg("posterior")
  backend <- if (inherits(model, "brmsfit")) "brms" else if (inherits(model, "stanreg")) "rstanarm" else
    repair_error("Supply a supported sampled brmsfit or rstanarm stanreg model.")
  parts <- if (backend == "brms") repair_bayes_brms(model) else repair_bayes_rstanarm(model)
  elements <- repair_bayes_elements(parts, model, backend)
  draws <- posterior::as_draws_array(parts$draws)
  if (posterior::nchains(draws) < 2L || posterior::niterations(draws) < 4L) {
    repair_error("Bayesian plug-in snapshots require at least two chains and four retained iterations per chain for diagnostics.")
  }
  diagnostics <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  sampler <- parts$sampler
  divergences <- sum(sampler$Value[sampler$Parameter == "divergent__"])
  depth <- sampler$Value[sampler$Parameter == "treedepth__"]
  depth_hits <- if (is.finite(parts$max_treedepth)) sum(depth >= parts$max_treedepth) else NA_integer_
  issues <- character()
  if (any(!is.finite(diagnostics$rhat)) || any(diagnostics$rhat > 1.01, na.rm = TRUE))
    issues <- c(issues, "Rhat is undefined or exceeds 1.01 for one or more parameters.")
  if (any(!is.finite(diagnostics$ess_bulk)) || any(!is.finite(diagnostics$ess_tail)) ||
      any(diagnostics$ess_bulk < 400 | diagnostics$ess_tail < 400, na.rm = TRUE))
    issues <- c(issues, "Bulk or tail ESS is undefined or below 400 for one or more parameters.")
  if (divergences > 0) issues <- c(issues, paste(divergences, "divergent transitions."))
  if (is.na(depth_hits)) issues <- c(issues, "Configured maximum tree depth is unavailable; inspect the original sampler fit.") else
    if (depth_hits > 0) issues <- c(issues, paste(depth_hits, "maximum tree-depth hits."))
  limitation <- paste("Experimental Bayesian plug-in REPAIR: conditional model-based inference;",
    "step-1 parameter uncertainty is not propagated; full-workflow coverage is unvalidated.")
  audit <- list(mlm_backend = backend, mlm_singular = NA,
    mlm_messages = c(limitation, issues), experimental = TRUE,
    inference = "conditional_model_based", stage1_uncertainty = "not_propagated",
    coverage = "full_refitted_workflow_unvalidated", summary = elements@diagnostics$summary,
    package_version = as.character(utils::packageVersion(backend)), priors = parts$priors,
    sampler = list(parameters = as.data.frame(diagnostics), divergences = divergences,
      max_observed_treedepth = max(depth), max_treedepth = parts$max_treedepth,
      max_treedepth_hits = depth_hits,
      chains = posterior::nchains(draws), iterations_per_chain = posterior::niterations(draws),
      issues = issues))
  repair_warn(paste(c(limitation, issues), collapse = " "))
  structure(list(fit = model, elements = elements, diagnostics = audit), class = "repair_bayes_plugin")
}

#' @exportS3Method base::print
print.repair_bayes_plugin <- function(x, ...) {
  cat("Experimental REPAIR Bayesian plug-in snapshot (", x$diagnostics$mlm_backend, ")\n", sep = "")
  cat("Summary:", x$diagnostics$summary, "\n")
  cat(paste(x$diagnostics$mlm_messages, collapse = "\n"), "\n")
  invisible(x)
}

