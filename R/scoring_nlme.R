# A linear lme object stores source data separately from its effective rows.
# getData()/getResponse() may reinsert exclusions or re-evaluate a subset; neither
# defines a safe common row order for score extraction.
repair_nlme_structure <- function(model) {
  repair_check_pkg("nlme")
  if (inherits(model, "nlme")) repair_error("Nonlinear `nlme` fits are not supported; supply a linear `nlme::lme()` fit.")
  if (is.null(model$groups) || ncol(model$groups) != 1L ||
      length(model$modelStruct$reStruct) != 1L) {
    repair_error("nlme scoring currently requires exactly one grouping factor.")
  }
  correlation <- model$modelStruct$corStruct
  variance <- model$modelStruct$varStruct
  # nlme converts AR(1) with nonconsecutive integer times to corARMA(p=1,q=0).
  ar1 <- inherits(correlation, "corAR1") ||
    (inherits(correlation, "corARMA") && isTRUE(attr(correlation, "p") == 1) &&
      isTRUE(attr(correlation, "q") == 0))
  if (!is.null(correlation) && !ar1 && !inherits(correlation, "corCAR1")) {
    repair_error("Supported nlme correlations are independent, corAR1, and corCAR1.")
  }
  if (!is.null(variance) && !inherits(variance, "varIdent") && !inherits(variance, "varPower")) {
    repair_error("Supported nlme variance structures are constant, varIdent, and varPower.")
  }
  invisible(NULL)
}

S7::method(repair_get_mlm_elements, class_lme) <- function(model) {
  repair_nlme_structure(model)
  source_data <- model$data
  if (!is.data.frame(source_data)) {
    repair_error("nlme extraction requires stored data. Refit with an explicit `data` frame and `keep.data = TRUE`.")
  }
  ids <- rownames(model$groups)
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) ||
      anyDuplicated(rownames(source_data)) || anyNA(match(ids, rownames(source_data)))) {
    repair_error("Cannot align the nlme fitted observation IDs with its stored data.")
  }
  effective_data <- source_data[match(ids, rownames(source_data)), , drop = FALSE]
  # nlme drops unused levels before fitting; restore the fitted contrasts rather
  # than the caller's current options(contrasts=...).
  for (name in names(effective_data)) {
    if (is.factor(effective_data[[name]])) effective_data[[name]] <- droplevels(effective_data[[name]])
  }
  fixed_frame <- stats::model.frame(stats::delete.response(model$terms),
    data = effective_data, na.action = stats::na.fail)
  fitted_contrasts <- model$contrasts %||% list()
  fixed_design <- stats::model.matrix(stats::delete.response(model$terms),
    data = fixed_frame,
    contrasts.arg = fitted_contrasts[intersect(names(fitted_contrasts), names(fixed_frame))])
  fixed_coefs <- nlme::fixef(model)
  random_vars <- unique(unlist(lapply(stats::formula(model$modelStruct$reStruct, asList = TRUE), all.vars)))
  random_design <- stats::model.matrix(model$modelStruct$reStruct,
    data = effective_data, contrast = fitted_contrasts[intersect(names(fitted_contrasts), random_vars)])
  ranef <- nlme::ranef(model)
  terms <- colnames(ranef)
  group <- names(model$groups)[[1L]]
  grouping <- model$groups[[1L]]
  levels <- rownames(ranef)
  if (!identical(colnames(fixed_design), names(fixed_coefs)) ||
      !identical(colnames(random_design), terms) || anyDuplicated(terms) ||
      nrow(fixed_design) != length(ids) || nrow(random_design) != length(ids) ||
      !identical(rownames(fixed_design), ids) || !identical(rownames(random_design), ids) ||
      any(!is.finite(fixed_design)) || any(!is.finite(random_design))) {
    repair_error("Could not reconstruct nlme design matrices with the fitted coefficient names and observation order.")
  }
  if (!identical(rownames(model$fitted), ids) || !identical(rownames(model$residuals), ids)) {
    repair_error("nlme fitted values, residuals, and grouping rows are not aligned.")
  }
  response <- as.numeric(model$fitted[, 1L] + model$residuals[, 1L])
  reconstructed_fixed <- drop(fixed_design %*% fixed_coefs)
  reconstructed_conditional <- reconstructed_fixed + rowSums(random_design *
    as.matrix(ranef[match(as.character(grouping), levels), , drop = FALSE]))
  close <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-7, check.attributes = FALSE))
  if (!close(reconstructed_fixed, model$fitted[, 1L]) ||
      !close(reconstructed_conditional, model$fitted[, ncol(model$fitted)])) {
    repair_error("Reconstructed nlme designs do not reproduce fitted values. Stored data or formula dependencies may have changed; refit using explicit data columns.")
  }
  if (any(!is.finite(response))) repair_error("The fitted nlme response contains nonfinite values.")
  group_rows <- repair_group_rows(grouping, levels)
  correlation <- model$modelStruct$corStruct
  # Correlation and variance structures are fitted in group-sorted order. With
  # one group and no finer correlation grouping, within-group order is unchanged.
  if (!is.null(correlation)) {
    correlation_groups <- nlme::getGroups(correlation)
    if (!identical(as.character(correlation_groups), as.character(grouping[order(grouping)]))) {
      repair_error("The nlme correlation groups must match the random-effects grouping factor; finer or shared correlation groups are not supported.")
    }
    correlation_blocks <- nlme::corMatrix(correlation)
    if (!is.list(correlation_blocks) || !setequal(names(correlation_blocks), levels)) {
      repair_error("Could not align nlme correlation blocks by grouping ID.")
    }
  }
  # These fitted SDs are stored by nlme in effective (original) row order. They
  # also preserve varPower(~ fitted(.)) without recomputing its covariate.
  residual_sd <- attr(model$residuals, "std")
  if (!is.numeric(residual_sd) || length(residual_sd) != length(ids) ||
      any(!is.finite(residual_sd)) || any(residual_sd <= 0)) {
    repair_error("The nlme fit must contain finite positive residual SDs for every effective observation.")
  }
  blocks <- stats::setNames(lapply(levels, function(level) {
    rows <- group_rows[[level]]
    corr <- if (is.null(correlation)) diag(length(rows)) else correlation_blocks[[level]]
    if (!identical(dim(corr), rep(length(rows), 2L)) || any(!is.finite(corr))) {
      repair_error("An nlme correlation block has incompatible dimensions or nonfinite values.")
    }
    covariance <- as.matrix(corr) * outer(residual_sd[rows], residual_sd[rows])
    dimnames(covariance) <- list(ids[rows], ids[rows])
    covariance
  }), levels)
  g <- as.matrix(nlme::getVarCov(model, type = "random.effects"))
  if (!identical(dim(g), rep(length(terms), 2L)) || any(!is.finite(g))) {
    repair_error("The nlme random-effect covariance does not match the fitted terms.")
  }
  operator <- list(type = "nlme_conditional", group = group, blocks = blocks,
    matrix = function(group = NULL, level, rows = NULL) {
      block <- blocks[[level]]
      if (is.null(block)) repair_error("Unknown nlme residual block: {.val {level}}.")
      if (!is.null(rows)) {
        positions <- match(rows, group_rows[[level]])
        if (anyNA(positions)) repair_error("Rows do not belong to the requested nlme residual block.")
        block <- block[positions, positions, drop = FALSE]
      }
      block
    })
  class_mlm_elements(
    backend = "nlme", model = model,
    ranef = stats::setNames(list(ranef), group),
    fixed_design = fixed_design, fixed_coefs = fixed_coefs, response = response,
    offset = numeric(), observation_weights = numeric(),
    random_design = stats::setNames(list(random_design), group),
    random_terms = stats::setNames(list(terms), group),
    grouping = stats::setNames(list(grouping), group),
    group_levels = stats::setNames(list(levels), group),
    group_rows = stats::setNames(list(group_rows), group), group_ptr = NULL,
    random_covariance = stats::setNames(list(kronecker(diag(length(levels)), g)), group),
    residual_covariance = operator,
    diagnostics = list(nobs = length(ids), nfixed = ncol(fixed_design), groups = group,
      effective_row_ids = ids, residual_group = group,
      correlation = if (is.null(correlation)) "independent" else class(correlation)[[1L]],
      variance = if (is.null(model$modelStruct$varStruct)) "constant" else class(model$modelStruct$varStruct)[[1L]],
      conditional_on_stage1 = TRUE)
  )
}
