get_ranef_design <- function(elements, group, rows, level = NULL) {
  design <- elements@random_design[[group]]
  terms <- elements@random_terms[[group]]
  levels <- elements@group_levels[[group]]
  grouping <- as.character(elements@grouping[[group]])
  nterms <- length(terms)
  nlevels <- length(levels)

  if (ncol(design) == nterms) {
    if (!is.null(level)) {
      if (!level %in% levels) {
        repair_error(
          "Grouping level {.val {level}} was not found for grouping ",
          "factor {.val {group}}."
        )
      }

      return(
        design[
          rows, ,
          drop = FALSE
        ]
      )
    }

    expanded_design <- do.call(
      cbind,
      lapply(
        levels,
        function(group_level) {
          indicator <- as.numeric(
            grouping[rows] == group_level
          )

          design[
            rows, ,
            drop = FALSE
          ] * indicator
        }
      )
    )

    design_names <- paste(
      rep(levels, each = nterms),
      rep(terms, times = nlevels),
      sep = "::"
    )

    colnames(expanded_design) <- design_names

    return(expanded_design)
  }

  expected_columns <- nterms * nlevels

  if (ncol(design) != expected_columns) {
    repair_error(
      "The random-effects design for grouping factor ",
      "{.val {group}} has incompatible dimensions."
    )
  }

  if (is.null(level)) {
    return(
      design[
        rows, ,
        drop = FALSE
      ]
    )
  }

  level_index <- match(level, levels)

  if (is.na(level_index)) {
    repair_error(
      "Grouping level {.val {level}} was not found for grouping ",
      "factor {.val {group}}."
    )
  }

  columns <- (
    (level_index - 1L) * nterms +
      seq_len(nterms)
  )

  design[
    rows,
    columns,
    drop = FALSE
  ]
}

get_ranef_base_covariance <- function(elements, rows) {
  operator <- elements@residual_covariance
  if (identical(elements@backend, "lme4")) {
    return(
      operator$matrix(group = NULL, level = NULL, rows = rows)
    )
  }

  residual_group <- operator$group
  if (is.null(residual_group)) {
    repair_error(
      "The residual covariance operator does not identify its ",
      "grouping factor."
    )
  }
  residual_values <- as.character(
    elements@grouping[[residual_group]][rows]
  )
  if (anyNA(residual_values)) {
    repair_error(
      "Target rows contain missing residual grouping values."
    )
  }
  covariance <- matrix(0, nrow = length(rows), ncol = length(rows))

  for (level in unique(residual_values)) {
    positions <- which(residual_values == level)
    block_rows <- elements@group_rows[[residual_group]][[level]]
    local_rows <- match(rows[positions], block_rows)
    if (anyNA(local_rows)) {
      repair_error(
        "Could not align target rows with the residual covariance ",
        "block for grouping level {.val {level}}."
      )
    }
    block <- operator$matrix(group = residual_group, level = level)
    covariance[positions, positions] <-
      block[local_rows, local_rows, drop = FALSE]
  }

  covariance
}

get_ranef_predictions <- function(elements, group, include_fixed = FALSE) {
  # group-level information
  terms <- elements@random_terms[[group]]
  levels <- elements@group_levels[[group]]
  rows_by_level <- elements@group_rows[[group]]

  # initialize output containers
  predictions <- matrix(
    NA_real_,
    nrow = length(levels),
    ncol = length(terms),
    dimnames = list(levels, terms)
  )
  error_covariance <- array(
    NA_real_,
    dim = c(
      length(terms),
      length(terms),
      length(levels)
    ),
    dimnames = list(terms, terms, levels)
  )
  diagnostics <- stats::setNames(
    vector("list", length(levels)),
    levels
  )
  nuisance_groups <- setdiff(
    names(elements@random_design),
    group
  )

  for (level_index in seq_along(levels)) {
    level <- levels[[level_index]]
    rows <- rows_by_level[[level]]
    target_design <- get_ranef_design(
      elements = elements,
      group = group,
      rows = rows,
      level = level
    )
    target_covariance <- get_ranef_base_covariance(
      elements = elements,
      rows = rows
    )

    for (nuisance_group in nuisance_groups) {
      nuisance_design <- get_ranef_design(
        elements = elements,
        group = nuisance_group,
        rows = rows
      )
      nuisance_covariance <- elements@random_covariance[[nuisance_group]]

      if (ncol(nuisance_design) != nrow(nuisance_covariance)) {
        repair_error(
          "The random-effects design and covariance dimensions do not ",
          "match for grouping factor {.val {nuisance_group}}."
        )
      }

      target_covariance <- target_covariance +
        nuisance_design %*%
          nuisance_covariance %*%
          t(nuisance_design)
    }

    residual <- elements@response[rows]
    if (length(elements@offset)) {
      residual <- residual - elements@offset[rows]
    }
    residual <- residual - as.vector(
      elements@fixed_design[rows, , drop = FALSE] %*%
        elements@fixed_coefs
    )

    weighted_design <- repair_solve(
      lhs = target_covariance,
      rhs = target_design,
      method = "auto",
      fail_on_error = FALSE
    )
    weighted_residual <- repair_solve(
      lhs = target_covariance,
      rhs = residual,
      method = "auto",
      fail_on_error = FALSE
    )

    diagnostics[[level]] <- list(
      covariance = target_covariance,
      design_solver = weighted_design@diagnostics,
      response_solver = weighted_residual@diagnostics
    )

    if (!isTRUE(weighted_design@success) ||
          !isTRUE(weighted_residual@success)) {
      diagnostics[[level]]$error <- paste(
        c(
          weighted_design@message,
          weighted_residual@message
        ),
        collapse = "; "
      )
      next
    }

    information <- crossprod(target_design, weighted_design@solution)
    prediction_rhs <- crossprod(target_design, weighted_residual@solution)
    prediction <- repair_solve(
      lhs = information,
      rhs = prediction_rhs,
      method = "auto",
      fail_on_error = FALSE
    )
    diagnostics[[level]]$information_solver <- prediction@diagnostics

    if (!isTRUE(prediction@success)) {
      diagnostics[[level]]$error <- prediction@message
      next
    }

    predictions[level_index, ] <- as.numeric(prediction@solution)
    error_covariance[, , level_index] <- prediction@inverse
  }

  if (isTRUE(include_fixed)) {
    fixed_indices <- match(
      terms,
      names(elements@fixed_coefs)
    )
    if (anyNA(fixed_indices)) {
      repair_error("`include_fixed = TRUE` requires matching fixed coefficients for every random coefficient.")
    }
    available <- !is.na(fixed_indices)

    if (any(available)) {
      predictions[, available] <- sweep(
        predictions[, available, drop = FALSE], 2L,
        elements@fixed_coefs[fixed_indices[available]], "+"
      )
    }
  }

  failed_levels <- names(diagnostics)[
    vapply(
      diagnostics,
      function(item) !is.null(item$error),
      logical(1)
    )
  ]

  if (length(failed_levels)) {
    repair_warn(
      "Likelihood-only predictions failed for ",
      "{length(failed_levels)} cluster(s); their scores and ",
      "covariances are missing. See {.field diagnostics}."
    )
  }

  list(
    success = !length(failed_levels),
    predictions = predictions,
    error_covariance = error_covariance,
    terms = terms,
    levels = levels,
    diagnostics = list(
      success = !length(failed_levels),
      step = "ranef_gls",
      solver = "auto",
      failed_levels = failed_levels,
      per_level = diagnostics
    )
  )
}

#' Extract likelihood-only random-coefficient scores
#'
#' Compute target-specific GLS scores after subtracting fitted fixed effects
#' and offsets. The target covariance excludes the target random effects and
#' includes residual and nuisance-group covariance contributions.
#' @param mlm A Gaussian `lme4::lmer()` or `nlme::lme()` fit.
#' @param target_group Grouping factor to score; NULL scores all groups.
#' @param include_fixed Add matching fixed coefficients to the scores. Default
#'   FALSE gives deviations. TRUE requires a fixed coefficient for every random
#'   coefficient; it does not reconstruct covariate-dependent individual effects.
#' @details Failed clusters have missing predictions and covariance blocks, issue
#'   a warning, and appear in `@diagnostics[[group]]$failed_levels`. A singular
#'   random-effect covariance is allowed when the target design is identifiable.
#'   Standard errors condition on fitted measurement-model parameters.
#'
#'   Single-group linear nlme fits must retain their data (`keep.data = TRUE`).
#'   Supported correlations are independent, `corAR1`, and `corCAR1`; an AR(1)
#'   represented internally as `corARMA(p = 1, q = 0)` is also accepted.
#'   Supported residual variance structures are constant, `varIdent`, and
#'   `varPower` (including fitted-value covariates). Correlation groups must
#'   match the random-effects grouping factor. Other structures, nonlinear
#'   nlme fits, and multilevel nlme fits are rejected explicitly.
#'
#'   nlme designs use stored effective observation IDs and fitted contrasts.
#'   Residual covariance uses fitted residual SDs and correlation blocks;
#'   `na.exclude` padding and re-evaluation of `subset` cannot add observations.
#'   See `vignette("nlme-workflow", package = "repair")` for reproduction examples.
#'
#'   Multi-group lme4 scoring remains experimental. Per-cluster blocks do not
#'   describe cross-cluster score errors from shared nuisance effects. The
#'   structural [repair()] workflow therefore requires one grouping factor.
#' @return An S7 score object with group-keyed `@predictions` matrices (ID by
#'   coefficient), `@error_covariance` arrays (coefficient by coefficient by ID),
#'   `@levels`, `@terms`, and `@diagnostics`; `@success` is true only if all selected
#'   clusters scored successfully. These are sampling-error covariances, not
#'   posterior BLUP covariances.
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' scores <- repair_ranef_gls(growth, target_group = "person")
#' head(scores@predictions$person)
#' scores@error_covariance$person[, , 1]
#' @export
repair_ranef_gls <- function(mlm, target_group = NULL, include_fixed = FALSE) {
  if (!is.null(target_group) &&
    (!is.character(target_group) ||
      length(target_group) != 1L || is.na(target_group) ||
      !nzchar(target_group))) {
    repair_error(
      "`target_group` must be NULL or one nonempty character string."
    )
  }
  if (!is.logical(include_fixed) || length(include_fixed) != 1L ||
    is.na(include_fixed)) {
    repair_error("`include_fixed` must be TRUE or FALSE.")
  }

  if (!inherits(mlm, "lmerMod") && !inherits(mlm, "lme")) {
    repair_error("`mlm` must be a Gaussian lmer or lme fit; GLMMs are not supported.")
  }
  elements <- repair_get_mlm_elements(mlm)
  group_names <- names(elements@random_terms)

  if (!length(group_names)) {
    repair_error(
      "The fitted model has no random-effects grouping factors."
    )
  }
  if (!is.null(target_group) && !target_group %in% group_names) {
    repair_error(
      "Grouping factor {.val {target_group}} was not found in the model."
    )
  }

  selected_groups <- if (is.null(target_group)) {
    group_names
  } else {
    target_group
  }

  group_results <- lapply(
    selected_groups,
    function(group) {
      get_ranef_predictions(
        elements = elements,
        group = group,
        include_fixed = include_fixed
      )
    }
  )
  names(group_results) <- selected_groups

  group_success <- vapply(
    group_results,
    function(result) isTRUE(result$success),
    logical(1)
  )

  class_repair_ranef(
    success = all(group_success),
    predictions = stats::setNames(
      lapply(
        group_results,
        function(result) result$predictions
      ),
      selected_groups
    ),
    error_covariance = stats::setNames(
      lapply(
        group_results,
        function(result) result$error_covariance
      ),
      selected_groups
    ),
    method = "gls",
    group = target_group,
    levels = stats::setNames(
      lapply(
        group_results,
        function(result) result$levels
      ),
      selected_groups
    ),
    terms = stats::setNames(
      lapply(
        group_results,
        function(result) result$terms
      ),
      selected_groups
    ),
    include_fixed = include_fixed,
    diagnostics = stats::setNames(
      lapply(
        group_results,
        function(result) result$diagnostics
      ),
      selected_groups
    )
  )
}
