# Centralized function for computing predictions from a mixed-model fit.
# repair_ranef <- S7::new_generic(
#   name = "repair_ranef",
#   dispatch_args = "model"
# )

# S7::method(repair_ranef, class_merMod) <- function(model, term, type = c("blup", "gls", "unshrink"), include_fixed = FALSE) {}
# S7::method(repair_ranef, class_lme) <- function(model, term, type = c("blup", "gls", "unshrink"), include_fixed = FALSE) {}


#' Extract likelihood-only random-effect scores.
#'
#' The direct method computes the closed-form GLS score for each cluster.
#' Currently only works for single term/simple residuals.
#'
#' @param model An `lme4` mixed-model fit.
#' @param term The term for which to extract random-effect scores.
#' @return A list containing cluster-level scores and covariance blocks.
#' @export
repair_ranef_gls <- function(mlm, method = c("direct", "unshrink")) {
  method <- match.arg(method)
  if (identical(method, "unshrink")) {
    repair_error("`method = \"unshrink\"` is not implemented yet.")
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    repair_error("The `lme4` package is required for direct scoring.")
  }
  if (!inherits(mlm, "lmerMod")) {
    repair_error("`mlm` must be a Gaussian lme4::lmer() fit.")
  }
  if (length(lme4::getME(mlm, "flist")) != 1L) {
    repair_error("Direct scoring currently requires exactly one grouping factor.")
  }

  # Extract fixed effects, random-effect designs, and cluster membership.
  random_effects <- lme4::ranef(mlm)
  group <- names(random_effects)[[1L]]
  grouping <- as.character(lme4::getME(mlm, "flist")[[group]])
  cluster_ids <- rownames(random_effects[[group]])
  cluster_rows <- split(seq_along(grouping), grouping, drop = TRUE)[cluster_ids]
  z_terms <- lme4::getME(mlm, "Ztlist")
  fixed_design <- lme4::getME(mlm, "X")
  fixed_coefs <- lme4::fixef(mlm)
  model_frame <- stats::model.frame(mlm)
  response <- stats::model.response(model_frame)
  offset <- stats::model.offset(model_frame) %||% rep(0, length(response))
  prior_weights <- stats::model.weights(model_frame) %||% rep(1, length(response))
  if (any(!is.finite(prior_weights)) || any(prior_weights <= 0)) {
    repair_error("Direct scoring requires finite, positive prior weights.")
  }
  resid_var <- stats::sigma(mlm)^2

  if (!is.finite(resid_var) || resid_var <= 0) {
    repair_error(
      "Direct scoring requires a finite, positive residual variance."
    )
  }

  random_names <- colnames(random_effects[[group]])
  nrandom <- length(random_names)

  scores <- matrix(
    NA_real_,
    length(cluster_ids),
    nrandom,
    dimnames = list(cluster_ids, random_names)
  )
  err_cov <- array(
    NA_real_,
    c(nrandom, nrandom, length(cluster_ids)),
    dimnames = list(random_names, random_names, cluster_ids)
  )
  score_diagnostics <- vector("list", length(cluster_ids))

  # Apply the closed-form GLS score independently within each cluster.
  for (cluster_pos in seq_along(cluster_ids)) {
    rows <- cluster_rows[[cluster_ids[[cluster_pos]]]]
    random_design <- do.call(
      cbind,
      lapply(z_terms, function(term) {
        as.numeric(term[cluster_pos, rows, drop = TRUE])
      })
    )
    residual <- response[rows] - offset[rows] - as.vector(
      fixed_design[rows, , drop = FALSE] %*% fixed_coefs
    )
    weighted_design <- random_design * sqrt(prior_weights[rows])
    design_rank <- qr(weighted_design, tol = 1e-10)$rank
    if (design_rank < nrandom) {
      score_diagnostics[[cluster_pos]] <- list(
        solver = "chol", rank = design_rank, dimension = nrandom,
        condition_number = Inf,
        error = "The within-cluster random-effect design was rank-deficient."
      )
      next
    }
    information <- crossprod(weighted_design) / resid_var
    score_rhs <- crossprod(
      weighted_design, residual * sqrt(prior_weights[rows])
    ) / resid_var
    score <- repair_solve(
      information, score_rhs,
      method = "chol"
    )
    score_diagnostics[[cluster_pos]] <- score@diagnostics
    if (!score@success) {
      score_diagnostics[[cluster_pos]]$error <- score@message
      next
    }
    scores[cluster_pos, ] <- as.vector(score@solution)
    err_cov[, , cluster_pos] <- score@inverse
  }
  failed_clusters <- cluster_ids[vapply(
    score_diagnostics, function(item) !is.null(item$error), logical(1)
  )]
  if (length(failed_clusters)) {
    repair_warn(paste0(
      "Likelihood-only scoring failed for {length(failed_clusters)} cluster(s); ",
      "their scores and covariances are missing. See `diagnostics$clusters`."
    ))
  }
  list(
    success = !length(failed_clusters),
    group = group,
    levels = cluster_ids,
    estimates = scores,
    err_cov = err_cov,
    score_cov = err_cov,
    method = method,
    diagnostics = list(
      step = "ranef_like",
      solver = "chol",
      failed_clusters = failed_clusters,
      clusters = stats::setNames(score_diagnostics, cluster_ids)
    )
  )
}


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
    covariance[positions, positions, drop = FALSE] <-
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
