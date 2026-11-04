#' Propagate observation error covariance between measurement models
#'
#' Convert explicitly assumed observation-level cross-covariances to the
#' likelihood-only prediction-error blocks accepted by [repair()]. No covariance
#' parameters are estimated and no measurement model is jointly refitted.
#' @param formula Structural formula selecting random coefficients from at least
#'   two models, as in [repair()].
#' @param models Named list of fitted models or [repair_mlm()] specifications.
#'   Specifications are fitted once in this call; reuse fitted models when
#'   constructing several scenarios.
#' @param observation_keys Named list, one entry per referenced model. Each entry
#'   is a character vector of occasion keys, named by the original data row names.
#'   Effective fitted rows are matched by name, including after subset/NA handling;
#'   extra rows are allowed. Keys must be unique within each target person.
#'   Equal keys across models identify the same occasion for the same person.
#'   Supply meaningful keys explicitly: row position and nearest times are never
#'   used to pair observations.
#' @param residual_cross_covariance List with one block for every unordered
#'   referenced model pair. Each block has `models = c("a", "b")` and
#'   `covariance`, a named list of observation covariance matrices, one per ID.
#'   Matrix row/column names must exactly cover that person's effective occasion
#'   keys in models a/b. Values are in observation units. Order does not matter.
#'   Include covariance of any effective nuisance contributions, but exclude the
#'   target random effects. Supply either this argument or `rho`.
#' @param rho Optional scalar correlation in `[-1, 1]`, restricted to exactly two
#'   models with iid homoscedastic errors and no nuisance random groups.
#'   Assumes covariance `rho * sigma_a * sigma_b` at equal occasion keys and zero
#'   at all other occasion pairs, including unmatched occasions. Marginal scales
#'   come from the fitted models. This is an assumption, not an estimate.
#' @param ids Character vector of unique target IDs to propagate. Default is the
#'   intersection of effective target levels across referenced models. Supply the
#'   analysis ID column to require coverage of your intended scope.
#' @details For each person, uses the same target design Z and effective marginal
#'   error covariance R as likelihood-only scoring, with extraction operator
#'   `L = (Z' R^-1 Z)^-1 Z' R^-1`. Returns selected entries of `L_a R_ab L_b'`.
#'   Fitted fixed effects and covariance parameters are treated as fixed.
#'   Covariance between true random coefficients is excluded from R_ab.
#'
#'   The complete joint observation covariance, using the fitted within-model
#'   blocks, must be positive semidefinite (relative tolerance 1e-8), including
#'   when there are more than two models. The assembled selected prediction
#'   covariance is also checked. Invalid inputs and singular extraction designs
#'   stop; there is no silent covariance repair or omission. For correlated or
#'   heterogeneous errors, supply compatible explicit observation blocks instead
#'   of `rho`. Separate marginal fits alone do not identify these blocks.
#' @return A plain list of `list(models, covariance)` blocks for
#'   `cross_model_error` in [repair()] or [repair_data()]. Arrays have dimensions
#'   selected coefficient in model a, selected coefficient in model b, and ID.
#' @seealso [repair_sensitivity()], `vignette("multiple-models")`
#' @examples
#' if (requireNamespace("nlme", quietly = TRUE)) {
#'   case <- readRDS(system.file("extdata", "multimodel-workflow-reference.rds",
#'                              package = "repair"))$cases$supplied
#'   models <- list(growth = nlme::lme(y ~ z, random = ~ z | id, data = case$a),
#'     response = nlme::lme(q ~ w, random = ~ w | subject, data = case$b))
#'   equation <- re(w, "response") ~ re(z, "growth")
#'   # These two time scales encode the same synthetic measurement occasions.
#'   keys <- list(growth = setNames(as.character(round(4 * case$a$z)), rownames(case$a)),
#'     response = setNames(as.character(round(5 * (case$b$w - .5))), rownames(case$b)))
#'   blocks <- repair_cross_covariance(equation, models, keys, rho = .25,
#'                                     ids = case$people$person)
#'   repair(equation, models, case$people, "person", cross_model_error = blocks)
#' }
#' @export
repair_cross_covariance <- function(formula, models, observation_keys,
                                    residual_cross_covariance = NULL,
                                    rho = NULL, ids = NULL) {
  models <- normalize_repair_models(models)
  parsed <- parse_struc_formula(formula, repair_effects(models))
  selected <- parsed$model
  if (length(selected) < 2L) repair_error("Select coefficients from at least two models.")
  if (is.null(rho) == is.null(residual_cross_covariance)) {
    repair_error("Supply exactly one of `rho` and `residual_cross_covariance`.")
  }
  if (!is.null(rho) && (!is.numeric(rho) || length(rho) != 1L ||
      !is.finite(rho) || abs(rho) > 1 || length(selected) != 2L)) {
    repair_error("`rho` must be one finite correlation in [-1, 1] for exactly two models.")
  }
  if (!is.list(observation_keys) || !repair_cov_names(names(observation_keys)) ||
      !setequal(names(observation_keys), selected)) {
    repair_error("`observation_keys` must be a uniquely named list covering the referenced models.")
  }
  info <- stats::setNames(lapply(selected, function(name) {
    effects <- Filter(function(e) e$model == name, parsed$effects)
    group <- effects[[1L]]$group
    nuisance <- repair_validate_target(models[[name]], group, name)
    elements <- repair_get_mlm_elements(models[[name]])
    rows <- repair_refit_rows(models[[name]])
    keys <- observation_keys[[name]]
    if (!is.character(keys) || !is.null(dim(keys)) ||
        !repair_cov_names(names(keys)) || anyNA(keys) || any(!nzchar(trimws(keys))) ||
        !all(rows %in% names(keys))) {
      repair_error("Occasion keys for model {.val {name}} must be nonempty character values named by all effective fitted row names.")
    }
    if (!is.null(rho) && length(nuisance)) {
      repair_error("`rho` requires iid homoscedastic errors without nuisance random groups; supply explicit observation covariance blocks.")
    }
    list(elements = elements, group = group, keys = keys[rows],
      terms = vapply(effects, `[[`, character(1), "term"))
  }), selected)
  available <- lapply(info, function(x) x$elements@group_levels[[x$group]])
  if (is.null(ids)) ids <- Reduce(intersect, available)
  if (!is.character(ids) || !repair_cov_names(ids) ||
      !all(vapply(available, function(x) all(ids %in% x), logical(1)))) {
    repair_error("`ids` must be nonempty unique character IDs present in every referenced model.")
  }
  pairs <- utils::combn(selected, 2L, simplify = FALSE)
  if (!is.null(residual_cross_covariance)) {
    supplied <- repair_observation_pairs(residual_cross_covariance, selected, ids)
  }
  result <- lapply(pairs, function(pair) list(models = pair,
    covariance = array(0, c(length(info[[pair[1]]]$terms), length(info[[pair[2]]]$terms), length(ids)),
      dimnames = list(info[[pair[1]]]$terms, info[[pair[2]]]$terms, ids))))
  marginal_scale <- stats::setNames(rep(NA_real_, length(selected)), selected)
  for (i in seq_along(ids)) {
    id <- ids[i]
    local <- lapply(selected, function(name) {
      x <- info[[name]]; e <- x$elements
      rows <- e@group_rows[[x$group]][[id]]
      keys <- unname(x$keys[rows])
      if (anyDuplicated(keys)) repair_error("Duplicate occasion keys for model {.val {name}}, ID {.val {id}}.")
      r <- get_ranef_target_covariance(e, x$group, rows)
      z <- get_ranef_design(e, x$group, rows, id)
      weighted <- repair_solve(r, z, fail_on_error = FALSE)
      if (!isTRUE(weighted@success)) repair_error("Observation covariance could not be solved for model {.val {name}}, ID {.val {id}}.")
      operator <- repair_solve(crossprod(z, weighted@solution), t(weighted@solution), fail_on_error = FALSE)
      if (!isTRUE(operator@success)) repair_error("Prediction extraction design could not be solved for model {.val {name}}, ID {.val {id}}.")
      if (!is.null(rho)) {
        variance <- r[1L, 1L]
        if (!is.finite(variance) || variance <= 0 ||
            max(abs(r - diag(variance, nrow(r)))) > 1e-10 * variance ||
            (!is.na(marginal_scale[name]) && abs(variance - marginal_scale[name]) > 1e-10 * variance)) {
          repair_error("`rho` requires iid homoscedastic errors; supply explicit observation covariance blocks for model {.val {name}}.")
        }
        marginal_scale[name] <<- variance
      }
      list(r = r, keys = keys, l = operator@solution[match(x$terms, e@random_terms[[x$group]]), , drop = FALSE])
    })
    names(local) <- selected
    positions <- split(seq_len(sum(vapply(local, function(x) nrow(x$r), integer(1)))),
      rep(selected, vapply(local, function(x) nrow(x$r), integer(1))))
    joint <- matrix(0, sum(lengths(positions)), sum(lengths(positions)))
    for (name in selected) joint[positions[[name]], positions[[name]]] <- local[[name]]$r
    for (j in seq_along(pairs)) {
      pair <- pairs[[j]]; a <- local[[pair[1]]]; b <- local[[pair[2]]]
      if (!is.null(rho)) {
        cross <- outer(a$keys, b$keys, `==`) * rho * sqrt(a$r[1, 1] * b$r[1, 1])
      } else {
        block <- supplied[[j]]
        cross <- block$covariance[[id]]
        if (block$models[1L] != pair[1L]) cross <- t(cross)
        if (!is.numeric(cross) || !is.matrix(cross) || any(!is.finite(cross)) ||
            !repair_cov_names(rownames(cross)) || !repair_cov_names(colnames(cross)) ||
            !setequal(rownames(cross), a$keys) || !setequal(colnames(cross), b$keys)) {
          repair_error("Observation covariance names/dimensions must exactly match effective occasion keys for {.val {pair}}, ID {.val {id}}.")
        }
        cross <- cross[a$keys, b$keys, drop = FALSE]
      }
      joint[positions[[pair[1]]], positions[[pair[2]]]] <- cross
      joint[positions[[pair[2]]], positions[[pair[1]]]] <- t(cross)
      result[[j]]$covariance[, , i] <- a$l %*% cross %*% t(b$l)
    }
    repair_cov_psd(joint, paste0("Joint observation covariance for ID ", id))
    # Validate the complete selected-effect block separately, using propagation
    # for the fitted marginal blocks as well as the returned cross blocks.
    effect_sizes <- vapply(local, function(x) nrow(x$l), integer(1))
    pos <- split(seq_len(sum(effect_sizes)), rep(selected, effect_sizes))
    prediction <- matrix(0, sum(effect_sizes), sum(effect_sizes))
    for (name in selected) {
      x <- local[[name]]
      prediction[pos[[name]], pos[[name]]] <- x$l %*% x$r %*% t(x$l)
    }
    for (j in seq_along(pairs)) {
      pair <- pairs[[j]]; cross <- result[[j]]$covariance[, , i]
      prediction[pos[[pair[1]]], pos[[pair[2]]]] <- cross
      prediction[pos[[pair[2]]], pos[[pair[1]]]] <- t(matrix(cross, effect_sizes[pair[1]], effect_sizes[pair[2]]))
    }
    repair_cov_psd(prediction, paste0("Joint prediction-error covariance for ID ", id))
  }
  result
}

repair_cov_names <- function(x) {
  is.character(x) && length(x) > 0L && !anyNA(x) &&
    all(nzchar(trimws(x))) && !anyDuplicated(x)
}

repair_cov_psd <- function(x, context) {
  scale <- max(abs(x))
  if (any(!is.finite(x)) || (scale > 0 &&
      min(eigen(x / scale, symmetric = TRUE, only.values = TRUE)$values) < -1e-8)) {
    repair_error("{context} is not positive semidefinite; check the assumed cross-model blocks.")
  }
}

repair_observation_pairs <- function(blocks, models, ids) {
  pairs <- utils::combn(models, 2L, simplify = FALSE)
  if (!is.list(blocks) || length(blocks) != length(pairs)) {
    repair_error("Supply one observation covariance block for every referenced model pair.")
  }
  result <- vector("list", length(pairs))
  for (block in blocks) {
    if (!is.list(block) || length(block) != 2L ||
        !setequal(names(block), c("models", "covariance")) ||
        !is.character(block$models) || length(block$models) != 2L ||
        anyNA(block$models) || anyDuplicated(block$models) || !all(block$models %in% models)) {
      repair_error("Each observation block must contain `models` naming two distinct referenced models and `covariance`.")
    }
    j <- which(vapply(pairs, function(pair) setequal(pair, block$models), logical(1)))
    if (!is.null(result[[j]])) repair_error("Duplicate observation covariance model pair.")
    if (!is.list(block$covariance) || !repair_cov_names(names(block$covariance)) ||
        !all(ids %in% names(block$covariance))) {
      repair_error("Observation covariance matrices must be in a uniquely ID-named list covering all requested IDs.")
    }
    result[[j]] <- block
  }
  result
}
