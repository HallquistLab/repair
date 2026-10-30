# Validate identities before prediction extraction; alignment to the retained scope happens
# afterwards. No positional fallback or partially specified zero blocks.
repair_cross_error_policy <- function(value, parsed) {
  models <- parsed$model
  if (length(models) == 1L) {
    if (!is.null(value)) repair_error("`cross_model_error` applies only when the formula references multiple models; leave it NULL for one model.")
    return(list(policy = "not_applicable", blocks = list()))
  }
  if (identical(value, "independent")) return(list(policy = "independent", blocks = list()))
  if (!is.list(value) || !length(value)) {
    repair_error("Cross-model prediction errors require an explicit `cross_model_error`: use \"independent\" or supply covariance blocks for every model pair.")
  }
  seen <- matrix(FALSE, length(models), length(models))
  for (block in value) {
    if (!is.list(block) || !setequal(names(block), c("models", "covariance")) || length(block) != 2L) {
      repair_error("Each cross-model block must contain exactly `models` and `covariance`.")
    }
    pair <- block$models
    if (!is.character(pair) || length(pair) != 2L || anyNA(pair) ||
        anyDuplicated(pair) || !all(pair %in% models)) {
      repair_error("Each block's `models` must name two distinct models referenced in the formula.")
    }
    index <- match(pair, models)
    if (seen[index[1], index[2]]) repair_error("Duplicate cross-model covariance pair: {.val {pair}}.")
    seen[index[1], index[2]] <- seen[index[2], index[1]] <- TRUE
    a <- block$covariance
    dn <- dimnames(a)
    valid_names <- function(x) !is.null(x) && !anyNA(x) &&
      all(nzchar(trimws(x))) && !anyDuplicated(x)
    if (!is.numeric(a) || length(dim(a)) != 3L || any(dim(a) == 0L) ||
        length(dn) != 3L || !all(vapply(dn, valid_names, logical(1))) || any(!is.finite(a))) {
      repair_error("Cross-model covariance arrays must be finite numeric 3D arrays with unique, nonempty coefficient and ID dimnames.")
    }
    for (j in 1:2) {
      terms <- vapply(Filter(function(e) e$model == pair[j], parsed$effects), `[[`, character(1), "term")
      if (!setequal(dn[[j]], terms)) {
        repair_error("Cross-model covariance term names for model {.val {pair[j]}} must exactly match its selected coefficients.")
      }
    }
  }
  if (!all(seen[upper.tri(seen)])) repair_error("Supply a covariance block for every pair of referenced models, including explicit zero blocks.")
  list(policy = "supplied", blocks = value)
}

repair_apply_cross_error <- function(covariances, policy, parsed, ids) {
  effect_models <- vapply(parsed$effects, `[[`, character(1), "model")
  effect_terms <- vapply(parsed$effects, `[[`, character(1), "term")
  for (block in policy$blocks) {
    a <- which(effect_models == block$models[1])
    b <- which(effect_models == block$models[2])
    missing <- setdiff(ids, dimnames(block$covariance)[[3L]])
    if (length(missing)) repair_error("Cross-model covariance for {.val {block$models}} is missing retained IDs: {.val {missing}}.")
    values <- block$covariance[effect_terms[a], effect_terms[b], ids, drop = FALSE]
    covariances[a, b, ] <- values
    covariances[b, a, ] <- aperm(values, c(2, 1, 3))
  }
  # Pairwise-valid covariances need not define a valid joint covariance when
  # three or more models are used. Check the entire selected-effect block.
  if (policy$policy == "supplied") {
    for (i in seq_along(ids)) {
      m <- covariances[, , i]
      scale <- max(abs(m))
      if (scale > 0 && min(eigen(m / scale, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
        repair_error("Assembled prediction-error covariance is not positive semidefinite for ID {.val {ids[i]}}. Check the supplied cross-model blocks.")
      }
    }
  }
  covariances
}
