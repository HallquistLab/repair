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
