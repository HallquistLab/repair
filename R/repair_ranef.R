#' Extract likelihood-only random-effect scores.
#'
#' The direct method computes the closed-form GLS score for each cluster.
#' Currently only works for single term/simple residuals.
#'
#' @param model An `lme4` mixed-model fit.
#' @param term The term for which to extract random-effect scores.
#' @return A list containing cluster-level scores and covariance blocks.
#' @export
repair_ranef_gls <- function(model, term, supported = c("merMod", )) {
  model_class <- class(model)
  if (!inherits(model, "merMod")) {
    repair_error("`model` must be an lme4 mixed-model fit.")
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    repair_error("The `lme4` package is required for direct scoring.")
  }

  # Extract fixed effects, random-effect designs, and cluster membership.
  random_effects <- lme4::ranef(model)
  group <- names(random_effects)[[1L]]
  grouping <- as.character(lme4::getME(model, "flist")[[group]])
  cluster_ids <- rownames(random_effects[[group]])
  cluster_rows <- split(seq_along(grouping), grouping, drop = TRUE)[cluster_ids]
  z_terms <- lapply(lme4::getME(model, "Ztlist"), as.matrix)
  fixed_design <- lme4::getME(model, "X")
  fixed_coefs <- lme4::fixef(model)
  response <- stats::model.response(stats::model.frame(model))
  resid_var <- stats::sigma(model)^2

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
    residual <- response[rows] - as.vector(
      fixed_design[rows, , drop = FALSE] %*% fixed_coefs
    )

    score <- tryCatch(
      {
        information <- crossprod(random_design) / resid_var
        score_rhs <- crossprod(random_design, residual) / resid_var
        chol_result <- repair_solve(information, score_rhs, method = "chol")
        list(
          score = as.vector(chol_result@solution),
          covariance = chol_result@inverse,
          diagnostics = chol_result@diagnostics
        )
      }, error = function(error) {
        list(
          score = rep(NA_real_, nrandom),
          covariance = matrix(NA_real_, nrandom, nrandom),
          diagnostics = list(
            solver = "chol",
            rank = NA_integer_,
            dimension = nrandom,
            condition_number = Inf,
            error = conditionMessage(error)
          )
        )
      }
    )
    scores[cluster_pos, ] <- score$score
    err_cov[, , cluster_pos] <- score$covariance
    score_diagnostics[[cluster_pos]] <- score$diagnostics
  }
  list(
    group = group,
    levels = cluster_ids,
    estimates = scores,
    # TODO: why both?
    err_cov = err_cov,
    score_cov = err_cov,
    diagnostics = list(
      step = "ranef_like",
      solver = "chol",
      clusters = stats::setNames(score_diagnostics, cluster_ids)
    )
  )
