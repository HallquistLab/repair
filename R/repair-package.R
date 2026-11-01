#' REPAIR: structural regression with estimated random coefficients
#'
#' Random-Effect Prediction Adjustment for Inference in Regression (REPAIR)
#' estimates relationships involving person- or cluster-specific processes
#' represented by random intercepts and slopes. Extracted predictions contain
#' sampling error, and shrunken BLUPs can distort relationships in a subsequent
#' regression. REPAIR uses likelihood-only predictions and their sampling-error
#' covariances in a modified method-of-moments regression.
#'
#' Discover effects with [repair_effects()], inspect matched inputs with
#' [repair_data()], fit with [repair()], and report results with
#' [repair_coefficients()]. [repair_example] provides simulated teaching data.
#'
#' The supported workflow uses Gaussian lmer models with one target group and
#' optional nuisance groups contained within targets, or supported single-group
#' linear lme models. Standard errors are conditional on fitted MLM quantities.
#' See `vignette("support-and-troubleshooting", package = "repair")` for the
#' support matrix, inference boundaries, diagnostics, and saved-fit guidance.
#' Use [repair_refit()] for verified single-model replay.
#' Run `citation("repair")` for software and manuscript references.
#' @author Zachary A. Vig and Michael N. Hallquist (co-maintainers).
#'   Zachary A. Vig is the primary R package maintainer contact.
#' @seealso [repair()], [repair_example]
#' @keywords internal
#' @useDynLib repair, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
