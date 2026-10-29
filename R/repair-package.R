#' REPAIR: structural regression with estimated random coefficients
#'
#' Discover effects with [repair_effects()], inspect matched inputs with
#' [repair_data()], fit with [repair()], and report results with
#' [repair_coefficients()]. [repair_example] provides simulated teaching data.
#'
#' The supported workflow uses Gaussian single-group lmer or supported linear
#' lme models. Standard errors are conditional on fitted measurement quantities.
#' See `vignette("support-and-troubleshooting", package = "repair")` for the
#' support matrix, inference boundaries, diagnostics, and saved-fit guidance.
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
