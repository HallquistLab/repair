#' Synthetic practice and learning data for a REPAIR analysis
#'
#' A completely simulated example of 60 people with 8--12 daily learning predictions
#' and a unit-level measure of prior practice. These are teaching data, not
#' measurements from real participants or a dataset used in the REPAIR paper.
#' @format A list with three elements:
#' \describe{
#'   \item{repeated}{A data frame with character `person` IDs, numeric `day`
#'     (days since baseline, starting at zero), and numeric `score`
#'     (arbitrary learning-score units).}
#'   \item{people}{A data frame with one row per `person` and numeric `practice`
#'     (simulated prior practice hours per week, between zero and ten).}
#'   \item{truth}{Generating constants and a `random_effects` data frame with
#'     unit IDs, intercept deviations, and slope deviations. For teaching
#'     comparisons only; these values must not enter an empirical REPAIR fit.}
#' }
#' @details The generating model is
#'   `score = 50 + b0 + (1 + b1) * day + error`, with residual SD 1.5,
#'   `b0 = 3 * z0`, and `b1 = 0.18 * (practice - 5) + 0.15 * z0 + 0.35 * z1`,
#'   where z0 and z1 are independent standard normal variables. Practice is
#'   generated independently from a uniform distribution and rounded to one
#'   decimal place. The population association of slope deviation with one
#'   extra hour of practice is 0.18 score units per day. One finite sample need
#'   not recover this constant exactly. Fitted predictions are deviations from the
#'   estimated fixed effects, so their intercept/centering can differ from the
#'   generating population coefficients.
#'
#'   IDs are shared between the two tables. Observation counts vary, and each
#'   unit has an identifiable intercept/slope design. Loading the saved data
#'   does not generate random numbers or require the companion research repository.
#' @source Generated with base R using seed 20261029, Mersenne-Twister/Inversion/
#'   Rejection RNG settings. The source repository contains the complete
#'   `dev/generate_example_data.R` script. These data are separate from the frozen
#'   numerical reference fixtures used to test the implementation.
#' @seealso [repair()], [repair_effects()], [repair_data()]
#' @examples
#' data(repair_example)
#' head(repair_example$repeated)
#' head(repair_example$people)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' repair_coefficients(fit)
"repair_example"
