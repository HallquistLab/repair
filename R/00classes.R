# S7 classes returned by the package.

#' Combined first-stage and structural repair model.
#' @param mlm_fit The first-stage multilevel-model fit.
#' @param struc_fit The structural M-MOM fit.
#' @param analysis_data Matched structural data from [repair_data()].
#' @param call The original call.
#' @details An advanced S7 class used by [repair()]. Obtain validated fitted
#'   objects through [repair()] rather than constructing them manually. The
#'   exported class remains available for type checks and advanced integration;
#'   its property layout is not the stable reporting interface. Use the public
#'   generics and table extractors for routine analysis.
#' @return An S7 object containing MLM fits, the structural fit,
#'   assembled analysis data, and the original call. Manually supplying these
#'   properties does not perform the statistical validation done by [repair()].
#' @examples
#' data(repair_example)
#' growth <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
#' fit <- repair(re("day") ~ practice, growth, repair_example$people, "person")
#' S7::S7_inherits(fit, class_repair_mod)
#' coef(fit)
#' @export
class_repair_mod <- S7::new_class(
  "class_repair_mod",
  properties = list(
    mlm_fit = S7::class_any,
    struc_fit = S7::class_any,
    analysis_data = S7::class_any,
    call = S7::class_any
  )
)

# The assembly boundary has explicit dimensions and IDs, including when inspected
# or modified by a developer. Matrix/covariance contents are validated by M-MOM.
class_repair_data <- S7::new_class(
  "class_repair_data",
  properties = list(
    formula = S7::class_any, ids = S7::class_character,
    outcome = S7::class_numeric, predictors = S7::class_any,
    err_x = S7::class_any, err_y = S7::class_numeric, err_xy = S7::class_any,
    effects = S7::class_data.frame, predictions = S7::class_any,
    diagnostics = S7::class_list
  ),
  validator = function(self) {
    n <- length(self@ids)
    q <- ncol(self@predictors)
    if (anyNA(self@ids) || anyDuplicated(self@ids)) return("IDs must be nonmissing and unique.")
    if (!is.matrix(self@predictors) || is.null(q) || q < 1L ||
        nrow(self@predictors) != n || length(self@outcome) != n ||
        length(self@err_y) != n || !identical(dim(self@err_x), c(q, q, n)) ||
        !identical(dim(self@err_xy), c(n, q))) {
      return("Outcome, predictor, covariance, and ID dimensions must agree.")
    }
    NULL
  }
)

# A fit and a numerical failure share one stable return type.
class_mmom_fit <- S7::new_class(
  "class_mmom_fit",
  properties = list(
    converged = S7::class_any,
    status_code = S7::class_any,
    message = S7::class_any,
    coefficients = S7::class_any,
    vcov = S7::class_any,
    preliminary_coef = S7::class_any,
    corrected_design_moment = S7::class_any,
    analysis_weights = S7::class_any,
    steps = S7::class_any,
    diagnostics = S7::class_any
  )
)

class_mmom_step <- S7::new_class(
  "class_mmom_step",
  properties = list(
    step = S7::class_any,
    converged = S7::class_any,
    status_code = S7::class_any,
    message = S7::class_any,
    diagnostics = S7::class_any
  )
)

class_mmom_step1 <- S7::new_class(
  "class_mmom_step1",
  parent = class_mmom_step,
  properties = list(
    err_x_sum = S7::class_any,
    err_cov_design = S7::class_any,
    lambda_prelim = S7::class_any,
    correction_prelim = S7::class_any,
    coef0 = S7::class_any,
    predictor_idx = S7::class_any
  )
)

class_mmom_step2 <- S7::new_class(
  "class_mmom_step2",
  parent = class_mmom_step,
  properties = list(
    err_resid_var = S7::class_any,
    resid_var = S7::class_any,
    case_var = S7::class_any,
    case_weight = S7::class_any
  )
)

class_mmom_step3 <- S7::new_class(
  "class_mmom_step3",
  parent = class_mmom_step,
  properties = list(
    weighted_design = S7::class_any,
    weighted_moment = S7::class_any,
    weighted_x_sum = S7::class_any,
    weighted_cov_design = S7::class_any,
    weighted_xy_sum = S7::class_any,
    weighted_err_moment = S7::class_any,
    lambda_final = S7::class_any,
    correction_final = S7::class_any,
    corrected_moment = S7::class_any,
    coef_moment = S7::class_any,
    coef = S7::class_any
  )
)

class_mmom_step_se <- S7::new_class(
  "class_mmom_step_se",
  parent = class_mmom_step,
  properties = list(
    standard_error = S7::class_any,
    coef_vcov = S7::class_any
  )
)


class_mat_solve <- S7::new_class(
  "class_mat_solve",
  properties = list(
    success = S7::class_any,
    message = S7::class_any,
    lhs = S7::class_any,
    rhs = S7::class_any,
    solution = S7::class_any,
    inverse = S7::class_any,
    solver = S7::class_any,
    diagnostics = S7::class_any
  )
)

class_mlm_elements <- S7::new_class(
  "class_mlm_elements",
  properties = list(
    backend = S7::class_character,
    model = S7::class_any,
    ranef = S7::class_any,
    fixed_design = S7::class_any,
    fixed_coefs = S7::class_numeric,
    response = S7::class_numeric,
    offset = S7::class_any,
    observation_weights = S7::class_any,
    random_design = S7::class_any,
    random_terms = S7::class_any,
    grouping = S7::class_any,
    group_levels = S7::class_any,
    group_rows = S7::class_any,
    group_ptr = S7::class_any,
    random_covariance = S7::class_any,
    residual_covariance = S7::class_any,
    diagnostics = S7::class_any
  )
)

class_repair_ranef <- S7::new_class(
  "class_repair_ranef",
  properties = list(
    success = S7::class_logical,
    predictions = S7::class_any,
    prediction_error_covariance = S7::class_any,
    method = S7::class_character,
    group = S7::class_any,
    levels = S7::class_any,
    terms = S7::class_any,
    include_fixed = S7::class_logical,
    diagnostics = S7::class_any
  )
)

# Classes from other packages
class_merMod <- methods::getClassDef("merMod", package = "lme4")
class_lme <- S7::new_S3_class("lme")
