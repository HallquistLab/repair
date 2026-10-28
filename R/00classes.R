# S7 classes returned by the package.

#' Combined first-stage and structural repair model.
#' @param mlm_fit The first-stage multilevel-model fit.
#' @param struc_fit The structural M-MOM fit.
#' @export
class_repair_mod <- S7::new_class(
  "class_repair_mod",
  properties = list(
    mlm_fit = S7::class_any,
    struc_fit = S7::class_any
  )
)

repair_model <- function(mlm_fit, struc_fit) {
  class_repair_mod(mlm_fit = mlm_fit, struc_fit = struc_fit)
}

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
    error_covariance = S7::class_any,
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
