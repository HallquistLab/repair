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
    diagnostics = S7::class_any
  )
)
