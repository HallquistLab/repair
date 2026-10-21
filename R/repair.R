# this is the workhorse function of the package will accepts formulas for level 1 and level 2
repair <- function(f_mlm, f_struc, data) {

  parsed_mlm_form <- parse_mlm_formula(f_mlm)
  parsed_struc_form <- parse_struc_formula(f_struc)
  validate_formulas(parsed_mlm_form, parsed_struc_form, data)

  mlm_fit <- run_mlm(parsed_mlm_form, data)

  struc_data <- build_struc_data(mlm_fit, parsed_struc_form, data)
  struc_fit <- run_struc(parsed_struc_form, struc_data)

  class_repair_mod(mlm_fit, struc_fit)

}