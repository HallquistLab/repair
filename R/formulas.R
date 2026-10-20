# this file will contain functions for parsing and validating MLM and structural formulas

# this function takes an MLM formula and parses:
#   - all observed variable names, including grouping variables
#   - the full fixed and random components
#   - the random intercepts and slopes identified by the group
# it also prevents the use of re()
parse_mlm_formula <- function(formula) {

}

# this function takes a structural formula and parses:
#   - all observed variable names
#   - uses of re() that define random slopes and/or intercepts
#     - it defines internal labels for these
parse_struc_formula <- function(formula) {

}

# this function ensures that:
#   - (a) the grouping variables are consistent between the MLM and structural formulas,
#   - (b) all required variables are present in the data and
#   - (c) that any uses of re() have corresponding defined REs in the MLM
# the function handles error messaging for these conditions
validate_formulas <- function(parsed_mlm_form, parsed_struc_form, data) {

}

# below are util functions for parsing formulas and working with the new re() implementation

# this is just a first pass at a group var extraction function
get_gr_terms <- function(formula) {
  bars <- reformulas::findbars(formula)

  lapply(bars, function(bar) {
    list(
      expression = bar[[3L]],
      variables = all.vars(bar[[3L]]),
      label = paste(deparse(bar[[3L]]), collapse = "")
    )
  })
}

# also a first pass at a logical for re()
has_re <- function(formula) {
  any(vapply(all.names(formula), identical, logical(1), "re"))
}
