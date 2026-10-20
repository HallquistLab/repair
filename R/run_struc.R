# this function fits the structural model based on the parsed formula and provided structural data
# it assumes inputs have been validated
# it will use the mmom backend implemented in R/mmom.R and src/mmom_kernels.cpp
run_struc <- function(parsed, struc_data) {
  
}

# this function will take the first stage model and formula, and any additional data to construct the structural data
#   - these data include the likelihood-only predictions and the error covariances
#   - multiple prediction methods will eventually be supported (e.g., explicit unshrinking)
build_struc_data <- function(mlm_fit, parsed_struc_form, data) {

}