multimodel_example <- function(kind = "supplied") {
  case <- readRDS(system.file("extdata", "multimodel-workflow-reference.rds", package = "repair"))$cases[[kind]]
  models <- list(growth = nlme::lme(y ~ z, random = ~ z | id, data = case$a),
    response = nlme::lme(q ~ w, random = ~ w | subject, data = case$b))
  list(case = case, models = models, formula = as.formula(case$formula))
}
