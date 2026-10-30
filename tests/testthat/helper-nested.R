nested_example <- function(kind = "nested_intercept") {
  case <- readRDS(system.file("extdata", "nested-workflow-reference.rds", package = "repair"))$cases[[kind]]
  model <- lme4::lmer(as.formula(case$measurement_formula), case$data,
    weights = prior, na.action = na.exclude,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  list(case = case, model = model, formula = as.formula(case$formula))
}
