workflow_example <- function() {
  model <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  people <- data.frame(person = levels(lme4::sleepstudy$Subject),
    x = seq_len(18), category = factor(rep(c("A", "B", "C"), 6)))
  list(models = list(growth = model), data = people)
}
