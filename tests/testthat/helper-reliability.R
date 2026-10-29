reliability_example <- function() {
  ex <- workflow_example()
  repair(re("Days") ~ re("(Intercept)") + x, ex$models, ex$data, "person", backend = "r")
}
