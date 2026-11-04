bayes_test_parts <- function() {
  set.seed(431)
  sizes <- rep(5:8, 5)
  data <- do.call(rbind, lapply(seq_along(sizes), function(i) {
    z <- seq(-1, 1, length.out = sizes[i]) + i / 40
    data.frame(id = paste0("unit", i), z = z, w = i / 20,
      y = .7 + .4 * z + i / 30 + rnorm(1, sd = .4) * z + rnorm(length(z), sd = .3))
  }))
  data <- data[sample(nrow(data)), ]
  levels <- rev(unique(data$id))
  x <- model.matrix(~ z + w, data)
  list(spec = list(group = "id"), x = x, beta = c("(Intercept)" = .7, z = .4, w = .2),
    y = data$y, z = model.matrix(~ z, data), grouping = data$id, levels = levels,
    terms = c("(Intercept)", "z"), g = matrix(c(.5, .12, .12, .3), 2L), sigma2 = .09)
}

bayes_test_snapshot <- function(parts = bayes_test_parts(), backend = "brms") {
  elements <- repair_bayes_elements(parts, NULL, backend)
  structure(list(fit = NULL, elements = elements,
    diagnostics = list(mlm_backend = backend, mlm_singular = NA,
      mlm_messages = "Experimental Bayesian plug-in REPAIR: conditional model-based inference.",
      experimental = TRUE, inference = "conditional_model_based", stage1_uncertainty = "not_propagated")),
    class = "repair_bayes_plugin")
}
