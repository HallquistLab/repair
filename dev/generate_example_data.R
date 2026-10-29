# Regenerate deliberately from the package root with Rscript dev/generate_example_data.R.
# Teaching data only: this is not a numerical acceptance fixture or real study.
local({
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20261029)
  n <- 60L
  person <- sprintf("P%03d", seq_len(n))
  practice <- round(runif(n, 0, 10), 1)
  z0 <- rnorm(n)
  z1 <- rnorm(n)
  intercept <- 3 * z0
  slope <- 0.18 * (practice - 5) + 0.15 * z0 + 0.35 * z1
  visits <- sample(8:12, n, replace = TRUE)
  repeated <- do.call(rbind, lapply(seq_len(n), function(i) {
    day <- seq_len(visits[i]) - 1L
    data.frame(person = person[i], day = as.numeric(day),
      score = 50 + intercept[i] + (1 + slope[i]) * day + rnorm(visits[i], sd = 1.5))
  }))
  rownames(repeated) <- NULL
  repair_example <- list(
    repeated = repeated,
    people = data.frame(person = person, practice = practice),
    truth = list(fixed_effects = c(intercept = 50, day = 1),
      slope_on_practice = 0.18, practice_center = 5, residual_sd = 1.5,
      random_effects = data.frame(person = person,
        intercept_deviation = intercept, slope_deviation = slope))
  )
  save(repair_example, file = "data/repair_example.rda", compress = "xz", version = 2)
})
