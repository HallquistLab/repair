bootstrap_example <- function() {
  d <- repair_example
  model <- lme4::lmer(score ~ day + (day | person), d$repeated,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  list(repeated = d$repeated, people = d$people,
    fit = repair(re("day") ~ practice, model, d$people, "person", backend = "r"))
}

test_that("bootstrap draws and full refits agree with an independent manual loop", {
  ex <- bootstrap_example()
  before <- list(coef = coef(ex$fit), vcov = vcov(ex$fit), scores = repair_scores(ex$fit))
  set.seed(572); caller <- .Random.seed; kind <- RNGkind()
  boot <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 4, seed = 41)
  expect_identical(.Random.seed, caller)
  expect_identical(RNGkind(), kind)
  expect_s3_class(boot, "repair_bootstrap")
  expect_identical(coef(boot), coef(ex$fit))
  expect_identical(boot$specification$backend, "r")
  expect_identical(boot$specification$structural_formula, 're("day") ~ practice')
  expect_identical(boot$specification$effects, ex$fit@analysis_data@effects)
  expect_identical(boot$specification$regularization,
    ex$fit@struc_fit@diagnostics[c("preliminary_moment_control", "final_moment_control")])
  expect_identical(nobs(boot), nrow(ex$people))
  expect_identical(before, list(coef = coef(ex$fit), vcov = vcov(ex$fit), scores = repair_scores(ex$fit)))
  expect_true(all(repair_bootstrap_diagnostics(boot)$usable))
  expect_true(all(repair_bootstrap_diagnostics(boot)$n_analyzed == nrow(ex$people)))
  on.exit({ do.call(RNGkind, as.list(kind)); assign(".Random.seed", caller, .GlobalEnv) }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(41)
  stream <- .Random.seed
  estimates <- matrix(NA_real_, 4, length(coef(ex$fit)), dimnames = dimnames(boot$replicates))
  for (i in 1:4) {
    stream <- parallel::nextRNGStream(stream); assign(".Random.seed", stream, .GlobalEnv)
    index <- sample.int(nrow(ex$people), nrow(ex$people), replace = TRUE)
    expect_identical(boot$draws[i, ], index)
    expect_identical(boot$rng$replicates[[i]], stream)
    ids <- sprintf("manual%07d", seq_along(index))
    blocks <- lapply(seq_along(index), function(j) {
      block <- ex$repeated[ex$repeated$person == ex$people$person[index[j]], ]
      block$person <- ids[j]; block
    })
    repeated <- do.call(rbind, blocks)
    people <- ex$people[index, ]; people$person <- ids
    model <- lme4::lmer(score ~ day + (day | person), repeated,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    manual <- repair(re("day") ~ practice, model, people, "person", backend = "r")
    estimates[i, ] <- coef(manual)
    expect_equal(boot$diagnostics$n_measurement[i], nrow(repeated))
  }
  expect_equal(boot$replicates, estimates, tolerance = 1e-6)
  expect_equal(vcov(boot)[, ], cov(estimates), tolerance = 1e-6)
  expect_identical(attr(vcov(boot), "usable"), 4L)
  intervals <- t(apply(estimates, 2, quantile, probs = c(.05, .95), names = FALSE, type = 7))
  actual <- confint(boot, level = .9)
  expect_equal(unname(actual[, ]), unname(intervals), tolerance = 1e-6)
  expect_identical(attr(actual, "method"), "percentile")
  expect_false(attr(actual, "conditional_on_success"))
  table <- repair_bootstrap_coefficients(boot, conf.level = .9)
  expect_identical(names(table), c("term", "estimate", "std.error", "conf.low", "conf.high",
    "conf.level", "method", "attempted", "usable", "failed", "conditional_on_success", "inference"))
  expect_equal(table$std.error, sqrt(diag(cov(estimates))), tolerance = 1e-6, ignore_attr = TRUE)
  expect_identical(table$estimate, unname(coef(ex$fit)))
  expect_false(any(c("p.value", "statistic") %in% names(table)))
  expect_output(print(boot), "experimental")
  expect_output(print(summary(boot)), "attempts: 4")
  expect_equal(confint(boot, parm = 2, level = .9)[1, ], actual[2, ])
  map <- repair_bootstrap_diagnostics(boot, "draws")
  expect_identical(map$source_id, boot$cohort[as.vector(t(boot$draws))])
  expect_identical(map$replicate, rep(1:4, each = nrow(ex$people)))
  expect_equal(anyDuplicated(map[c("replicate", "refit_id")]), 0L)
})

test_that("attempt streams are reproducible, extendible, and insulated from fitter RNG use", {
  ex <- bootstrap_example()
  short <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 2, seed = 704)
  long <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 3, seed = 704)
  expect_identical(long$draws[1:2, ], short$draws)
  expect_equal(long$replicates[1:2, ], short$replicates, tolerance = 1e-10)
  recipe <- function(data) {
    runif(19)
    lme4::lmer(score ~ day + (day | person), data,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  }
  extra_rng <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 3, seed = 704, refitter = recipe)
  expect_identical(extra_rng$draws, long$draws)
  expect_equal(extra_rng$replicates, long$replicates, tolerance = 1e-10)
  expect_identical(extra_rng$rng, long$rng)
})

test_that("source IDs retain leading zeros with different structural and measurement column names", {
  d <- repair_example$repeated; people <- repair_example$people
  ids <- sprintf("%03d", seq_len(nrow(people)))
  d$person <- ids[match(as.character(d$person), people$person)]
  people$participant <- ids; people$person <- NULL
  model <- lme4::lmer(score ~ day + (day | person), d)
  fit <- repair(re("day") ~ practice, model, people, "participant")
  boot <- repair_bootstrap(fit, d, people, B = 2, seed = 41)
  map <- repair_bootstrap_diagnostics(boot, "draws")
  expect_true(all(grepl("^[0-9]{3}$", map$source_id)))
  expect_identical(boot$cohort, ids)
  expect_true(all(boot$diagnostics$usable))
})

test_that("factor loss records failed attempts without replacements or hidden deletion", {
  ex <- bootstrap_example(); people <- ex$people
  people$rare <- factor(seq_len(nrow(people)) == 1L)
  fit <- repair(re("day") ~ practice + rare, ex$fit@mlm_fit, people, "person")
  expect_warning(boot <- repair_bootstrap(fit, ex$repeated, people, B = 10, seed = 41), "attempts failed")
  expected_failures <- apply(boot$draws, 1, function(index) !1L %in% index)
  expect_true(any(expected_failures))
  expect_true(sum(!expected_failures) >= 2)
  expect_identical(!boot$diagnostics$usable, expected_failures)
  expect_equal(nrow(boot$diagnostics), 10L)
  expect_equal(nrow(boot$replicates), 10L)
  expect_true(all(is.na(boot$replicates[expected_failures, ])))
  expect_true(all(boot$diagnostics$stage[expected_failures] == "replacement_data"))
  expect_true(all(grepl("new or absent levels", boot$diagnostics$reason[expected_failures])))
  expect_error(vcov(boot), "withheld")
  expect_error(confint(boot), "withheld")
  expect_error(repair_bootstrap_coefficients(boot), "withheld")
  expect_false(summary(boot)$inference_available)
  expect_null(summary(boot)$coefficients)
  expect_identical(summary(boot)$failed, sum(expected_failures))
  vectors <- boot$replicates[!expected_failures, , drop = FALSE]
  expect_equal(vcov(boot, allow_failed = TRUE)[, ], cov(vectors))
  expect_identical(attr(vcov(boot, allow_failed = TRUE), "failed"), sum(expected_failures))
  table <- repair_bootstrap_coefficients(boot, allow_failed = TRUE)
  expect_true(all(table$conditional_on_success))
  expect_true(all(grepl("conditional on successful", table$inference)))
  expect_true(all(table$usable == nrow(vectors)))
  expect_true(summary(boot, allow_failed = TRUE)$inference_available)
  expect_true(summary(boot, allow_failed = TRUE)$conditional_on_success)
  expect_true(attr(confint(boot, allow_failed = TRUE), "conditional_on_success"))
})

test_that("measurement errors preserve warning messages and stage for every attempted draw", {
  ex <- bootstrap_example()
  recipe <- function(data) {
    if (all(startsWith(as.character(data$person), "r"))) {
      warning("specific fitting diagnostic")
      message("specific fitting message")
      stop("deliberate replacement failure")
    }
    lme4::lmer(score ~ day + (day | person), data,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  }
  set.seed(810); before <- .Random.seed
  expect_warning(boot <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 2,
    refitter = recipe), "2 of 2")
  expect_identical(.Random.seed, before)
  expect_false(any(boot$diagnostics$usable))
  expect_true(all(is.na(boot$replicates)))
  expect_identical(boot$diagnostics$stage, rep("replacement_measurement", 2))
  expect_true(all(grepl("deliberate replacement failure", boot$diagnostics$reason)))
  expect_identical(boot$diagnostics$warnings[[1]], "specific fitting diagnostic")
  expect_identical(boot$diagnostics$messages[[1]], "specific fitting message\n")
  expect_error(vcov(boot, allow_failed = TRUE), "at least two")
  expect_false(summary(boot, allow_failed = TRUE)$inference_available)
})

test_that("structural failures retain the draw and successful measurement warnings", {
  ex <- bootstrap_example(); people <- ex$people
  people$spike <- as.numeric(seq_len(nrow(people)) == 1L)
  fit <- repair(re("day") ~ practice + spike, ex$fit@mlm_fit, people, "person")
  recipe <- function(data) {
    warning("measurement information")
    lme4::lmer(score ~ day + (day | person), data,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  }
  expect_warning(boot <- repair_bootstrap(fit, ex$repeated, people, B = 10,
    seed = 41, refitter = recipe), "attempts failed")
  missing_spike <- apply(boot$draws, 1, function(index) !1L %in% index)
  expect_identical(!boot$diagnostics$usable, missing_spike)
  expect_true(all(boot$diagnostics$stage[missing_spike] == "replacement_structural"))
  expect_true(all(vapply(boot$diagnostics$warnings, function(x) "measurement information" %in% x, logical(1))))
  expect_true(all(is.na(boot$replicates[missing_spike, ])))
})

test_that("singular but usable sampled models remain eligible", {
  set.seed(541)
  d <- expand.grid(day = 0:7, person = factor(1:30))
  people <- data.frame(person = levels(d$person), practice = seq_len(30))
  b1 <- .1 * (people$practice - 15) + rnorm(30, sd = .3)
  d$score <- 50 + (1 + b1[d$person]) * d$day + rep(c(.3, -.3, -.3, .3, .3, -.3, -.3, .3), times = 30)
  model <- suppressMessages(lme4::lmer(score ~ day + (day | person), d))
  fit <- repair(re("day") ~ practice, model, people, "person")
  boot <- repair_bootstrap(fit, d, people, B = 2, seed = 41)
  expect_true(all(boot$diagnostics$usable))
  expect_true(all(boot$diagnostics$measurement_singular))
  expect_true(all(vapply(boot$diagnostics$messages, length, integer(1)) > 0L))
  expect_true(all(is.finite(boot$replicates)))
})

test_that("nested bootstrap duplicates keep complete histories and independent nuisance levels", {
  ex <- nested_example("nested_slope"); case <- ex$case
  fit <- repair(ex$formula, list(growth = ex$model), case$people, "person")
  boot <- repair_bootstrap(fit, case$data, case$people, B = 2, seed = 304)
  map <- repair_bootstrap_diagnostics(boot, "draws")
  expect_true(all(boot$diagnostics$usable))
  expect_true(anyDuplicated(map$source_id[map$replicate == 1]) > 0L)
  for (i in 1:2) {
    ids <- map$source_id[map$replicate == i]
    raw <- case$data[match(rownames(model.frame(ex$model)), rownames(case$data)), ]
    copies <- lapply(seq_along(ids), function(j) {
      block <- raw[raw$person == ids[j], ]; block$person <- sprintf("manual%07d", j); block
    })
    repeated <- do.call(rbind, copies)
    people <- case$people[match(ids, case$people$person), ]; people$person <- sprintf("manual%07d", seq_along(ids))
    model <- lme4::lmer(y ~ x + offset(off) + (x | person) + (within | person:day) + (1 | person:block),
      repeated, weights = prior, control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    manual <- repair(ex$formula, list(growth = model), people, "person")
    expect_equal(boot$replicates[i, ], coef(manual), tolerance = 1e-5)
    expect_equal(boot$diagnostics$n_measurement[i], nrow(repeated))
  }
})

test_that("nlme bootstrap rebuilds correlation and variance structures for sampled people", {
  case <- readRDS(system.file("extdata", "nlme-workflow-reference.rds", package = "repair"))$cases$both_ar1_ident
  model <- nlme::lme(y ~ z, random = ~ z | id, data = case$data,
    correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varIdent(form = ~ 1 | phase),
    control = nlme::lmeControl(msMaxIter = 300))
  fit <- repair(as.formula(case$formula), list(growth = model), case$people, "id")
  boot <- repair_bootstrap(fit, case$data, case$people, B = 2, seed = 143)
  expect_true(all(boot$diagnostics$usable))
  expect_true(all(is.na(boot$diagnostics$measurement_singular)))
  ids <- boot$cohort[boot$draws[1, ]]
  copies <- lapply(seq_along(ids), function(i) {
    block <- case$data[case$data$id == ids[i], ]; block$id <- sprintf("manual%07d", i); block
  })
  repeated <- do.call(rbind, copies)
  people <- case$people[match(ids, case$people$id), ]; people$id <- sprintf("manual%07d", seq_along(ids))
  manual_model <- nlme::lme(y ~ z, random = ~ z | id, data = repeated,
    correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varIdent(form = ~ 1 | phase),
    control = nlme::lmeControl(msMaxIter = 300))
  manual <- repair(as.formula(case$formula), list(growth = manual_model), people, "id")
  expect_equal(boot$replicates[1, ], coef(manual), tolerance = 1e-6)
})

test_that("preflight errors stop before sampling and preserve caller RNG", {
  ex <- bootstrap_example(); wrong <- ex$repeated; wrong$score[1] <- 100
  set.seed(730); before <- .Random.seed; kind <- RNGkind()
  expect_error(repair_bootstrap(ex$fit, wrong, ex$people, B = 2), "Reference replay does not reproduce")
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  people <- ex$people; people$practice[1] <- NA
  expect_warning(fit <- repair(re("day") ~ practice, ex$fit@mlm_fit, people, "person", na_action = "omit"), "Omitting")
  expect_error(repair_bootstrap(fit, ex$repeated, people, B = 2), "cohorts to coincide")
  expect_identical(.Random.seed, before)
  multi <- multimodel_example("independent")
  fit <- repair(multi$formula, multi$models, multi$case$people, "person", cross_model_error = "independent")
  expect_error(repair_bootstrap(fit, multi$case$data, multi$case$people, B = 2), "one measurement model")
  expect_identical(.Random.seed, before)
})

test_that("arguments, interval selection, and numerical eligibility are validated", {
  ex <- bootstrap_example()
  for (B in list(0, -1, 1.5, NA, Inf, c(1, 2), "2"))
    expect_error(repair_bootstrap(ex$fit, ex$repeated, ex$people, B = B), "`B`")
  for (seed in list(-1, 1.5, NA, Inf, c(1, 2), "2"))
    expect_error(repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 1, seed = seed), "`seed`")
  boot <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 1)
  expect_false(summary(boot)$inference_available)
  expect_error(confint(boot), "at least two")
  expect_error(vcov(boot), "at least two")
  boot <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 2)
  for (level in list(0, 1, NA, Inf, c(.9, .95))) expect_error(confint(boot, level = level), "`level`")
  for (parm in list(0, -1, 1.5, 3, NA, Inf, "bad")) expect_error(confint(boot, parm = parm), "coefficient")
  expect_error(confint(boot, type = "bca"), "percentile")
  expect_error(vcov(boot, allow_failed = NA), "allow_failed")
  expect_error(summary(boot, allow_failed = 1), "allow_failed")
  expect_error(repair_bootstrap_diagnostics(boot, "draw"), "table")
  expect_error(repair_bootstrap_coefficients(ex$fit), "bootstrap")
  kind <- RNGkind(); seed <- .Random.seed
  on.exit({ do.call(RNGkind, as.list(kind)); assign(".Random.seed", seed, .GlobalEnv) }, add = TRUE)
  RNGkind(normal.kind = "Box-Muller"); rnorm(1)
  before <- .Random.seed
  expect_error(repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 1), "hidden Box-Muller")
  expect_identical(.Random.seed, before)
})

test_that("saved bootstrap results dispatch and retain absent RNG state in a fresh R process", {
  ex <- bootstrap_example()
  boot <- repair_bootstrap(ex$fit, ex$repeated, ex$people, B = 2, seed = 45)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(boot, path)
  package_path <- system.file(package = "repair")
  source <- !file.exists(file.path(package_path, "Meta", "package.rds"))
  outputs <- callr::r(function(path, package_path, source) {
    if (source) pkgload::load_all(package_path, quiet = TRUE) else library(repair)
    b <- readRDS(path)
    # No reference model/data/global bindings are needed for saved extraction.
    a <- list(coef = stats::coef(b), vcov = stats::vcov(b), intervals = stats::confint(b),
      table = repair::repair_bootstrap_coefficients(b), summary = summary(b), specification = b$specification,
      draws = repair::repair_bootstrap_diagnostics(b, "draws"),
      printed = capture.output(print(b)))
    d <- repair::repair_example
    model <- lme4::lmer(score ~ day + (day | person), d$repeated)
    f <- repair::repair(as.formula('re("day") ~ practice'), model, d$people, "person")
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    before <- RNGkind()
    repair::repair_bootstrap(f, d$repeated, d$people, B = 1)
    a$absent_seed <- !exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    a$kind_preserved <- identical(RNGkind(), before)
    a
  }, args = list(path = path, package_path = package_path, source = source), libpath = .libPaths(), timeout = 120)
  expect_identical(outputs$coef, coef(boot))
  expect_identical(outputs$vcov, vcov(boot))
  expect_identical(outputs$intervals, confint(boot))
  expect_identical(outputs$table, repair_bootstrap_coefficients(boot))
  expect_identical(outputs$summary, summary(boot))
  expect_identical(outputs$specification, boot$specification)
  expect_identical(outputs$draws, repair_bootstrap_diagnostics(boot, "draws"))
  expect_identical(outputs$printed, capture.output(print(boot)))
  expect_true(outputs$absent_seed)
  expect_true(outputs$kind_preserved)
})
