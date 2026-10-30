refit_teaching <- function() {
  d <- repair_example$repeated; people <- repair_example$people
  model <- lme4::lmer(score ~ day + (day | person), d,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  list(data = d, people = people, model = model,
    fit = repair(re("day") ~ practice, model, people, "person", backend = "r"))
}

test_that("native lmer identity replay preserves coefficients, inference, and controls", {
  ex <- refit_teaching()
  snapshot <- function(fit) list(coefficients = coef(fit), vcov = vcov(fit),
    scores = repair_scores(fit), diagnostics = fit@analysis_data@diagnostics,
    weights = model.weights(model.frame(fit@mlm_fit[[1L]])))
  before <- snapshot(ex$fit)
  replay <- repair_refit(ex$fit, ex$data, ex$people)
  expect_true(S7::S7_inherits(replay, class_repair_mod))
  expect_equal(coef(replay), coef(ex$fit), tolerance = 1e-6)
  expect_equal(vcov(replay), vcov(ex$fit), tolerance = 1e-6)
  expect_equal(repair_scores(replay), repair_scores(ex$fit), tolerance = 1e-6)
  expect_equal(repair_score_covariance(replay), repair_score_covariance(ex$fit), tolerance = 1e-6)
  expect_equal(repair_reliability(replay), repair_reliability(ex$fit), tolerance = 1e-6)
  expect_identical(replay@mlm_fit$measurement@optinfo$control, ex$model@optinfo$control)
  expect_identical(replay@mlm_fit$measurement@optinfo$optimizer, "bobyqa")
  expect_identical(replay@struc_fit@diagnostics$backend, "r")
  expect_true(summary(replay)$refit$checks$verified)
  expect_identical(summary(replay)$refit$mode, "identity")
  expect_output(print(replay), "original-data replay verified")
  expect_identical(snapshot(ex$fit), before)
  # A replayed native fit remains usable as the reference for a later replay.
  again <- repair_refit(replay, ex$data, ex$people)
  expect_equal(coef(again), coef(replay), tolerance = 1e-6)
})

test_that("replacement fits agree with independent manual reconstruction", {
  ex <- refit_teaching(); d <- ex$data
  d$score <- d$score + .1 * d$day^2
  p <- ex$people; p$practice <- p$practice + .1
  manual_model <- lme4::lmer(score ~ day + (day | person), d,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  manual <- repair(re("day") ~ practice, manual_model, p, "person", backend = "r")
  r <- repair_refit(ex$fit, ex$data, ex$people, new_measurement_data = d, new_data = p)
  expect_equal(coef(r), coef(manual), tolerance = 1e-6)
  expect_equal(vcov(r), vcov(manual), tolerance = 1e-6)
  expect_equal(repair_scores(r), repair_scores(manual), tolerance = 1e-6)
  expect_identical(summary(r)$refit$checks$stage, c("reference", "replacement"))
  expect_equal(summary(r)$refit$checks$n_measurement, rep(nrow(d), 2))
  wrong <- ex$data; wrong$score[1] <- wrong$score[1] + 2
  expect_error(repair_refit(ex$fit, wrong, ex$people, new_measurement_data = d), "Reference replay does not reproduce")
  wrong <- ex$people; wrong$practice[1] <- 100
  expect_error(repair_refit(ex$fit, ex$data, wrong), "Reference replay does not reproduce")
  expect_error(repair_refit(ex$fit, ex$data, ex$people[60:1, ]), "requested/analyzed cohort")
})

test_that("nested deterministic draws duplicate entire people with independent identities", {
  for (kind in c("nested_intercept", "nested_slope")) {
    ex <- nested_example(kind); case <- ex$case
    fit <- repair(ex$formula, list(growth = ex$model), case$people, "person")
    draw <- c(case$people$person[-1], case$people$person[2])
    set.seed(703); rng <- .Random.seed
    result <- repair_refit(fit, list(growth = case$data), case$people, draw = draw)
    expect_identical(.Random.seed, rng)
    map <- summary(result)$refit$draw
    expect_identical(map$source_id, draw)
    expect_equal(anyDuplicated(map$refit_id), 0L)
    expect_equal(nobs(result), length(draw))
    model <- result@mlm_fit$growth; frame <- model.frame(model)
    expect_setequal(levels(frame$person), map$refit_id)
    expect_equal(nrow(frame), sum(vapply(draw, function(id) sum(case$data$person == id & !is.na(case$data$y)), integer(1))))
    # Separate occurrences have disjoint nuisance levels, while their data agree.
    flist <- lme4::getME(model, "flist")
    duplicates <- map$refit_id[map$source_id == case$people$person[2]]
    for (nuisance in setdiff(names(flist), "person")) {
      a <- as.character(flist[[nuisance]][flist$person == duplicates[1]])
      b <- as.character(flist[[nuisance]][flist$person == duplicates[2]])
      expect_length(intersect(a, b), 0L)
    }
    expect_equal(frame$y[frame$person == duplicates[1]], frame$y[frame$person == duplicates[2]])
    # Independent manual reconstruction, using different fresh ID labels.
    raw <- case$data[match(rownames(model.frame(ex$model)), rownames(case$data)), ]
    manual_ids <- sprintf("manual%07d", seq_along(draw))
    blocks <- lapply(seq_along(draw), function(i) {
      rows <- raw[as.character(raw$person) == draw[i], , drop = FALSE]
      rows$person <- manual_ids[i]; rows
    })
    repeated <- do.call(rbind, blocks)
    people <- case$people[match(draw, case$people$person), ]; people$person <- manual_ids
    m <- lme4::lmer(as.formula(case$measurement_formula), repeated, weights = prior,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
    manual <- repair(ex$formula, list(growth = m), people, "person")
    expect_equal(coef(result), coef(manual), tolerance = 1e-6)
    expect_equal(vcov(result), vcov(manual), tolerance = 1e-6)
    expect_identical(unique(repair_scores(result)$group), "person")
    expect_true(all(repair_diagnostics(result, "status")$converged))
  }
})

test_that("standalone nuisance ID columns are re-keyed without merging copies", {
  ex <- nested_example(); d <- ex$case$data
  d$visit <- interaction(d$person, d$day, drop = TRUE)
  m <- lme4::lmer(y ~ x + offset(off) + (x | person) + (1 | visit), d,
    weights = prior, na.action = na.exclude)
  f <- repair(ex$formula, list(growth = m), ex$case$people, "person")
  r <- repair_refit(f, d, ex$case$people, draw = rep(ex$case$people$person[1:24], 2))
  expect_equal(nobs(r), 48L)
  groups <- lme4::getME(r@mlm_fit$growth, "flist")
  owners <- split(as.character(groups$person), groups$visit, drop = TRUE)
  expect_true(all(vapply(owners, function(x) length(unique(x)) == 1L, logical(1))))
})

test_that("fixed factor contrasts and the structural schema remain stable", {
  ex <- refit_teaching(); p <- ex$people
  p$category <- factor(rep(c("A", "B", "C"), 20))
  contrasts(p$category) <- contr.sum(3)
  fit <- repair(re("day") ~ practice + category, ex$model, p, "person", regularization = "custom", alpha = 9)
  original_options <- options(contrasts = c("contr.helmert", "contr.poly"))
  on.exit(options(original_options), add = TRUE)
  r <- repair_refit(fit, ex$data, p, new_data = p)
  expect_equal(coef(r), coef(fit), tolerance = 1e-6)
  expect_equal(vcov(r), vcov(fit), tolerance = 1e-6)
  expect_identical(r@struc_fit@diagnostics$final_moment_control$alpha, 9)
  missing <- p; missing$category[missing$category == "C"] <- "A"
  expect_error(repair_refit(fit, ex$data, p, new_data = missing), "absent levels")
  new <- p; new$category <- as.character(new$category); new$category[1] <- "new"
  expect_error(repair_refit(fit, ex$data, p, new_data = new), "new or absent levels")
})

test_that("measurement contrasts, ML, weights and offsets survive changed global defaults", {
  ex <- nested_example(); d <- ex$case$data
  d$phase <- factor(rep(c("a", "b"), length.out = nrow(d)))
  m <- lme4::lmer(y ~ x + phase + offset(off) + (x | person) + (1 | person:day), d,
    weights = prior, na.action = na.exclude, REML = FALSE,
    contrasts = list(phase = contr.sum(2)),
    control = lme4::lmerControl(optimizer = "bobyqa"))
  fit <- repair(ex$formula, list(growth = m), ex$case$people, "person")
  old <- options(contrasts = c("contr.helmert", "contr.poly")); on.exit(options(old), add = TRUE)
  r <- repair_refit(fit, d, ex$case$people)
  expect_false(lme4::isREML(r@mlm_fit$growth))
  expect_equal(coef(r), coef(fit), tolerance = 1e-6)
  expect_equal(model.weights(model.frame(r@mlm_fit$growth)), model.weights(model.frame(m)))
  expect_equal(model.offset(model.frame(r@mlm_fit$growth)), model.offset(model.frame(m)))
})

test_that("nlme residual structures and fitting settings replay from observations", {
  reference <- readRDS(system.file("extdata", "nlme-workflow-reference.rds", package = "repair"))
  for (case in reference$cases) {
    correlated <- case$kind == "both_ar1_ident"
    # These explicit constructors are deliberately visible in the original call.
    m <- if (correlated) nlme::lme(y ~ z, random = ~ z | id, data = case$data,
      correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varIdent(form = ~ 1 | phase),
      control = nlme::lmeControl(msMaxIter = 300)) else
        nlme::lme(y ~ z, random = ~ z | id, data = case$data, control = nlme::lmeControl(msMaxIter = 300))
    f <- repair(as.formula(case$formula), list(growth = m), case$people, "id", backend = "r")
    r <- repair_refit(f, case$data, case$people)
    expect_equal(coef(r), coef(f), tolerance = 1e-6)
    expect_equal(vcov(r), vcov(f), tolerance = 1e-6)
    expect_equal(repair_scores(r), repair_scores(f), tolerance = 1e-6)
    expect_equal(r@mlm_fit$growth$call$control$msMaxIter, 300)
    expect_identical(r@mlm_fit$growth$method, m$method)
  }
})

test_that("callbacks are explicit and must preserve rows, method and schema", {
  ex <- refit_teaching()
  recipe <- function(data) lme4::lmer(score ~ day + (day | person), data,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  replay <- repair_refit(ex$fit, ex$data, ex$people, refitter = recipe)
  expect_equal(coef(replay), coef(ex$fit), tolerance = 1e-6)
  expect_identical(summary(replay)$refit$adapter, "callback")
  wrong <- function(data) lme4::lmer(score ~ day + (day | person), data, REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  expect_error(repair_refit(ex$fit, ex$data, ex$people, refitter = wrong), "specification")
  dropping <- function(data) recipe(data[-1, ])
  expect_error(repair_refit(ex$fit, ex$data, ex$people, refitter = dropping), "effective observations")
  expect_error(repair_refit(ex$fit, ex$data, ex$people, refitter = function(data) stop("failure")), "Reference measurement refit failed: failure")
  expect_error(repair_refit(ex$fit, ex$data, ex$people, refitter = function(data) list()), "must return")
})

test_that("nonliteral recipes do not read mutable global bindings", {
  ex <- refit_teaching(); ctrl <- lme4::lmerControl(optimizer = "bobyqa")
  m <- lme4::lmer(score ~ day + (day | person), ex$data, control = ctrl)
  f <- repair(re("day") ~ practice, m, ex$people, "person")
  ctrl <- lme4::lmerControl(optimizer = "Nelder_Mead")
  expect_error(repair_refit(f, ex$data, ex$people), "global bindings are not replayed")
  d <- ex$data; d$external <- NULL
  external <- d$day
  m2 <- lme4::lmer(score ~ external + (external | person), d)
  f2 <- repair(re("external") ~ practice, m2, ex$people, "person")
  expect_error(repair_refit(f2, d, ex$people), "columns absent")
})

test_that("replacement omission, row loss, and unsupported scope fail explicitly", {
  ex <- refit_teaching(); d <- ex$data; d$score[1] <- NA_real_
  expect_error(repair_refit(ex$fit, ex$data, ex$people, new_measurement_data = d), "Replacement measurement refit failed")
  absent <- ex$data; absent$day <- NULL
  expect_error(repair_refit(ex$fit, ex$data, ex$people, new_measurement_data = absent), "columns absent")
  p <- ex$people; p$practice[1] <- NA_real_
  expect_error(repair_refit(ex$fit, ex$data, ex$people, new_data = p), "incomplete")
  expect_error(repair_refit(ex$fit, ex$data, ex$people, id = "wrong"), "original structural ID")
  expect_error(repair_refit(ex$fit, ex$data, ex$people, draw = "absent"), "analyzed source IDs")
  expect_error(repair_refit(ex$fit, ex$data, ex$people, draw = ex$people$person, new_data = p), "either")
  for (tol in list(NA_real_, 0, -1, Inf, "x", c(.1, .2)))
    expect_error(repair_refit(ex$fit, ex$data, ex$people, tolerance = tol), "positive finite")
  source <- multimodel_example()
  f <- repair(source$formula, source$models, source$case$people, "person", cross_model_error = source$case$cross_model_error)
  expect_error(repair_refit(f, list(), data.frame()), "one measurement model")
})

test_that("original exclusions are audited and measurement-only people cannot silently disappear in draws", {
  ex <- refit_teaching(); p <- ex$people; p$practice[1] <- NA_real_
  expect_warning(fit <- repair(re("day") ~ practice, ex$model, p, "person", na_action = "omit"), "Omitting 1")
  replay <- repair_refit(fit, ex$data, p)
  expect_equal(nobs(replay), 59L)
  expect_equal(summary(replay)$refit$checks$n_requested, 60L)
  expect_true(any(grepl("Omitting 1", summary(replay)$refit$warnings$reference)))
  expect_error(repair_refit(fit, ex$data, p, draw = fit@analysis_data@ids), "cohorts to coincide")
})

test_that("warnings are visible in the audit and optimizer failure is fatal", {
  ex <- refit_teaching()
  recipe <- function(data) {
    warning("Informational fitting warning")
    lme4::lmer(score ~ day + (day | person), data,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
  }
  r <- repair_refit(ex$fit, ex$data, ex$people, refitter = recipe)
  expect_identical(summary(r)$refit$warnings$reference, "Informational fitting warning")
  broken <- function(data) {
    model <- suppressWarnings(recipe(data)); model@optinfo$conv$opt <- 1L; model
  }
  expect_error(repair_refit(ex$fit, ex$data, ex$people, refitter = broken), "optimizer did not converge")
})

test_that("nlme power and continuous-time correlation recipes are rebuilt with explicit data", {
  d <- nlme_example_data()
  models <- list(
    power = nlme::lme(y ~ z + phase, random = ~ z | id, data = d, na.action = na.exclude,
      correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varPower(form = ~ scale_var),
      control = nlme::lmeControl(msMaxIter = 300)),
    fitted_power = nlme::lme(y ~ z + phase, random = ~ z | id, data = d, na.action = na.exclude,
      weights = nlme::varPower(form = ~ fitted(.)), control = nlme::lmeControl(msMaxIter = 300)),
    car1 = nlme::lme(y ~ z + phase, random = ~ z | id, data = d, na.action = na.exclude,
      correlation = nlme::corCAR1(form = ~ z | id), control = nlme::lmeControl(msMaxIter = 300)),
    diag_ml = nlme::lme(y ~ z + phase, random = list(id = nlme::pdDiag(~ z)),
      data = d, method = "ML", na.action = na.exclude, control = nlme::lmeControl(msMaxIter = 300))
  )
  for (model in models) {
    people <- data.frame(person = rownames(nlme::ranef(model)), observed = seq_len(35))
    fit <- repair(re("z") ~ observed, model, people, "person")
    replay <- repair_refit(fit, d, people)
    expect_equal(coef(replay), coef(fit), tolerance = 1e-6)
    expect_equal(vcov(replay), vcov(fit), tolerance = 1e-6)
    expect_equal(repair_score_covariance(replay), repair_score_covariance(fit), tolerance = 1e-6)
    expect_identical(summary(replay)$refit$checks$controls, "recorded_recipe")
  }
})

test_that("nlme control changes require an explicit matching recipe", {
  ex <- repair_example
  model <- nlme::lme(score ~ day, random = ~ day | person, data = ex$repeated)
  f <- repair(re("day") ~ practice, model, ex$people, "person")
  changed <- function(data) nlme::lme(score ~ day, random = ~ day | person, data = data,
    control = nlme::lmeControl(msMaxIter = 1000))
  expect_error(repair_refit(f, ex$repeated, ex$people, refitter = changed), "optimizer control recipe")
  settings <- nlme::lmeControl(msMaxIter = 300)
  model <- nlme::lme(score ~ day, random = ~ day | person, data = ex$repeated, control = settings)
  f <- repair(re("day") ~ practice, model, ex$people, "person")
  expect_error(repair_refit(f, ex$repeated, ex$people), "global bindings are not replayed")
  explicit <- function(data) nlme::lme(score ~ day, random = ~ day | person, data = data,
    control = nlme::lmeControl(msMaxIter = 300))
  r <- repair_refit(f, ex$repeated, ex$people, refitter = explicit)
  expect_identical(summary(r)$refit$checks$controls, "callback_declared")
  expect_equal(coef(r), coef(f), tolerance = 1e-6)
})

test_that("a singular but usable measurement fit is retained", {
  set.seed(541)
  d <- expand.grid(day = 0:7, person = factor(1:30))
  people <- data.frame(person = levels(d$person), practice = seq_len(30))
  b1 <- .1 * (people$practice - 15) + rnorm(30, sd = .3)
  d$score <- 50 + (1 + b1[d$person]) * d$day + rep(c(.3, -.3, -.3, .3, .3, -.3, -.3, .3), times = 30)
  model <- suppressMessages(lme4::lmer(score ~ day + (day | person), d))
  expect_true(lme4::isSingular(model))
  fit <- repair(re("day") ~ practice, model, people, "person")
  replay <- repair_refit(fit, d, people)
  expect_equal(coef(replay), coef(fit), tolerance = 1e-6)
  expect_true(lme4::isSingular(replay@mlm_fit$measurement))
  expect_true(length(summary(replay)$refit$messages$reference) > 0L)
})

test_that("native replay succeeds in a fresh R process without original call bindings", {
  ex <- refit_teaching()
  # Remove formula environments that could otherwise retain helper/global data.
  model <- ex$fit@mlm_fit$measurement
  formula <- stats::formula(model); environment(formula) <- baseenv()
  model@call$formula <- formula; model@call$data <- as.name("missing_original_binding")
  terms <- attr(model@frame, "terms"); attr(terms, ".Environment") <- baseenv()
  attr(model@frame, "terms") <- terms
  ex$fit@mlm_fit <- list(measurement = model)
  structural <- ex$fit@analysis_data@formula; environment(structural) <- baseenv()
  ex$fit@analysis_data@formula <- structural
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(list(fit = ex$fit, repeated = ex$data, people = ex$people), path)
  package_path <- system.file(package = "repair")
  source <- !file.exists(file.path(package_path, "Meta", "package.rds"))
  actual <- callr::r(function(path, package_path, source) {
    if (source) pkgload::load_all(package_path, quiet = TRUE) else library(repair)
    payload <- readRDS(path)
    result <- repair::repair_refit(payload$fit, payload$repeated, payload$people)
    list(coefficients = stats::coef(result), covariance = stats::vcov(result),
      checks = summary(result)$refit$checks)
  }, args = list(path = path, package_path = package_path, source = source), libpath = .libPaths(), timeout = 120)
  expect_equal(actual$coefficients, coef(ex$fit), tolerance = 1e-6)
  expect_equal(actual$covariance, vcov(ex$fit), tolerance = 1e-6)
  expect_true(actual$checks$verified)
})

test_that("logical predictor coding is preserved alongside factor coding", {
  ex <- refit_teaching(); d <- ex$data; p <- ex$people
  d$batch <- seq_len(nrow(d)) %% 3 == 0
  p$flag <- rep(c(TRUE, FALSE), length.out = nrow(p))
  model <- lme4::lmer(score ~ day + batch + (day | person), d)
  fit <- repair(re("day") ~ practice + flag, model, p, "person")
  old <- options(contrasts = c("contr.helmert", "contr.poly")); on.exit(options(old), add = TRUE)
  replay <- repair_refit(fit, d, p, new_data = p)
  expect_identical(names(coef(replay)), names(coef(fit)))
  expect_equal(coef(replay), coef(fit), tolerance = 1e-6)
  expect_identical(names(lme4::fixef(replay@mlm_fit$measurement)), names(lme4::fixef(model)))
})

test_that("nlme person draws preserve paired observations and agree with a manual refit", {
  case <- readRDS(system.file("extdata", "nlme-workflow-reference.rds", package = "repair"))$cases$both_ar1_ident
  model <- nlme::lme(y ~ z, random = ~ z | id, data = case$data,
    correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varIdent(form = ~ 1 | phase),
    control = nlme::lmeControl(msMaxIter = 300))
  fit <- repair(as.formula(case$formula), list(growth = model), case$people, "id")
  ids <- c(case$people$id[-1], case$people$id[2])
  result <- repair_refit(fit, case$data, case$people, draw = ids)
  copies <- lapply(seq_along(ids), function(i) {
    block <- case$data[case$data$id == ids[i], ]; block$id <- sprintf("manual%07d", i); block
  })
  repeated <- do.call(rbind, copies)
  people <- case$people[match(ids, case$people$id), ]; people$id <- sprintf("manual%07d", seq_along(ids))
  manual_model <- nlme::lme(y ~ z, random = ~ z | id, data = repeated,
    correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varIdent(form = ~ 1 | phase),
    control = nlme::lmeControl(msMaxIter = 300))
  manual <- repair(as.formula(case$formula), list(growth = manual_model), people, "id")
  expect_equal(coef(result), coef(manual), tolerance = 1e-6)
  expect_equal(vcov(result), vcov(manual), tolerance = 1e-6)
  expect_identical(summary(result)$refit$draw$source_id, ids)
  expect_equal(nobs(result), length(ids))
})
