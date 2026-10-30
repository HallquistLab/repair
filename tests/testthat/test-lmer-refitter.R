test_that("robust recipe replays original fits and whole-person draws without RNG use", {
  ex <- repair_example
  fitter <- repair_lmer_refitter(score ~ day + (day | person))
  set.seed(730); before <- .Random.seed
  model <- fitter(ex$repeated)
  expect_identical(.Random.seed, before)
  audit <- attr(model, "repair_lmer_optimization")
  expect_true(audit$usable)
  expect_identical(model@optinfo$optimizer, "bobyqa")
  expect_true(nrow(audit$attempts) >= 2L)
  expect_lte(audit$attempts$objective[audit$selected], min(audit$attempts$objective, na.rm = TRUE) + 1e-6)
  fit <- repair(re("day") ~ practice, model, ex$people, "person")
  replay <- repair_refit(fit, ex$repeated, ex$people)
  expect_equal(coef(replay), coef(fit), tolerance = 1e-8)
  expect_equal(vcov(replay), vcov(fit), tolerance = 1e-8)
  expect_true(summary(replay)$refit$optimization$reference$usable)
  draw <- rep(ex$people$person[1:30], 2)
  changed <- repair_refit(fit, ex$repeated, ex$people, draw = draw)
  expect_identical(.Random.seed, before)
  expect_equal(nobs(changed), length(draw))
  expect_true(summary(changed)$refit$optimization$replacement$usable)
  boot <- repair_bootstrap(fit, ex$repeated, ex$people, B = 3, seed = 730)
  expect_length(boot$diagnostics$optimization, 3L)
  expect_true(all(vapply(boot$diagnostics$optimization, function(x) x$usable, logical(1))))
  expect_identical(.Random.seed, before)
  # Restoring a saved model needs neither the original closure nor global data.
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path); restored <- readRDS(path)
  expect_equal(coef(repair_refit(restored, ex$repeated, ex$people)), coef(fit), tolerance = 1e-8)
  wrong <- repair_lmer_refitter(score ~ day + (day | person), maxfun = 300000)
  expect_error(repair_refit(fit, ex$repeated, ex$people, refitter = wrong), "recipe|specification")
})

test_that("robust fitting preserves ML, weights, offsets, contrasts and nested groups", {
  ex <- nested_example(); d <- ex$case$data[match(rownames(model.frame(ex$model)), rownames(ex$case$data)), ]
  d$phase <- factor(rep(c("a", "b"), length.out = nrow(d)))
  fitter <- repair_lmer_refitter(y ~ x + phase + offset(off) + (x | person) + (1 | person:day),
    REML = FALSE, weights = "prior", contrasts = list(phase = contr.sum(2)))
  model <- fitter(d)
  fit <- repair(ex$formula, list(growth = model), ex$case$people, "person")
  old <- options(contrasts = c("contr.helmert", "contr.poly")); on.exit(options(old), add = TRUE)
  replay <- repair_refit(fit, d, ex$case$people)
  expect_false(lme4::isREML(replay@mlm_fit$growth))
  expect_equal(model.weights(model.frame(replay@mlm_fit$growth)), d$prior)
  expect_equal(model.offset(model.frame(replay@mlm_fit$growth)), d$off)
  expect_equal(coef(replay), coef(fit), tolerance = 1e-7)
  expect_equal(vcov(replay), vcov(fit), tolerance = 1e-7)
})

test_that("unresolved convergence and inferior objectives remain failures with audit trails", {
  ex <- repair_example
  real_attempt <- repair_lmer_attempt
  # Exercise genuine safety branches without depending on platform-specific warnings.
  local_mocked_bindings(repair_lmer_attempt = function(args, label) {
    x <- real_attempt(args, label)
    if (label != "nlopt_start") {
      x$record$usable <- FALSE
      x$record$reason <- "Injected optimizer nonconvergence"
      x$model@optinfo$conv$opt <- 1L
    }
    x
  })
  fitter <- repair_lmer_refitter(score ~ day + (day | person))
  expect_warning(model <- fitter(ex$repeated), "remains unresolved")
  audit <- attr(model, "repair_lmer_optimization")
  expect_false(audit$usable); expect_equal(nrow(audit$attempts), 4L)
  fit <- repair(re("day") ~ practice, model, ex$people, "person")
  expect_error(repair_refit(fit, ex$repeated, ex$people), "optimizer did not converge")
})

test_that("likelihood gate rejects converged but inferior retry fits", {
  ex <- repair_example; real_attempt <- repair_lmer_attempt
  local_mocked_bindings(repair_lmer_attempt = function(args, label) {
    x <- real_attempt(args, label)
    if (label != "nlopt_start") x$record$objective <- x$record$objective + 100
    x
  })
  fitter <- repair_lmer_refitter(score ~ day + (day | person))
  expect_warning(model <- fitter(ex$repeated), "remains unresolved")
  expect_false(attr(model, "repair_lmer_optimization")$usable)
  fit <- repair(re("day") ~ practice, model, ex$people, "person")
  expect_error(repair_refit(fit, ex$repeated, ex$people), "optimizer did not converge")
})

test_that("recipe inputs and implicit omissions fail explicitly", {
  expect_error(repair_lmer_refitter("y ~ x"), "formula")
  expect_error(repair_lmer_refitter(y ~ x, REML = NA), "REML")
  expect_error(repair_lmer_refitter(y ~ x, weights = 1:3), "column")
  expect_error(repair_lmer_refitter(y ~ x, maxfun = 0), "positive integer")
  expect_error(repair_lmer_refitter(y ~ x, objective_tolerance = Inf), "finite")
  fitter <- repair_lmer_refitter(score ~ day + (day | person))
  d <- repair_example$repeated; d$score[1] <- NA_real_
  expect_error(fitter(d), "no finite BOBYQA")
})

test_that("a failed warm attempt retries deterministically and retains the discarded warning", {
  ex <- repair_example; real_attempt <- repair_lmer_attempt
  local_mocked_bindings(repair_lmer_attempt = function(args, label) {
    x <- real_attempt(args, label)
    if (label == "bobyqa_warm") {
      x$record$usable <- FALSE; x$record$reason <- "Warm-start convergence failure"
      x$record$warnings[[1]] <- "Warm-start convergence failure"
    }
    x
  })
  set.seed(402); before <- .Random.seed
  m <- repair_lmer_refitter(score ~ day + (day | person))(ex$repeated)
  a <- attr(m, "repair_lmer_optimization")
  expect_true(a$usable)
  expect_equal(a$attempts$attempt[a$selected], "bobyqa_cold")
  expect_equal(a$attempts$warnings[[2]], "Warm-start convergence failure")
  expect_identical(.Random.seed, before)
})

test_that("a recorded simulation failure is recovered without changing the likelihood target", {
  case <- readRDS(system.file("extdata", "lmer-convergence-reference.rds", package = "repair"))
  baseline <- suppressWarnings(suppressMessages(lme4::lmer(response ~ time + (time | person), case$data)))
  fitter <- repair_lmer_refitter(response ~ time + (time | person))
  recovered <- fitter(case$data)
  expect_true(attr(recovered, "repair_lmer_optimization")$usable)
  expect_identical(rownames(model.frame(recovered)), rownames(case$data))
  expect_true(lme4::isREML(recovered))
  expect_lte(-2 * as.numeric(logLik(recovered)), -2 * as.numeric(logLik(baseline)) + 1e-6)
  # Independent optimizer fit from the reported starting theta.
  manual <- suppressMessages(lme4::lmer(response ~ time + (time | person), case$data,
    start = lme4::getME(baseline, "theta"),
    control = lme4::lmerControl(optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000, rhobeg = .002, rhoend = 1e-8))))
  expect_equal(lme4::fixef(recovered), lme4::fixef(manual), tolerance = 1e-5)
  expect_equal(sigma(recovered), sigma(manual), tolerance = 1e-5)
})
