expect_spec_fit_equal <- function(actual, expected) {
  expect_identical(actual@analysis_data@ids, expected@analysis_data@ids)
  for (field in c('outcome', 'predictors', 'err_x', 'err_y', 'err_xy'))
    expect_equal(S7::prop(actual@analysis_data, field), S7::prop(expected@analysis_data, field), tolerance=1e-12)
  for (field in c('coefficients', 'vcov', 'preliminary_coef', 'corrected_design_moment',
                  'analysis_weights', 'steps', 'diagnostics'))
    expect_equal(S7::prop(actual@struc_fit, field), S7::prop(expected@struc_fit, field), tolerance=1e-12)
}

test_that('MLM specifications defer fitting and are consumed only once per public call', {
  count <- 0L; fitter <- repair_fit_mlm_spec
  local_mocked_bindings(repair_fit_mlm_spec = function(...) {
    count <<- count + 1L; fitter(...)
  })
  specification <- repair_mlm(score ~ day + (day | person), repair_example$repeated)
  expect_s3_class(specification, 'repair_mlm_spec')
  expect_equal(count, 0L)
  expect_output(print(specification), 'first-stage MLM specification \\(lmer\\)')
  expect_equal(count, 0L)
  fit <- repair(re('day') ~ practice, specification, repair_example$people, 'person')
  expect_equal(count, 1L)
  expect_identical(names(fit@mlm_fit), 'mlm')
  repair_data(re('day') ~ practice, specification, repair_example$people, 'person')
  expect_equal(count, 2L)
  repair_effects(specification)
  expect_equal(count, 3L)
  expect_identical(normalize_repair_models(fit@mlm_fit), fit@mlm_fit)
  expect_equal(count, 3L)
})

test_that('inline lmer fitting preserves native controls, rows, weights, offsets and contrasts', {
  d <- repair_example$repeated
  d$phase <- factor(ifelse(d$day %% 2, 'a', 'b'))
  d$precision <- 1 + d$day / 20; d$off <- .03 * d$day
  d$score[2L] <- NA_real_
  control <- lme4::lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000))
  contrast <- list(phase=contr.sum(2))
  specification <- repair_mlm(score ~ day + phase + (day | person), d,
    control=control, REML=FALSE, weights=precision, offset=off,
    subset=day > 0, na.action=na.exclude, contrasts=contrast)
  manual <- lme4::lmer(score ~ day + phase + (day | person), d,
    control=control, REML=FALSE, weights=precision, offset=off,
    subset=day > 0, na.action=na.exclude, contrasts=contrast)
  # Specifications snapshot their data and evaluated arguments.
  original <- d; d$score <- 999; control$optimizer <- 'not_an_optimizer'
  for (backend in c('r','rcpp')) {
    a <- repair(re('day') ~ practice, specification, repair_example$people, 'person',
      backend=backend, regularization='step_up')
    b <- repair(re('day') ~ practice, manual, repair_example$people, 'person',
      backend=backend, regularization='step_up')
    expect_spec_fit_equal(a, b)
    expect_identical(rownames(model.frame(a@mlm_fit$mlm)), rownames(model.frame(manual)))
    expect_equal(model.weights(model.frame(a@mlm_fit$mlm)), model.weights(model.frame(manual)))
    expect_equal(model.offset(model.frame(a@mlm_fit$mlm)), model.offset(model.frame(manual)))
    replay <- repair_refit(a, original, repair_example$people)
    expect_equal(coef(replay), coef(a), tolerance=1e-6)
    expect_equal(vcov(replay), vcov(a), tolerance=1e-6)
    expect_identical(replay@mlm_fit$mlm@optinfo$optimizer, 'bobyqa')
  }
})

test_that('inline lme fitting preserves residual structures and supports native replay', {
  skip_if_not_installed('nlme')
  d <- nlme_example_data()
  control <- nlme::lmeControl(msMaxIter=300)
  specification <- repair_mlm(y ~ z + phase, d, engine='lme', random=~z|id,
    correlation=nlme::corAR1(form=~time|id), weights=nlme::varIdent(form=~1|phase),
    method='ML', subset=time > 0, na.action=na.exclude, control=control,
    contrasts=list(phase=contr.sum(2)))
  manual <- nlme::lme(y ~ z + phase, data=d, random=~z|id,
    correlation=nlme::corAR1(form=~time|id), weights=nlme::varIdent(form=~1|phase),
    method='ML', subset=time > 0, na.action=na.exclude, control=control,
    contrasts=list(phase=contr.sum(2)))
  people <- data.frame(id=rownames(nlme::ranef(manual)), x=seq_len(35))
  for (backend in c('r','rcpp')) {
    a <- repair(re('z') ~ x, specification, people, 'id', backend=backend)
    b <- repair(re('z') ~ x, manual, people, 'id', backend=backend)
    expect_spec_fit_equal(a, b)
    expect_equal(as.numeric(logLik(a@mlm_fit$mlm)), as.numeric(logLik(manual)), tolerance=1e-12)
    replay <- repair_refit(a, d, people)
    expect_equal(coef(replay), coef(a), tolerance=1e-6)
    expect_equal(vcov(replay), vcov(a), tolerance=1e-6)
    expect_identical(replay@mlm_fit$mlm$method, 'ML')
  }
})

test_that('named model lists may mix specifications and fitted models without implicit covariance assumptions', {
  m <- lme4::lmer(score ~ day + (day | person), repair_example$repeated)
  specification <- repair_mlm(score ~ day + (day | person), repair_example$repeated)
  models <- list(a=m, b=specification)
  resolved <- normalize_repair_models(models)
  expect_identical(resolved$a, m)
  expect_equal(repair_effects(models), repair_effects(resolved))
  form <- practice ~ re('day', 'a') + re('(Intercept)', 'b')
  expect_error(repair_data(form, models, repair_example$people, 'person'), 'Cross-model')
  a <- repair(form, models, repair_example$people, 'person', cross_model_error='independent')
  b <- repair(form, resolved, repair_example$people, 'person', cross_model_error='independent')
  expect_spec_fit_equal(a, b)
})

test_that('specification fits retain bootstrap draws and numerical results', {
  specification <- repair_mlm(score ~ day + (day | person), repair_example$repeated,
    control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000)))
  native <- lme4::lmer(score ~ day + (day | person), repair_example$repeated,
    control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000)))
  a <- repair(re('day') ~ practice, specification, repair_example$people, 'person', regularization='step_up')
  b <- repair(re('day') ~ practice, native, repair_example$people, 'person', regularization='step_up')
  first <- repair_bootstrap(a, repair_example$repeated, repair_example$people, B=3, seed=41)
  second <- repair_bootstrap(b, repair_example$repeated, repair_example$people, B=3, seed=41)
  expect_identical(first$draws, second$draws)
  expect_equal(first$replicates, second$replicates, tolerance=1e-8)
  columns <- c('usable','preliminary_alpha','final_alpha','preliminary_fallback','final_fallback')
  expect_equal(first$diagnostics[columns], second$diagnostics[columns], tolerance=1e-8)
  expect_equal(vcov(first), vcov(second), tolerance=1e-8)
  expect_equal(confint(first), confint(second), tolerance=1e-8)
})

test_that('specifications and their fitted results survive a fresh process', {
  skip_if_not_installed('callr')
  specification <- repair_mlm(score ~ day + (day | person), repair_example$repeated,
    control=lme4::lmerControl(optimizer='bobyqa'), REML=FALSE)
  fit <- repair(re('day') ~ practice, specification, repair_example$people, 'person')
  path <- tempfile(fileext='.rds'); on.exit(unlink(path), add=TRUE)
  saveRDS(list(specification=specification, fit=fit, repeated=repair_example$repeated,
               people=repair_example$people), path)
  result <- callr::r(function(path, source) {
    if (nzchar(source)) pkgload::load_all(source, quiet=TRUE) else library(repair)
    x <- readRDS(path)
    rebuilt <- repair(re('day') ~ practice, x$specification, x$people, 'person')
    replay <- repair_refit(x$fit, x$repeated, x$people)
    list(coef=coef(rebuilt), vcov=vcov(rebuilt), replay=coef(replay))
  }, list(path=path, source=if (pkgload::is_dev_package('repair')) normalizePath('.') else ''),
  libpath=.libPaths())
  expect_equal(result$coef, coef(fit), tolerance=1e-8)
  expect_equal(result$vcov, vcov(fit), tolerance=1e-8)
  expect_equal(result$replay, coef(fit), tolerance=1e-6)
})

test_that('invalid specifications and fitting failures have explicit boundaries', {
  d <- repair_example$repeated
  expect_error(repair_mlm(~day,d), 'two-sided')
  expect_error(repair_mlm(score~day,d[FALSE,]), 'nonempty')
  expect_error(repair_mlm(score~day,d, eng='bad'), 'arg')
  expect_error(repair_mlm(score~day,d, 1), 'arg')
  expect_error(repair_mlm(score~day,d, REML=TRUE,REML=FALSE), 'uniquely')
  expect_error(repair_mlm(score~day,d, engine='lmer', TRUE), 'uniquely')
  expect_error(repair_mlm(score~day,d, reml=TRUE), 'Unsupported')
  for (value in list(TRUE, 1, NA))
    expect_error(repair_mlm(score~day,d, devFunOnly=value), 'fitted MLM')
  expect_error(repair_mlm(score~day,d, engine='lme'), 'random')
  expect_error(repair_mlm(score~day,d, engine='lme',random=~day|person,keep.data=FALSE), 'keep.data')
  bad <- repair_mlm(score ~ absent + (day | person),d)
  expect_error(repair_effects(list(growth=bad)), 'growth.*lmer.*failed')
  d$duplicate <- d$day
  bad <- repair_mlm(score ~ day + duplicate + (day | person),d,
    control=lme4::lmerControl(check.rankX='warn+drop.cols'))
  expect_warning(repair_effects(bad), 'rank deficient')
})

test_that('lme specifications, fitted controls and bootstrap results persist', {
  skip_if_not_installed('nlme'); skip_if_not_installed('callr')
  specification <- repair_mlm(score ~ day, repair_example$repeated, engine='lme',
    random=~day|person, control=nlme::lmeControl(msMaxIter=300),
    correlation=NULL, weights=NULL, method='ML')
  fit <- repair(re('day') ~ practice, specification, repair_example$people, 'person')
  native <- nlme::lme(score ~ day, repair_example$repeated, random=~day|person,
    control=nlme::lmeControl(msMaxIter=300),method='ML')
  manual <- repair(re('day') ~ practice, native, repair_example$people, 'person')
  expect_spec_fit_equal(fit,manual)
  a <- repair_bootstrap(fit, repair_example$repeated, repair_example$people, B=2,seed=41)
  b <- repair_bootstrap(manual, repair_example$repeated, repair_example$people,B=2,seed=41)
  expect_identical(a$draws,b$draws)
  expect_equal(a$replicates,b$replicates,tolerance=1e-8)
  expect_identical(a$diagnostics$usable,b$diagnostics$usable)
  path <- tempfile(fileext='.rds'); on.exit(unlink(path), add=TRUE)
  saveRDS(list(specification=specification, fit=fit, repeated=repair_example$repeated,
               people=repair_example$people),path)
  package_path <- system.file(package='repair')
  from_source <- !file.exists(file.path(package_path,'Meta','package.rds'))
  result <- callr::r(function(path, source) {
    if (nzchar(source)) pkgload::load_all(source,quiet=TRUE) else library(repair)
    x <- readRDS(path)
    rebuilt <- repair(re('day')~practice,x$specification,x$people,'person')
    replay <- repair_refit(x$fit,x$repeated,x$people)
    list(coef=coef(rebuilt),vcov=vcov(rebuilt),replay=coef(replay))
  },list(path=path,source=if(from_source) package_path else ''),libpath=.libPaths())
  expect_equal(result$coef,coef(fit),tolerance=1e-8)
  expect_equal(result$vcov,vcov(fit),tolerance=1e-8)
  expect_equal(result$replay,coef(fit),tolerance=1e-6)
})

test_that('lme list random specifications preserve pdMat structure through replay', {
  skip_if_not_installed('nlme')
  random <- list(person=nlme::pdDiag(~day))
  spec <- repair_mlm(score ~ day, repair_example$repeated, engine='lme',
    random=random, control=nlme::lmeControl(msMaxIter=300))
  native <- nlme::lme(score ~ day, data=repair_example$repeated,random=random,
    control=nlme::lmeControl(msMaxIter=300))
  a <- repair(re('day') ~ practice,spec,repair_example$people,'person')
  b <- repair(re('day') ~ practice,native,repair_example$people,'person')
  expect_spec_fit_equal(a,b)
  replay <- repair_refit(a,repair_example$repeated,repair_example$people)
  expect_equal(coef(replay),coef(a),tolerance=1e-6)
  expect_identical(class(replay@mlm_fit$mlm$modelStruct$reStruct$person),class(random$person))
})
