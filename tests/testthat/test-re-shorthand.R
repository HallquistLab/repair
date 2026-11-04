test_that('unqualified quoted and bare terms agree for a single fitted MLM', {
  ex <- workflow_example()
  quoted <- repair(re('Days') ~ x, ex$models, ex$data, 'person', regularization='step_up')
  bare <- repair(re(Days) ~ x, ex$models, ex$data, 'person', regularization='step_up')
  qualified <- repair(re('Days', model='growth') ~ x, ex$models, ex$data, 'person', regularization='step_up')
  expect_identical(bare@struc_fit,quoted@struc_fit)
  expect_identical(qualified@struc_fit,quoted@struc_fit)
  expect_identical(bare@analysis_data@effects,quoted@analysis_data@effects)
  expect_identical(repair_predictions(bare),repair_predictions(quoted))
  expect_identical(repair_prediction_covariance(bare),repair_prediction_covariance(quoted))
  expect_identical(repair_diagnostics(bare,'regularization'),repair_diagnostics(quoted,'regularization'))
  expect_output(print(bare),'REPAIR structural regression')
})

test_that('bare names are literal syntax and backticks support scalar coefficient names', {
  catalog <- data.frame(model='one',group='id',term=c('x','(Intercept)','x:y','non syntactic'),n_levels=30L)
  env <- new.env(parent=baseenv())
  makeActiveBinding('x',function(...) stop('x must never be evaluated'),env)
  bare <- as.formula(quote(response ~ re(x)),env=env)
  expected <- parse_struc_formula(response ~ re('x'),catalog)
  expect_identical(parse_struc_formula(bare,catalog)$effects,expected$effects)
  for (equation in list(response ~ re(`x:y`),response ~ re(`non syntactic`),
                         response ~ re(Intercept),response ~ re(`(Intercept)`))) {
    parsed <- parse_struc_formula(equation,catalog)
    term <- parsed$effects[[1L]]$term
    quoted <- as.formula(call('~',as.name('response'),call('re',term)))
    expect_identical(parsed$effects,parse_struc_formula(quoted,catalog)$effects)
  }
  expect_identical(parse_struc_formula(response ~ re(model='one',term=x),catalog)$effects,
                   expected$effects)
  env$selected <- 'x'
  other <- as.formula(quote(response ~ re(selected)),env=env)
  expect_error(parse_struc_formula(other,catalog),'no match for.*selected')
  expect_error(re(x),'formula marker')
})

test_that('model and group ambiguity is preserved and error messages identify matches', {
  catalog <- data.frame(model=c('a','a','b','a'),group=c('id','id:visit','id','id'),
                         term=c('x','x','x','unique'),n_levels=30L)
  expect_error(parse_struc_formula(y ~ re(x),catalog),'Specify.*model.*group')
  expect_error(parse_struc_formula(y ~ re(x,model='a'),catalog),'Specify.*group')
  a <- parse_struc_formula(y ~ re(x,model='a',group='id'),catalog)
  b <- parse_struc_formula(y ~ re('x',model='a',group='id'),catalog)
  expect_identical(a$effects,b$effects)
  expect_identical(parse_struc_formula(y ~ re(unique),catalog)$effects[[1]]$model,'a')
  # One MLM can still need a group qualifier; a model qualifier alone is insufficient.
  one <- catalog[catalog$model=='a',]
  expect_error(parse_struc_formula(y ~ re('x'),one),'Specify.*group')
  expect_error(parse_struc_formula(y ~ re(x),one),'Specify.*group')
  expect_identical(parse_struc_formula(y ~ re(x,group='id'),one)$effects,
                   parse_struc_formula(y ~ re('x',group='id'),one)$effects)
})

test_that('marker expressions and qualifier variables are not evaluated', {
  catalog <- data.frame(model='one',group='id',term='x',n_levels=30L)
  for (equation in list(y ~ re(get('x')), y ~ re(paste0('x')),y ~ re(x+1),
                        y ~ re({stop('must not execute'); 'x'}), y ~ re(1), y ~ re(.),
                        y ~ re(x,model=model_name),y ~ re(x,group=group_name),y ~ re()))
    expect_error(parse_struc_formula(equation,catalog),'literal string')
})

test_that('unquoted terms work in integrated calls for both native engines', {
  data <- repair_example$repeated; people <- repair_example$people
  specs <- list(repair_mlm(score ~ day + (day | person),data))
  if (requireNamespace('nlme',quietly=TRUE)) specs[[2]] <- repair_mlm(score ~ day,data,
    engine='lme',random=~day|person,control=nlme::lmeControl(msMaxIter=300))
  for (spec in specs) for (backend in c('r','rcpp')) {
    bare <- repair(practice ~ re(day),spec,people,'person',backend=backend)
    quoted <- repair(practice ~ re('day'),spec,people,'person',backend=backend)
    expect_identical(bare@struc_fit,quoted@struc_fit)
    expect_identical(bare@analysis_data@effects,quoted@analysis_data@effects)
  }
})

test_that('refitting and bootstrap preserve the bare-name analysis', {
  d <- repair_example$repeated; people <- repair_example$people
  model <- lme4::lmer(score ~ day + (day | person),d,
    control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000)))
  bare <- repair(re(day) ~ practice,model,people,'person',regularization='step_up')
  quoted <- repair(re('day') ~ practice,model,people,'person',regularization='step_up')
  replay <- repair_refit(bare,d,people)
  expect_equal(coef(replay),coef(bare),tolerance=1e-6)
  expect_equal(vcov(replay),vcov(bare),tolerance=1e-6)
  a <- repair_bootstrap(bare,d,people,B=2,seed=41)
  b <- repair_bootstrap(quoted,d,people,B=2,seed=41)
  expect_identical(a$draws,b$draws)
  expect_identical(a$replicates,b$replicates)
  expect_identical(a$diagnostics$preliminary_alpha,b$diagnostics$preliminary_alpha)
  expect_identical(a$diagnostics$final_alpha,b$diagnostics$final_alpha)
  expect_identical(vcov(a),vcov(b))
  expect_identical(confint(a),confint(b))
})
