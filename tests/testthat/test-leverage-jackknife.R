lev_example <- function(REML = TRUE) {
  d <- repair_example
  m <- repair_lmer_refitter(score ~ day + (day | person), REML = REML)(d$repeated)
  list(model = m, repeated = d$repeated, people = d$people,
    fit = repair(re("day") ~ practice, m, d$people, "person", backend = "r"))
}

test_that("deleted ML and REML scores retain the exact global likelihood", {
  for (reml in c(FALSE, TRUE)) {
    ex <- lev_example(reml); ctx <- repair_if_context(ex$fit)
    out <- repair_if_likelihood(ctx, ctx$psi, deletions = TRUE)
    for (i in c(2L, 21L)) {
      f <- rep(1, length(ctx$ids)); f[i] <- 0
      manual <- repair_if_likelihood(ctx, ctx$psi, frequency = f)
      expect_equal(unname(out$deletion_scores[i, ]), unname(manual$score), tolerance = 1e-9)
    }
  }
})

test_that("augmented equations reproduce the complete canonical estimator", {
  ex <- lev_example()
  for (reml in c(TRUE, FALSE)) {
    if (!reml) ex <- lev_example(FALSE)
    for (policy in c("minimum", "none", "custom")) {
      f <- repair(practice ~ re('(Intercept)') + re('day'), ex$model, ex$people, 'person',
        regularization = policy, alpha = if (policy == "custom") 7 else NULL, backend = "r")
      ctx <- repair_if_context(f); inputs <- repair_if_inputs(f, ctx, ctx$psi)
      state <- repair_lev_state(inputs); eq <- repair_lev_equations(inputs, state$value, state$rows)
      expect_equal(tail(state$value, 3), unname(coef(f)), ignore_attr = TRUE, tolerance = 1e-8)
      expect_equal(colMeans(eq$value), rep(0, ncol(eq$value)), tolerance = 1e-8, ignore_attr = TRUE)
      expect_identical(state$signature, eq$signature)
    }
  }
})

test_that("full-system correction has an independent OLS HC3 anchor", {
  set.seed(8251); n <- 40; time <- seq(-1, 1, length.out = 8)
  x <- rnorm(n); b0 <- rnorm(n); b1 <- .5*x + rnorm(n)
  repeated <- do.call(rbind, lapply(1:n, function(i) data.frame(person=sprintf('p%03d',i),time=time,
    y=.4*time+b0[i]+b1[i]*time+rnorm(8,sd=.5))))
  people <- data.frame(person=sprintf('p%03d',1:n), x=x)
  m <- repair_lmer_refitter(y~time+(time|person))(repeated)
  f <- repair(re('time')~x,m,people,'person',backend='r')
  value <- repair_leverage(f)
  expect_true(value$diagnostics$full_available)
  y <- f@analysis_data@outcome; X <- cbind(1,x); inverse <- solve(crossprod(X))
  beta <- drop(inverse %*% crossprod(X,y)); residual <- y-drop(X%*%beta)
  h <- rowSums((X%*%inverse)*X)
  manual <- inverse %*% crossprod(X * (residual/(1-h))) %*% inverse
  expect_equal(unname(vcov(value,component='stage2')),unname(manual),tolerance=2e-6)
  expect_equal(vcov(value)[2,2],manual[2,2],tolerance=2e-6)
  expect_lt(value$diagnostics$hc0_relative_gap,1e-5)
})

test_that("leverage keeps cross-stage terms, IDs, original estimates and RNG", {
  ex <- lev_example(); f <- repair(practice~re('(Intercept)')+re('day'),ex$model,ex$people,'person',backend='r')
  before <- list(coef=coef(f),vcov=vcov(f)); set.seed(24); rng <- .Random.seed
  value <- repair_leverage(f)
  expect_true(value$diagnostics$full_available); expect_identical(.Random.seed,rng)
  expect_equal(coef(value),coef(f));expect_identical(rownames(value$deletions),ex$people$person)
  expect_equal(vcov(value),crossprod(value$deletions),tolerance=1e-12)
  expect_equal(vcov(value),vcov(value,component='stage2')+vcov(value,component='measurement')+vcov(value,component='cross'),tolerance=1e-12)
  expect_equal(before,list(coef=coef(f),vcov=vcov(f)))
  expect_match(repair_sandwich_coefficients(value)$inference[1],'HC3')
  expect_output(print(value),'HC3')
  saved <- tempfile();saveRDS(value,saved);expect_equal(vcov(readRDS(saved)),vcov(value));unlink(saved)
  expect_error(repair_lev_solve(matrix(0,2,2),c(1,1),c(1,1)),'Zero/invalid')
})

test_that("jackknife refits complete people and uses the centered delete-one formula", {
  ex <- lev_example(); before <- coef(ex$fit);set.seed(921); rng <- .Random.seed
  value <- repair_jackknife(ex$fit,ex$repeated,ex$people)
  expect_true(all(value$diagnostics$usable)); expect_identical(.Random.seed,rng)
  expect_identical(value$ids,ex$people$person);expect_equal(nrow(value$replicates),nrow(ex$people))
  expect_equal(coef(value),before);expect_equal(coef(ex$fit),before)
  n <- nobs(value); centered <- sweep(value$replicates,2,colMeans(value$replicates))
  expect_equal(vcov(value),(n-1)/n*crossprod(centered),tolerance=1e-12)
  expect_equal(value$bias,(n-1)*(colMeans(value$replicates)-before))
  for (i in c(2L,17L)) {
    actual <- repair_refit(ex$fit,ex$repeated,ex$people,draw=ex$people$person[-i])
    expect_equal(unname(value$replicates[i,]),unname(coef(actual)),tolerance=1e-10)
  }
  table <- repair_jackknife_coefficients(value)
  expect_equal(table$std.error,sqrt(diag(vcov(value))),ignore_attr=TRUE)
  expect_equal(unname(confint(value)),unname(cbind(table$conf.low,table$conf.high)))
  expect_equal(repair_jackknife_coefficients(value,reference='t')$df,rep(n-1,ncol(value$replicates)))
  expect_output(print(value),'delete-one')
  value$diagnostics$usable[1] <- FALSE;value$replicates[1,] <- NA
  expect_error(vcov(value),'all N vectors');expect_null(summary(value)$coefficients)
  expect_error(repair_jackknife_coefficients(value),'all N vectors')
})

test_that("leverage failure never exposes uncorrected covariance as corrected", {
  ex <- lev_example(); ctx <- repair_if_context(ex$fit); ctx$singular <- TRUE
  base <- repair_sandwich_compute(ex$fit,1e-4,.03,context=ctx)
  expect_false(base$diagnostics$full_available)
  out <- repair_leverage_compute(ex$fit,1e-4,.03,ctx,base)
  expect_false(out$diagnostics$full_available)
  expect_true(all(is.na(out$covariance$full)))
  expect_error(vcov(out),'singular')
  expect_error(repair_leverage(ex$fit,step=0),'positive finite')
})

test_that("actual deletion failures retain every ID and restore callback RNG use", {
  ex <- lev_example(); set.seed(42); seed <- .Random.seed
  callback <- function(data) {
    runif(1)
    if (length(unique(data$person)) < nrow(ex$people)) stop('deliberate deletion failure')
    ex$model
  }
  out <- repair_jackknife(ex$fit,ex$repeated,ex$people,refitter=callback)
  expect_identical(.Random.seed,seed)
  expect_equal(nrow(out$diagnostics),nrow(ex$people))
  expect_true(all(!out$diagnostics$usable));expect_true(all(is.na(out$replicates)))
  expect_true(all(grepl('deliberate deletion failure',out$diagnostics$reason)))
  expect_error(confint(out),'all N vectors')
})
