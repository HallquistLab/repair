multimodel_keys <- function(ex) list(
  growth = setNames(as.character(round(4 * ex$case$a$z)), rownames(ex$case$a)),
  response = setNames(as.character(round(5 * (ex$case$b$w - .5))), rownames(ex$case$b)))

observation_blocks <- function(ex, keys, rho = .25) {
  ids <- ex$case$people$person
  blocks <- lapply(ids, function(id) {
    a <- keys$growth[rownames(ex$case$a)[ex$case$a$id == id]]
    b <- keys$response[rownames(ex$case$b)[ex$case$b$subject == id]]
    value <- outer(unname(a), unname(b), `==`) * rho *
      sigma(ex$models$growth) * sigma(ex$models$response)
    dimnames(value) <- list(a, b)
    value
  })
  list(list(models = c("growth", "response"), covariance = setNames(blocks, ids)))
}

test_that("iid propagation reproduces independently generated prediction covariances and fits", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); keys <- multimodel_keys(ex)
  blocks <- repair_cross_covariance(ex$formula, ex$models, keys, rho = .25,
    ids = ex$case$people$person)
  expect_identical(blocks[[1]]$models, c("response", "growth"))
  expect_equal(aperm(blocks[[1]]$covariance, c(2, 1, 3)),
    ex$case$cross_model_error[[1]]$covariance, tolerance = 1e-12)
  explicit <- repair_cross_covariance(ex$formula, ex$models, keys,
    observation_blocks(ex, keys), ids = ex$case$people$person)
  expect_equal(explicit, blocks, tolerance = 1e-12)
  for (backend in c("r", "rcpp")) {
    actual <- repair(ex$formula, ex$models, ex$case$people, "person",
      cross_model_error = blocks, backend = backend)
    expect_equal(unname(coef(actual)), ex$case$coefficients$estimate, tolerance = 1e-5)
    expect_equal(unname(vcov(actual)), unname(ex$case$vcov), tolerance = 1e-5)
    expect_equal(unname(confint(actual)),
      unname(as.matrix(ex$case$coefficients[, c("ci_low", "ci_high")])), tolerance = 1e-5)
    zero <- repair_cross_covariance(ex$formula, ex$models, keys, rho = 0)
    direct <- repair(ex$formula, ex$models, ex$case$people, "person",
      cross_model_error = "independent", backend = backend)
    from_zero <- repair(ex$formula, ex$models, ex$case$people, "person",
      cross_model_error = zero, backend = backend)
    expect_identical(from_zero@struc_fit, direct@struc_fit)
  }
})

test_that("pair axes, matrix order, ID order, model order and input row order are keyed", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); keys <- multimodel_keys(ex)
  ids <- rev(ex$case$people$person)[1:12]
  raw <- observation_blocks(ex, keys)
  expected <- repair_cross_covariance(ex$formula, ex$models, keys, raw, ids = ids)
  raw[[1]]$models <- rev(raw[[1]]$models)
  raw[[1]]$covariance <- rev(lapply(raw[[1]]$covariance, function(x) {
    x <- t(x); x[rev(seq_len(nrow(x))), rev(seq_len(ncol(x))), drop = FALSE]
  }))
  actual <- repair_cross_covariance(ex$formula, rev(ex$models), lapply(keys, rev), raw, ids = ids)
  expect_identical(actual, expected)
  expect_identical(dimnames(actual[[1]]$covariance)[[3]], ids)
  # One coefficient per model preserves all three dimensions.
  f <- re(w, "response") ~ re(z, "growth")
  one <- repair_cross_covariance(f, ex$models, keys, raw, ids = ids[1])
  expect_identical(dim(one[[1]]$covariance), c(1L, 1L, 1L))
})

test_that("iid matched occasions allow explicitly unmatched observations without row pairing", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); keys <- multimodel_keys(ex)
  id <- ex$case$people$person[1]
  arows <- which(ex$case$a$id == id); brows <- which(ex$case$b$subject == id)
  keys$response[rownames(ex$case$b)[brows[1]]] <- "unmatched"
  actual <- repair_cross_covariance(ex$formula, ex$models, keys, rho = -.3, ids = id)
  za <- cbind(1, ex$case$a$z[arows]); zb <- cbind(1, ex$case$b$w[brows])
  la <- solve(crossprod(za), t(za)); lb <- solve(crossprod(zb), t(zb))
  rab <- diag(-.3 * sigma(ex$models$growth) * sigma(ex$models$response), nrow(za))
  rab[, 1] <- 0
  expect_equal(unname(actual[[1]]$covariance[, , 1]),
    unname((lb %*% t(rab) %*% t(la))[2, , drop = TRUE]), tolerance = 1e-12)
  for (rho in c(-1, 1)) expect_no_error(
    repair_cross_covariance(ex$formula, ex$models, multimodel_keys(ex), rho = rho, ids = id))
})

test_that("rectangular observation blocks follow native subset rows and explicit keys", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); keys <- multimodel_keys(ex)
  id <- ex$case$people$person[1]
  dat <- ex$case$b
  removed <- which(dat$subject == id)[1]
  dat$keep <- seq_len(nrow(dat)) != removed
  # Original keys include the omitted row; the fitted model does not.
  ex$models$response <- nlme::lme(q ~ w, random = ~ w | subject, data = dat,
    subset = keep)
  ex$case$b <- dat[-removed, ]
  raw <- observation_blocks(ex, keys, rho = -.2)
  expect_identical(dim(raw[[1]]$covariance[[id]]), c(7L, 6L))
  actual <- repair_cross_covariance(ex$formula, ex$models, keys, raw, ids = id)
  shortcut <- repair_cross_covariance(ex$formula, ex$models, keys, rho = -.2, ids = id)
  expect_equal(actual, shortcut, tolerance = 1e-12)
  da <- ex$case$a[ex$case$a$id == id, ]; db <- ex$case$b[ex$case$b$subject == id, ]
  za <- cbind(1, da$z); zb <- cbind(1, db$w)
  la <- solve(crossprod(za), t(za)); lb <- solve(crossprod(zb), t(zb))
  oracle <- lb %*% t(raw[[1]]$covariance[[id]]) %*% t(la)
  expect_equal(unname(actual[[1]]$covariance[, , id]), unname(oracle[2, ]), tolerance = 1e-12)
})

test_that("explicit propagation uses GLS for correlated heterogeneous errors and omitted rows", {
  skip_if_not_installed("nlme")
  dat <- nlme_example_data()
  a <- nlme::lme(y ~ z + phase, random = ~ z | id, data = dat,
    correlation = nlme::corAR1(form = ~ time | id), weights = nlme::varIdent(form = ~ 1 | phase),
    na.action = na.exclude)
  dat$b <- dat$y + .1 * dat$z^2
  b <- lme4::lmer(b ~ z + phase + (z | id), dat, weights = scale_var, na.action = na.exclude)
  models <- list(a = a, b = b); f <- re(z, "b") ~ re(z, "a")
  keys <- list(a = setNames(as.character(dat$time), rownames(dat)),
    b = setNames(as.character(dat$time), rownames(dat)))
  oracle <- independent_nlme_gls(a, dat, "ar1_ident")
  ids <- rownames(oracle$estimates)[1:4]
  raw <- lapply(ids, function(id) {
    rows <- match(rownames(a$groups), rownames(dat))
    d <- dat[rows, ]; d <- d[as.character(d$id) == id, ]
    ra <- oracle$residual_blocks[[id]]
    rb <- diag(sigma(b)^2 / d$scale_var)
    # A rank-one contraction of whitened errors guarantees joint PSD.
    cross <- .4 * outer(t(chol(ra))[, 1], t(chol(rb))[, 1])
    dimnames(cross) <- list(as.character(d$time), as.character(d$time))
    cross
  })
  blocks <- list(list(models = c("a", "b"), covariance = setNames(raw, ids)))
  actual <- repair_cross_covariance(f, models, keys, blocks, ids = ids)
  for (id in ids) {
    d <- dat[match(rownames(a$groups), rownames(dat)), ]
    d <- d[as.character(d$id) == id, ]
    z <- cbind(1, d$z); ra <- oracle$residual_blocks[[id]]
    rb <- diag(sigma(b)^2 / d$scale_var)
    la <- solve(crossprod(z, solve(ra, z)), t(solve(ra, z)))
    lb <- solve(crossprod(z, solve(rb, z)), t(solve(rb, z)))
    expect_equal(unname(actual[[1]]$covariance[, , id]),
      unname((lb %*% t(blocks[[1]]$covariance[[id]]) %*% t(la))[2, 2]), tolerance = 1e-10)
  }
  expect_error(repair_cross_covariance(f, models, keys, rho = .1, ids = ids), "iid homoscedastic")
})

test_that("joint observation PSD is checked even when propagated covariances look valid", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); keys <- multimodel_keys(ex)
  id <- ex$case$people$person[1]
  raw <- observation_blocks(ex, keys)
  # This error covariance lies in directions removed by the extraction operator.
  d <- ex$case$a[ex$case$a$id == id, ]; z <- cbind(1, d$z)
  u <- qr.resid(qr(z), seq_len(nrow(d))^2); u <- u / sqrt(sum(u^2))
  v <- u
  x <- 2 * sigma(ex$models$growth) * sigma(ex$models$response) * outer(u, v)
  dimnames(x) <- dimnames(raw[[1]]$covariance[[id]])
  raw[[1]]$covariance[[id]] <- x
  expect_error(repair_cross_covariance(ex$formula, ex$models, keys, raw, ids = id),
    "Joint observation covariance.*positive semidefinite")
  models <- c(ex$models, list(third = ex$models$growth))
  keys$third <- keys$growth
  f <- re(w, "response") ~ re(z, "growth") + re(z, "third")
  pairs <- list(c("growth", "response"), c("growth", "third"), c("response", "third"))
  scales <- c(growth = sigma(models$growth), response = sigma(models$response), third = sigma(models$third))
  raw <- lapply(pairs, function(pair) {
    x <- diag(-.6 * prod(scales[pair]), nrow(d)); dimnames(x) <- list(as.character(4*d$z), as.character(4*d$z))
    list(models = pair, covariance = setNames(list(x), id))
  })
  # Each pair with rho=-.6 is valid, but the joint three-model correlation is not.
  expect_error(repair_cross_covariance(f, models, keys, raw, ids = id), "Joint observation covariance.*positive semidefinite")
  for (j in seq_along(raw)) raw[[j]]$covariance[[id]][] <- 0
  expect_length(repair_cross_covariance(f, models, keys, raw, ids = id), 3)
  expect_error(repair_cross_covariance(f, models, keys, raw[c(1, 1, 3)], ids = id), "Duplicate")
  expect_error(repair_cross_covariance(f, models, keys, raw[1:2], ids = id), "every.*pair")
})

test_that("invalid identities and covariance specifications are rejected", {
  skip_if_not_installed("nlme")
  ex <- multimodel_example(); keys <- multimodel_keys(ex); id <- ex$case$people$person[1]
  raw <- observation_blocks(ex, keys)
  expect_error(repair_cross_covariance(ex$formula, ex$models, keys), "exactly one")
  expect_error(repair_cross_covariance(ex$formula, ex$models, keys, raw, rho = .2), "exactly one")
  for (bad in list(NA_real_, Inf, 1.1, c(0, 1), "0")) expect_error(
    repair_cross_covariance(ex$formula, ex$models, keys, rho = bad), "finite correlation")
  for (bad in list(character(), c(id, id), "unknown", NA_character_, 1)) expect_error(
    repair_cross_covariance(ex$formula, ex$models, keys, raw, ids = bad), "unique character IDs")
  bad <- keys; bad$growth <- unname(bad$growth)
  expect_error(repair_cross_covariance(ex$formula, ex$models, bad, raw), "fitted row names")
  bad <- keys; bad$growth[ex$case$a$id == id] <- "duplicate"
  expect_error(repair_cross_covariance(ex$formula, ex$models, bad, raw, ids = id), "Duplicate occasion")
  bad <- raw; dimnames(bad[[1]]$covariance[[id]]) <- NULL
  expect_error(repair_cross_covariance(ex$formula, ex$models, keys, bad, ids = id), "names/dimensions")
  bad <- raw; bad[[1]]$covariance[[id]][1] <- NA
  expect_error(repair_cross_covariance(ex$formula, ex$models, keys, bad, ids = id), "names/dimensions")
  bad <- raw; bad[[1]]$covariance[[id]] <- NULL
  expect_error(repair_cross_covariance(ex$formula, ex$models, keys, bad, ids = id), "covering all")
  expect_error(repair_cross_covariance(y ~ re(z), ex$models, keys, raw), "at least two")
})

test_that("effective nuisance covariance is propagated and iid shortcut rejects it", {
  ex <- nested_example(); case <- ex$case
  models <- list(a = ex$model, b = ex$model)
  rows <- rownames(model.frame(ex$model))
  # Row identity itself is a valid occasion key for these identical data sets.
  keys <- list(a = setNames(rows, rows), b = setNames(rows, rows))
  ids <- case$people$person[1:3]
  raw <- lapply(ids, function(id) {
    obs <- which(as.character(model.frame(ex$model)$person) == id)
    x <- .2 * case$covariance[[id]]
    dimnames(x) <- list(rows[obs], rows[obs])
    x
  })
  f <- re(x, "b", "person") ~ re(x, "a", "person")
  blocks <- list(list(models = c("a", "b"), covariance = setNames(raw, ids)))
  actual <- repair_cross_covariance(f, models, keys, blocks, ids = ids)
  expect_equal(unname(actual[[1]]$covariance),
    unname(.2 * case$error["x", "x", ids, drop = FALSE]), tolerance = 1e-5)
  expect_error(repair_cross_covariance(f, models, keys, rho = .2, ids = ids),
    "without nuisance random groups")
})

test_that("failed likelihood extraction stops instead of dropping an ID", {
  ex <- workflow_example(); models <- list(a = ex$models$growth, b = ex$models$growth)
  rows <- rownames(model.frame(models$a))
  keys <- list(a = setNames(rows, rows), b = setNames(rows, rows))
  extract <- repair_get_mlm_elements
  local_mocked_bindings(repair_get_mlm_elements = function(model) {
    elements <- extract(model)
    design <- elements@random_design
    # The scorer also fails a target design whose columns are identical.
    design$Subject[, 2] <- design$Subject[, 1]
    elements@random_design <- design
    elements
  })
  expect_error(repair_cross_covariance(re(Days, "b") ~ re(Days, "a"),
    models, keys, rho = 0), "extraction design could not be solved")
})
