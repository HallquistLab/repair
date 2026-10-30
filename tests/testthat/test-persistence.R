# Keep this collector independent of the test environment: the child receives no
# original model/data bindings, and must dispatch from a reloaded fitted object.
persistence_outputs <- function(fit) {
  report <- base::summary(fit)
  report$formula <- base::paste(base::deparse(report$formula), collapse = " ")
  list(coefficients = stats::coef(fit), vcov = stats::vcov(fit),
    intervals = stats::confint(fit), n = stats::nobs(fit), summary = report,
    printed_summary = utils::capture.output(base::print(base::summary(fit))),
    coefficient_table = repair::repair_coefficients(fit),
    scores = repair::repair_scores(fit),
    measurement_scores = repair::repair_scores(fit, cohort = "measurement"),
    covariance = repair::repair_score_covariance(fit),
    reliability = repair::repair_reliability(fit, level = "person"),
    diagnostics = repair::repair_diagnostics(fit))
}
environment(persistence_outputs) <- baseenv()

expect_fresh_session_roundtrip <- function(fit) {
  expected <- persistence_outputs(fit)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  package_path <- system.file(package = "repair")
  # R CMD check uses the installed package. test_local() uses source loading in
  # the new process so it never accidentally checks a stale installed version.
  from_source <- !file.exists(file.path(package_path, "Meta", "package.rds"))
  # pkgload's synthetic package environment can be retained by nlme formulas.
  # R warns when saving that development-only environment; the fresh child loads
  # the same source explicitly. Never silence serialization warnings in checks
  # of the real installed package, or any other warning in either mode.
  withCallingHandlers(saveRDS(fit, path), warning = function(w) {
    if (from_source && identical(conditionMessage(w),
                                "'package:repair' may not be available when loading")) {
      invokeRestart("muffleWarning")
    }
  })
  restored <- callr::r(function(path, package_path, from_source, collect) {
    if (from_source) pkgload::load_all(package_path, quiet = TRUE) else library(repair)
    collect(readRDS(path))
  }, args = list(path = path, package_path = package_path,
                 from_source = from_source, collect = persistence_outputs),
     libpath = .libPaths(), timeout = 120)
  expect_identical(names(restored), names(expected))
  for (name in names(expected)) {
    expect_equal(restored[[name]], expected[[name]], info = name, tolerance = 1e-12)
  }
  invisible(restored)
}

test_that("a saved teaching fit retains the omitted cohort in a new R session", {
  example <- repair_example
  growth <- lme4::lmer(score ~ day + (day | person), example$repeated)
  people <- example$people[nrow(example$people):1, ]
  people$practice[2] <- NA_real_
  expect_warning(fit <- repair(re("day", "growth") ~ practice,
    list(growth = growth), people, "person", na_action = "omit"), "Omitting 1")
  expect_true(repair_diagnostics(fit, "status")$converged)
  restored <- expect_fresh_session_roundtrip(fit)
  expect_equal(restored$n, nrow(people) - 1L)
  expect_identical(restored$scores$id, people$person[-2])
  expect_identical(restored$diagnostics$excluded$id, people$person[2])
  expect_equal(sum(!restored$measurement_scores$analyzed), 1L)
})

test_that("a saved multi-model fit retains supplied cross-error covariance in a new R session", {
  ex <- multimodel_example("supplied")
  fit <- repair(ex$formula, ex$models, ex$case$people, "person",
    cross_model_error = ex$case$cross_model_error)
  restored <- expect_fresh_session_roundtrip(fit)
  cross <- restored$covariance[restored$covariance$model1 != restored$covariance$model2, ]
  expect_true(nrow(cross) > 0L)
  expect_true(all(cross$source == "supplied"))
  expect_true(any(cross$covariance != 0))
  # Retain agreement with an independent production-code reference, not just
  # agreement with this package's own pre-serialization output.
  expect_equal(restored$coefficient_table$estimate, ex$case$coefficients$estimate,
    tolerance = 1e-5)
  expect_equal(restored$coefficient_table$std.error, ex$case$coefficients$se,
    tolerance = 1e-5)
})

test_that("a saved nested-nuisance fit retains group identities and reliability in a new R session", {
  ex <- nested_example("nested_slope")
  fit <- repair(ex$formula, list(growth = ex$model), ex$case$people, "person")
  restored <- expect_fresh_session_roundtrip(fit)
  expect_identical(unique(restored$scores$group), "person")
  expect_setequal(restored$summary$measurement_models$growth$nuisance_groups,
    c("person:day", "person:block"))
})
