# Run from the repair source root. Never load repair to generate expectations.
# Rscript dev/generate_reliability_fixtures.R /path/to/blup_unshrinking
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else "../blup_unshrinking")
path <- "R/reliability_calibration.R"
ref <- new.env(parent = globalenv()); sys.source(file.path(root, path), ref)
source_cases <- readRDS("inst/extdata/nlme-workflow-reference.rds")$cases
cases <- lapply(source_cases, function(case) {
  correlated <- case$kind == "both_ar1_ident"
  model <- nlme::lme(y ~ z, random = ~ z | id, data = case$data,
    correlation = if (correlated) nlme::corAR1(form = ~ time | id) else NULL,
    weights = if (correlated) nlme::varIdent(form = ~ 1 | phase) else NULL,
    control = nlme::lmeControl(msMaxIter = 300))
  g <- as.matrix(nlme::getVarCov(model, type = "random.effects"))
  ids <- rownames(nlme::ranef(model))
  z <- lapply(ids, function(id) model.matrix(~ z, case$data[case$data$id == id, ]))
  r <- lapply(ids, function(id) as.matrix(nlme::getVarCov(model, individual = id, type = "conditional")[[1]]))
  v <- lapply(seq_along(ids), function(i) ref$posterior_random_effect_covariance(g, z[[i]], r[[i]]))
  a <- c(-g[1, 2] / g[1, 1], 1)
  per_person <- data.frame(id = ids, n_observations = vapply(z, nrow, integer(1)),
    intercept = vapply(v, function(x) 1 - x[1, 1] / g[1, 1], numeric(1)),
    slope = vapply(v, function(x) 1 - x[2, 2] / g[2, 2], numeric(1)),
    residualized = vapply(v, function(x) 1 - drop(t(a) %*% x %*% a) / drop(t(a) %*% g %*% a), numeric(1)))
  summary <- c(intercept = ref$expected_contrast_reliability(g, z, r, c(1, 0)),
    slope = ref$expected_slope_reliability(g, z, r),
    residualized = ref$expected_residualized_slope_reliability(g, z, r))
  stopifnot(isTRUE(all.equal(unname(summary), unname(colMeans(per_person[3:5])), tolerance = 1e-12)))
  list(data = case$data, kind = case$kind, g = g, posterior = v,
       per_person = per_person, summary = summary)
})
fixture <- list(provenance = list(source_commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
  source_path = path, source_md5 = unname(tools::md5sum(file.path(root, path))),
  source_fixture_md5 = unname(tools::md5sum("inst/extdata/nlme-workflow-reference.rds")),
  R_version = R.version.string, nlme_version = as.character(packageVersion("nlme")),
  definition = "Manuscript empirical marginal and intercept-residualized posterior reliability; equal person weights"), cases = cases)
saveRDS(fixture, "inst/extdata/reliability-reference.rds", version = 2)
cat("Saved", length(cases), "independent reliability references.\n")
