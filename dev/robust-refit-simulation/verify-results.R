#!/usr/bin/env Rscript
# Independently recompute the published rates from retained interval endpoints.
file_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
dir <- dirname(normalizePath(sub('^--file=', '', file_arg[[1]])))
manifest <- read.csv(file.path(dir, 'bundle-manifest.csv'))
stopifnot(identical(manifest$md5, unname(tools::md5sum(file.path(dir, manifest$path)))))
x <- readRDS(file.path(dir, 'dataset-intervals.rds'))
p <- read.csv(file.path(dir, 'performance.csv'))
a <- read.csv(file.path(dir, 'old_robust_performance.csv'))
stopifnot(nrow(x) == 6400L, !anyDuplicated(x[c('profile', 'rep', 'gamma', 'method')]))
rate <- function(z) if (length(z)) mean(z) else NA_real_
close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-11, check.attributes = FALSE))
for (d in split(x, interaction(x$profile, x$gamma, x$method, drop = TRUE))) {
  stopifnot(nrow(d) == 100L)
  perf <- p[p$profile == d$profile[1] & p$gamma == d$gamma[1] & p$method == d$method[1], ]
  cmp <- a[a$profile == d$profile[1] & a$gamma == d$gamma[1] & a$method == d$method[1], ]
  stopifnot(nrow(perf) == 1L, nrow(cmp) == 1L)
  for (side in c('old', 'robust')) {
    available <- d[[paste0('available_', side)]]
    low <- d[[paste0('low_', side)]]; high <- d[[paste0('high_', side)]]
    stopifnot(identical(available, is.finite(low) & is.finite(high)))
    covered <- low[available] <= d$gamma[available] & high[available] >= d$gamma[available]
    rejected <- low[available] > 0 | high[available] < 0
    stopifnot(sum(available) == cmp[[paste0('available_', side)]],
      close(rate(covered), cmp[[paste0('coverage_', side)]]),
      close(rate(rejected), cmp[[paste0('rejection_', side)]]))
    if (side == 'robust') stopifnot(sum(available) == perf$available,
      close(rate(covered), perf$coverage), close(rate(rejected), perf$rejection))
  }
}
full <- a[a$method == 'Full person percentile' & a$gamma == .4, ]
stopifnot(sum(full$available_old) == 121L, sum(full$available_robust) == 584L, sum(full$lost) == 0L)
r <- read.csv(file.path(dir, 'optimizer_recovery.csv'))
stopifnot(sum(r$measurement_calls) == 349300L, sum(r$unresolved) == 27L,
  sum(r$recovered_old_failures, na.rm = TRUE) == 2349L,
  sum(r$retained_old_failures, na.rm = TRUE) == 14L, sum(r$new_failures, na.rm = TRUE) == 13L)
t <- read.csv(file.path(dir, 'tail_sensitivity.csv')); t <- t[t$inference == 'default', ]
stopifnot(nrow(t) == 42L, all(t$matched), all(t$covers_499 == t$covers_1999))
cat('PASS bundle checksums and all 64 method/profile/effect summaries against dataset-level endpoints.\n')
cat('PASS 121 -> 584 default ordinary intervals, 2349 recovered failures, 27 unresolved.\n')
cat('PASS 42 B499/B1999 paired default coverage decisions.\n')
