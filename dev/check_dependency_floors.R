#!/usr/bin/env Rscript
# Check a built source archive with the supported older dependency stack.
# Installs only into a new temporary library; does not alter the user's library.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !file.exists(args[[1L]]))
  stop('Usage: Rscript dev/check_dependency_floors.R /path/to/repair_version.tar.gz')
archive <- normalizePath(args[[1L]])
work <- tempfile('repair-dependency-floors-', tmpdir = dirname(tempdir()))
dir.create(work)
lib <- file.path(work, 'library')
dir.create(lib)
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
versions <- c(Rcpp = '1.0.12', S7 = '0.2.0',
              RcppArmadillo = '0.12.8.4.0', lme4 = '1.1-35.5')
cat('Check directory: ', work, '\n', sep = '')
for (package in names(versions)) {
  filename <- paste0(package, '_', versions[[package]], '.tar.gz')
  url <- paste0('https://cran.r-project.org/src/contrib/Archive/',
                package, '/', filename)
  source <- file.path(work, filename)
  utils::download.file(url, source, mode = 'wb')
  utils::install.packages(source, repos = NULL, type = 'source', lib = lib)
  stopifnot(as.character(utils::packageVersion(package, lib.loc = lib)) ==
              versions[[package]])
}
metadata <- installed.packages()
write.table(metadata[, c('Package', 'Version', 'LibPath')],
            file.path(work, 'dependencies.tsv'), sep = '\t',
            row.names = FALSE, quote = FALSE)
setwd(work)
# R CMD check installs/recompiles the archive against these older headers;
# it does not reuse a DLL compiled with the current RcppArmadillo headers.
status <- system2(file.path(R.home('bin'), 'R'),
                  c('CMD', 'check', '--as-cran', '--no-manual', shQuote(archive)),
                  stdout = 'check.log', stderr = 'check.log')
cat('R CMD check exit status: ', status, '\nLogs: ', work, '\n', sep = '')
check_log <- file.path('repair.Rcheck', '00check.log')
problems <- if (file.exists(check_log))
  any(grepl('ERROR|WARNING', readLines(check_log, warn = FALSE))) else TRUE
if (status != 0L || problems)
  stop('Dependency-floor check failed or reported warnings; inspect check.log.')
