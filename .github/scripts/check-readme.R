# Run after check-r-package, using the installed package inside its check directory.
# Rendering in a temporary directory keeps generated documentation out of the tree.
library_path <- normalizePath("check/repair.Rcheck", mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
stopifnot(normalizePath(find.package("repair")) == normalizePath(file.path(library_path, "repair")))
output <- tempfile("repair-readme-")
dir.create(output)
source_file <- file.path(output, "README.Rmd")
if (!file.copy("README.Rmd", source_file)) stop("Could not copy README source.")
rmarkdown::render(source_file, quiet = TRUE, envir = new.env(parent = globalenv()))
unlink(output, recursive = TRUE)
cat("README executed against the checked package.\n")
