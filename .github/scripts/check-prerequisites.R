# Fail before R CMD check if optional test/vignette dependencies are missing.
# In particular, nlme tests must not quietly skip on an incomplete CI library.
description <- read.dcf("DESCRIPTION")
fields <- intersect(c("Depends", "Imports", "LinkingTo", "Suggests"), colnames(description))
packages <- trimws(unlist(strsplit(paste(description[1, fields], collapse = ","), ",")))
packages <- unique(trimws(sub("\\s*\\(.*\\)$", "", packages)))
packages <- setdiff(packages[nzchar(packages)], "R")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing check dependencies: ", paste(missing, collapse = ", "))
if (!rmarkdown::pandoc_available()) stop("Pandoc is required to execute package vignettes.")
cat("All declared dependencies available; Pandoc ", as.character(rmarkdown::pandoc_version()), "\n", sep = "")
print(sessionInfo())
