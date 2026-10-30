#!/usr/bin/env Rscript
# Copy a compact review archive; never alter original study checkpoints.
args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
file_arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
here <- dirname(normalizePath(sub('^--file=','',file_arg[[1]])))
package_root <- normalizePath(args[1])
out <- file.path(package_root,'dev/influence-simulation');dir.create(out,showWarnings=FALSE)
summary_dir <- file.path(here,'results/summary')
copy <- function(from,to) stopifnot(all(file.copy(from,to,overwrite=TRUE)))
files <- list.files(summary_dir,'[.](csv|png|pdf|rds|txt)$',full.names=TRUE)
files <- files[!basename(files)%in%c('dataset_metrics.csv','raw_manifest.csv','analysis_source_manifest.csv')]
copy(files,out)
saveRDS(read.csv(file.path(summary_dir,'dataset_metrics.csv')),file.path(out,'dataset-intervals.rds'),compress='xz')
copy(file.path(here,c('report.Rmd','report.html')),out)
copy(file.path(here,'design.md'),file.path(out,'protocol.md'))
copy(file.path(here,'README.md'),file.path(out,'study-README.md'))
drivers <- file.path(out,'drivers');dir.create(drivers,showWarnings=FALSE)
sources <- c('engine.R','run.R','summarize.R','validate.R','validate_guarded.R','make_bundle.R')
copy(file.path(here,sources),drivers)
copy(file.path(here,'results/main/sessionInfo.txt'),out)
normalize_paths <- function(paths) {
 paths <- sub('^.*/repair_bootstrap_pilot/','',paths)
 paths <- sub('^.*/blup_unshrinking/','research/',paths)
 sub('^.*/repair-if-lib/','repair-if-lib/',paths)
}
provenance <- lapply(c('main','sensitivity'),function(phase){
 p <- readRDS(file.path(here,'results',phase,'provenance.rds'))
 for(field in c('source_md5','comparator_md5','library_md5','raw_md5'))
  names(p[[field]]) <- normalize_paths(names(p[[field]]))
 p
})
names(provenance) <- c('main','sensitivity')
saveRDS(provenance,file.path(out,'provenance.rds'))
source_rows <- do.call(rbind,lapply(names(provenance),function(phase){
 p <- provenance[[phase]]
 do.call(rbind,lapply(c('source_md5','comparator_md5','library_md5'),function(field)
  data.frame(phase=phase,kind=field,path=names(p[[field]]),md5=unname(p[[field]]))))
}))
write.csv(source_rows,file.path(out,'source-manifest.csv'),row.names=FALSE)
raw <- read.csv(file.path(summary_dir,'raw_manifest.csv'));raw$path<-normalize_paths(raw$path)
write.csv(raw,file.path(out,'raw-manifest.csv'),row.names=FALSE)
final_sources <- list.files(file.path(package_root,'R'),'[.]R$',full.names=TRUE)
write.csv(data.frame(path=paste0('R/',basename(final_sources)),md5=unname(tools::md5sum(final_sources))),
 file.path(out,'final-source-manifest.csv'),row.names=FALSE)
# Normalize display whitespace in copied logs only; original check/session
# artifacts remain unchanged in the research archive.
for(name in c('sessionInfo.txt','testthat.Rout','00check.log','00install.out')) {
 path <- file.path(out,name)
 if(file.exists(path))writeLines(sub('[[:blank:]]+$','',readLines(path)),path)
}
files <- list.files(out,recursive=TRUE,full.names=TRUE)
files <- files[basename(files)!='bundle-manifest.csv']
write.csv(data.frame(file=substring(files,nchar(out)+2L),md5=unname(tools::md5sum(files))),
 file.path(out,'bundle-manifest.csv'),row.names=FALSE)
cat('Copied',length(files),'review files. Run the base-R verifier from the package.\n')
