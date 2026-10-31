#!/usr/bin/env Rscript
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE);here<-dirname(normalizePath(sub('^--file=','',arg[[1]])))
package<-normalizePath(args[1]);out<-file.path(package,'dev/leverage-simulation');dir.create(out,showWarnings=FALSE)
copy<-function(from,to)stopifnot(all(file.copy(from,to,overwrite=TRUE)))
summary<-file.path(here,'results/summary')
copy(list.files(summary,'[.](csv|rds|png|pdf|txt|json)$',full.names=TRUE),out)
copy(file.path(here,c('report.Rmd','report.html')),out)
copy(file.path(here,'design.md'),file.path(out,'protocol.md'))
copy(file.path(here,'README.md'),file.path(out,'study-README.md'))
dir.create(file.path(out,'drivers'),showWarnings=FALSE)
copy(file.path(here,c('engine.R','run.R','summarize.R','plots.R','check-plots.R','validate.R','audit_extremes.R','make_bundle.R')),file.path(out,'drivers'))
files<-list.files(file.path(here,'results/main'),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE);stopifnot(length(files)==600L)
shifts<-lapply(files,function(f){
 r<-readRDS(f);stopifnot(is.null(r$error))
 lev<-lapply(r$analytic,function(a)if(a$leverage$diagnostics$full_available)a$leverage$deletions else NULL)
 list(profile=r$profile$profile,rep=r$job$rep,n=r$profile$n,tau=r$profile$tau,ids=r$analytic[[1]]$influence$ids,
  estimate=lapply(r$analytic,function(a)a$influence$estimate),jackknife=r$deletions$estimates,usable=r$deletions$usable,leverage=lev)
})
saveRDS(shifts,file.path(out,'person-shifts.rds'),compress='xz')
provenance<-lapply(c('main','sensitivity'),function(phase){
 p<-readRDS(file.path(here,'results',phase,'provenance.rds'))
 for(key in c('source_md5','reference_md5','library_md5')){
  names(p[[key]])<-sub('^.*/repair_bootstrap_pilot/','',names(p[[key]]))
  names(p[[key]])<-sub('^.*/blup_unshrinking/','research/',names(p[[key]]))
  names(p[[key]])<-sub('^.*/repair-leverage-lib/','repair-leverage-lib/',names(p[[key]]))
 }
 p
})
names(provenance)<-c('main','sensitivity');saveRDS(provenance,file.path(out,'provenance.rds'))
copy(file.path(here,'results/main/sessionInfo.txt'),out)
for(name in c('sessionInfo.txt','initial-00check.log','initial-00install.out','testthat.Rout')){
 path<-file.path(out,name);if(file.exists(path))writeLines(sub('[[:blank:]]+$','',readLines(path)),path)
}
sources<-list.files(file.path(package,'R'),'[.]R$',full.names=TRUE)
write.csv(data.frame(path=paste0('R/',basename(sources)),md5=unname(tools::md5sum(sources))),file.path(out,'final-source-manifest.csv'),row.names=FALSE)
files<-list.files(out,recursive=TRUE,full.names=TRUE);files<-files[basename(files)!='bundle-manifest.csv']
write.csv(data.frame(file=substring(files,nchar(out)+2L),md5=unname(tools::md5sum(files))),file.path(out,'bundle-manifest.csv'),row.names=FALSE)
cat('Copied',length(files),'review files. Run the base-R verifier.\n')
