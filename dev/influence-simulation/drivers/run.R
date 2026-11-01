#!/usr/bin/env Rscript
args<-commandArgs(TRUE);mode<-if(length(args))args[1] else 'main'
file_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
study_dir<-dirname(normalizePath(sub('^--file=','',file_arg[[1]])))
research_root<-normalizePath(file.path(study_dir,'../../../..'))
lib<-Sys.getenv('REPAIR_STUDY_LIB','/private/tmp/repair-if-lib');.libPaths(c(lib,.libPaths()))
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
suppressPackageStartupMessages(library(repair));source(file.path(study_dir,'engine.R'));legacy<-load_comparison_sources(research_root)
stopifnot(mode%in%c('benchmark','main','sensitivity'))
source_dir<-file.path(study_dir,'../robust_refit/results',mode)
paths<-list.files(source_dir,'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
expected<-switch(mode,benchmark=6L,main=600L,sensitivity=18L);stopifnot(length(paths)==expected)
out_dir<-file.path(study_dir,'results',mode);dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
provenance<-list(time=Sys.time(),mode=mode,
 source_md5=tools::md5sum(file.path(study_dir,c('run.R','engine.R','design.md'))),
 comparator_md5=tools::md5sum(file.path(research_root,c('R/stacked_sandwich_helpers.R','R/derivative_backends.R'))),
 library_md5=tools::md5sum(file.path(find.package('repair'),c('DESCRIPTION','R/repair.rdb','libs/repair.so'))),
 raw_md5=tools::md5sum(paths),session=sessionInfo())
if(file.exists(file.path(out_dir,'provenance.rds'))){
 old<-readRDS(file.path(out_dir,'provenance.rds'))
 for(key in c('source_md5','comparator_md5','library_md5','raw_md5'))stopifnot(identical(old[[key]],provenance[[key]]))
}else{
 saveRDS(provenance,file.path(out_dir,'provenance.rds'))
 writeLines(capture.output(sessionInfo()),file.path(out_dir,'sessionInfo.txt'))
}
workers<-as.integer(Sys.getenv('REPAIR_STUDY_WORKERS',if(mode=='benchmark')'1' else '12'))
cl<-parallel::makePSOCKcluster(workers,outfile=file.path(out_dir,'workers.log'))
tryCatch({
 parallel::clusterExport(cl,c('study_dir','research_root','lib','out_dir','mode'),envir=.GlobalEnv)
 parallel::clusterEvalQ(cl,{
  .libPaths(c(lib,.libPaths()));Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  suppressPackageStartupMessages(library(repair));source(file.path(study_dir,'engine.R'));legacy<-load_comparison_sources(research_root);NULL
 })
 invisible(parallel::parLapplyLB(cl,paths,function(path){
  rep<-as.integer(sub('^p[0-9]+_r([0-9]+)[.]rds$','\\1',basename(path)))
  run_comparison(path,mode,out_dir,legacy,check_public=mode=='benchmark'||(mode=='main'&&rep<=2L),
    step_check=mode=='sensitivity')
 },chunk.size=1L))
},finally=parallel::stopCluster(cl))
