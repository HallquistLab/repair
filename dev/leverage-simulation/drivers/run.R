#!/usr/bin/env Rscript
args<-commandArgs(TRUE);phase<-if(length(args))args[1] else 'main'
file_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
study_dir<-dirname(normalizePath(sub('^--file=','',file_arg[[1]])))
root<-normalizePath(file.path(study_dir,'../../../..'))
lib<-Sys.getenv('REPAIR_STUDY_LIB','/private/tmp/repair-leverage-lib');.libPaths(c(lib,.libPaths()))
suppressPackageStartupMessages(library(repair));source(file.path(study_dir,'engine.R'));load_sources(root)
stopifnot(phase%in%c('benchmark','main','sensitivity'))
profiles<-make_profiles(root)
RNGkind("L'Ecuyer-CMRG",'Inversion','Rejection');set.seed(202610319L)
streams<-vector('list',1000L);streams[[1]]<-.Random.seed
for(i in 2:1000)streams[[i]]<-parallel::nextRNGStream(streams[[i-1]])
reps<-switch(phase,benchmark=801L,main=1:100,sensitivity=1:3)
B<-switch(phase,benchmark=9L,main=499L,sensitivity=1999L)
jobs<-unlist(lapply(reps,function(rep)lapply(1:6,function(p)list(profile=p,rep=rep,stream=streams[[rep]],bootstrap_seed=340000000L+rep))),recursive=FALSE)
out_dir<-file.path(study_dir,'results',phase);dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
sources<-file.path(study_dir,c('run.R','engine.R','design.md'))
reference<-file.path(root,c('vig_hallquist_2026/analysis/repair_bootstrap_pilot/robust_refit/engine.R','R/reliability_calibration.R'))
provenance<-list(time=Sys.time(),phase=phase,B=B,master_seed=202610319L,profiles=profiles,
 source_md5=tools::md5sum(sources),reference_md5=tools::md5sum(reference),
 library_md5=tools::md5sum(file.path(find.package('repair'),c('DESCRIPTION','R/repair.rdb','libs/repair.so'))),session=sessionInfo())
manifest<-file.path(out_dir,'provenance.rds')
if(file.exists(manifest)){
 old<-readRDS(manifest)
 for(key in c('source_md5','reference_md5','library_md5'))stopifnot(identical(old[[key]],provenance[[key]]))
}else{saveRDS(provenance,manifest);writeLines(capture.output(sessionInfo()),file.path(out_dir,'sessionInfo.txt'))}
workers<-as.integer(Sys.getenv('REPAIR_STUDY_WORKERS',if(phase=='benchmark')'2' else '12'))
cl<-parallel::makePSOCKcluster(workers,outfile=file.path(out_dir,'workers.log'))
tryCatch({
 parallel::clusterExport(cl,c('lib','study_dir','root','profiles','B','phase','out_dir'),envir=.GlobalEnv)
 parallel::clusterEvalQ(cl,{
  .libPaths(c(lib,.libPaths()));Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
  suppressPackageStartupMessages(library(repair));source(file.path(study_dir,'engine.R'));load_sources(root);NULL
 })
 invisible(parallel::parLapplyLB(cl,jobs,function(job)run_fresh(job,profiles,B,phase,out_dir),chunk.size=1L))
},finally=parallel::stopCluster(cl))
