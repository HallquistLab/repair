#!/usr/bin/env Rscript
file_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
study_dir<-dirname(normalizePath(sub('^--file=','',file_arg[[1]])))
.libPaths(c(Sys.getenv('REPAIR_STUDY_LIB','/private/tmp/repair-if-lib'),.libPaths()))
suppressPackageStartupMessages(library(repair));out<-file.path(study_dir,'results/summary')
read_phase<-function(phase,N){
 f<-list.files(file.path(study_dir,'results',phase),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
 stopifnot(length(f)==N);x<-lapply(f,readRDS);names(x)<-basename(f)
 stopifnot(all(vapply(x,function(r)is.null(r$error),logical(1))));list(files=f,runs=x)
}
main<-read_phase('main',600L);tail<-read_phase('sensitivity',18L)
for(phase in c('main','sensitivity')){
 pr<-readRDS(file.path(study_dir,'results',phase,'provenance.rds'))
 for(field in c('source_md5','comparator_md5','library_md5','raw_md5'))
  stopifnot(identical(unname(tools::md5sum(names(pr[[field]]))),unname(pr[[field]])))
}
records<-read.csv(file.path(out,'dataset_metrics.csv'));performance<-read.csv(file.path(out,'performance.csv'))
stopifnot(nrow(records)==17200L,!anyDuplicated(records[c('profile','rep','arm','gamma','method')]))
for(x in split(records,interaction(records$profile,records$arm,records$gamma,records$method,drop=TRUE))){
 stopifnot(nrow(x)==100L)
 p<-performance[performance$profile==x$profile[1] & performance$arm==x$arm[1] & performance$gamma==x$gamma[1] & performance$method==x$method[1],]
 stopifnot(nrow(p)==1L,sum(x$available)==p$available)
 use<-x$available
 if(any(use))stopifnot(abs(mean(x$low[use]<=x$gamma[use] & x$high[use]>=x$gamma[use])-p$coverage)<1e-12,
   abs(mean(x$low[use]>0 | x$high[use]<0)-p$rejection)<1e-12)
 stopifnot(all(is.na(x$low[!use])),all(is.na(x$covers[!use])),all(is.na(x$rejects_zero[!use])))
}
checks<-read.csv(file.path(out,'public_bootstrap_validation.csv'))
stopifnot(nrow(checks)==24L,sum(checks$B)==2376L,max(checks$error,na.rm=TRUE)<1e-7)
maximum_error<-maximum_mean<-0;full_count<-stage2_count<-0L;anchor<-list()
for(r in main$runs){
 old<-readRDS(file.path(study_dir,'../robust_refit/results/main',sprintf('p%02d_r%04d.rds',r$job$profile,r$job$rep)))
 stopifnot(identical(r$job,old$job))
 for(arm in c('predictor','outcome'))for(g in 1:2){
  inf<-r$influences[[paste(arm,g,sep='_')]]
  stopifnot(identical(inf$ids,old$data$people$person))
  if(grepl('reconstruct|do not reproduce',inf$diagnostics$full_reason,ignore.case=TRUE))stop(inf$diagnostics$full_reason)
  if(inf$diagnostics$stage2_available){
    stage2_count<-stage2_count+1L
    stopifnot(identical(rownames(inf$influence$stage2),inf$ids))
    maximum_mean<-max(maximum_mean,max(abs(colMeans(inf$influence$stage2))))
  }
  if(inf$diagnostics$full_available){
   full_count<-full_count+1L
   scale<-1+max(abs(inf$covariance$full))
   error<-max(abs(inf$covariance$full-inf$covariance$stage2-inf$covariance$measurement-inf$covariance$cross))/scale
   maximum_error<-max(maximum_error,error)
   stopifnot(error<1e-10,max(abs(inf$covariance$full-crossprod(inf$influence$full)/length(inf$ids)^2))/scale<1e-10,
    identical(rownames(inf$influence$full),inf$ids),max(abs(colMeans(inf$influence$full)))<1e-6)
  }
  if(arm=='predictor'){
   for(mode in c('ordinary','stratified')){
    b<-old[[mode]];if(is.null(b))next
    row<-r$rows[r$rows$arm=='predictor' & r$rows$gamma==c(0,.4)[g] &
      r$rows$method==if(mode=='ordinary')'Full person percentile' else 'Stratified person percentile',]
    available<-old$original[[g]]$preflight_ok && all(b$usable[,g])
    stopifnot(nrow(row)==1L,identical(row$available,available))
    if(available){ci<-quantile(b$estimates[,3,g],c(.025,.975))*r$profile$tau
      stopifnot(max(abs(c(row$low,row$high)-unname(ci)))<1e-12)}
   }
  }else{
   for(mode in c('ordinary','stratified')){
    b<-r$outcome_bootstrap[[paste(g,mode,sep='_')]];if(is.null(b))next
    stopifnot(identical(b$draws,old[[mode]]$draws))
    # Every eligible cached estimate is a complete finite structural vector in
    # the original driver; focal outputs deliberately do not claim intercept CIs.
    stopifnot(all(is.finite(b$estimates[b$usable])))
   }
   legacy<-r$checks[[paste(g,'legacy',sep='_')]]
   if(!is.null(legacy)){
    stopifnot(is.finite(legacy$hc0_gap))
    if(r$profile$small==r$profile$large){
     stopifnot(legacy$balanced_point_gap<1e-7)
     i<-r$rows[r$rows$arm=='outcome' & r$rows$gamma==c(0,.4)[g],]
     a<-i[i$method=='Stage2 influence Wald',];b<-i[i$method=='Legacy OLS stacked HC0 (t)',]
     gap<-if(a$available && b$available)abs(a$se^2-b$se^2) else NA_real_
     if(is.finite(gap))stopifnot(gap<1e-6)
     anchor[[length(anchor)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
       point_gap=legacy$balanced_point_gap,hc0_source_gap=legacy$hc0_gap,if_legacy_variance_gap=gap)
    }
   }
  }
 }
}
anchor<-do.call(rbind,anchor);write.csv(anchor,file.path(out,'balanced_ols_validation.csv'),row.names=FALSE)
for(r in tail$runs){
 a<-main$runs[[sprintf('p%02d_r%04d.rds',r$job$profile,r$job$rep)]]
 for(g in 1:2)for(mode in c('ordinary','stratified')){
  b<-r$outcome_bootstrap[[paste(g,mode,sep='_')]];old<-a$outcome_bootstrap[[paste(g,mode,sep='_')]]
  if(is.null(b))next
  stopifnot(identical(b$draws[1:499,,drop=FALSE],old$draws),identical(b$usable[1:499],old$usable),
   isTRUE(all.equal(b$estimates[1:499],old$estimates,tolerance=1e-12)),
   identical(b$frozen_usable[1:499],old$frozen_usable),isTRUE(all.equal(b$frozen[1:499],old$frozen,tolerance=1e-12)))
 }
}
steps<-read.csv(file.path(out,'derivative_sensitivity.csv'))
# Step sensitivity is evidence, not a condition used to cherry-pick a default.
# Changes or relatively large differences must be reported, not silently removed.
text<-c('PASS 600 main and 18 sensitivity comparison checkpoints, no driver/public-API errors.',
 'PASS frozen analysis/comparator/library hashes and all original native-bootstrap file hashes.',
 'PASS 172 method/profile/arm/effect summaries independently recomputed from 17,200 interval endpoints.',
 'PASS unavailable intervals remain NA; not counted as valid non-rejections.',
 'PASS all original noisy-predictor bootstrap intervals and person draws reproduced.',
 sprintf('PASS 24 complete public outcome-bootstrap calls / %d attempts; maximum focal coefficient error %.3g.',sum(checks$B),max(checks$error,na.rm=TRUE)),
 sprintf('PASS %d full / %d stage2 influence results; covariance decomposition max relative error %.3g.',full_count,stage2_count,maximum_error),
 sprintf('PASS centered, correctly aligned person influences; maximum stage2 mean %.3g.',maximum_mean),
 sprintf('PASS %d balanced OLS comparisons; max point gap %.3g / IF-HC0 variance gap %.3g.',nrow(anchor),max(anchor$point_gap),max(anchor$if_legacy_variance_gap,na.rm=TRUE)),
 'PASS B499 prefixes of B1999 outcome bootstraps: draws, finite-vector decisions, focal estimates and frozen estimates.',
 sprintf('Derivative-step sensitivity: full availability changes %d, stage2 changes %d, max matched full relative covariance change %.3g / stage2 %.3g.',
 sum(!steps$full_same),sum(!steps$stage2_same),max(steps$full_relative,na.rm=TRUE),max(steps$stage2_relative,na.rm=TRUE)))
writeLines(text,file.path(out,'validation.txt'));cat(paste(text,collapse='\n'),'\n')
manifest<-function(paths,name)write.csv(data.frame(path=paths,md5=unname(tools::md5sum(paths))),file.path(out,name),row.names=FALSE)
manifest(c(main$files,tail$files),'raw_manifest.csv')
manifest(file.path(study_dir,c('design.md','engine.R','run.R','summarize.R','validate.R')),'analysis_source_manifest.csv')
