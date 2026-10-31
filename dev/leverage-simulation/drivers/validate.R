#!/usr/bin/env Rscript
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE);here<-dirname(normalizePath(sub('^--file=','',arg[[1]])))
out<-file.path(here,'results/summary')
read_phase<-function(phase,N){
 f<-list.files(file.path(here,'results',phase),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
 stopifnot(length(f)==N);x<-lapply(f,readRDS);names(x)<-basename(f)
 stopifnot(all(vapply(x,function(r)is.null(r$error),logical(1))));list(files=f,runs=x)
}
main<-read_phase('main',600L);sens<-read_phase('sensitivity',18L)
for(phase in c('main','sensitivity')){
 p<-readRDS(file.path(here,'results',phase,'provenance.rds'))
 for(field in c('source_md5','reference_md5','library_md5'))stopifnot(identical(unname(tools::md5sum(names(p[[field]]))),unname(p[[field]])))
}
d<-readRDS(file.path(out,'dataset-intervals.rds'));p<-read.csv(file.path(out,'performance.csv'))
stopifnot(nrow(d)==12200L,nrow(p)==122L,!anyDuplicated(d[c('profile','rep','gamma','method')]))
for(x in split(d,interaction(d$profile,d$gamma,d$method,drop=TRUE))){
 r<-p[p$profile==x$profile[1] & p$gamma==x$gamma[1] & p$method==x$method[1],]
 stopifnot(nrow(x)==100L,nrow(r)==1L,sum(x$available)==r$available,
  all(is.na(x$low[!x$available])),all(is.na(x$high[!x$available])),all(is.na(x$covers[!x$available])))
 e<-x[x$available,]
 if(nrow(e))stopifnot(abs(mean(e$low<=e$gamma & e$high>=e$gamma)-r$coverage)<1e-12,
  abs(mean(e$low>0 | e$high<0)-r$rejection)<1e-12)
}
checks<-read.csv(file.path(out,'public_api_validation.csv'))
stopifnot(sum(checks$method=='bootstrap')==24L,sum(checks$attempts[checks$method=='bootstrap'])==2376L,
 sum(checks$method=='jackknife')==12L,sum(checks$attempts[checks$method=='jackknife'])==760L,
 max(checks$error)<1e-7)
max_cov_error<-max_gap<-0;failures<-list()
for(r in main$runs){
 raw<-readRDS(file.path(here,'results/main/native',r$raw));stopifnot(identical(raw$job,r$job))
 for(g in 1:2){
  inf<-r$analytic[[g]]$influence;lev<-r$analytic[[g]]$leverage
  stopifnot(identical(inf$ids,raw$data$people$person),identical(lev$ids,inf$ids))
  if(lev$diagnostics$full_available){
   stopifnot(inf$diagnostics$full_available,all(lev$diagnostics$deletions$usable),
    identical(rownames(lev$deletions),inf$ids))
   V<-lev$covariance;err<-max(abs(V$full-V$stage2-V$measurement-V$cross))/(1+max(abs(V$full)))
   stopifnot(err<1e-10,max(abs(V$full-crossprod(lev$deletions)))<1e-10)
   max_cov_error<-max(max_cov_error,err);max_gap<-max(max_gap,lev$diagnostics$hc0_relative_gap)
  } else {
   stopifnot(all(is.na(lev$covariance$full)))
   failures[[length(failures)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],reason=lev$diagnostics$full_reason)
  }
  usable<-r$deletions$usable[,g];values<-r$deletions$estimates[,,g]
  stopifnot(all(is.na(values[!usable,,drop=FALSE])),all(is.finite(values[usable,,drop=FALSE])))
  row<-r$rows[r$rows$gamma==c(0,.4)[g] & r$rows$method=='Delete-one Wald',]
  stopifnot(identical(row$available,all(usable)))
  if(all(usable)){
   expected<-(r$profile$n-1)/r$profile$n*crossprod(scale(values,scale=FALSE))
   stopifnot(abs(row$se^2-expected[3,3]*r$profile$tau^2)<1e-9)
  }
  row<-r$rows[r$rows$gamma==c(0,.4)[g] & r$rows$method=='Person percentile',]
  available<-raw$original[[g]]$preflight_ok && all(raw$ordinary$usable[,g])
  stopifnot(identical(row$available,available))
  if(available)stopifnot(max(abs(c(row$low,row$high)-quantile(raw$ordinary$estimates[,3,g],c(.025,.975))*r$profile$tau))<1e-12)
 }
}
failures<-if(length(failures))do.call(rbind,failures) else data.frame()
write.csv(failures,file.path(out,'leverage_refusals.csv'),row.names=FALSE)
for(r in sens$runs){
 name<-sprintf('p%02d_r%04d.rds',r$job$profile,r$job$rep)
 a<-readRDS(file.path(here,'results/main/native',name));b<-readRDS(file.path(here,'results/sensitivity/native',name))
 for(mode in c('ordinary','stratified')){
  if(is.null(a[[mode]]))next
  for(field in c('draws','usable','frozen_usable','measurement_ok')){
   x<-a[[mode]][[field]];y<-b[[mode]][[field]]
   prefix<-if(is.matrix(y))y[1:499,,drop=FALSE] else y[1:499]
   stopifnot(identical(x,prefix))
  }
  stopifnot(isTRUE(all.equal(a[[mode]]$estimates,b[[mode]]$estimates[1:499,,,drop=FALSE],tolerance=1e-12)),
   isTRUE(all.equal(a[[mode]]$frozen,b[[mode]]$frozen[1:499,,,drop=FALSE],tolerance=1e-12)))
 }
}
# Independent-stream / no reused-data verification, not just a seed label.
hash_input<-function(r){
 temp<-tempfile();on.exit(unlink(temp));writeBin(serialize(r$data[c('repeated','people','outcomes')],NULL,version=2),temp)
 unname(tools::md5sum(temp))
}
old_files<-list.files(file.path(here,'../robust_refit/results/main'),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
stopifnot(length(old_files)==600L)
old_hashes<-vapply(old_files,function(f)hash_input(readRDS(f)),character(1))
write.csv(data.frame(dataset=basename(old_files),md5=old_hashes),file.path(out,'previous-input-hashes.csv'),row.names=FALSE)
fresh_files<-list.files(file.path(here,'results/main/native'),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
hashes<-vapply(fresh_files,function(f)hash_input(readRDS(f)),character(1))
stopifnot(length(hashes)==600L,!anyDuplicated(hashes),!any(hashes%in%old_hashes))
write.csv(data.frame(dataset=basename(fresh_files),md5=hashes,reused_old=hashes%in%old_hashes),file.path(out,'freshness.csv'),row.names=FALSE)
steps<-read.csv(file.path(out,'derivative_sensitivity.csv'))
text<-c('PASS 600 fresh main and 18 sensitivity checkpoints; no driver/API errors.',
 'PASS frozen source, reference and installed-library hashes.',
 'PASS 12,200 endpoints / 122 performance cells independently recomputed; missing intervals remain NA.',
 sprintf('PASS 24 public B99 bootstrap calls / 2376 attempts; 12 complete public jackknives / 760 deletions; max coefficient gap %.3g.',max(checks$error)),
 sprintf('PASS HC3 uncentered shift covariance, signed decomposition and IDs; max decomposition error %.3g / augmented HC0 relative gap %.3g.',max_cov_error,max_gap),
 'PASS exact centered delete-one covariance, strict complete-vector gates and recorded percentile endpoints.',
 'PASS B499 arrays are exact B1999 prefixes for draws, eligibility and full/frozen coefficients.',
 'PASS 600 distinct fresh input hashes; none equals any of the 600 old study inputs.',
 sprintf('Half-step sensitivity: %d changed availability decisions; max matched relative covariance change %.3g.',sum(!steps$same),max(steps$relative,na.rm=TRUE)))
writeLines(text,file.path(out,'validation.txt'));cat(paste(text,collapse='\n'),'\n')
files<-c(main$files,sens$files,fresh_files,old_files)
write.csv(data.frame(path=sub('^.*/repair_bootstrap_pilot/','',files),md5=unname(tools::md5sum(files))),file.path(out,'raw-manifest.csv'),row.names=FALSE)
sources<-file.path(here,c('design.md','engine.R','run.R','summarize.R','plots.R','check-plots.R','validate.R'))
write.csv(data.frame(path=basename(sources),md5=unname(tools::md5sum(sources))),file.path(out,'analysis-source-manifest.csv'),row.names=FALSE)
