#!/usr/bin/env Rscript
# Base R verification from retained endpoints and actual person-shift arrays.
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-dirname(normalizePath(sub('^--file=','',arg[[1]])))
read<-function(name)read.csv(file.path(root,name),stringsAsFactors=FALSE)
same<-function(a,b)isTRUE(all.equal(a,b,tolerance=1e-10,check.attributes=FALSE))
manifest<-read('bundle-manifest.csv')
stopifnot(identical(unname(tools::md5sum(file.path(root,manifest$file))),manifest$md5))
d<-readRDS(file.path(root,'dataset-intervals.rds'));p<-read('performance.csv')
stopifnot(nrow(d)==12200L,nrow(p)==122L,!anyDuplicated(d[c('profile','rep','gamma','method')]))
for(x in split(d,interaction(d$profile,d$gamma,d$method,drop=TRUE))){
 row<-p[p$profile==x$profile[1] & p$gamma==x$gamma[1] & p$method==x$method[1],]
 stopifnot(nrow(x)==100L,nrow(row)==1L,sum(x$available)==row$available,
  all(is.na(x$low[!x$available])),all(is.na(x$high[!x$available])),
  all(is.na(x$covers[!x$available])),all(is.na(x$rejects_zero[!x$available])))
 e<-x[x$available,]
 if(nrow(e))stopifnot(same(row$coverage,mean(e$low<=e$gamma & e$high>=e$gamma)),
  same(row$rejection,mean(e$low>0 | e$high<0)),same(row$median_width,median(e$high-e$low)),
  same(row$mean_se,mean(e$se)))
 stopifnot(same(row$bias,mean(x$estimate-x$gamma)),
  same(row$interval_and_cover,sum(x$covers,na.rm=TRUE)/100),
  same(row$interval_and_reject,sum(x$rejects_zero,na.rm=TRUE)/100))
}
paired<-read('paired_comparisons.csv')
for(i in seq_len(nrow(paired))){
 row<-paired[i,]
 a<-d[d$profile==row$profile & d$gamma==row$gamma & d$method==row$method_a,]
 b<-d[d$profile==row$profile & d$gamma==row$gamma & d$method==row$method_b,]
 e<-merge(a,b,by='rep',suffixes=c('_a','_b'));e<-e[e$available_a & e$available_b,]
 delta<-as.numeric(e[[paste0(row$metric,'_a')]])-as.numeric(e[[paste0(row$metric,'_b')]])
 stopifnot(nrow(e)==row$n)
 if(nrow(e))stopifnot(same(row$difference,mean(delta)),row$a_only==sum(delta==1),row$b_only==sum(delta==-1))
 if(nrow(e)>1)stopifnot(same(row$mc_se,sd(delta)/sqrt(length(delta))))
}
diag<-read('analytic_diagnostics.csv');cv<-read('covariance_decomposition.csv')
stopifnot(nrow(diag)==1200L,!any(diag$leverage & !diag$full_if),!any(diag$full_if & diag$singular),
 nrow(cv)==sum(diag$leverage),max(abs(cv$full-cv$direct-cv$measurement-cv$cross)/(1+abs(cv$full)))<1e-10)
shifts<-readRDS(file.path(root,'person-shifts.rds'));stopifnot(length(shifts)==600L)
for(r in shifts)for(g in 1:2){
 rows<-d[d$profile==r$profile & d$rep==r$rep & d$gamma==c(0,.4)[g],]
 stopifnot(all(abs(rows$estimate-r$estimate[[g]][3]*r$tau)<1e-10),length(r$ids)==r$n)
 jack<-rows[rows$method=='Delete-one Wald',];lev<-rows[rows$method=='Full-system HC3 Wald',]
 values<-r$jackknife[,,g];usable<-r$usable[,g]
 stopifnot(identical(jack$available,all(usable)),all(is.na(values[!usable,,drop=FALSE])))
 if(all(usable))stopifnot(same(jack$se^2,(r$n-1)/r$n*sum((values[,3]-mean(values[,3]))^2)*r$tau^2))
 if(lev$available){
  x<-r$leverage[[g]]
  stopifnot(identical(rownames(x),r$ids),all(dim(x)==c(r$n,3L)),same(lev$se^2,sum(x[,3]^2)*r$tau^2))
 }else stopifnot(is.null(r$leverage[[g]]))
}
api<-read('public_api_validation.csv')
stopifnot(sum(api$method=='bootstrap')==24L,sum(api$attempts[api$method=='bootstrap'])==2376L,
 sum(api$method=='jackknife')==12L,sum(api$attempts[api$method=='jackknife'])==760L,max(api$error)<1e-7)
fresh<-read('freshness.csv');stopifnot(nrow(fresh)==600L,!anyDuplicated(fresh$md5),!any(fresh$reused_old))
old<-read('previous-input-hashes.csv');stopifnot(nrow(old)==600L,!any(fresh$md5%in%old$md5))
steps<-read('derivative_sensitivity.csv');stopifnot(nrow(steps)==36L)
extreme<-read('posthoc_extreme_audit.csv')
stopifnot(nrow(extreme)==24L,all(extreme$available),all(is.finite(extreme$relative)),
 all(extreme$step%in%c(5e-5,2.5e-5)))
cat('PASS checksums; 12,200 endpoints / 122 rates; paired decisions; actual centered jackknife\n',
 'and uncentered HC3 person-shift covariances; strict failure gates; signed decomposition;\n',
 '24 public bootstraps / 12 jackknives; 600 fresh hashes; 36 sensitivity records.\n',sep='')
source(file.path(root,'drivers/check-plots.R'))
