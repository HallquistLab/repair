#!/usr/bin/env Rscript
# No repair, lme4 or research helper dependency: independently check endpoints.
arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
root <- dirname(normalizePath(sub('^--file=','',arg[[1]])))
read <- function(name) read.csv(file.path(root,name),stringsAsFactors=FALSE)
manifest <- read('bundle-manifest.csv')
stopifnot(all(file.exists(file.path(root,manifest$file))),
 identical(unname(tools::md5sum(file.path(root,manifest$file))),manifest$md5))
d <- readRDS(file.path(root,'dataset-intervals.rds'))
p <- read('performance.csv')
stopifnot(nrow(d)==17200L,nrow(p)==172L,
 !anyDuplicated(d[c('profile','rep','arm','gamma','method')]))
same <- function(a,b) isTRUE(all.equal(a,b,tolerance=1e-10,check.attributes=FALSE))
for(x in split(d,interaction(d$profile,d$arm,d$gamma,d$method,drop=TRUE))){
 stopifnot(nrow(x)==100L)
 row <- p[p$profile==x$profile[1] & p$arm==x$arm[1] &
          p$gamma==x$gamma[1] & p$method==x$method[1],]
 stopifnot(nrow(row)==1L, row$planned==100L,sum(x$available)==row$available,
  all(is.na(x$low[!x$available])),all(is.na(x$high[!x$available])),
  all(is.na(x$covers[!x$available])),all(is.na(x$rejects_zero[!x$available])))
 e <- x[x$available,]
 stopifnot(all(is.finite(as.matrix(e[,c('estimate','se','low','high')]))),
  all(e$low<=e$high),all(e$se>=0),
  same(e$covers,e$low<=e$gamma & e$high>=e$gamma),
  same(e$rejects_zero,e$low>0 | e$high<0),
  same(row$coverage,mean(e$low<=e$gamma & e$high>=e$gamma)),
  same(row$rejection,mean(e$low>0 | e$high<0)),
  same(row$bias,mean(x$estimate-x$gamma,na.rm=TRUE)),
  same(row$interval_and_cover,sum(x$covers,na.rm=TRUE)/100),
  same(row$interval_and_reject,sum(x$rejects_zero,na.rm=TRUE)/100),
  same(row$median_width,median(e$high-e$low)))
}
paired <- read('paired_comparisons.csv')
for(j in seq_len(nrow(paired))){
 row <- paired[j,]
 a <- d[d$profile==row$profile & d$arm==row$arm & d$gamma==row$gamma & d$method==row$method_a,]
 b <- d[d$profile==row$profile & d$arm==row$arm & d$gamma==row$gamma & d$method==row$method_b,]
 e <- merge(a,b,by='rep',suffixes=c('_a','_b')); e <- e[e$available_a & e$available_b,]
 delta <- as.numeric(e[[paste0(row$metric,'_a')]])-as.numeric(e[[paste0(row$metric,'_b')]])
 stopifnot(nrow(e)==row$n,same(row$difference,mean(delta)),
  same(row$mc_se,sd(delta)/sqrt(length(delta))),row$a_only==sum(delta==1),row$b_only==sum(delta==-1))
}
dg <- read('analytic_diagnostics.csv'); cv <- read('covariance_decomposition.csv')
stopifnot(nrow(dg)==2400L,sum(dg$full)==1351L,sum(dg$stage2)==2400L,
 sum(dg$measurement_singular)==1048L,sum(!dg$full & !dg$measurement_singular)==1L,
 !any(dg$full & dg$measurement_singular),nrow(cv)==1351L,
 max(abs(cv$full-cv$stage2-cv$measurement-cv$cross)/(1+abs(cv$full)))<1e-10)
a <- read('balanced_ols_validation.csv'); v <- read('public_bootstrap_validation.csv')
stopifnot(nrow(a)==1000L,max(a$point_gap)<1e-7,max(a$if_legacy_variance_gap)<1e-6,
 nrow(v)==24L,sum(v$B)==2376L,max(v$error)<1e-7)
st <- read('derivative_sensitivity.csv'); bs <- read('bootstrap_sensitivity.csv')
stopifnot(nrow(st)==72L,all(st$full_same),all(st$stage2_same),
 max(st$full_relative,na.rm=TRUE)<2e-7,max(st$stage2_relative,na.rm=TRUE)<2e-7,
 nrow(bs)==156L,all(bs$available_499),all(bs$available_1999),
 sum(!bs$covers_499 & bs$covers_1999)==1L,sum(bs$covers_499 & !bs$covers_1999)==0L)
g <- read('guarded_api_validation.csv')
stopifnot(nrow(g)==72L,max(g$stage2_relative)<1e-10,max(g$full_relative,na.rm=TRUE)<1e-10)
cat('PASS bundle checksums; 17,200 endpoints / 172 rates; all matched decision contrasts;\n',
 'availability and signed decomposition; 1,000 OLS anchors; 24 public bootstraps;\n',
 '72 derivative and guarded-API checks; 156 B-sensitivity records.\n',sep='')
