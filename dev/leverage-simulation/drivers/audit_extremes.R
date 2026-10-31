#!/usr/bin/env Rscript
# POST HOC numerical audit, specified after seeing the wide interval tails.
# Select the 12 largest reported HC3/HC0 focal variance ratios. Compare the
# existing result with h/2 and h/4. Never replace or select the main estimates.
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE);here<-dirname(normalizePath(sub('^--file=','',arg[[1]])))
.libPaths(c('/private/tmp/repair-leverage-lib',.libPaths()));library(repair)
source(file.path(here,'engine.R'));load_sources(normalizePath(file.path(here,'../../../..')))
ratios<-read.csv(file.path(here,'results/summary/covariance_decomposition.csv'))
selected<-head(ratios[order(ratios$to_hc0,decreasing=TRUE),],12L)
profiles<-make_profiles(normalizePath(file.path(here,'../../../..')))
rows<-list()
for(i in seq_len(nrow(selected))){
 x<-selected[i,];index<-match(x$profile,profiles$profile);name<-sprintf('p%02d_r%04d.rds',index,x$rep)
 native<-readRDS(file.path(here,'results/main/native',name));record<-readRDS(file.path(here,'results/main',name))
 sim<-native$data;g<-match(x$gamma,c(0,.4));people<-sim$people;people$outcome<-sim$outcomes[,g]
 m<-native_measurement(sim$repeated)$model
 f<-repair(outcome~re('(Intercept)','growth')+re('time','growth'),list(growth=m),people,'person',backend='rcpp')
 baseline<-record$analytic[[g]]$leverage
 for(step in c(5e-5,2.5e-5)){
  smaller<-repair_leverage(f,step=step)
  relative<-if(smaller$diagnostics$full_available)
    max(abs(smaller$covariance$full-baseline$covariance$full))/(1e-8+max(abs(baseline$covariance$full))) else NA_real_
  rows[[length(rows)+1L]]<-data.frame(profile=x$profile,rep=x$rep,gamma=x$gamma,variance_ratio=x$to_hc0,
   step=step,available=smaller$diagnostics$full_available,relative=relative,reason=smaller$diagnostics$full_reason)
 }
}
rows<-do.call(rbind,rows)
write.csv(rows,file.path(here,'results/summary/posthoc_extreme_audit.csv'),row.names=FALSE)
cat(nrow(rows),'post hoc comparisons; unavailable:',sum(!rows$available),'; max matched relative covariance change:',max(rows$relative,na.rm=TRUE),'\n')
