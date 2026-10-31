#!/usr/bin/env Rscript
# Check plotting inputs against the retained endpoints, using base R only.
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE);here<-dirname(normalizePath(sub('^--file=','',arg[[1]])))
study_output<-file.exists(file.path(here,'results/summary/dataset-intervals.rds'))
out<-if(study_output)file.path(here,'results/summary') else
 if(file.exists(file.path(here,'dataset-intervals.rds')))here else normalizePath(file.path(here,'..'))
read<-function(name)read.csv(file.path(out,name),stringsAsFactors=FALSE)
d<-readRDS(file.path(out,'dataset-intervals.rds'));p<-read('performance.csv')
methods<-c('Conditional Wald','Full IF Wald','Full-system HC3 Wald','Delete-one Wald','Person percentile')
same<-function(a,b)isTRUE(all.equal(a,b,tolerance=1e-10,check.attributes=FALSE))
bias<-read('bias_plot_summary.csv');points<-read('bias_plot_points.csv')
stopifnot(nrow(points)==1200L,nrow(bias)==12L,all(bias$n==100L),same(points$error,points$estimate-points$gamma))
for(i in seq_len(nrow(bias))){
 b<-bias[i,];x<-subset(d,profile==b$profile & gamma==b$gamma & method=='Conditional Wald')
 err<-x$estimate-x$gamma;half<-qt(.975,99)*sd(err)/10
 stopifnot(same(b$bias,mean(err)),same(c(b$mc_low,b$mc_high),mean(err)+c(-half,half)))
}
cohorts<-read('matched_plot_cohorts.csv');matched<-read('matched_coverage_plot_data.csv')
stopifnot(!anyDuplicated(cohorts[c('profile','gamma','rep')]),nrow(matched)==60L)
for(x in split(d[d$method%in%methods,],interaction(d$profile[d$method%in%methods],d$gamma[d$method%in%methods],drop=TRUE))){
 ids<-Reduce(intersect,lapply(methods,function(method)x$rep[x$method==method & x$available]))
 c<-cohorts[cohorts$profile==x$profile[1] & cohorts$gamma==x$gamma[1],]
 stopifnot(setequal(ids,c$rep))
 for(method in methods){
  e<-x[x$rep%in%ids & x$method==method,];r<-matched[matched$profile==x$profile[1] & matched$gamma==x$gamma[1] & matched$method==method,]
  stopifnot(nrow(r)==1L,r$available==length(ids),same(r$coverage,mean(e$low<=e$gamma & e$high>=e$gamma)))
  # R's score interval is an independent check of the plotted Wilson bounds.
  # Only the interval is used; edge-rate chi-square p-value warnings are irrelevant.
  bounds<-suppressWarnings(stats::prop.test(sum(e$covers),nrow(e),correct=FALSE)$conf.int)
  stopifnot(same(c(r$coverage_mc_low,r$coverage_mc_high),unname(bounds)))
 }
}
# Every coverage interval fits the 50-100% view; the nominal line is 95%.
stopifnot(all(subset(p,method%in%methods)$coverage_mc_low>=.5),all(matched$coverage_mc_low>=.5),
 all(subset(p,method%in%methods)$coverage_mc_high<=1+1e-12),all(matched$coverage_mc_high<=1+1e-12))
tails<-read('width_tail_plot_data.csv');curves<-read('width_exceedance_plot_data.csv')
stopifnot(nrow(tails)==30L)
for(i in seq_len(nrow(tails))){
 t<-tails[i,];x<-d[d$profile==t$profile & d$gamma==.4 & d$method==t$method & d$available,]
 stopifnot(t$n==nrow(x),same(unlist(t[c('median','p90','p95','p99')]),quantile(x$high-x$low,c(.5,.9,.95,.99),names=FALSE)),
  same(t$maximum,max(x$width)))
 c<-curves[curves$profile==t$profile & curves$method==t$method,]
 thresholds<-c(c$value[1],sort(unique(x$width)))
 stopifnot(all(c$n==nrow(x)),same(max(c$value),max(x$width)),tail(c$exceedance,1)==0,
  same(c$value,thresholds),
  all(vapply(seq_len(nrow(c)),function(j)same(c$exceedance[j],mean(x$width>thresholds[j])),logical(1))))
}
wald<-read('studentized_plot_points.csv');zcurves<-read('studentized_exceedance_plot_data.csv')
wald_methods<-c(methods[-5],'Bootstrap-SD Wald')
stopifnot(setequal(unique(wald$method),wald_methods))
for(x in split(wald,interaction(wald$profile,wald$method,drop=TRUE))){
 ids<-subset(cohorts,profile==x$profile[1] & gamma==.4)$rep
 stopifnot(setequal(x$rep,ids))
 e<-d[d$profile==x$profile[1] & d$gamma==.4 & d$method==x$method[1] & d$rep%in%ids,]
 e<-e[match(x$rep,e$rep),]
 real_z<-abs(e$estimate-e$gamma)/e$se
 stopifnot(same(x$absolute_z,real_z),same(mean(real_z>qnorm(.975)),mean(!e$covers)))
 c<-zcurves[zcurves$profile==x$profile[1] & zcurves$method==x$method[1],]
 thresholds<-c(c$value[1],sort(unique(real_z)))
 stopifnot(all(c$n==nrow(x)),same(max(c$value),max(real_z)),tail(c$exceedance,1)==0,same(c$value,thresholds),
  all(vapply(seq_len(nrow(c)),function(j)same(c$exceedance[j],mean(real_z>thresholds[j])),logical(1))))
}
text<-c('PASS 1,200 original point estimates / 12 mean-bias MC intervals; one shared point estimator.',
 'PASS all-method intersection IDs / 60 matched coverage cells; 50-100% coverage views contain every MC bound.',
 'PASS 30 full returned-set width quantile/max summaries and exact exceedance curves; no extreme intervals omitted.',
 'PASS matched-cohort studentized tails, including exact 1.96 noncoverage identities for normal Wald intervals.')
if(study_output)writeLines(text,file.path(out,'plot-validation.txt'))
cat(paste(text,collapse='\n'),'\n')
