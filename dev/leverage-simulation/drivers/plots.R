#!/usr/bin/env Rscript
# Presentation only: reuse all frozen endpoints; do not refit or select estimates.
if (!exists('here', inherits=FALSE)) {
  arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
  here<-dirname(normalizePath(sub('^--file=','',arg[[1]])))
}
out<-file.path(here,'results/summary')
d<-readRDS(file.path(out,'dataset-intervals.rds'))
s<-readRDS(file.path(out,'summary.rds'));p<-s$performance
suppressPackageStartupMessages(library(ggplot2))
methods<-c('Conditional Wald','Full IF Wald','Full-system HC3 Wald','Delete-one Wald','Person percentile')
palette<-setNames(c('#355F8D','#8956A0','#C08327','#B14F55','#187B65'),methods)
method_labels<-setNames(c('Conditional','Full IF','Full-system HC3','Delete-one','Bootstrap percentile'),methods)
profile_order<-s$profiles$profile
profile_labels<-c(low_balanced='Low, balanced',low_negative='Negative correlation',small_short='Small, short',
 low_larger='Low, larger N',high_control='High reliability',low_unequal='Unequal information')
decorate<-function(x){
 x$profile<-factor(x$profile,levels=profile_order)
 if('method'%in%names(x))x$method<-factor(x$method,levels=rev(methods))
 x
}
theme_report<-theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
 strip.text=element_text(face='bold',size=11),plot.title=element_text(face='bold',size=15),
 plot.subtitle=element_text(size=10.5),plot.caption=element_text(hjust=0,size=9),
 legend.position='bottom',legend.title=element_blank(),legend.box.margin=margin(t=4,b=4),
 axis.title.x=element_text(margin=margin(t=8,b=8)),panel.spacing=grid::unit(1.2,'lines'))
save_plot<-function(g,name,width=13,height=7.5){
 ggsave(file.path(out,paste0(name,'.png')),g,width=width,height=height,dpi=180)
 ggsave(file.path(out,paste0(name,'.pdf')),g,width=width,height=height)
}
write_table<-function(x,name)write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE)
wilson<-function(success,n){
 if(!n)return(c(NA_real_,NA_real_));z<-qnorm(.975);fraction<-success/n;den<-1+z^2/n
 center<-(fraction+z^2/(2*n))/den;half<-z*sqrt(fraction*(1-fraction)/n+z^2/(4*n^2))/den
 pmax(0,pmin(1,c(center-half,center+half)))
}

# Point estimates are identical across inference methods: show them once.
points<-d[d$method=='Conditional Wald',c('profile','rep','gamma','estimate')]
points$error<-points$estimate-points$gamma
stopifnot(nrow(points)==1200L,!anyDuplicated(points[c('profile','rep','gamma')]))
bias<-do.call(rbind,lapply(split(points,interaction(points$profile,points$gamma,drop=TRUE)),function(x){
 n<-nrow(x);mean<-mean(x$error);mc_se<-sd(x$error)/sqrt(n);half<-qt(.975,n-1)*mc_se
 data.frame(profile=x$profile[1],gamma=x$gamma[1],n=n,bias=mean,mc_se=mc_se,mc_low=mean-half,mc_high=mean+half)
}))
write_table(points,'bias_plot_points');write_table(bias,'bias_plot_summary')
points$profile<-factor(points$profile,levels=rev(profile_order));bias$profile<-factor(bias$profile,levels=rev(profile_order))
truth_labels<-c('0'='Null partial effect: truth = 0','0.4'='Nonzero partial effect: truth = 0.40')
g<-ggplot(points,aes(error,profile))+geom_vline(xintercept=0,linetype=2,color='#555555')+
 geom_point(position=position_jitter(width=0,height=.13,seed=20261031),size=1.2,alpha=.35,color='#73818B')+
 geom_segment(data=bias,aes(x=mc_low,xend=mc_high,y=profile,yend=profile),linewidth=1.15,color='#174F76')+
 geom_point(data=bias,aes(x=bias),size=3,shape=21,fill='#174F76',color='white',stroke=.4)+
 facet_wrap(~gamma,ncol=2,labeller=labeller(gamma=truth_labels))+
 scale_y_discrete(labels=profile_labels)+labs(title='Bias and sampling errors of the REPAIR point estimate',
 subtitle='All 100 fresh datasets per profile and truth; inference methods share the same point estimate',
 x='Estimation error (estimate minus true coefficient)',y=NULL,
 caption='Gray dots: individual dataset errors. Blue marker and bar: mean bias and its 95% Monte Carlo interval.\nThese bars describe uncertainty in mean bias across simulations, not a method confidence interval.')+theme_report
save_plot(g,'bias',12.5,6)

# Common eligibility across the five primary methods; retain the exact IDs.
primary<-d[d$method%in%methods,]
cohorts<-do.call(rbind,lapply(split(primary,interaction(primary$profile,primary$gamma,drop=TRUE)),function(x){
 ids<-Reduce(intersect,lapply(methods,function(method)x$rep[x$method==method & x$available]))
 data.frame(profile=rep(x$profile[1],length(ids)),gamma=rep(x$gamma[1],length(ids)),rep=ids)
}))
matched<-merge(primary,cohorts,by=c('profile','gamma','rep'))
matched_coverage<-do.call(rbind,lapply(split(matched,interaction(matched$profile,matched$gamma,matched$method,drop=TRUE)),function(x){
 n<-nrow(x);success<-sum(x$low<=x$gamma & x$high>=x$gamma);bounds<-wilson(success,n)
 data.frame(profile=x$profile[1],gamma=x$gamma[1],method=x$method[1],available=n,coverage=success/n,
  coverage_mc_low=bounds[1],coverage_mc_high=bounds[2])
}))
write_table(cohorts,'matched_plot_cohorts');write_table(matched_coverage,'matched_coverage_plot_data')
coverage_plot<-function(data,truth,title,subtitle){
 z<-decorate(data[data$gamma==truth & data$method%in%methods,])
 stopifnot(all(z$coverage_mc_low>=.5),all(z$coverage<=1))
 ggplot(z,aes(coverage,method,color=method))+
  geom_vline(xintercept=.95,linetype=2,color='#505050',linewidth=.6)+
  geom_segment(aes(x=coverage_mc_low,xend=coverage_mc_high,yend=method),linewidth=.8)+
  geom_point(size=2.8)+geom_text(aes(label=paste0('n=',available)),nudge_y=.27,hjust=1,size=2.8,color='#555555')+
  facet_wrap(~profile,ncol=3,labeller=labeller(profile=profile_labels))+
  scale_color_manual(values=palette,guide='none')+scale_y_discrete(labels=method_labels,expand=expansion(add=c(.45,.7)))+
  scale_x_continuous(limits=c(.5,1.015),breaks=seq(.5,1,.1),labels=scales::label_percent(accuracy=1))+
  labs(title=title,subtitle=subtitle,x='Coverage (displayed range: 50-100%)',y=NULL,
   caption='Dashed line: nominal 95% coverage. Bars: 95% Wilson Monte Carlo intervals. Labels give the number of intervals.\nAll plotted points and Monte Carlo bounds are shown; no uncertainty bars are clipped.')+theme_report
}
save_plot(coverage_plot(p,.4,'Coverage of the nonzero effect (truth = 0.40)',
 'Returned intervals for each method; denominators differ and are shown above the points'),'coverage')
save_plot(coverage_plot(p,0,'Coverage of the null partial effect (truth = 0)',
 'Returned intervals for each method; one minus coverage is the null rejection rate'),'coverage_null')
save_plot(coverage_plot(matched_coverage,.4,'Coverage on the same datasets (truth = 0.40)',
 'Only datasets where all five primary methods returned an interval; the selected cohort is identical within each panel'),
 'coverage_matched')

# Upper-width summaries retain every returned interval, including all maxima.
widths<-primary[primary$available & primary$gamma==.4,]
stopifnot(all(is.finite(widths$width)),all(widths$width>0))
tails<-do.call(rbind,lapply(split(widths,interaction(widths$profile,widths$method,drop=TRUE)),function(x){
 q<-unname(quantile(x$width,c(.5,.9,.95,.99),type=7))
 data.frame(profile=x$profile[1],gamma=.4,method=x$method[1],n=nrow(x),median=q[1],p90=q[2],p95=q[3],p99=q[4],maximum=max(x$width))
}))
write_table(tails,'width_tail_plot_data')
z<-decorate(tails)
long<-do.call(rbind,lapply(c('median','p95','maximum'),function(key){x<-z;x$value<-x[[key]];x$statistic<-key;x}))
g<-ggplot(z,aes(y=method,color=method))+
 geom_segment(aes(x=median,xend=maximum,yend=method),linewidth=.55,alpha=.6)+
 geom_segment(aes(x=median,xend=p95,yend=method),linewidth=1.4)+
 geom_point(data=long,aes(x=value,shape=statistic),size=2.5)+
 facet_wrap(~profile,ncol=3,labeller=labeller(profile=profile_labels))+
 scale_y_discrete(labels=method_labels)+scale_color_manual(values=palette,guide='none')+
 scale_shape_manual(values=c(median=16,p95=18,maximum=17),breaks=c('median','p95','maximum'),
  labels=c('Median','95th percentile','Largest observed interval'))+
 scale_x_log10(breaks=c(.1,.3,1,3,10,30,100),labels=scales::label_number())+
 labs(title='Interval-width tails: typical width versus the largest corrections',
  subtitle='Truth = 0.40; each method uses all of its returned intervals; logarithmic width axis',
  x='Interval width (logarithmic scale)',y=NULL,
  caption='Thick segment: median to empirical 95th percentile. Thin segment ends at the observed maximum.\nNo winsorization, clipping, or omission of wide intervals. Empirical upper quantiles use 45-100 returned intervals and are imprecise.')+
 theme_report
save_plot(g,'width_tails',13,8)

exceedance<-function(x){
 unique_values<-sort(unique(x));counts<-tabulate(match(x,unique_values),nbins=length(unique_values));n<-length(x)
 data.frame(value=c(min(x)/1.05,unique_values),exceedance=c(1,(n-cumsum(counts))/n))
}
width_curves<-do.call(rbind,lapply(split(widths,interaction(widths$profile,widths$method,drop=TRUE)),function(x){
 curve<-exceedance(x$width)
 data.frame(profile=x$profile[1],gamma=.4,method=x$method[1],n=nrow(x),curve)
}))
write_table(width_curves,'width_exceedance_plot_data')
z<-decorate(width_curves)
g<-ggplot(z,aes(value,exceedance,color=method))+
 geom_hline(yintercept=.05,linetype=2,color='#666666')+geom_step(linewidth=.8,direction='hv')+
 facet_wrap(~profile,ncol=3,labeller=labeller(profile=profile_labels))+
 scale_color_manual(values=palette,labels=method_labels,breaks=methods)+
 scale_x_log10(breaks=c(.1,.3,1,3,10,30,100),labels=scales::label_number())+
 scale_y_continuous(breaks=seq(0,.2,.05),labels=scales::label_percent(accuracy=1))+
 coord_cartesian(ylim=c(0,.2))+
 labs(title='How often are the returned intervals exceptionally wide?',
  subtitle='Truth = 0.40; upper 20% of the width distribution; logarithmic width axis',
  x='Width threshold',y='Fraction of returned intervals wider than threshold',
  caption='Steps use every returned interval and reach zero at each method\'s observed maximum. Dashed line: 5% exceedance.\nOnly the upper-width tail is displayed vertically; all observed width thresholds, including the largest, are retained.')+
 theme_report+guides(color=guide_legend(nrow=2,byrow=TRUE))
save_plot(g,'width_exceedance',13,8)

# A fair Wald-tail diagnostic uses the same selected cohort for all methods.
# Percentile CIs are not studentized Wald intervals; substitute their resample SD
# comparator and label it explicitly, rather than mislabeling percentile tails.
wald_methods<-c(methods[-5],'Bootstrap-SD Wald')
wald<-merge(d[d$method%in%wald_methods & d$gamma==.4,],cohorts,by=c('profile','gamma','rep'))
stopifnot(all(wald$available),all(is.finite(wald$se)),all(wald$se>0))
wald$absolute_z<-abs(wald$estimate-wald$gamma)/wald$se
write_table(wald[c('profile','rep','gamma','method','absolute_z')],'studentized_plot_points')
zcurves<-do.call(rbind,lapply(split(wald,interaction(wald$profile,wald$method,drop=TRUE)),function(x){
 curve<-exceedance(x$absolute_z)
 data.frame(profile=x$profile[1],gamma=.4,method=x$method[1],n=nrow(x),curve)
}))
write_table(zcurves,'studentized_exceedance_plot_data')
zcurves$profile<-factor(zcurves$profile,levels=profile_order)
normal<-data.frame(value=seq(0,max(wald$absolute_z),length.out=800))
normal$exceedance<-2*pnorm(-normal$value)
wald_palette<-setNames(unname(palette),wald_methods)
wald_labels<-setNames(c(unname(method_labels[methods[-5]]),'Bootstrap SD (Wald diagnostic)'),wald_methods)
g<-ggplot(zcurves,aes(value,exceedance,color=method))+
 geom_hline(yintercept=.05,color='#AAAAAA')+
 geom_vline(xintercept=qnorm(.975),linetype=3,color='#555555')+
 geom_step(linewidth=.8,direction='hv')+
 geom_line(data=normal,aes(value,exceedance),inherit.aes=FALSE,color='#222222',linetype=2,linewidth=.8)+
 facet_wrap(~profile,ncol=3,labeller=labeller(profile=profile_labels))+
 scale_color_manual(values=wald_palette,labels=wald_labels,breaks=wald_methods)+
 scale_x_continuous(trans=scales::transform_log1p(),breaks=c(0,1,qnorm(.975),3,5,10,15),
  labels=c('0','1','1.96','3','5','10','15'))+
 scale_y_continuous(breaks=seq(0,.25,.05),labels=scales::label_percent(accuracy=1))+
 coord_cartesian(ylim=c(0,.25))+
 labs(title='Tail calibration of symmetric Wald intervals',
  subtitle='Truth = 0.40; identical five-method eligibility cohort; upper 25% of absolute studentized errors',
  x='Absolute studentized error |estimate - truth| / SE\n(log(1 + x) axis spacing)',y='Fraction exceeding the error threshold',
  caption='Black dashed curve: standard-normal reference. Vertical line: 1.96; horizontal line: nominal 5% noncoverage.\nBootstrap SD is a Wald diagnostic, not the percentile interval. All observed studentized-error thresholds are retained.')+
 theme_report+guides(color=guide_legend(nrow=2,byrow=TRUE))
save_plot(g,'studentized_tails',13,8.5)
cat('Saved bias, available/matched coverage, and untrimmed width/studentized-tail figures (PNG/PDF) plus plotting tables.\n')
