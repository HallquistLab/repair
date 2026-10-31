#!/usr/bin/env Rscript
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE);here<-dirname(normalizePath(sub('^--file=','',arg[[1]])))
read_phase<-function(phase,N){
 files<-list.files(file.path(here,'results',phase),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
 stopifnot(length(files)==N);x<-lapply(files,readRDS);names(x)<-basename(files)
 errors<-Filter(function(r)!is.null(r$error),x)
 if(length(errors))stop(paste(vapply(errors,function(r)paste(r$job$profile,r$job$rep,r$stage,r$error),character(1)),collapse='\n'))
 x
}
main<-read_phase('main',600L);sens<-read_phase('sensitivity',18L)
out<-file.path(here,'results/summary');dir.create(out,showWarnings=FALSE)
d<-do.call(rbind,lapply(main,`[[`,'rows'));rownames(d)<-NULL
d$covers<-ifelse(d$available,d$low<=d$gamma & d$high>=d$gamma,NA)
d$rejects_zero<-ifelse(d$available,d$low>0 | d$high<0,NA);d$width<-d$high-d$low
saveRDS(d,file.path(out,'dataset-intervals.rds'),compress='xz')
profiles<-do.call(rbind,lapply(main,`[[`,'profile'));profiles<-profiles[!duplicated(profiles$profile),]
write.csv(profiles,file.path(out,'profiles.csv'),row.names=FALSE)
wilson<-function(p,n){if(!n)return(c(NA_real_,NA_real_));z<-qnorm(.975);den<-1+z*z/n
 c((p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/den,(p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/den)}
performance<-do.call(rbind,lapply(split(d,interaction(d$profile,d$gamma,d$method,drop=TRUE)),function(x){
 e<-x[x$available,];n<-nrow(e);coverage<-if(n)mean(e$covers) else NA_real_;rejection<-if(n)mean(e$rejects_zero) else NA_real_
 ci<-wilson(coverage,n);ri<-wilson(rejection,n)
 data.frame(profile=x$profile[1],gamma=x$gamma[1],method=x$method[1],planned=nrow(x),available=n,availability=n/nrow(x),
 coverage=coverage,coverage_mc_low=ci[1],coverage_mc_high=ci[2],rejection=rejection,rejection_mc_low=ri[1],rejection_mc_high=ri[2],
 interval_and_cover=sum(x$covers,na.rm=TRUE)/nrow(x),interval_and_reject=sum(x$rejects_zero,na.rm=TRUE)/nrow(x),
 bias=mean(x$estimate-x$gamma[1],na.rm=TRUE),empirical_sd=sd(x$estimate,na.rm=TRUE),
 eligible_bias=if(n)mean(e$estimate-e$gamma) else NA_real_,eligible_sd=if(n>1)sd(e$estimate) else NA_real_,
 mean_se=if(n)mean(e$se) else NA_real_,median_width=if(n)median(e$width) else NA_real_,
 p95_width=if(n)unname(quantile(e$width,.95)) else NA_real_)
}))
rownames(performance)<-NULL;write.csv(performance,file.path(out,'performance.csv'),row.names=FALSE)
pairs<-list(c('Full-system HC3 Wald','Full IF Wald'),c('Full-system HC3 Wald','Person percentile'),
 c('Delete-one Wald','Person percentile'),c('Full-system HC3 Wald','Delete-one Wald'),
 c('Full-system HC3 t','Full-system HC3 Wald'),c('Delete-one t','Delete-one Wald'),
 c('Person percentile','Bootstrap-SD Wald'),c('Full IF Wald','Conditional Wald'))
paired<-list()
for(profile in profiles$profile)for(gamma in c(0,.4))for(pair in pairs){
 a<-d[d$profile==profile & d$gamma==gamma & d$method==pair[1],];b<-d[d$profile==profile & d$gamma==gamma & d$method==pair[2],]
 e<-merge(a,b,by='rep',suffixes=c('_a','_b'));e<-e[e$available_a & e$available_b,]
 for(metric in c('covers','rejects_zero')){
 delta<-as.numeric(e[[paste0(metric,'_a')]])-as.numeric(e[[paste0(metric,'_b')]]);n<-length(delta)
 difference<-if(n)mean(delta) else NA_real_;mc<-if(n>1)sd(delta)/sqrt(n) else NA_real_;zero<-n>0&&all(delta==0)
 paired[[length(paired)+1L]]<-data.frame(profile=profile,gamma=gamma,method_a=pair[1],method_b=pair[2],metric=metric,n=n,
 difference=difference,mc_se=mc,mc_low=if(!zero)difference-qnorm(.975)*mc else NA_real_,mc_high=if(!zero)difference+qnorm(.975)*mc else NA_real_,
 a_only=sum(delta==1),b_only=sum(delta==-1),zero_discordance=zero)
 }
}
paired<-do.call(rbind,paired);write.csv(paired,file.path(out,'paired_comparisons.csv'),row.names=FALSE)
diagnostics<-components<-deletions<-api<-branch_rows<-deletion_statistics<-list()
for(r in main)for(g in 1:2){
 inf<-r$analytic[[g]]$influence;lev<-r$analytic[[g]]$leverage;dg<-lev$diagnostics
 diagnostic<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
  singular=r$measurement$singular,branch=inf$diagnostics$branch,full_if=inf$diagnostics$full_available,
  if_reason=inf$diagnostics$full_reason,leverage=dg$full_available,leverage_reason=dg$full_reason,
  delete_one=all(r$deletions$usable[,g]),failed_deletions=sum(!r$deletions$usable[,g]),
  hc0_gap=if(is.null(dg$hc0_relative_gap))NA_real_ else dg$hc0_relative_gap,
  derivative_error=if(is.null(dg$leverage_derivative_error))NA_real_ else dg$leverage_derivative_error,
  min_rcond=if(is.null(dg$deletions)||!any(dg$deletions$usable))NA_real_ else min(dg$deletions$rcond[dg$deletions$usable]))
 diagnostics[[length(diagnostics)+1L]]<-diagnostic
 original_floor<-grepl('forced_zero|clipped_zero',inf$diagnostics$branch)
 original_final<-tail(strsplit(inf$diagnostics$branch,':',fixed=TRUE)[[1]],1)
 usable<-r$deletions$usable[,g];n<-r$profile$n
 delete_alpha<-abs(r$deletions$correction[,g]-(1-4/(n-1)))<1e-10
 deletion_statistics[[length(deletion_statistics)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
  n=n,usable=sum(usable),original_floor=original_floor,singular_deletions=sum(r$deletions$singular,na.rm=TRUE),
  floor_switches=sum((r$deletions$residual[,g]==0)!=original_floor,na.rm=TRUE),
  final_branch_switches=if(original_final%in%c('alpha','root'))sum(delete_alpha!=(original_final=='alpha'),na.rm=TRUE) else NA_integer_,
  jackknife_bias=if(all(usable))(n-1)*(mean(r$deletions$estimates[,3,g])-inf$estimate[3])*r$profile$tau else NA_real_)
 branch_rows[[length(branch_rows)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
  id=inf$ids,usable=r$deletions$usable[,g],singular=r$deletions$singular,
  floor_zero=ifelse(r$deletions$usable[,g],r$deletions$residual[,g]==0,NA),
  correction=r$deletions$correction[,g],reason=r$deletions$reason[,g])
 if(dg$full_available){
  V<-vapply(lev$covariance,function(x)x[3,3]*r$profile$tau^2,numeric(1))
  components[[length(components)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
   full=V['full'],direct=V['stage2'],measurement=V['measurement'],cross=V['cross'],
   to_hc0=V['full']/(inf$covariance$full[3,3]*r$profile$tau^2))
  point<-lev$estimate;exact<-sweep(r$deletions$estimates[,,g],2,unname(point))
  if(all(r$deletions$usable[,g])){
   delta<-lev$deletions-exact
   deletions[[length(deletions)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
    singular_deletion=any(r$deletions$singular,na.rm=TRUE),
    floor_switch=any((r$deletions$residual[,g]==0)!=original_floor,na.rm=TRUE),
    final_branch_switch=if(original_final%in%c('alpha','root'))any(delete_alpha!=(original_final=='alpha'),na.rm=TRUE) else NA,
    relative_shift_rmse=sqrt(sum(delta^2)/(1e-8+sum(exact^2))),
    focal_relative_shift_rmse=sqrt(sum(delta[,3]^2)/(1e-8+sum(exact[,3]^2))),
    leverage_variance=V['full'],jackknife_variance=(r$profile$n-1)/r$profile$n*sum((exact[,3]-mean(exact[,3]))^2)*r$profile$tau^2)
  }
 }
 b<-if(length(r$bootstrap_checks)>=g)r$bootstrap_checks[[g]] else NULL
 if(!is.null(b))api[[length(api)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
  method='bootstrap',attempts=b$B,error=b$max_error)
 j<-r$checks[[paste0('jackknife_',g)]]
 if(!is.null(j))api[[length(api)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
  method='jackknife',attempts=j$n,error=j$error)
}
bind<-function(x)if(length(x))do.call(rbind,x) else data.frame()
diagnostics<-bind(diagnostics);components<-bind(components);deletions<-bind(deletions);api<-bind(api);branch_rows<-bind(branch_rows)
deletion_statistics<-bind(deletion_statistics)
write.csv(diagnostics,file.path(out,'analytic_diagnostics.csv'),row.names=FALSE)
write.csv(components,file.path(out,'covariance_decomposition.csv'),row.names=FALSE)
write.csv(deletions,file.path(out,'leverage_exact_comparison.csv'),row.names=FALSE)
write.csv(api,file.path(out,'public_api_validation.csv'),row.names=FALSE)
saveRDS(branch_rows,file.path(out,'deletion_diagnostics.rds'),compress='xz')
write.csv(deletion_statistics,file.path(out,'deletion_statistics.csv'),row.names=FALSE)
strat<-merge(d,diagnostics[c('profile','rep','gamma','branch','singular')],by=c('profile','rep','gamma'))
strat$state<-ifelse(strat$singular,'Singular measurement',ifelse(strat$branch=='alpha:positive:alpha','Interior structural branch','Active correction/floor'))
branch_performance<-bind(lapply(split(strat,interaction(strat$profile,strat$gamma,strat$method,strat$state,drop=TRUE)),function(x){
 e<-x[x$available,];data.frame(profile=x$profile[1],gamma=x$gamma[1],method=x$method[1],state=x$state[1],planned=nrow(x),available=nrow(e),
 coverage=if(nrow(e))mean(e$covers) else NA_real_,rejection=if(nrow(e))mean(e$rejects_zero) else NA_real_)
}))
write.csv(branch_performance,file.path(out,'branch_performance.csv'),row.names=FALSE)
steps<-tails<-list()
for(r in sens){
 old<-main[[sprintf('p%02d_r%04d.rds',r$job$profile,r$job$rep)]]
 for(g in 1:2){a<-r$checks[[paste0('step_',g)]];steps[[length(steps)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],same=a$same,relative=a$relative)}
 e<-merge(r$rows,old$rows,by=c('profile','rep','gamma','method'),suffixes=c('_1999','_499'))
 tails[[length(tails)+1L]]<-e[e$method%in%c('Person percentile','Frozen percentile','Stratified percentile'),
 c('profile','rep','gamma','method','available_499','available_1999','low_499','low_1999','high_499','high_1999')]
}
steps<-bind(steps);tails<-bind(tails)
write.csv(steps,file.path(out,'derivative_sensitivity.csv'),row.names=FALSE)
write.csv(tails,file.path(out,'bootstrap_sensitivity.csv'),row.names=FALSE)

suppressPackageStartupMessages(library(ggplot2))
methods<-c('Conditional Wald','Full IF Wald','Full-system HC3 Wald','Delete-one Wald','Person percentile')
palette<-setNames(c('#355F8D','#8956A0','#C08327','#B14F55','#187B65'),methods)
z<-performance[performance$method%in%methods,];z$profile<-factor(z$profile,levels=profiles$profile);z$method<-factor(z$method,levels=methods)
labels<-c(low_balanced='Low, balanced',low_negative='Negative correlation',small_short='Small, short',low_larger='Low, larger N',high_control='High reliability',low_unequal='Unequal information')
style<-theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position='bottom',legend.title=element_blank(),axis.title.x=element_blank(),axis.text.x=element_text(angle=20,hjust=1))
g<-ggplot(subset(z,gamma==.4),aes(profile,coverage,color=method))+geom_hline(yintercept=.95,linetype=2,color='grey50')+
 geom_errorbar(aes(ymin=coverage_mc_low,ymax=coverage_mc_high),position=position_dodge(.7),width=.08)+geom_point(position=position_dodge(.7),size=2.2)+
 scale_color_manual(values=palette)+scale_x_discrete(labels=labels)+scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(100*x),'%'))+
 labs(title='Fresh noisy-predictor simulations: nonzero-effect coverage',subtitle='Available intervals only; bars are 95% Monte Carlo intervals',y='Coverage')+style+guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave(file.path(out,'coverage.png'),g,width=12,height=6,dpi=160);ggsave(file.path(out,'coverage.pdf'),g,width=12,height=6)
g<-ggplot(z,aes(profile,availability,fill=method))+geom_col(position=position_dodge(.8),width=.7)+facet_wrap(~gamma)+
 scale_fill_manual(values=palette)+scale_x_discrete(labels=labels)+scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(100*x),'%'))+
 labs(title='Returned intervals / 100 planned fresh datasets',subtitle='Failed analytic sensitivities or complete-vector fits remain unavailable',y='Availability')+style+guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(file.path(out,'availability.png'),g,width=14,height=6,dpi=160);ggsave(file.path(out,'availability.pdf'),g,width=14,height=6)
saveRDS(list(performance=performance,paired=paired,diagnostics=diagnostics,components=components,exact=deletions,api=api,
 deletion_statistics=deletion_statistics,branch_performance=branch_performance,steps=steps,tails=tails,profiles=profiles),file.path(out,'summary.rds'))
cat('Wrote',nrow(d),'endpoints /',nrow(performance),'performance cells.\n')
source(file.path(here,'plots.R'))
