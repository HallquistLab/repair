#!/usr/bin/env Rscript
file_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
study_dir<-dirname(normalizePath(sub('^--file=','',file_arg[[1]])))
read_phase<-function(phase,expected){
 f<-list.files(file.path(study_dir,'results',phase),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
 if(length(f)!=expected)stop('Expected ',expected,' ',phase,' files, found ',length(f))
 x<-lapply(f,readRDS);names(x)<-basename(f)
 errors<-Filter(function(r)!is.null(r$error),x)
 if(length(errors))stop('Driver/validation errors: ',paste(vapply(errors,function(r)paste(r$source,r$stage,r$error),character(1)),collapse=' | '))
 x
}
main<-read_phase('main',600L);sensitivity<-read_phase('sensitivity',18L)
out<-file.path(study_dir,'results/summary');dir.create(out,showWarnings=FALSE)
d<-do.call(rbind,lapply(main,`[[`,'rows'));rownames(d)<-NULL
write.csv(d,file.path(out,'dataset_metrics.csv'),row.names=FALSE)
profiles<-do.call(rbind,lapply(main,function(r)r$profile));profiles<-profiles[!duplicated(profiles$profile),]
write.csv(profiles,file.path(out,'profiles.csv'),row.names=FALSE)
wilson<-function(p,n){if(!n)return(c(NA_real_,NA_real_));z<-qnorm(.975);den<-1+z^2/n
 c((p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/den,(p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/den)}
performance<-do.call(rbind,lapply(split(d,interaction(d$profile,d$arm,d$gamma,d$method,drop=TRUE)),function(x){
 e<-x[x$available,];k<-nrow(e);N<-nrow(x);cov<-if(k)mean(e$covers) else NA_real_;rej<-if(k)mean(e$rejects_zero) else NA_real_
 ci<-wilson(cov,k);ri<-wilson(rej,k)
 data.frame(profile=x$profile[1],arm=x$arm[1],gamma=x$gamma[1],method=x$method[1],planned=N,
  point_n=sum(is.finite(x$estimate)),available=k,availability=k/N,coverage=cov,coverage_mc_low=ci[1],coverage_mc_high=ci[2],
  rejection=rej,rejection_mc_low=ri[1],rejection_mc_high=ri[2],
  interval_and_cover=sum(x$available & x$covers,na.rm=TRUE)/N,interval_and_reject=sum(x$available & x$rejects_zero,na.rm=TRUE)/N,
  bias=mean(x$estimate-x$gamma,na.rm=TRUE),empirical_sd=sd(x$estimate,na.rm=TRUE),
  eligible_bias=if(k)mean(e$estimate-e$gamma) else NA_real_,eligible_sd=if(k>1)sd(e$estimate) else NA_real_,
  mean_se=if(k)mean(e$se) else NA_real_,se_to_empirical_sd=if(k>1)mean(e$se)/sd(e$estimate) else NA_real_,
  median_width=if(k)median(e$width) else NA_real_,p95_width=if(k)unname(quantile(e$width,.95)) else NA_real_)
}))
rownames(performance)<-NULL;write.csv(performance,file.path(out,'performance.csv'),row.names=FALSE)
pairs<-list(c('Full influence Wald','Conditional M-MOM Wald'),c('Full influence Wald','Stage2 influence Wald'),
 c('Full influence Wald','Full person percentile'),c('Stage2 influence Wald','Frozen measurement percentile'),
 c('Full influence Wald','Legacy OLS stacked HC0 (t)'),c('Full influence Wald','Legacy OLS stacked HC3 (t)'),
 c('Full person percentile','Legacy OLS stacked HC3 (t)'))
paired<-list()
for(profile in profiles$profile)for(arm in c('predictor','outcome'))for(gamma in c(0,.4))for(pair in pairs){
 a<-d[d$profile==profile & d$arm==arm & d$gamma==gamma & d$method==pair[1],]
 b<-d[d$profile==profile & d$arm==arm & d$gamma==gamma & d$method==pair[2],]
 if(!nrow(a)||!nrow(b))next
 e<-merge(a,b,by='rep',suffixes=c('_a','_b'));e<-e[e$available_a & e$available_b,]
 for(metric in c('covers','rejects_zero')){
  delta<-as.numeric(e[[paste0(metric,'_a')]])-as.numeric(e[[paste0(metric,'_b')]])
  n<-length(delta);change<-if(n)mean(delta) else NA_real_;se<-if(n>1)sd(delta)/sqrt(n) else NA_real_
  zero<-n>0&&all(delta==0)
  paired[[length(paired)+1L]]<-data.frame(profile=profile,arm=arm,gamma=gamma,method_a=pair[1],method_b=pair[2],metric=metric,n=n,
   difference=change,mc_se=se,mc_low=if(!zero)change-qnorm(.975)*se else NA_real_,mc_high=if(!zero)change+qnorm(.975)*se else NA_real_,
   a_only=sum(delta==1),b_only=sum(delta==-1),zero_discordance=zero,
   mean_abs_point_difference=if(n)mean(abs(e$estimate_a-e$estimate_b)) else NA_real_,
   mean_se_ratio=if(n)mean(e$se_a/e$se_b) else NA_real_)
 }
}
paired<-do.call(rbind,paired);write.csv(paired,file.path(out,'paired_comparisons.csv'),row.names=FALSE)

diag_rows<-components<-checks<-list()
for(r in main)for(arm in c('predictor','outcome'))for(g in 1:2){
 inf<-r$influences[[paste(arm,g,sep='_')]];dg<-inf$diagnostics;idx<-if(arm=='predictor')3 else 2
 diag_rows[[length(diag_rows)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,arm=arm,gamma=c(0,.4)[g],
   measurement_singular=r$measurement$singular,full=dg$full_available,stage2=dg$stage2_available,
   reason=dg$full_reason,stage2_reason=dg$stage2_reason,branch=dg$branch,
   frequency_error=if(is.null(dg$frequency_derivative_error))NA_real_ else dg$frequency_derivative_error,
   parameter_error=if(is.null(dg$measurement_derivative_error))NA_real_ else dg$measurement_derivative_error)
 if(dg$full_available){
  f<-if(arm=='predictor')r$profile$tau^2 else 1/r$profile$tau^2
  V<-vapply(inf$covariance,function(x)x[idx,idx]*f,numeric(1))
  components[[length(components)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,arm=arm,gamma=c(0,.4)[g],
    full=V['full'],stage2=V['stage2'],measurement=V['measurement'],cross=V['cross'],
    full_to_stage2=V['full']/V['stage2'],measurement_fraction=V['measurement']/V['full'],cross_fraction=V['cross']/V['full'],
    theta_sensitivity=max(abs(inf$measurement$structural_sensitivity[idx,grep('^theta:',names(inf$measurement$psi))])))
 }
 if(arm=='outcome'){
  legacy<-r$checks[[paste(g,'legacy',sep='_')]];v<-r$checks[[paste(g,'public',sep='_')]]
  if(!is.null(v))checks[[length(checks)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,gamma=c(0,.4)[g],
    B=v$B,error=if(is.null(v$error))NA_real_ else v$error)
 }
}
diag_rows<-do.call(rbind,diag_rows);write.csv(diag_rows,file.path(out,'analytic_diagnostics.csv'),row.names=FALSE)
components<-do.call(rbind,components);rownames(components)<-NULL;write.csv(components,file.path(out,'covariance_decomposition.csv'),row.names=FALSE)
checks<-do.call(rbind,checks);write.csv(checks,file.path(out,'public_bootstrap_validation.csv'),row.names=FALSE)

# Separately report legacy t-vs-normal construction without calling that an SE
# improvement. This secondary diagnostic uses the same old estimates and SEs.
legacy_normal<-d[grepl('^Legacy',d$method),]
legacy_normal$low<-legacy_normal$estimate-qnorm(.975)*legacy_normal$se
legacy_normal$high<-legacy_normal$estimate+qnorm(.975)*legacy_normal$se
legacy_normal$covers<-ifelse(legacy_normal$available,legacy_normal$low<=legacy_normal$gamma & legacy_normal$high>=legacy_normal$gamma,NA)
legacy_normal$rejects_zero<-ifelse(legacy_normal$available,legacy_normal$low>0|legacy_normal$high<0,NA)
normal<-do.call(rbind,lapply(split(legacy_normal,interaction(legacy_normal$profile,legacy_normal$gamma,legacy_normal$method,drop=TRUE)),function(x){
 e<-x[x$available,];data.frame(profile=x$profile[1],gamma=x$gamma[1],method=x$method[1],n=nrow(e),
  normal_coverage=if(nrow(e))mean(e$covers) else NA_real_,normal_rejection=if(nrow(e))mean(e$rejects_zero) else NA_real_)
}))
write.csv(normal,file.path(out,'legacy_normal_diagnostic.csv'),row.names=FALSE)

sensitivity_rows<-step_rows<-list()
for(r in sensitivity){
 old<-main[[sprintf('p%02d_r%04d.rds',r$job$profile,r$job$rep)]]
 # Same original fit/IF calculation; sensitivity bootstraps extend existing draws.
 for(arm in c('predictor','outcome'))for(g in 1:2){
  key<-paste(arm,g,'step',sep='_');st<-r$checks[[key]]
  step_rows[[length(step_rows)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,arm=arm,gamma=c(0,.4)[g],
    full_same=st$full_same,stage2_same=st$stage2_same,full_relative=st$full_relative,stage2_relative=st$stage2_relative)
 }
 a<-r$rows;b<-old$rows;e<-merge(a,b,by=c('arm','gamma','method'),suffixes=c('_1999','_499'))
 e<-e[grepl('percentile',e$method),]
 sensitivity_rows[[length(sensitivity_rows)+1L]]<-data.frame(profile=r$profile$profile,rep=r$job$rep,e[,c('arm','gamma','method',
  'available_499','available_1999','covers_499','covers_1999','low_499','low_1999','high_499','high_1999')])
}
sensitivity_rows<-do.call(rbind,sensitivity_rows);step_rows<-do.call(rbind,step_rows)
write.csv(sensitivity_rows,file.path(out,'bootstrap_sensitivity.csv'),row.names=FALSE)
write.csv(step_rows,file.path(out,'derivative_sensitivity.csv'),row.names=FALSE)

suppressPackageStartupMessages(library(ggplot2))
labels<-c(low_balanced='Low, balanced',low_negative='Negative correlation',small_short='Small, short',low_larger='Low, larger N',high_control='High reliability',low_unequal='Unequal information')
methods<-c('Conditional M-MOM Wald','Stage2 influence Wald','Full influence Wald','Full person percentile','Legacy OLS stacked HC0 (t)','Legacy OLS stacked HC3 (t)')
palette<-setNames(c('#355F8D','#C18B2D','#8B5AA4','#167D65','#788895','#BB5A60'),methods)
plot_data<-performance[performance$method%in%methods,]
plot_data$profile<-factor(plot_data$profile,levels=profiles$profile);plot_data$method<-factor(plot_data$method,levels=methods)
plot_data$arm<-factor(plot_data$arm,levels=c('predictor','outcome'),labels=c('A: Noisy predictors','B: Random-slope outcome'))
plot_data$effect<-ifelse(plot_data$gamma==0,'Partial/null effect','Nonzero effect (.40)')
theme<-theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position='bottom',legend.title=element_blank(),
 axis.title.x=element_blank(),plot.title=element_text(face='bold'))
for(arm in levels(plot_data$arm)){
 z<-plot_data[plot_data$arm==arm & plot_data$gamma==.4,]
 g<-ggplot(z,aes(profile,coverage,color=method))+geom_hline(yintercept=.95,linetype=2,color='grey50')+
  geom_errorbar(aes(ymin=coverage_mc_low,ymax=coverage_mc_high),position=position_dodge(width=.75),width=.1)+
  geom_point(position=position_dodge(width=.75),size=2.2)+scale_color_manual(values=palette,drop=TRUE)+
  scale_x_discrete(labels=labels)+scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(100*x),'%'))+
  labs(title=paste(arm,'— nonzero-effect coverage'),subtitle='Available intervals only; bars show 95% Monte Carlo intervals',y='Coverage')+
  theme+guides(color=guide_legend(nrow=3,byrow=TRUE))
 stem<-if(grepl('^A:',arm))'predictor_coverage' else 'outcome_coverage'
 ggsave(file.path(out,paste0(stem,'.png')),g,width=12,height=6,dpi=160)
 ggsave(file.path(out,paste0(stem,'.pdf')),g,width=12,height=6)
}
g<-ggplot(plot_data,aes(profile,availability,fill=method))+geom_col(position=position_dodge(width=.8),width=.75)+
 facet_grid(arm~effect)+scale_fill_manual(values=palette,drop=FALSE)+scale_x_discrete(labels=labels)+
 scale_y_continuous(limits=c(0,1),labels=function(x)paste0(round(100*x),'%'))+
 labs(title='Interval availability is part of the comparison',subtitle='Full influences refuse singular/unstable measurement linearizations; bootstrap retains its failure gate',y='Returned intervals / 100')+
 theme+theme(axis.text.x=element_text(angle=25,hjust=1))+guides(fill=guide_legend(nrow=3,byrow=TRUE))
ggsave(file.path(out,'availability.png'),g,width=14,height=8,dpi=160)
ggsave(file.path(out,'availability.pdf'),g,width=14,height=8)
saveRDS(list(performance=performance,paired=paired,diagnostics=diag_rows,components=components,validation=checks,
 sensitivity=sensitivity_rows,steps=step_rows,profiles=profiles,legacy_normal=normal),file.path(out,'summary.rds'))
cat('Primary comparisons:\n');print(performance[performance$method%in%c('Conditional M-MOM Wald','Full influence Wald','Full person percentile','Legacy OLS stacked HC3 (t)'),
 c('profile','arm','gamma','method','available','coverage','rejection','bias')],row.names=FALSE)
