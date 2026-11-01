load_comparison_sources <- function(research_root) {
  source(file.path(research_root,'vig_hallquist_2026/analysis/repair_bootstrap_pilot/robust_refit/engine.R'),local=.GlobalEnv)
  legacy <- new.env(parent=globalenv())
  sys.source(file.path(research_root,'R/derivative_backends.R'),legacy)
  sys.source(file.path(research_root,'R/stacked_sandwich_helpers.R'),legacy)
  legacy
}

ci_from_values <- function(values,usable,eligible=TRUE) {
  available <- eligible && length(values)>=2L && all(usable) && all(is.finite(values))
  if(!available)return(c(se=NA_real_,low=NA_real_,high=NA_real_))
  c(se=sd(values),low=unname(quantile(values,.025,type=7)),high=unname(quantile(values,.975,type=7)))
}

cached_outcome_bootstrap <- function(r,model,people,g,mode='ordinary') {
  b<-r[[mode]];B<-nrow(b$draws);n<-nrow(people)
  estimates<-frozen<-rep(NA_real_,B);usable<-frozen_usable<-rep(FALSE,B)
  fixed<-lme4::fixef(model);sigma<-stats::sigma(model)
  h<-r$data$individual_coefficients[,2]-fixed['time']
  error_shape<-r$data$design_inverse[2,2,]
  # Focal slope only: common score centering changes the intercept, not this
  # slope. This is NOT an approximation to the bootstrap intercept covariance.
  for(i in 1:B){
    ids<-b$draws[i,];data<-list(outcome=h[ids],predictors=matrix(people$x[ids],n,1,
      dimnames=list(NULL,'x')),err_x=array(0,c(1L,1L,n)),err_y=sigma^2*error_shape[ids],err_xy=matrix(0,n,1L))
    fit<-fit_mmom(data,backend='rcpp')
    if(isTRUE(fit@converged) && all(is.finite(fit@coefficients$estimate))){
      frozen[i]<-fit@coefficients$estimate[2];frozen_usable[i]<-TRUE
    }
    if(!b$measurement_ok[i])next
    data$err_y<-b$sigma_hat[i]^2*error_shape[ids]
    fit<-fit_mmom(data,backend='rcpp')
    if(isTRUE(fit@converged) && all(is.finite(fit@coefficients$estimate))){
      estimates[i]<-fit@coefficients$estimate[2];usable[i]<-TRUE
    }
  }
  list(estimates=estimates,usable=usable,frozen=frozen,frozen_usable=frozen_usable,draws=b$draws)
}

run_comparison <- function(path,phase,out_dir,legacy,check_public=FALSE,step_check=FALSE){
  destination<-file.path(out_dir,basename(path));if(file.exists(destination))return(destination)
  started<-proc.time()[[3]];stage<-'read';warnings<-messages<-character()
  result<-tryCatch(withCallingHandlers({
    r<-readRDS(path);stopifnot(is.null(r$error));sim<-r$data;p<-r$profile
    stage<-'original_measurement';measured<-native_measurement(sim$repeated)
    stopifnot(!is.null(measured$model));model<-measured$model
    rows<-list();influences<-list();checks<-list();bootB<-list();legacy_results<-list()
    add<-function(arm,g,method,estimate,se,low,high,reason='',eligible=TRUE){
      scale<-if(arm=='predictor')p$tau else 1/p$tau
      available<-eligible && all(is.finite(c(estimate,se,low,high))) && se>=0
      rows[[length(rows)+1L]]<<-data.frame(profile=p$profile,rep=r$job$rep,arm=arm,gamma=sim$truth[g],method=method,
        estimate=estimate*scale,se=if(available)se*scale else NA_real_,low=if(available)low*scale else NA_real_,
        high=if(available)high*scale else NA_real_,available=available,
        covers=if(available)low*scale<=sim$truth[g] && high*scale>=sim$truth[g] else NA,
        rejects_zero=if(available)low>0||high<0 else NA,width=if(available)(high-low)*scale else NA_real_,
        reason=if(available)'' else reason,measurement_singular=lme4::isSingular(model),stringsAsFactors=FALSE)
    }
    # One marginal measurement fit serves both arms and both structural outcomes.
    context<-measurement_if<-NULL
    z<-(sim$outcomes[,1]-1.5-.4*sim$people$latent_intercept)/sqrt(1-(.4*p$tau)^2)
    legacy_data<-data.frame(id=sim$repeated$person,z=sim$repeated$time,y=sim$repeated$response)
    stage<-'legacy_measurement'
    legacy_model<-repair_lmer_refitter(y~z+(z|id),REML=FALSE)(legacy_data)
    legacy_split<-split(legacy_data,legacy_data$id)
    legacy_psi<-legacy$pack_psi(legacy_model)
    for(arm in c('predictor','outcome'))for(g in 1:2){
      people<-sim$people
      if(arm=='predictor'){
        people$outcome<-sim$outcomes[,g]
        fit<-repair(outcome~re('(Intercept)','growth')+re('time','growth'),list(growth=model),people,'person',backend='rcpp')
        focal<-3L
        stopifnot(max(abs(repair_coefficients(fit)$estimate-r$original[[g]]$coefficients$estimate))<1e-7,
          max(abs(repair_coefficients(fit)$std.error-r$original[[g]]$coefficients$std.error))<1e-7)
      }else{
        people$x<-sim$truth[g]*sim$people$latent_slope/p$tau+sqrt(1-sim$truth[g]^2)*z
        fit<-repair(re('time','growth')~x,list(growth=model),people,'person',backend='rcpp');focal<-2L
      }
      stopifnot(isTRUE(fit@struc_fit@converged));table<-repair_coefficients(fit)
      add(arm,g,'Conditional M-MOM Wald',table$estimate[focal],table$std.error[focal],table$conf.low[focal],table$conf.high[focal])
      stage<-'influence'
      if(is.null(context)){
        context<-repair:::repair_if_context(fit)
        measurement_if<-tryCatch(repair:::repair_if_measurement(context,1e-4,.03),error=function(e)e)
      }
      inf<-repair:::repair_sandwich_compute(fit,1e-4,.03,context,measurement_if)
      influences[[paste(arm,g,sep='_')]]<-inf
      for(component in c('stage2','full')){
        ok<-inf$diagnostics[[paste0(component,'_available')]]
        se<-if(ok)sqrt(inf$covariance[[component]][focal,focal]) else NA_real_
        add(arm,g,if(component=='full')'Full influence Wald' else 'Stage2 influence Wald',
          table$estimate[focal],se,table$estimate[focal]-qnorm(.975)*se,table$estimate[focal]+qnorm(.975)*se,
          inf$diagnostics[[paste0(component,'_reason')]],ok)
      }
      if(step_check){
        smaller_measurement<-tryCatch(repair:::repair_if_measurement(context,5e-5,.03),error=function(e)e)
        smaller<-repair:::repair_sandwich_compute(fit,5e-5,.03,context,smaller_measurement)
        checks[[paste(arm,g,'step',sep='_')]]<-list(full_same=identical(smaller$diagnostics$full_available,inf$diagnostics$full_available),
          stage2_same=identical(smaller$diagnostics$stage2_available,inf$diagnostics$stage2_available),
          full_relative=if(inf$diagnostics$full_available && smaller$diagnostics$full_available)
            max(abs(smaller$covariance$full-inf$covariance$full))/(1e-8+max(abs(inf$covariance$full))) else NA_real_,
          stage2_relative=if(inf$diagnostics$stage2_available && smaller$diagnostics$stage2_available)
            max(abs(smaller$covariance$stage2-inf$covariance$stage2))/(1e-8+max(abs(inf$covariance$stage2))) else NA_real_)
      }
      stage<-'bootstrap_comparison'
      if(arm=='predictor'){
        for(mode in c('ordinary','stratified')){
          b<-r[[mode]];if(is.null(b))next
          ci<-ci_from_values(b$estimates[,3,g],b$usable[,g],r$original[[g]]$preflight_ok)
          add(arm,g,if(mode=='ordinary')'Full person percentile' else 'Stratified person percentile',
            table$estimate[focal],ci['se'],ci['low'],ci['high'],'recorded bootstrap failed-vector gate')
          if(mode=='ordinary'){
            ci<-ci_from_values(b$frozen[,3,g],b$frozen_usable[,g],r$original[[g]]$point_ok)
            add(arm,g,'Frozen measurement percentile',table$estimate[focal],ci['se'],ci['low'],ci['high'],'recorded frozen failed-vector gate')
          }
        }
      }else{
        replay<-tryCatch(repair_refit(fit,sim$repeated,people),error=function(e)e)
        eligible<-!inherits(replay,'error')
        for(mode in c('ordinary','stratified')){
          if(is.null(r[[mode]]))next
          b<-cached_outcome_bootstrap(r,model,people,g,mode)
          bootB[[paste(g,mode,sep='_')]]<-b
          ci<-ci_from_values(b$estimates,b$usable,eligible)
          add(arm,g,if(mode=='ordinary')'Full person percentile' else 'Stratified person percentile',
            table$estimate[focal],ci['se'],ci['low'],ci['high'],'cached bootstrap failed-vector gate',eligible)
          if(mode=='ordinary'){
            ci<-ci_from_values(b$frozen,b$frozen_usable)
            add(arm,g,'Frozen measurement percentile',table$estimate[focal],ci['se'],ci['low'],ci['high'],'frozen failed-vector gate')
          }
        }
        stage<-'legacy_comparator'
        old<-tryCatch(legacy$stacked_sandwich_for_corrected_scores(legacy_split,
          data.frame(id=people$person,x=people$x),legacy_model,legacy_psi,legacy$make_derivative_backend('analytical')),error=function(e)e)
        legacy_results[[as.character(g)]]<-old
        for(hc in 0:3){
          ok<-!inherits(old,'error') && is.finite(old[[paste0('vcov_hc',hc)]][2,2]) && old[[paste0('vcov_hc',hc)]][2,2]>=0
          point<-if(!inherits(old,'error'))old$alpha_hat[2] else NA_real_
          se<-if(ok)sqrt(old[[paste0('vcov_hc',hc)]][2,2]) else NA_real_;crit<-qt(.975,nrow(people)-2L)
          add(arm,g,paste0('Legacy OLS stacked HC',hc,' (t)'),point,se,point-crit*se,point+crit*se,
            if(inherits(old,'error'))conditionMessage(old) else 'legacy covariance nonfinite',ok)
        }
        if(!inherits(old,'error')){
          X<-cbind(1,people$x);ols<-drop(solve(crossprod(X),crossprod(X,old$corrected_scores)))
          resid<-old$corrected_scores-drop(X%*%ols);inv<-solve(crossprod(X));hc0<-inv%*%crossprod(resid*X)%*%inv
          checks[[paste(g,'legacy',sep='_')]]<-list(criterion='ML',singular=lme4::isSingular(legacy_model),
            hc0_gap=abs(hc0[2,2]-old$vcov_hc0[2,2]),
            balanced_point_gap=if(p$small==p$large)abs(point-table$estimate[focal]) else NA_real_)
          if(p$small==p$large)stopifnot(abs(point-table$estimate[focal])<1e-7)
        }
        if(check_public){
          stage<-'public_validation';k<-min(99L,r$B)
          public<-tryCatch(repair_bootstrap(fit,sim$repeated,people,B=k,seed=r$job$bootstrap_seed),error=function(e)e)
          if(inherits(public,'error')){stopifnot(!eligible);checks[[paste(g,'public',sep='_')]]<-list(B=0L,preflight_agreed=TRUE)} else {
            b<-bootB[[paste(g,'ordinary',sep='_')]]
            stopifnot(identical(public$draws,b$draws[1:k,,drop=FALSE]),identical(public$diagnostics$usable,b$usable[1:k]))
            error<-max(abs(public$replicates[,focal]-b$estimates[1:k]),na.rm=TRUE)
            stopifnot(error<1e-7)
            checks[[paste(g,'public',sep='_')]]<-list(B=k,error=error,usable=public$diagnostics$usable)
          }
        }
      }
    }
    list(profile=p,job=r$job,B=r$B,phase=phase,rows=do.call(rbind,rows),influences=influences,checks=checks,
      outcome_bootstrap=bootB,legacy=legacy_results,
      measurement=list(singular=context$singular,converged=context$converged,criterion='REML',
        psi=context$psi,measurement_if=if(inherits(measurement_if,'error'))conditionMessage(measurement_if) else measurement_if[c('stationarity','rcond','derivative_error')]))
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},
    message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart('muffleMessage')}),
    error=function(e)list(error=conditionMessage(e),stage=stage,source=path))
  result$elapsed<-proc.time()[[3]]-started;result$warnings<-warnings;result$messages<-messages
  temp<-paste0(destination,'.tmp');saveRDS(result,temp,compress='xz');file.rename(temp,destination)
  cat(basename(path),phase,sprintf('%.1fs',result$elapsed),if(is.null(result$error))'ok' else paste(result$stage,result$error),'\n');flush.console()
  destination
}
