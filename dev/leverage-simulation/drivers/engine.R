load_sources <- function(root) source(file.path(root,'vig_hallquist_2026/analysis/repair_bootstrap_pilot/robust_refit/engine.R'),local=.GlobalEnv)

delete_pair <- function(sim,original) {
 n <- nrow(sim$people); estimates <- array(NA_real_,c(n,3L,2L)); usable <- matrix(FALSE,n,2L)
 reason <- matrix('',n,2L); residual <- correction <- matrix(NA_real_,n,2L)
 optimization <- vector('list',n); singular <- rep(NA,n); measurement_ok <- rep(FALSE,n)
 for(i in seq_len(n)) {
  keep <- seq_len(n)[-i]; measured <- native_measurement(resample_data(sim,keep))
  optimization[[i]] <- measured$optimization; singular[i] <- measured$singular; measurement_ok[i] <- measured$usable
  for(g in 1:2) {
   if(!original[[g]]$preflight_ok) { reason[i,g]<-'Original reference replay unavailable';next }
   if(!measured$usable){reason[i,g]<-measured$reason;next}
   fit <- tryCatch(structural_from_cache(sim,keep,lme4::fixef(measured$model),stats::sigma(measured$model),
     sim$outcomes[,g],original[[g]]$controls),error=function(e)e)
   result <- extract_structural(fit); usable[i,g]<-result$usable;reason[i,g]<-result$reason
   if(result$usable){estimates[i,,g]<-result$beta;residual[i,g]<-fit@diagnostics$resid_var;correction[i,g]<-fit@diagnostics$correction_final}
  }
 }
 list(estimates=estimates,usable=usable,reason=reason,measurement_ok=measurement_ok,singular=singular,
      residual=residual,correction=correction,optimization=optimization)
}

run_fresh <- function(job,profiles,B,phase,out_dir) {
 path <- file.path(out_dir,sprintf('p%02d_r%04d.rds',job$profile,job$rep))
 if(file.exists(path))return(path)
 raw_dir <- file.path(out_dir,'native');dir.create(raw_dir,showWarnings=FALSE)
 started <- proc.time()[[3]]; warnings <- messages <- character();stage <- 'native'
 result <- tryCatch(withCallingHandlers({
  raw <- run_dataset(job,profiles,B,raw_dir,api_check=phase=='benchmark'||(phase=='main'&&job$rep<=2L))
  r <- readRDS(raw);if(!is.null(r$error))stop(paste(r$error_stage,r$error))
  sim <- r$data;p <- r$profile;measured <- native_measurement(sim$repeated);model <- measured$model
  stopifnot(measured$usable)
  stage <- 'delete-one';deleted <- delete_pair(sim,r$original)
  rows <- analytic <- checks <- list()
  add <- function(g,method,point,se=NA_real_,low=NA_real_,high=NA_real_,reason='',eligible=TRUE) {
   available <- eligible && all(is.finite(c(point,se,low,high))) && se>=0
   rows[[length(rows)+1L]] <<- data.frame(profile=p$profile,rep=job$rep,gamma=sim$truth[g],method=method,
    estimate=point*p$tau,se=if(available)se*p$tau else NA_real_,low=if(available)low*p$tau else NA_real_,
    high=if(available)high*p$tau else NA_real_,available=available,
    reason=if(available)'' else reason,measurement_singular=measured$singular)
  }
  add_wald <- function(g,method,point,V,eligible,reason='',reference='normal') {
   se <- if(eligible)sqrt(V[3,3]) else NA_real_
   critical <- if(reference=='normal')qnorm(.975) else qt(.975,p$n-1L)
   add(g,method,point,se,point-critical*se,point+critical*se,reason,eligible)
  }
  for(g in 1:2) {
   people <- sim$people;people$outcome<-sim$outcomes[,g]
   fit <- repair(outcome~re('(Intercept)','growth')+re('time','growth'),list(growth=model),people,'person',backend='rcpp')
   point <- coef(fit)[3];table<-repair_coefficients(fit)
   stopifnot(max(abs(table$estimate-r$original[[g]]$coefficients$estimate))<1e-7)
   add(g,'Conditional Wald',point,table$std.error[3],table$conf.low[3],table$conf.high[3])
   stage <- 'analytic';influence <- repair_sandwich(fit);leverage <- repair_leverage(fit)
   analytic[[g]] <- list(influence=influence,leverage=leverage)
   for(key in c('full','stage2'))add_wald(g,if(key=='full')'Full IF Wald' else 'Stage2 IF Wald',point,
    influence$covariance[[key]],influence$diagnostics[[paste0(key,'_available')]],influence$diagnostics[[paste0(key,'_reason')]])
   add_wald(g,'Full-system HC3 Wald',point,leverage$covariance$full,leverage$diagnostics$full_available,leverage$diagnostics$full_reason)
   add_wald(g,'Full-system HC3 t',point,leverage$covariance$full,leverage$diagnostics$full_available,leverage$diagnostics$full_reason,'t')
   if(phase=='sensitivity') {
    half <- repair_leverage(fit,step=5e-5)
    checks[[paste0('step_',g)]] <- list(same=identical(half$diagnostics$full_available,leverage$diagnostics$full_available),
     relative=if(half$diagnostics$full_available && leverage$diagnostics$full_available)
      max(abs(half$covariance$full-leverage$covariance$full))/(1e-8+max(abs(leverage$covariance$full))) else NA_real_)
   }
   stage <- 'intervals';ok <- all(deleted$usable[,g])
   V <- matrix(NA_real_,3,3)
   if(ok){centered<-scale(deleted$estimates[,,g],scale=FALSE);V<-(p$n-1)/p$n*crossprod(centered)}
   add_wald(g,'Delete-one Wald',point,V,ok,'At least one complete deletion vector failed')
   add_wald(g,'Delete-one t',point,V,ok,'At least one complete deletion vector failed','t')
   b<-r$ordinary;ok<-r$original[[g]]$preflight_ok && all(b$usable[,g])
   values<-b$estimates[,3,g];se<-if(ok)sd(values) else NA_real_
   ci<-if(ok)quantile(values,c(.025,.975)) else c(NA_real_,NA_real_)
   add(g,'Person percentile',point,se,ci[1],ci[2],'Bootstrap failed-vector gate',ok)
   add(g,'Bootstrap-SD Wald',point,se,point-qnorm(.975)*se,point+qnorm(.975)*se,'Bootstrap failed-vector gate',ok)
   ok<-all(b$frozen_usable[,g]);values<-b$frozen[,3,g]
   ci<-if(ok)quantile(values,c(.025,.975)) else c(NA_real_,NA_real_)
   add(g,'Frozen percentile',point,if(ok)sd(values) else NA_real_,ci[1],ci[2],'Frozen failed-vector gate',ok)
   if(!is.null(r$stratified)){
    b<-r$stratified;ok<-r$original[[g]]$preflight_ok && all(b$usable[,g]);values<-b$estimates[,3,g]
    ci<-if(ok)quantile(values,c(.025,.975)) else c(NA_real_,NA_real_)
    add(g,'Stratified percentile',point,if(ok)sd(values) else NA_real_,ci[1],ci[2],'Stratified failed-vector gate',ok)
   }
   if(phase=='benchmark'||(phase=='main'&&job$rep==1L)) {
    stage <- 'public_delete_validation'
    public <- repair_jackknife(fit,sim$repeated,people)
    stopifnot(identical(public$diagnostics$usable,deleted$usable[,g]))
    error<-max(abs(public$replicates-deleted$estimates[,,g]),na.rm=TRUE)
    stopifnot(error<1e-7)
    if(all(public$diagnostics$usable))stopifnot(max(abs(vcov(public)-V))<1e-7)
    checks[[paste0('jackknife_',g)]]<-list(n=p$n,error=error,public=public)
   }
  }
  list(job=job,profile=p,B=B,rows=do.call(rbind,rows),analytic=analytic,deletions=deleted,checks=checks,
   bootstrap_checks=r$validation,raw=basename(raw),measurement=measured[c('singular','usable','optimization')])
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},
   message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart('muffleMessage')}),
   error=function(e)list(error=conditionMessage(e),stage=stage,job=job))
 result$elapsed<-proc.time()[[3]]-started;result$warnings<-warnings;result$messages<-messages
 temp<-paste0(path,'.tmp');saveRDS(result,temp,compress='xz');file.rename(temp,path)
 cat(basename(path),sprintf('%.1fs',result$elapsed),if(is.null(result$error))'ok' else paste(stage,result$error),'\n');flush.console()
 path
}
