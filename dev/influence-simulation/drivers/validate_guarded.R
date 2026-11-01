#!/usr/bin/env Rscript
# Final public API at the fixed first-three/profile subset, including the later
# pointwise-recipe scope guard. Do not overwrite the frozen simulation library.
args <- commandArgs(TRUE)
lib <- if(length(args)) args[1] else '/private/tmp/repair-if-pr-lib'
file_arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
study_dir <- dirname(normalizePath(sub('^--file=','',file_arg[[1]])))
research_root <- normalizePath(file.path(study_dir,'../../../..'))
.libPaths(c(lib,.libPaths())); library(repair)
source(file.path(study_dir,'engine.R')); legacy <- load_comparison_sources(research_root)
paths <- list.files(file.path(study_dir,'results/sensitivity'),'^p[0-9]+_r[0-9]+[.]rds$',full.names=TRUE)
stopifnot(length(paths)==18L)
cl <- parallel::makePSOCKcluster(6L)
tryCatch({
 parallel::clusterExport(cl,c('lib','study_dir','research_root'),envir=.GlobalEnv)
 parallel::clusterEvalQ(cl,{
  .libPaths(c(lib,.libPaths())); library(repair)
  source(file.path(study_dir,'engine.R')); legacy <- load_comparison_sources(research_root);NULL
 })
 rows <- parallel::parLapply(cl,paths,function(path){
  reference <- readRDS(path)
  r <- readRDS(file.path(study_dir,'../robust_refit/results/main',basename(path)))
  sim <- r$data; model <- native_measurement(sim$repeated)$model
  z <- (sim$outcomes[,1]-1.5-.4*sim$people$latent_intercept)/sqrt(1-(.4*r$profile$tau)^2)
  rows <- list()
  for(arm in c('predictor','outcome'))for(g in 1:2){
   people <- sim$people
   if(arm=='predictor'){
    people$outcome <- sim$outcomes[,g]
    fit <- repair(outcome~re('(Intercept)','growth')+re('time','growth'),list(growth=model),people,'person',backend='rcpp')
   }else{
    people$x <- sim$truth[g]*sim$people$latent_slope/r$profile$tau+sqrt(1-sim$truth[g]^2)*z
    fit <- repair(re('time','growth')~x,list(growth=model),people,'person',backend='rcpp')
   }
   inf <- repair_sandwich(fit)
   old <- reference$influences[[paste(arm,g,sep='_')]]
   stopifnot(identical(inf$diagnostics$full_available,old$diagnostics$full_available),
    identical(inf$diagnostics$stage2_available,old$diagnostics$stage2_available),identical(inf$ids,old$ids))
   errors <- vapply(c('stage2','full'),function(key){
    if(!inf$diagnostics[[paste0(key,'_available')]])return(NA_real_)
    delta <- max(abs(inf$covariance[[key]]-old$covariance[[key]]))/(1e-8+max(abs(old$covariance[[key]])))
    stopifnot(delta<1e-10);delta
   },numeric(1))
   rows[[length(rows)+1L]] <- data.frame(profile=r$profile$profile,rep=r$job$rep,arm=arm,gamma=sim$truth[g],
     full_available=inf$diagnostics$full_available,stage2_relative=errors['stage2'],full_relative=errors['full'])
  }
  do.call(rbind,rows)
 })
 rows <- do.call(rbind,rows)
 write.csv(rows,file.path(study_dir,'results/summary/guarded_api_validation.csv'),row.names=FALSE)
 cat('PASS 72 guarded public API comparisons; identical availability, IDs, and covariance (relative tolerance 1e-10).\n')
},finally=parallel::stopCluster(cl))
