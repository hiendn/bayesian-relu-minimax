#!/usr/bin/env Rscript
# Reproduce original reported fits, optionally starting every fit at stage1.
args <- commandArgs(trailingOnly=TRUE)
if(!length(args))stop('Usage: Rscript Refit.R OUTPUT [--workers=4] [--rows=1,2,...] [--all-stages] [--check-only]')
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1L])
root <- dirname(normalizePath(script,mustWork=TRUE))
source(file.path(root,'source/study.R'))
source(file.path(root,'source/credible.R'))
native <- read.csv(file.path(root,'records/selected_native_fits.csv'),stringsAsFactors=FALSE)
expected <- read.csv(file.path(root,'records/selected_risks.csv'),stringsAsFactors=FALSE)
study_validate(native,expected)
get_flag <- function(prefix,default) {
  z <- args[startsWith(args,prefix)];stopifnot(length(z)<=1L)
  if(length(z))substring(z,nchar(prefix)+1L)else default
}
workers <- as.integer(get_flag('--workers=','4'))
rows <- as.integer(strsplit(get_flag('--rows=',paste(seq_len(nrow(native)),collapse=',')),',',fixed=TRUE)[[1L]])
stopifnot(workers>=1L,!anyNA(rows),all(rows%in%seq_len(nrow(native))),!anyDuplicated(rows),
  all(args[-1L]%in%c('--check-only','--all-stages')|startsWith(args[-1L],'--workers=')|startsWith(args[-1L],'--rows=')))
if('--check-only'%in%args) {
  cat('Validated final study settings:',length(rows),'of 2100 original native cases; no sampling.\n')
  quit(save='no')
}
out <- args[1L]
if(file.exists(out))stop('Output already exists; choose a new directory.')
dir.create(out,recursive=TRUE);out <- normalizePath(out,mustWork=TRUE)
native_dir <- file.path(root,'source/native')
preload <- new.env(parent=globalenv());preload$native_dir <- native_dir
sys.source(file.path(native_dir,'native_sampler.R'),envir=preload)
recorded <- new.env(parent=baseenv())
compile_recorded <- function() {
  old <- Sys.getenv(c('PKG_CPPFLAGS','R_MAKEVARS_USER'),unset=NA_character_)
  on.exit(for(key in names(old))if(is.na(old[[key]]))Sys.unsetenv(key)else
    do.call(Sys.setenv,setNames(list(old[[key]]),key)))
  if(Sys.info()[['sysname']]=='Darwin') {
    headers <- '/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1'
    if(file.exists(file.path(headers,'cmath')))Sys.setenv(
      PKG_CPPFLAGS=paste(Sys.getenv('PKG_CPPFLAGS'),paste0('-I',headers)))
    if(!dir.exists('/opt/gfortran/lib')&&file.exists(file.path(R.home('lib'),'libgfortran.5.dylib')))
      Sys.setenv(R_MAKEVARS_USER=file.path(native_dir,'local_Makevars'))
  }
  Rcpp::sourceCpp(file.path(native_dir,'direct_kernel_recording.cpp'),env=recorded,
    cacheDir=file.path(out,'recording_cpp_cache'),showOutput=FALSE,rebuild=FALSE)
}
compile_recorded()
reference <- new.env(parent=globalenv())
sys.source(file.path(root,'source/gaussian.R'),envir=reference)
atomic_rds <- function(value,path) {
  temporary <- paste0(path,'.tmp-',Sys.getpid());saveRDS(value,temporary)
  if(!file.rename(temporary,path))stop('Could not publish fit file.')
}
log_width <- function(k,law)switch(law,G=-k*log(2),P=dpois(k-1L,1,log=TRUE),H=-2.5*log(k))
run_fit <- function(i,stage=NULL) {
  settings <- native[i,];published_stage <- settings$attempt
  stopifnot(study_stage_allowed(settings,published_stage))
  if(!is.null(stage)) {
    stopifnot(study_stage_allowed(settings,stage))
    settings$attempt <- as.integer(stage);settings$iterations <- study_iterations[stage]
    settings$protocol_amendment_sha256 <- if(stage>8L)study_amendment_sha else ''
    settings$chains <- study_chains[stage];settings$burn <- settings$iterations%/%2L
    settings$chain_seeds <- paste(settings$data_seed+1000000L*stage+
      100000L*match(settings$n,study_n_seed)+10000L*match(settings$family,study_families)+
      500L+seq_len(settings$chains),collapse=';')
  }
  target <- settings$target_id;n <- settings$n
  family <- settings$family;seeds <- study_seeds(settings$chain_seeds)
  dat <- study_data(settings$data_seed,n,target);x <- dat$x;y <- dat$y
  cut <- study_cutoff(x,y);stopifnot(abs(cut$Bhat-settings$Bhat)<1e-12)
  grid <- study_grid(1009L);f0 <- study_truth(grid,target)
  e <- new.env(parent=globalenv());e$eta <- .9;e$native_dir <- native_dir
  sys.source(file.path(native_dir,'native_sampler.R'),envir=e);e$Bclip <- cut$Bhat
  iterations <- settings$iterations;burn <- iterations%/%2L;chains <- settings$chains
  is_uniform_panel <- family=='uniform'&&target%in%study_targets&&n%in%c(1024L,4096L)&&settings$replicate==1L
  is_gaussian_panel <- family=='gaussian'&&n%in%c(1024L,4096L)&&settings$replicate==1L
  if(family=='gaussian') {
    cs <- lapply(seq_len(chains),function(j)e$run_tempered_chain(x,y,grid,seeds[j],
      Kmax=settings$Kmax,iterations=iterations,burn=burn,thin=settings$thin,
      lambdas=study_ladder,init_k=j,swap_every=settings$swap_every))
    widths <- lapply(cs,`[[`,'kept_width');probes <- lapply(cs,`[[`,'kept_probe')
    cap <- max(vapply(cs,`[[`,0.,'at_cap'));means <- lapply(cs,`[[`,'predictive_mean')
    pred <- lapply(cs,function(ch) {
      if(!is.null(ch$kept_predictions))return(ch$kept_predictions)
      vapply(ch$kept_states,function(st)reference$conditional_predict(
        reference$state_eval(x,y,st$u,st$b,.9),grid,cut$Bhat),numeric(length(grid)))
    })
    stopifnot(all(vapply(seq_along(cs),function(j)
      max(abs(rowMeans(pred[[j]])-means[[j]]))<1e-10,TRUE)))
    pooled <- do.call(cbind,pred);allwidth <- unlist(widths,use.names=FALSE)
    estimates <- list()
    risk_rows <- lapply(c('G','P','H'),function(law) {
      lw <- log_width(allwidth,law)-log_width(allwidth,'G')
      w <- exp(lw-max(lw));w <- w/sum(w)
      estimate <- drop(pooled%*%w);variant <- c(G='N/G',P='N/P',H='N/H')[[law]]
      estimates[[variant]] <<- estimate
      chain_risks <- vapply(seq_len(chains),function(j) {
        l <- log_width(widths[[j]],law)-log_width(widths[[j]],'G')
        ww <- exp(l-max(l));ww <- ww/sum(ww)
        mean((drop(pred[[j]]%*%ww)-f0)^2)
      },0.)
      data.frame(variant=variant,risk=mean((estimate-f0)^2),
        chain_risk_span=diff(range(chain_risks)),overlap_fraction=1/(length(w)*sum(w^2)))
    })
  } else {
    cs <- lapply(seq_len(chains),function(j) {
      original <- e$run_chain(x,y,grid,f0,family,seeds[j],iterations=iterations,
        burn=burn,thin=settings$thin,Kmax=settings$Kmax,lambdas=study_ladder,
        swap_every=settings$swap_every)
      original_rng <- .Random.seed
      if(!is_uniform_panel)return(original)
      set.seed(seeds[j])
      saved <- recorded$native_direct_chain(x,y,grid,f0,family,as.integer(iterations),
        as.integer(burn),as.integer(settings$thin),as.integer(settings$Kmax),
        study_ladder,as.integer(settings$swap_every),cut$Bhat)
      recorded_rng <- .Random.seed
      fields <- setdiff(names(original),'seconds')
      stopifnot(identical(original_rng,recorded_rng),
        all(vapply(fields,function(field)identical(original[[field]],saved[[field]]),TRUE)))
      saved
    })
    widths <- lapply(cs,`[[`,'width');probes <- lapply(cs,`[[`,'pred_probes')
    cap <- max(vapply(cs,`[[`,0.,'cap_fraction'));means <- lapply(cs,`[[`,'posterior_mean')
    estimate <- Reduce('+',means)/chains
    chain_risks <- vapply(means,function(m)mean((m-f0)^2),0.)
    variant <- if(family=='uniform')'U/G'else't10/G'
    estimates <- setNames(list(estimate),variant)
    risk_rows <- list(data.frame(variant=variant,risk=mean((estimate-f0)^2),
      chain_risk_span=diff(range(chain_risks)),overlap_fraction=NA_real_))
  }
  raw <- vapply(seq_len(5L),function(j)study_diagnostic(lapply(probes,function(p)p[,j])),c(rhat=0.,ess=0.))
  if(family!='gaussian')probes <- lapply(probes,function(p){p[]<-pmax(-cut$Bhat,pmin(cut$Bhat,p));p})
  clipped <- vapply(seq_len(5L),function(j)study_diagnostic(lapply(probes,function(p)p[,j])),c(rhat=0.,ess=0.))
  wd <- study_diagnostic(widths)
  observed <- settings
  observed$width_rhat <- unname(wd['rhat']);observed$width_ess <- unname(wd['ess'])
  observed$max_probe_rhat <- max(clipped['rhat',]);observed$min_probe_ess <- min(clipped['ess',])
  observed$max_original_probe_rhat <- max(raw['rhat',]);observed$min_original_probe_ess <- min(raw['ess',])
  observed$cap_fraction <- cap;observed$max_width_seen <- max(unlist(widths,use.names=FALSE))
  observed$Bhat <- cut$Bhat
  observed$elapsed_seconds <- sum(vapply(cs,`[[`,0.,if(family=='gaussian')'elapsed_seconds'else'seconds'))
  observed$accepted <- study_pass(observed)
  result <- do.call(rbind,risk_rows)
  observed$max_chain_risk_span <- max(result$chain_risk_span)
  if(!is.null(stage)&&stage!=published_stage) {
    stopifnot(stage<published_stage,!observed$accepted)
    atomic_rds(list(settings=observed,risks=result,estimates=estimates,
      diagnostics=list(width=wd,raw=raw,clipped=clipped)),
      file.path(out,sprintf('%s_%s_n%d_rep%03d_stage%d.rds',family,target,n,settings$replicate,stage)))
    cat('Prescribed stage',stage,'did not pass:',family,target,n,settings$replicate,'\n');flush(stdout())
    return(list(native=observed,risk=result,band=NULL))
  }
  stopifnot(observed$accepted)
  published <- expected[study_key(expected)==study_key(settings),]
  published <- published[match(result$variant,published$variant),]
  stopifnot(max(abs(result$risk-published$risk))<1e-10,
    max(abs(result$chain_risk_span-published$chain_risk_span))<1e-10)
  for(column in c('width_rhat','width_ess','max_probe_rhat','min_probe_ess',
    'max_original_probe_rhat','min_original_probe_ess'))
    stopifnot(abs(observed[[column]]-settings[[column]])<1e-7)
  observed$max_chain_risk_span <- max(result$chain_risk_span)
  new_metrics <- result
  result <- published
  for(column in c('risk','chain_risk_span','overlap_fraction'))
    result[[column]] <- new_metrics[[column]]
  result$q <- study_q(n,settings$beta)
  result$risk_over_q <- result$risk/result$q
  for(column in c('width_rhat','width_ess','max_probe_rhat','min_probe_ess',
    'max_original_probe_rhat','min_original_probe_ess','cap_fraction','Bhat'))
    result[[column]] <- observed[[column]]
  mean_matrix <- do.call(cbind,means)
  result$mean_function_disagreement <- sqrt(mean(apply(mean_matrix,1L,function(v)diff(range(v)))^2))
  result$cells <- cut$cells;result$min_count <- cut$min_count
  band <- NULL
  if(is_gaussian_panel)band <- gaussian_band(cs,reference,x,y,settings)
  if(is_uniform_panel)band <- uniform_band(cs,settings)
  name <- sprintf('%s_%s_n%d_rep%03d',family,target,n,settings$replicate)
  atomic_rds(list(settings=observed,risks=result,estimates=estimates,
    diagnostics=list(width=wd,raw=raw,clipped=clipped),band=band),file.path(out,paste0(name,'.rds')))
  cat('Reproduced',name,'with recorded seeds/settings.\n');flush(stdout())
  list(native=observed,risk=result,band=band)
}
dispatch <- function(indices,stage=NULL) {
  count <- if(is.null(stage))workers else min(workers,study_workers[stage])
  task <- function(i)run_fit(i,stage)
  result <- if(count==1L||.Platform$OS.type=='windows')lapply(indices,task)else
    parallel::mclapply(indices,task,mc.cores=count,mc.preschedule=FALSE)
  stopifnot(!any(vapply(result,inherits,TRUE,'try-error')))
  result
}
if('--all-stages'%in%args) {
  pending <- rows;results <- vector('list',length(rows))
  for(stage in seq_along(study_iterations)) {
    if(!length(pending))break
    stopifnot(all(study_stage_allowed(native[pending,],stage)))
    current <- dispatch(pending,stage)
    pass <- vapply(current,function(z)isTRUE(z$native$accepted),TRUE)
    for(j in which(pass))results[[match(pending[j],rows)]] <- current[[j]]
    pending <- pending[!pass]
  }
  stopifnot(!length(pending),all(vapply(results,Negate(is.null),TRUE)))
} else results <- dispatch(rows)
write.csv(do.call(rbind,lapply(results,`[[`,'native')),file.path(out,'selected_native_fits.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(results,`[[`,'risk')),file.path(out,'selected_risks.csv'),row.names=FALSE)
bands <- Filter(Negate(is.null),lapply(results,`[[`,'band'))
if(length(bands)) {
  all_bands <- do.call(rbind,bands)
  write.csv(all_bands[all_bands$family=='uniform',],file.path(out,'uniform_credible_intervals_data.csv'),row.names=FALSE)
  write.csv(all_bands[all_bands$family=='gaussian',],file.path(out,'gaussian_credible_intervals_data.csv'),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(out,'runtime_session.txt'))
cat('Completed',length(rows),'original configurations. Replicate indices identify the original data streams.\n')
