# Full original study: 50 replications per case, eta=0.9, n=64..4096.
study_ns <- c(64L,128L,256L,512L,1024L,2048L,4096L)
study_reps <- 50L
# Stages 1–8 were fixed before sampling. Stages 9–12 are a later
# user-authorized computation amendment for original job621 only.
study_iterations <- c(24000L,48000L,96000L,192000L,384000L,768000L,1536000L,3072000L,
  6144000L,12288000L,24576000L,49152000L)
study_chains <- c(3L,4L,6L,8L,8L,8L,8L,8L,8L,8L,8L,8L)
study_workers <- c(10L,10L,10L,8L,4L,2L,1L,1L,1L,1L,1L,1L)
study_n_seed <- c(64L,128L,256L,512L,1024L,4096L,8192L,2048L)
study_n_jobs <- c(64L,128L,256L,512L,1024L,4096L,2048L)
study_targets <- c('half','one')
study_families <- c('gaussian','uniform','student10')
study_variants <- c('N/G','N/P','N/H','U/G','t10/G')
study_beta <- c(half=.5,one=1)
study_ladder <- c(0,.0005,.001,.002,.004,.008,.016,.032,
  .064,.128,.192,.256,.32,.384,.448,.5)*(.9/.5)
study_truth <- function(x,target) {
  a <- abs(sin(2*pi*x)); if(target=='half')sqrt(a)else a
}
study_plot_truth <- function(x,target) {
  value <- study_truth(x,target);value[x%in%c(0,.5,1)] <- 0;value
}
study_grid <- function(G) (seq_len(G)-.5+sqrt(2)-1)/G
study_q <- function(n,beta) (log(n)/n)^(2*beta/(2*beta+1))
study_data <- function(seed,n,target) {
  set.seed(as.integer(seed))
  x <- runif(1024L); errors <- rnorm(1024L)
  if(n>1024L) {
    x <- c(x,runif(3072L)); errors <- c(errors,rnorm(3072L))
  }
  x <- x[seq_len(n)]
  list(x=x,y=study_truth(x,target)+errors[seq_len(n)])
}
study_cutoff <- function(x,y) {
  m <- ceiling(log(length(y))); cells <- pmin(m,floor(m*x)+1L)
  counts <- tabulate(cells,nbins=m)
  averages <- vapply(seq_len(m),function(j)
    if(counts[j])mean(y[cells==j])else 0,0.)
  list(Bhat=2+max(abs(averages)),cells=m,min_count=min(counts))
}
study_seeds <- function(text) as.integer(strsplit(as.character(text),';',fixed=TRUE)[[1]])
study_key <- function(x) paste(x$target_id,x$family,x$n,x$replicate,sep=':')
study_amendment_sha <- 'f3341cb4aa8e57c787fcc5df60038c2dc004c76b92bca01247bc02a174f603b0'
study_extension_fit <- function(x) {
  required <- c('job','target_id','family','n','replicate','data_seed')
  if(!all(required%in%names(x)))return(rep(FALSE,nrow(x)))
  x$job==621L & x$target_id=='half' & x$family=='uniform' &
    x$n==4096L & x$replicate==15L & x$data_seed==20675929L
}
study_stage_allowed <- function(x,stage) {
  valid <- is.finite(stage) & stage==as.integer(stage) & stage%in%seq_along(study_iterations)
  valid & (stage<=8L | study_extension_fit(x))
}
study_protocol_binding <- function(x) {
  if(!all(c('attempt','protocol_amendment_sha256')%in%names(x)))return(rep(FALSE,nrow(x)))
  hash <- as.character(x$protocol_amendment_sha256)
  supplemental <- x$attempt>8L
  (supplemental & !is.na(hash) & hash==study_amendment_sha) |
    (!supplemental & (is.na(hash) | hash==''))
}
study_diagnostic <- function(xs) {
  if(any(vapply(xs,function(x)!all(is.finite(x))||sd(x)==0,TRUE)))
    return(c(rhat=NA_real_,ess=NA_real_))
  m <- coda::mcmc.list(lapply(xs,coda::mcmc))
  c(rhat=tryCatch(unname(coda::gelman.diag(m,autoburnin=FALSE)$psrf[1L,1L]),
    error=function(e)NA_real_),ess=tryCatch(sum(vapply(xs,function(x)
      as.numeric(coda::effectiveSize(coda::mcmc(x))),0.)),error=function(e)NA_real_))
}
study_pass <- function(x) {
  values <- unlist(x[c('width_rhat','width_ess','max_probe_rhat','min_probe_ess',
    'max_original_probe_rhat','min_original_probe_ess','cap_fraction')],use.names=FALSE)
  all(is.finite(values)) && x$width_rhat<=1.10 && x$width_ess>=50 &&
    x$max_probe_rhat<=1.10 && x$min_probe_ess>=50 &&
    x$max_original_probe_rhat<=1.10 && x$min_original_probe_ess>=50 && x$cap_fraction==0
}
study_validate <- function(native,risk) {
  stopifnot(all(c('job','target_id','family','n','replicate','data_seed','attempt')%in%names(native)),
    all(study_protocol_binding(native)),all(study_protocol_binding(risk)),
    nrow(native)==2100L,nrow(risk)==3500L,
    all(native$target_id%in%study_targets),all(risk$target_id%in%study_targets),
    all(native$n%in%study_ns),all(risk$n%in%study_ns),
    all(native$family%in%study_families),all(risk$variant%in%study_variants),
    all(native$eta==.9),all(risk$eta==.9),all(native$accepted),all(risk$accepted),
    !anyDuplicated(study_key(native)),all(is.finite(risk$risk)),all(risk$risk>0))
  for(t in study_targets)for(n in study_ns)for(f in study_families) {
    z <- native[native$target_id==t&native$n==n&native$family==f,]
    stopifnot(nrow(z)==50L,setequal(z$replicate,1:50))
    stopifnot(all(vapply(seq_len(nrow(z)),function(i)study_pass(z[i,]),TRUE)))
  }
  for(t in study_targets)for(n in study_ns)for(v in study_variants) {
    z <- risk[risk$target_id==t&risk$n==n&risk$variant==v,]
    stopifnot(nrow(z)==50L,setequal(z$replicate,1:50))
  }
  for(i in seq_len(nrow(native))) {
    z <- native[i,]; seeds <- study_seeds(z$chain_seeds)
    stopifnot(length(seeds)==z$chains,all(is.finite(seeds)),!anyDuplicated(seeds),
      z$iterations>0L,z$iterations%%40L==0L,z$thin==20L,z$Kmax==64L,
      z$swap_every==5L,z$ladder_levels==16L,z$Bhat>0)
    stage <- as.integer(z$attempt);stopifnot(z$attempt==stage,study_stage_allowed(z,stage),
      z$iterations==study_iterations[stage],z$chains==study_chains[stage],
      z$burn==z$iterations%/%2L,z$original_replication_id==z$replicate)
    dseed <- (if(z$target_id=='half')20260929L+400000L else 20261001L+500000L)+1000L*z$replicate
    cseed <- dseed+1000000L*stage+100000L*match(z$n,study_n_seed)+
      10000L*match(z$family,study_families)+500L+seq_len(z$chains)
    job <- (z$replicate-1L)*42L+(match(z$n,study_n_jobs)-1L)*6L+
      (match(z$family,study_families)-1L)*2L+match(z$target_id,study_targets)
    stopifnot(z$data_seed==dseed,z$job==job,identical(as.integer(seeds),as.integer(cseed)),
      z$max_width_seen<64L,z$Bhat>=2)
    ladder <- as.numeric(strsplit(z$ladder_values,',',fixed=TRUE)[[1]])
    stopifnot(length(ladder)==16L,max(abs(ladder-study_ladder))<1e-15)
    p <- risk[study_key(risk)==study_key(z),]
    expected <- if(z$family=='gaussian')study_variants[1:3]else
      if(z$family=='uniform')'U/G'else't10/G'
    stopifnot(setequal(p$variant,expected),all(p$data_seed==z$data_seed),
      all(p$attempt==z$attempt),all(p$chains==z$chains),
      all(p$iterations==z$iterations),max(abs(p$Bhat-z$Bhat))<1e-12)
  }
  invisible(TRUE)
}
study_cells <- function(risk) {
  output <- list()
  for(t in study_targets)for(n in study_ns)for(v in study_variants) {
    z <- risk[risk$target_id==t&risk$n==n&risk$variant==v,]
    stopifnot(nrow(z)==50L,setequal(z$replicate,1:50))
    meanrisk <- mean(z$risk); se <- sd(z$risk)/sqrt(50L)
    q <- study_q(n,study_beta[[t]])
    output[[length(output)+1L]] <- data.frame(target=t,beta=study_beta[[t]],n=n,
      variant=v,family=z$family[1L],planned_count=50L,accepted_count=50L,
      failed_count=0L,pending_count=0L,complete=TRUE,datasets=50L,
      mean_Bhat=mean(z$Bhat),mean_risk=meanrisk,se_risk=se,ci_level=.95,
      ci_half=qnorm(.975)*se,ci_lower=meanrisk-qnorm(.975)*se,
      ci_upper=meanrisk+qnorm(.975)*se,q=q,mean_risk_over_q=meanrisk/q,
      se_risk_over_q=se/q,normalized_ci_lower=(meanrisk-qnorm(.975)*se)/q,
      normalized_ci_upper=(meanrisk+qnorm(.975)*se)/q,
      ci_below_zero=meanrisk-qnorm(.975)*se<=0)
  }
  do.call(rbind,output)
}
study_rates <- function(risk,B=2000L) {
  cells<-study_cells(risk)
  pairs<-rbind(data.frame(n_from=study_ns[-length(study_ns)],n_to=study_ns[-1L],type='adjacent'),
    data.frame(n_from=min(study_ns),n_to=max(study_ns),type='endpoint'))
  out<-list()
  for(t in study_targets) {
    set.seed(20261007L+if(t=='half')500L else 1000L)
    indices<-replicate(B,sample.int(study_reps,study_reps,replace=TRUE))
    for(v in study_variants)for(i in seq_len(nrow(pairs))) {
      pair<-pairs[i,];c<-cells[cells$target==t&cells$variant==v&cells$n%in%c(pair$n_from,pair$n_to),]
      complete<-nrow(c)==2L&&all(c$complete)
      estimate<-lower<-upper<-NA_real_
      if(complete) {
        p<-risk[risk$target_id==t&risk$variant==v,]
        a<-p[p$n==pair$n_from,];b<-p[p$n==pair$n_to,]
        a<-a$risk[match(seq_len(study_reps),a$replicate)];b<-b$risk[match(seq_len(study_reps),b$replicate)]
        stopifnot(!anyNA(a),!anyNA(b))
        denominator<-log(pair$n_to/pair$n_from)
        estimate<-log(mean(a)/mean(b))/denominator
        boot<-log(colMeans(matrix(a[indices],study_reps,B))/colMeans(matrix(b[indices],study_reps,B)))/denominator
        bounds<-quantile(boot,c(.025,.975),names=FALSE,type=7)
        lower<-bounds[1L];upper<-bounds[2L]
      }
      out[[length(out)+1L]]<-data.frame(target=t,beta=study_beta[[t]],variant=v,
        type=pair$type,n_from=pair$n_from,n_to=pair$n_to,complete=complete,
        paired_replications=if(complete)study_reps else NA_integer_,decay_slope=estimate,
        ci_lower=lower,ci_upper=upper,ci_level=.95,bootstrap_draws=B,
        ci_method='paired percentile bootstrap of log ratio of mean errors',
        illustrative_q_slope=log(study_q(pair$n_from,study_beta[[t]])/study_q(pair$n_to,study_beta[[t]]))/
          log(pair$n_to/pair$n_from),row.names=NULL)
    }
  }
  do.call(rbind,out)
}
