# Clipped latent-function summaries; no observation-noise variance is added.
credible_moments <- function(states, x, reference) {
  mu <- matrix(0, nrow = length(x), ncol = length(states))
  variance <- mu
  for (j in seq_along(states)) {
    state <- states[[j]]
    Phi <- reference$design_matrix(x, state$u, state$b)
    mu[, j] <- drop(Phi %*% state$m)
    variance[, j] <- colSums(forwardsolve(t(state$R), t(Phi))^2)
  }
  list(mu = mu, sd = sqrt(pmax(variance, 0)))
}

# Generalised inverse of an equally weighted Gaussian mixture CDF. Bracketing
# on [-B,B] accounts for clipping atoms, including zero conditional variances.
credible_mixture_quantile <- function(mu, sd, probability, B,
                                      tolerance = 1e-11) {
  point_count <- nrow(mu)
  zero <- sd == 0
  mixture_cdf <- function(value) {
    difference <- value - mu
    zz <- difference / sd
    zz[zero] <- ifelse(difference[zero] >= 0, Inf, -Inf)
    rowMeans(pnorm(zz))
  }
  lower <- rep(-B, point_count)
  upper <- rep(B, point_count)
  at_lower <- mixture_cdf(lower) >= probability
  at_upper <- mixture_cdf(upper) < probability
  active <- !at_lower & !at_upper
  for (iteration in seq_len(50L)) {
    midpoint <- (lower + upper) / 2
    cdf <- mixture_cdf(midpoint)
    move_upper <- active & cdf >= probability
    move_lower <- active & cdf < probability
    upper[move_upper] <- midpoint[move_upper]
    lower[move_lower] <- midpoint[move_lower]
    if (max(upper[active] - lower[active], 0) <= tolerance) break
  }
  answer <- upper
  answer[at_lower] <- -B
  answer[at_upper] <- B
  answer
}

credible_plot <- function(records,output_pdf,targets) {
  stopifnot(setequal(targets,study_targets))
  uniform_plot(records,sub('[.]pdf$','',output_pdf))
  invisible(output_pdf)
}

gaussian_band <- function(chains,reference,x,y,settings) {
  grid <- seq.int(0L,1200L)/1200
  states <- unlist(lapply(chains,function(chain)lapply(chain$kept_states,function(st)
    reference$state_eval(x,y,as.numeric(st$u),as.numeric(st$b),.9))),recursive=FALSE)
  moments <- credible_moments(states,grid,reference)
  mean <- rowMeans(matrix(reference$clipped_normal_mean(as.vector(moments$mu),
    as.vector(moments$sd)^2,settings$Bhat),nrow=length(grid)))
  lower <- credible_mixture_quantile(moments$mu,moments$sd,.025,settings$Bhat)
  upper <- credible_mixture_quantile(moments$mu,moments$sd,.975,settings$Bhat)
  data.frame(target=settings$target_id,beta=settings$beta,n=settings$n,
    replicate=settings$replicate,eta=.9,family='gaussian',variant='N/G',
    x=grid,truth=study_plot_truth(grid,settings$target_id),posterior_mean=mean,
    credible_lower=lower,credible_upper=upper)
}
uniform_band <- function(chains,settings) {
  grid <- seq.int(0L,1200L)/1200
  states <- unlist(lapply(chains,`[[`,'kept_states'),recursive=FALSE,use.names=FALSE)
  network <- function(state,grid) {
    value <- numeric(length(grid))
    for(j in seq_along(state$a))
      value <- value+state$a[j]*pmax(state$u[j]*grid+state$b[j],0)
    value
  }
  mean <- lower <- upper <- numeric(length(grid))
  for(indices in split(seq_along(grid),ceiling(seq_along(grid)/32L))) {
    values <- vapply(states,network,numeric(length(indices)),grid=grid[indices])
    values[] <- pmax(-settings$Bhat,pmin(settings$Bhat,values))
    mean[indices] <- rowMeans(values)
    quantiles <- t(apply(values,1L,quantile,probs=c(.025,.975),type=1,names=FALSE))
    lower[indices] <- quantiles[,1L];upper[indices] <- quantiles[,2L]
  }
  truth <- study_truth(grid,settings$target_id);truth[grid%in%c(0,.5,1)] <- 0
  data.frame(target=settings$target_id,beta=settings$beta,n=settings$n,replicate=settings$replicate,
    eta=.9,family='uniform',variant='U/G',x=grid,truth=truth,
    posterior_mean=mean,credible_lower=lower,credible_upper=upper)
}
