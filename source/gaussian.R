# Gaussian state evaluation and conditional clipped prediction.
# The native parallel-tempered transition kernel is in native/gaussian_kernel.cpp.
eta <- 0.9
Bclip <- as.numeric(Sys.getenv('SIM_BCLIP', '32'))
stopifnot(is.finite(Bclip), Bclip > 0)
quadrature_grid <- function(G) (seq_len(G)-0.5+sqrt(2)-1)/G

design_matrix <- function(x, u, b) {
  k <- length(u)
  if (k == 1L) return(matrix(pmax(x * u[1L] + b[1L], 0), ncol = 1L))
  pmax(outer(x, u) + matrix(rep(b, each = length(x)), nrow = length(x)), 0)
}

# The data-constant -eta*y'y/2 is omitted from log_marginal because it
# cancels in every Metropolis ratio. All K-dependent terms are retained.
state_eval <- function(x, y, u, b, eta = 0.9) {
  k <- length(u)
  Phi <- design_matrix(x, u, b)
  tau2 <- k^3
  G <- crossprod(Phi)
  rhs <- eta * drop(crossprod(Phi, y))
  Q <- eta * G + diag(1 / tau2, k)
  R <- chol(Q)
  m <- drop(backsolve(R, forwardsolve(t(R), rhs)))
  log_det <- k * log(tau2) + 2 * sum(log(diag(R)))
  log_marginal <- -0.5 * log_det + 0.5 * sum(rhs * m)
  list(u = u, b = b, k = k, R = R, m = m,
       log_marginal = log_marginal)
}

birth_probability <- function(k, Kmax) {
  if (k >= Kmax) 0 else if (k <= 1L) 1 else 0.5
}
death_probability <- function(k, Kmax) 1 - birth_probability(k, Kmax)

insert_at <- function(z, value, index) append(z, value, after = index - 1L)

clipped_normal_mean <- function(mu, variance, B = Bclip) {
  s <- sqrt(pmax(variance, 0))
  out <- pmax(-B, pmin(B, mu))
  active <- which(s > 1e-10)
  if (length(active)) {
    a <- (-B - mu[active]) / s[active]
    b <- ( B - mu[active]) / s[active]
    out[active] <- -B * pnorm(a) + B * (1 - pnorm(b)) +
      mu[active] * (pnorm(b) - pnorm(a)) +
      s[active] * (dnorm(a) - dnorm(b))
  }
  out
}

conditional_predict <- function(state, grid, B = Bclip) {
  P <- design_matrix(grid, state$u, state$b)
  mu <- drop(P %*% state$m)
  # Q=R'R, so P Q^{-1} P' has diagonal colSums((R'\P')^2).
  v <- colSums(forwardsolve(t(state$R), t(P))^2)
  clipped_normal_mean(mu, v, B)
}
