library(mvtnorm)

h <- function(params, x)
{
  V <- params[1]; sigma2 <- params[2]; p <- params[3]
  1/(1 + (1-p)/p * sqrt(sigma2/sigma2+V) * exp(x^2*V/(2*sigma2*(sigma2+V))))
}

log_prior_p <- function(p, alpha)
{
  log(alpha + 1) + alpha * log(p)
}


log_prior_sv <- function(sigma2, V)
{
  -2*log(sigma2 + V)
}

log_pi <- function(V, sigma2, p, x, alpha = 10)
{
  sum(log(p/sqrt(2*pi*sigma2) * exp(-x*x/(2*sigma2)) + (1-p)/sqrt(2*pi*(sigma2+V)) * exp(-x*x/(2*(sigma2+V))))) +
    log_prior_sv(sigma2, V) +
    log_prior_p(p, alpha)
}

log_pi_star <- function(params, data)
{
  xi <- params[1]; eta <- params[2]; lambda <- params[3]
  log_pi(exp(xi), exp(eta), 1/(1+exp(-lambda)), data) + xi + eta + lambda - 2*log(1 + exp(lambda))
}

calc_posterior <- function(samples, weights, data, n)
{
  samps_trans <- cbind(exp(samples[ ,1]),
                       exp(samples[ ,2]),
                       1/(1+exp(-samples[ ,3])))
  post_p <- 1 - apply(samps_trans, 1, function(s) h(s, data)) %*% weights / sum(weights)
  return(post_p[1:n])
}

# inputs:
#     n         number of "noise" points to generate in addition to the 10 signals
#     signals   vector of the 10 signal points
#     a         scale tuning parameter
run <- function(n, signals, a = 5)
{
  # generate noise points
  x <- c(signals, rnorm(n, 0, 1))
  M <- length(x)
  # numerically find the posterior mode (xi, eta, lambda)
  o <- optim(c(0, 0, 0), log_pi_star, data = x, hessian = T, control = list(fnscale = -1))
  post_mode <- o$par
  # also keep track of the inverse Hessian H^{-1}
  hess_inv <- solve(-o$hessian)
  # generate M samples from t3 centered at the mode
  samps <- rmvt(M, sigma = a*hess_inv/3, df = 3, delta = post_mode, type = 'shifted')
  # calculate the log weights
  w <- apply(samps, 1, function(s) log_pi_star(s, x)) -
    apply(samps, 1, function(s) dmvt(s, delta = post_mode, sigma = a*hess_inv/3, log = TRUE, type = 'shifted'))
  # shift the log weights to prevent underflow, and exponentiate them
  w2 <- exp(w - max(w))
  # given the samples and the weights, calculate posterior inclusion probabilities
  calc_posterior(samps, w2, x, length(signals))
}

set.seed(597)
a <- 5
n <- c(25, 100, 500, 5000)
signalobs <- c(-5.65, -5.56, -2.62, -1.20, -1.01, -0.90, -0.15, 1.65, 1.94, 3.57)
results <- round(t(sapply(n, function(nn) run(nn, signalobs, a))), 2)
results <- as.data.frame(results)
names(results) <- signalobs
rownames(results) <- n
knitr::kable(results)