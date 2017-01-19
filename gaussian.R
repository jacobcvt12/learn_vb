# input
set.seed(1)
n <- 100000
mu.true <- 12
tau.true <- 0.3
x <- rnorm(n, mu.true, sqrt(1/tau.true))

# priors
mu.mu <- 0
tau.mu <- 1e-5
A.0 <- 0.001
B.0 <- 0.001

# initialize
m <- 0
s.2 <- 1000
A.1 <- 1
B.1 <- 1

# function to calculate ELBO
elbo <- function(s.2, m, tau.mu, mu.mu, n) {
  #0.5 + 0.5 * log(s.2 * tau.mu) - n / 2 * log(2*pi) - 
  #tau.mu / 2 * ((s.2 + m^2)-2 * m * mu.mu + mu.mu ^ 2)
    0.5 - (n / 2) * log(2 * pi) + 0.5 * log(s.2 * tau.mu) -
        tau.mu * ((m - mu.mu) ^ 2 + s.2) / 2
}

# while the ELBO has not converged
i <- 1
elbo.vec <- rep(NA, 10000)
epsilon <- 1e-5
elbo.old <- elbo(s.2, m, tau.mu, mu.mu, n)
elbo.new <- elbo.old + epsilon * 1000
elbo.vec[1] <- elbo.old

while (abs(elbo.new - elbo.old) >= epsilon) {
  # update mean variational density
  denom <- n * A.1 / B.1 + tau.mu
  m <- (n * mean(x) * A.1 / B.1 + tau.mu * mu.mu) / denom
  s.2 <- 1 / denom
  
  # update variance variational density
  A.1 <- A.0 + n / 2
  B.1 <- B.0 + 0.5 * sum(x ^ 2) - m * sum(x) + 0.5 * n * (m^2 + s.2)
  
  # update elbo
  elbo.old <- elbo.new
  elbo.new <- elbo(s.2, m, tau.mu, mu.mu, n)
  i <- i + 1
  elbo.vec[i] <- elbo.new
}

m
A.1 / B.1

elbo.vec <- elbo.vec[!is.na(elbo.vec)]
ts.plot(elbo.vec)