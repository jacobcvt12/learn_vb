# input
set.seed(1)
n <- 100000
K <- 3
prop.true <- c(0.2, 0.17, 0.63)
z.true <- sample(1:3, n, TRUE, prop.true)
mu.true <- c(-5, 0, 5)
x <- rnorm(n, mu.true[z.true], 1)

# priors
sigma.2 <- 10

# initialize
m.k <- rnorm(K, 0, sigma.2)
s2.k <- c(0.1, 0.1, 0.1)
phi <- t(mapply(function(y) exp(m.k * y - (s2.k + m.k ^ 2) / 2), x))
phi <- t(apply(phi, 1, function(x) x / sum(x)))

# function to calculate ELBO
elbo <- function(m.k, s.2, phi, sigma.2, n) {
  K <- length(m.k)
  
  sum(-0.5 * log(2 * sigma.2 * pi) - (s2.k + m.k ^ 2) / (2 * sigma.2)) +
  n * log(1 / K) +
  sum(rowSums(-0.5 * phi * log(2 * pi) - 
              phi * (x ^ 2) / 2 + 
              t(t(phi * x) * as.vector(m.k)) - 
              t(t(phi) * as.vector(s2.k + m.k ^ 2)) / 2)) -
  sum(rowSums(phi * log(phi))) -
  sum(-0.5 - 0.5 * log(2 * s2.k * pi))
}

# while the ELBO has not converged
i <- 1
elbo.vec <- rep(NA, 10000)
epsilon <- 0.01
elbo.old <- elbo(m.k, s.2, phi, sigma.2, n)
elbo.new <- elbo.old + epsilon * 1000
elbo.vec[1] <- elbo.old

while ((elbo.new - elbo.old) >= epsilon) {
  # update probability of latent class
  phi <- t(mapply(function(y) exp(m.k * y - (s2.k + m.k ^ 2) / 2), x))
  phi <- t(apply(phi, 1, function(x) x / sum(x)))
  
  # update mixture components
  denom <- 1 / sigma.2 + colSums(phi)
  m.k <- (x %*% phi) / denom
  s2.k <- 1 / denom
  
  # update elbo
  elbo.old <- elbo.new
  elbo.new <- elbo(m.k, s.2, phi, sigma.2, n)
  i <- i + 1
  elbo.vec[i] <- elbo.new
}

m.k
s2.k
colSums(phi) / sum(colSums(phi))

elbo.vec <- elbo.vec[!is.na(elbo.vec)]
ts.plot(elbo.vec)
