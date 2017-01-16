# input
set.seed(1)
n <- 1000000
K <- 3
prop.true <- c(0.2, 0.17, 0.63)
z.true <- sample(1:3, n, TRUE, prop.true)
mu.true <- c(-5, 0, 5)
x <- rnorm(n, mu.true[z.true], 1)

# priors
sigma.2 <- 10 ^ 2

# initialize
m.k <- c(0, 3, 5)
s2.k <- c(0.1, 0.1, 0.1)
phi <- t(mapply(function(y) exp(m.k * y - (s2.k + m.k ^ 2) / 2), x))
phi <- t(apply(phi, 1, function(x) x / sum(x)))

# function to calculate ELBO
elbo <- function(m, s.2, phi) {
  
}

# while the ELBO has not converged
elbo.diff <- 20
epsilon <- 0.01
#while (elbo.diff >= epsilon) {
for (i in 1:10) {
  # update probability of latent class
  phi <- t(mapply(function(y) exp(m.k * y - (s2.k + m.k ^ 2) / 2), x))
  phi <- t(apply(phi, 1, function(x) x / sum(x)))
  
  # update mixture components
  denom <- 1 / sigma.2 + colSums(phi)
  m.k <- (x %*% phi) / denom
  s2.k <- 1 / denom
}

m.k
s2.k
colSums(phi) / sum(colSums(phi))
