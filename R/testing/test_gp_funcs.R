#!/usr/bin/Rscript
#  test_gp_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.04.2018

## Some sanity testing on the basic GP functions

source('../lib/some_gp_funcs.R')
require(geometry)

##A simple 1D example
set.seed(123)
n <- 500
p <- 1
nugget <- 0.01

##First test that we can generate from a GP
X <- matrix(runif(n*p), ncol = p)
kern <- kernel_factory(lengthscale=0.1, covariance = 1)
y <- gen_gp(X, kern, nugget)

plot(X, y)

## And test that we can make reasonable looking predictions given the truth
XX <- as.matrix(seq(0,1,length.out=200), ncol = p)
mu <- gp_post_mean_factory(X, y, kern, nugget)
yy <- sapply(1:nrow(XX), function(xx) mu(XX[xx,]))
points(XX, yy, col = 'red', type = 'l')

## Now try generating from a local-global GP
set.seed(123)
n <- 500
p <- 1
nugget <- 0.01

X <- matrix(runif(n*p), ncol = p)

kern1 <- kernel_factory(lengthscale=0.1, covariance = 1)
kern2 <- kernel_factory(lengthscale=0.01, covariance = 1)
kern12 <- kernel_factory(lengthscale=1, covariance = 1)

A1 <- c(0,0.5)
A2 <- c(0.5, 1.0)

y <- gen_mult_gp(X, A1, A2, kern1, kern2, kern12) 

plot(X, y)
abline(v=0.5, col = 'red')

#Man I really wish the two met in the middle. Maybe we could try:
## See if we can get a smooth local GP by looking at little islands, and interpolating
# in between
set.seed(123)
n <- 200
p <- 1
nugget <- 0.01

X <- matrix(runif(n*p), ncol = p)
X <- matrix(seq(0,1, length.out=n), ncol = p)

a1 <- 0.4
a2 <- 0.6

kern1 <- kernel_factory(lengthscale=0.1, covariance = 1)
kern2 <- kernel_factory(lengthscale=0.01, covariance = 1)
par(mfrow=c(1,2))
y <- gen_nomansland(X, a1, a2, kern1, kern2)

par(mfrow=c(1,1))
plot(X, y)
abline(v=c(a1, a2), col = 'red')

## See if we can get a smooth local GP by partitioning the space
set.seed(1234)

a <- sample(1:100000, 1)
set.seed(a)
n <- 200
p <- 1
nugget <- 0.01

#X <- matrix(runif(n*p), ncol = p)
X <- matrix(seq(0,1, length.out=n), ncol = p)

cp <- 0.5

kern1 <- kernel_factory(lengthscale=0.1, covariance = 1)
kern2 <- kernel_factory(lengthscale=0.01, covariance = 1)
#par(mfrow=c(1,2))
y <- gen_partition_gp(X, cp, kern1, kern2)

plot(X, y)
abline(v=cp, col = 'red')

kerns <- list(kern1, kern2)
kern_to_use <- (X > cp) + 1

#Create the covariance matrix
n <- nrow(X)
SIGMA <- matrix(NA, nrow = n, ncol = n)
for (i in 1:n) {
    for (j in 1:n) {
        SIGMA[i, j] <- sqrt(kerns[[kern_to_use[i]]](X[i,], X[j,]) * 
        kerns[[kern_to_use[j]]](X[i,], X[j,]))
        SIGMA[i, j] <- SIGMA[i, j] + (i==j) * nugget
    }
}

#Sometimes this bad boi won't be PD
PSD_APPROX <- nearPD(SIGMA)$mat
YB_PSD <- yabois_psd(SIGMA)
GAMMA <- chol(PSD_APPROX)

#Plot the result sometimes
image(SIGMA[nrow(SIGMA):1,])
image(as.matrix(PSD_APPROX[nrow(PSD_APPROX):1,]))

print(norm(SIGMA - PSD_APPROX))

#Sample a bunch of draws from this GP
ones <- rep(NA, iters)
twos <- rep(NA, iters)
iters <- 10000
for (i in 1:iters) {
    z <- rnorm(n)
    y <- t(GAMMA) %*% z
    ones[i] <- y[1]
    twos[i] <- y[2]
}
mean(twos * ones)

plot(X, y)
abline(v=cp, col = 'red')
