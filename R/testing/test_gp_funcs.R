#!/usr/bin/Rscript
#  test_gp_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.04.2018

## Some sanity testing on the basic GP functions

source('../lib/some_gp_funcs.R')
require(geometry)

##A simple 1D example
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

## See if we can get a smooth local GP
n <- 200
p <- 1
nugget <- 0.01

X <- matrix(runif(n*p), ncol = p)

a1 <- 0.45
a2 <- 0.55

kern1 <- kernel_factory(lengthscale=-0.1, covariance = 1)
kern2 <- kernel_factory(lengthscale=0.01, covariance = 1)
y <- gen_mult_gp_smooth(X, a1, a2, kern1, kern2)

par(mfrow=c(1,1))
plot(X, y)
abline(v=c(a1, a2), col = 'red')
