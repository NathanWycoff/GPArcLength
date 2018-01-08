#!/usr/bin/Rscript
#  test_gp_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.04.2018

## Some sanity testing on the basic GP functions

source('some_funcs.R')
require(geometry)

##A simple 1D example
n <- 500
p <- 5
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
X <- matrix(runif(n*p), ncol = p)

kern1 <- kernel_factory(lengthscale=0.1, covariance = 1)
kern2 <- kernel_factory(lengthscale=0.01, covariance = 1)
kern12 <- kernel_factory(lengthscale=0.1, covariance = 1)

A1 <- c(0,0.4)
A2 <- c(0.6, 1.0)

y <- gen_mult_gp(X, A1, A2, kern1, kern2, kern12) 

plot(X, y)
