#!/usr/bin/Rscript
#  test_gp_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.04.2018

## Some sanity testing on the basic GP functions

source('some_funcs.R')

##A simple 1D example
n <- 500
p <- 1
nugget <- 0.01

X <- matrix(rnorm(n*p), ncol = p)
kern <- kernel_factory(lengthscale=0.1, covariance = 1)
y <- gen_gp(X, kern, nugget)

plot(X, y)
