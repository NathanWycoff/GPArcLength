#!/usr/bin/Rscript
#  playground.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.31.2017

require(mds.methods)
source('some_funcs.R')

#Generate the points where we're going to observe the response
n <- 20
p <- 10
X <- matrix(rnorm(n*p), ncol = p)

kern <- kernel_factory()
nugget <- 1
y <- gen_GP(X, kern, nugget)
