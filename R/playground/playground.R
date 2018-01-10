#!/usr/bin/Rscript
#  playground.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.31.2017

require(mds.methods)
source('../lib/some_gp_funcs.R')

#Generate the points where we're going to observe the response
n <- 200
p <- 1
X <- matrix(runif(n), ncol = p)

kern <- kernel_factory('rational quadratic', lengthscale = 0.1, alpha = 1e-1)
nugget <- 0.01
y <- gen_gp(X, kern, nugget)

plot(X, y)
