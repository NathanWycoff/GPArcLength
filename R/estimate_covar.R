#!/usr/bin/Rscript
#  estimate_covar.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.04.2018

#One ought to start with the simple things

##Randomly generate a covariance matrix by creating a gaussian upper triangular matrix
n <- 5
GAMMA <- matrix(0, ncol = 5, nrow = 5)
for (i in 1:n) {
    for (j in i:n) {
        GAMMA[i,j] <- rnorm(1)
    }
}
SIGMA <- t(GAMMA) %*% GAMMA

##Simulate a ton of data 
iters <- 1000
ys <- lapply(1:iters, function(i) t(GAMMA) %*% rnorm(n))

##Try to estimate the covariance matrix
Reduce(function(x1, x2) x1 + x2 %*% t(x2), ys, 0) / iters
