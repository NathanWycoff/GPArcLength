#!/usr/bin/Rscript
#  expected_arclen.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.27.2018

## When looking at arc lengths between points, it may be useful context to have 
## the expected arc length averaged over all response realizations, conditioned
## on sampling locations (X) and GP parameters (theta)
## We're going to try to calculate that here.

source('../lib/some_gp_funcs.R')

######### Prepare to generate some data with specified weirdness
set.seed(seed)
n <- 50
p <- 1
#X <- matrix(runif(n*p), ncol = p)
X <- matrix(seq(0,1,length.out=n), ncol = p)

## Denote the top s which have the largest inner product with the p-vector of 1's as being in the different kernel space
#Compute the line separating the two
sp <- 0.4#proportion of things in different kernel space
s <- ceiling(sp * n)
r <- rank(-X %*% rep(1,p))
o <- r > s
cp <- (min(X[!o,]) + max(X[o,])) / 2

#Pick some kernels
kernn <- kernel_factory(lengthscale=0.1)
kernd <- kernel_factory(lengthscale=0.01)
nugget <- 0.01

## We will now choose 3 pairs of points: two in the correct nugget, two in the
## incorrect nugget, and one from each, and compare expected arc lens from
## the correct as well as incorrect model

## Under the assumed model (a single kernel), what is the expected arclen for each?
iters <- 100
bci <- 1
bcj <- 2
bc_lens_sin <- rep(NA, iters)
bwi <- n-1
bwj <- n
bw_lens_sin <- rep(NA, iters)
oei <- 1
oej <- n
oe_lens_sin <- rep(NA, iters)
for (iter in 1:iters) {
    #Generate data from assumed model (single kernel gp).
    y <- gen_gp(X, kernn, nugget)
    post_mean <- gp_post_mean_factory(X, y, kernn, nugget)
    gp.dist <- function(a, b) gp_arclen(post_mean, a, b)

    #Both in correct area
    bc_lens_sin[iter] <- gp.dist(X[bci,], X[bcj])
    #Both in incorrect area
    bw_lens_sin[iter] <- gp.dist(X[bwi,], X[bwj])
    #One in each
    oe_lens_sin[iter] <- gp.dist(X[oei,], X[oej])
}

## Under the true model (multi kernel), what is the expected arclen for each?
iters <- 100
bci <- 1
bcj <- 2
bc_lens_mult <- rep(NA, iters)
bwi <- n-1
bwj <- n
bw_lens_mult <- rep(NA, iters)
oei <- 1
oej <- n
oe_lens_mult <- rep(NA, iters)
for (iter in 1:iters) {
    #Generate data from assumed model (multiple kernel gp).
    y <- gen_partition_gp(X, cp, kernn, kernd, nugget)
    post_mean <- gp_post_mean_factory(X, y, kernn, nugget)
    gp.dist <- function(a, b) gp_arclen(post_mean, a, b)

    #Both in correct area
    bc_lens_mult[iter] <- gp.dist(X[bci,], X[bcj])
    #Both in incorrect area
    bw_lens_mult[iter] <- gp.dist(X[bwi,], X[bwj])
    #One in each
    oe_lens_mult[iter] <- gp.dist(X[oei,], X[oej])
}

## See how they vary
mean(bc_lens_sin) / mean(bc_lens_mult)
mean(bw_lens_sin) / mean(bw_lens_mult)
mean(oe_lens_sin) / mean(oe_lens_mult)
