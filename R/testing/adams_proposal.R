#!/usr/bin/Rscript
#  adams_proposal.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.16.2018

## Adam's proposal on how to sample from a GP with changing lengthscale.
## First, sample independently from 2 GP's
## Then, get the posterior predictive distribution and weigh by predictive accuracy

source('../lib/some_gp_funcs.R')

######### Generate some data with specified weirdness
seed <- sample(1:10000, 1)
seed <- 7027
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

y <- gen_partition_gp(X, cp, kernn, kernd)

##For 1D only, plot the points as well as the normal GP fit.
if (p == 1) {
    quartz()
    cols <- c('red', 'blue')
    plot(X, y, lwd=0)
    text(X, y, 1:n, col = cols[o + 1])
    XX <- as.matrix(seq(0,1,length.out=200), ncol = p)
    mu <- gp_post_mean_factory(X, y, kernn, nugget)
    yy <- sapply(1:nrow(XX), function(xx) mu(XX[xx,]))
    points(XX, yy, col = 'red', type = 'l')
    abline(v=cp, col = 'red')
}

