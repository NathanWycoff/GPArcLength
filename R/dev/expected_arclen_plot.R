#!/usr/bin/Rscript
#  expected_arclen_plot.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.27.2018

## For one multi GP realization, calculate the difference between expected and 
## observed arclens.
## This file expects to be run in /dev/
source('../lib/some_gp_funcs.R')
require(mds.methods)

######### Prepare to generate some data with specified weirdness
seed <- 123
set.seed(seed)
n <- 3
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

# Generate the actual observed data
y_obs <- gen_partition_gp(X, cp, kernn, kernd, nugget)

## Under the assumed model (a single kernel), what is the expected arclen for pair of points?
iters <- 100
exp_al <- matrix(NA, nrow = n, ncol = n)
diag(exp_al) <- 0
for (i in 2:n) {
    for (j in 1:(i-1)) {
        dists <- rep(NA, iters)

        for (iter in 1:iters) {
            #Generate data from assumed model (single kernel gp).
            y <- gen_gp(X, kernn, nugget)
            post_mean <- gp_post_mean_factory(X, y, kernn, nugget)
            gp.dist <- function(a, b) gp_arclen(post_mean, a, b)

            dists[iter] <- gp.dist(X[i,], X[j,])
        }
        exp_al[i,j] <- exp_al[j,i] <- mean(dists)
    }
}

#Calculate the observed arclens
post_mean <- gp_post_mean_factory(X, y_obs, kernn, nugget)
gp.dist <- function(a, b) gp_arclen(post_mean, a, b)
obs_al <- good.dist(X, gp.dist, symm = TRUE)

#Calculate the mean arc len for both expected and observed
mean_obs <- sapply(1:nrow(X), function(i) mean(obs_al[-i,i]))
mean_exp <- sapply(1:nrow(X), function(i) mean(exp_al[-i,i]))
mean_ration <- sapply(1:nrow(X), function(i) mean(exp_al[-i,i] / obs_al[-i,i]))

#Plot the expected vs observed arc lengths.
png('../images/obs_min_expected.png')
plot(mean_exp, mean_obs, xlab = "Mean Expected Arc Length", 
     ylab  = "Mean Observed Arc Length", main = "Arc Length Plot")
abline(0,1,col='red')
dev.off()
