#!/usr/bin/Rscript
#  gen_gp_2_kerns.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.03.2018

## This file generates from a mixture of Gaussian Processes and then visualizes the output.
## In one part of the space, the lenthscale is a tenth of what it is in the rest of the space
## It is different from old_gen_ in that it generates from 1 GP model

require(mds.methods)
source('../lib/some_gp_funcs.R')

#Pick a seed
seed <- 1234

#Quickload
load(paste('./data/partition_save_', seed, '.RData', sep = ''))

######### Generate some data with specified weirdness
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

######## Visualize them using MDS with arc length distance, do it twice to ensure
# convergence.
post_mean <- gp_post_mean_factory(X, y, kernn, nugget)
gp.dist <- function(a, b) gp_arclen(post_mean, a, b)

mds_ret1 <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist, 
                               n.inits = 1)
mds_ret2 <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist,
                               n.inits = 1)

# Compare the output
low_d1 <- mds_ret1$par
low_d2 <- mds_ret2$par

cols <- c('red', 'blue')

par(mfrow=c(2,1))
plot(low_d1, lwd = 0, main = paste('MDS baby', 'stress:',mds_ret1$value), lty=0)
text(low_d1[,1], low_d1[,2], 1:n, col = cols[(r > s) + 1])

plot(low_d2, lwd = 0, main = paste('MDS baby', 'stress:',mds_ret2$value), lty=0)
text(low_d2[,1], low_d2[,2], 1:n, col = cols[(r > s) + 1])

#Get the two distance matrices and compare them
low_1 <- as.matrix(dist(low_d1))
low_2 <- as.matrix(dist(low_d2))
print(norm(low_1 - low_2))

##### Compare visualizations using arc length and euclidean 2 distance
post_mean <- gp_post_mean_factory(X, y, kernn, nugget)
gp.dist <- function(a, b) gp_arclen(post_mean, a, b, rel.tol = 1e1)
euclidean.dist <- function(a, b) norm(b-a, '2')

system.time(E <- good.dist(X, gp.dist))

Rprof()
gparc_mds <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist, 
                               n.inits = 1)
Rprof(NULL)
euc_mds <- smacof_forward_mds(high_d = X, weights = rep(1,p), 
                               dist.func = euclidean.dist,
                               n.inits = 1)

##Plot the GParclength one in polar coords
cart2polar <- function(X) {
    Y <- matrix(NA, ncol = ncol(X), nrow = nrow(X))
    Y[,1] <- sqrt(X[,1]^2 + X[,2]^2)
    Y[,2] <- atan(X[,2] / X[,1])
    return(Y)
}

# Compare the output
low_d1 <- gparc_mds$par
low_d2 <- euc_mds$par

cols <- c('red', 'blue')

par(mfrow=c(1,2))
plot(low_d1, lwd = 0, main = paste('GPArclength MDS', 'stress:',gparc_mds$value), 
     lty=0)
text(low_d1[,1], low_d1[,2], rank(X), col = cols[(r > s) + 1])

plot(low_d2, lwd = 0, main = paste('Euclidean MDS', 'stress:',euc_mds$value), lty=0)
text(low_d2[,1], low_d2[,2], rank(X), col = cols[(r > s) + 1])

#Get the two distance matrices and compare them
low_1 <- as.matrix(dist(low_d1))
low_2 <- as.matrix(dist(low_d2))
print(norm(low_1 - low_2))

save.image(paste('partition_save_', seed, '.RData', sep = ''))
