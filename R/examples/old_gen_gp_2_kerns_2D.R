#!/usr/bin/Rscript
#  gen_gp_2_kerns.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.03.2018

## This file generates from two seperate Gaussian Processes and then visualizes the output.
## In one part of the space, the lenthscale is a tenth of what it is in the rest of the space

require(mds.methods)
require(rgl)
source('../lib/some_gp_funcs.R')

#Pick a seed
seed <- 123

#Quickload
load(paste('HD_2gps_save_', seed, '.RData', sep = ''))

######### Generate some data with specified weirdness
set.seed(seed)
n <- 75
p <- 2
X <- matrix(runif(n*p), ncol = p)

## Denote the top s which have the largest inner product with the p-vector of 1's as being in the different kernel space
sp <- 0.25#proportion of things in different kernel space
s <- ceiling(sp * n)
r <- rank(-X %*% rep(1,p))

##Xn -- x normal, obeys the kernel for most of the space
##Xd -- x different, obeys a kernel with a tenth of the lengthscale
Xn <- as.matrix(X[r > s,], ncol = p)
Xd <- as.matrix(X[r <= s,], ncol = p)

#Create the kernels and the response
kernn <- kernel_factory(lengthscale=0.5)
kernd <- kernel_factory(lengthscale=0.1, covariance = 10)
nugget <- 0.01
yn <- gen_gp(Xn, kernn, nugget)
yd <- gen_gp(Xd, kernd, nugget)

#Store them in one vector
y <- rep(NA, n)
y[r > s] <- yn
y[r <= s] <- yd

##For 2D, make a 3D plot of the 2D X's plus the y, and add the GP mean surface to it.
if (p == 2) {
    cols <- c('red', 'blue')
    plot3d(X[,1], X[,2], y, col = cols[(r > s) + 1])

    XX1 <- seq(0,1,length.out=10)
    XX2 <- seq(0,1,length.out=10)
    mu <- gp_post_mean_factory(X, y, kernn, nugget)
    #yy <- sapply(1:nrow(XX), function(xx) mu(XX[xx,]))
    yy <- outer(XX1, XX2, Vectorize(function(i, j) as.numeric(mu(cbind(i, j)))))
    persp3d(XX1, XX2, yy, col = 'green', add = TRUE)
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

##### Compare visualizations using arc length and euclidean distance
post_mean <- gp_post_mean_factory(X, y, kernn, nugget)
gp.dist <- function(a, b) gp_arclen(post_mean, a, b)
euclidean.dist <- function(a, b) norm(b-a, '2')

gparc_mds <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist, 
                               n.inits = 1)
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
text(low_d1[,1], low_d1[,2], 1:n, col = cols[(r > s) + 1])

plot(low_d2, lwd = 0, main = paste('Euclidean MDS', 'stress:',euc_mds$value), lty=0)
text(low_d2[,1], low_d2[,2], 1:n, col = cols[(r > s) + 1])

#Get the two distance matrices and compare them
low_1 <- as.matrix(dist(low_d1))
low_2 <- as.matrix(dist(low_d2))
print(norm(low_1 - low_2))
save.image(paste('2gps_save_', seed, '.RData', sep = ''))
