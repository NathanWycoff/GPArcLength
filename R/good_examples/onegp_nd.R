#!/usr/bin/Rscript
#  gen_gp_2_kerns.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.03.2018

## This file generates from a standard gaussian process

require(mds.methods)
source('../lib/some_gp_funcs.R')

#Pick a seed
seed <- 1234

#Quickload
load(paste('./data/onegpnd_', seed, '.RData', sep = ''))

######### Generate some data with specified weirdness
set.seed(seed)
n <- 50
p <- 5
X <- matrix(runif(n*p), ncol = p)

#Create the kernels and the response
kern <- kernel_factory(lengthscale=0.1)
nugget <- 0.01
y <- gen_gp(X, kern, nugget)

##For 1D only, plot the points as well as the normal GP fit.
if (p == 1) {
    quartz()
    cols <- c('red', 'blue')
    plot(X, y, lwd=0)
    text(X, y, 1:n, col = cols[(r > s) + 1])
    XX <- as.matrix(seq(0,1,length.out=200), ncol = p)
    mu <- gp_post_mean_factory(X, y, kernn, nugget)
    yy <- sapply(1:nrow(XX), function(xx) mu(XX[xx,]))
    points(XX, yy, col = 'red', type = 'l')
}

##### Compare visualizations using arc length and euclidean distance
post_mean <- gp_post_mean_factory(X, y, kern, nugget)
gp.dist <- function(a, b) gp_arclen(post_mean, a, b)
euclidean.dist <- function(a, b) norm(b-a, '2')

gparc_mds <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist, 
                               n.inits = 1, std = F)
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
text(low_d1[,1], low_d1[,2], 1:n)

plot(low_d2, lwd = 0, main = paste('Euclidean MDS', 'stress:',euc_mds$value), lty=0)
text(low_d2[,1], low_d2[,2], 1:n)

#Get the two distance matrices and compare them
low_1 <- as.matrix(dist(low_d1))
low_2 <- as.matrix(dist(low_d2))
print(norm(low_1 - low_2))

## Save the output
save.image(paste('./data/onegpnd_', seed, '.RData', sep = ''))
