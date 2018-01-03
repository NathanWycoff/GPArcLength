#!/usr/bin/Rscript
#  gen_gp_2_kerns.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.03.2018

## This file generates from a mixture of Gaussian Processes and then visualizes the output.
## In one part of the space, the lenthscale is a tenth of what it is in the rest of the space

require(mds.methods)
source('some_funcs.R')

######### Generate some data with specified weirdness
set.seed(123)
n <- 20
p <- 10
X <- matrix(rnorm(n*p), ncol = p)

## Denote the top s which have the largest inner product with the p-vector of 1's as being in the different kernel space
s <- 5
r <- rank(-X %*% rep(1,p))

##Xn -- x normal, obeys the kernel for most of the space
##Xd -- x different, obeys a kernel with a tenth of the lengthscale
Xn <- as.matrix(X[r > s,], ncol = 1)
Xd <- as.matrix(X[r <= s,], ncol = 1)

#Create the kernels and the response
kernn <- kernel_factory(lengthscale=1)
kernd <- kernel_factory(lengthscale=0.1)
nugget <- 1
yn <- gen_gp(Xn, kernn, nugget)
yd <- gen_gp(Xd, kernd, nugget)

#Store them in one vector
y <- rep(NA, n)
y[r > s] <- yn
y[r <= s] <- yd

######## Visualize them using MDS with arc length distance.
post_mean <- gp_post_mean_factory(X, y, kernn, nugget)
gp.dist <- function(a, b) gp_arclen(post_mean, a, b)

mds_ret1 <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist, 
                               n.inits = 5e3)
mds_ret2 <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist,
                               n.inits = 5e3)

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
