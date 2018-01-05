#!/usr/bin/Rscript
#  viz_testing.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.02.2018

#Try visualizing some very basic data using arc length

require(mds.methods)
source('some_funcs.R')

######### Generate some example data
set.seed(123)
n <- 5
p <- 10
X <- matrix(rnorm(n*p), ncol = p)

kern <- kernel_factory()
nugget <- 1
y <- gen_gp(X, kern, nugget)

######## Visualize them using MDS with arc length distance.
post_mean <- gp_post_mean_factory(X, y, kern, nugget)
gp.dist <- function(a, b) gp_arclen(post_mean, a, b)

mds_ret1 <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist, 
                               n.inits = 1e1)
mds_ret2 <- smacof_forward_mds(high_d = X, weights = rep(1,p), dist.func = gp.dist,
                               n.inits = 1e1)

# Compare the output
low_d1 <- mds_ret1$par
low_d2 <- mds_ret2$par

par(mfrow=c(2,1))
plot(low_d1, lwd = 0, main = paste('MDS baby', 'stress:',mds_ret1$value), lty=0)
text(low_d1[,1], low_d1[,2], 1:n)

plot(low_d2, lwd = 0, main = paste('MDS baby', 'stress:',mds_ret2$value), lty=0)
text(low_d2[,1], low_d2[,2], 1:n)

#Get the two distance matrices and compare them
low_1 <- as.matrix(dist(low_d1))
low_2 <- as.matrix(dist(low_d2))
print(norm(low_1 - low_2))

