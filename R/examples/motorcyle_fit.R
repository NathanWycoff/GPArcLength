#!/usr/bin/Rscript
#  motorcyle_fit.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.13.2018

## Apply the method to the motorcycle data from:
###        Silverman,  B.W.  (1985)  Some  aspects  of   the   spline
###        smoothing   approach   to  non-parametric  curve  fitting.
###        Journal of the Royal Statistical Society, B, 47, 1-52.
load('../data/motorscyle.RData')

source('../lib/some_gp_funcs.R')

require(tgp)
require(mds.methods)

## Make the predictor data lie in [0,1]
motorsim$times <- (motorsim$times - min(motorsim$times)) / max(motorsim$times)

#Fit a GP to the data and get predictions
nn <- 10
XX <- with(motorsim, seq(min(times), max(times), length.out = nn))
fit <- with(motorsim, bgp(times, accel, XX, meanfn = 'constant', corr = 'exp', 
          trace = T)) 
res <- fit$trace$XX[[1]]

#Extract posterior quantities
lengthscale <- mean(res$d)
nugget <- mean(res$nug)
covar <- mean(res$s2)

#Create the GP with the point estimates
kern <- kernel_factory(lengthscale = lengthscale, covariance = covar)

X <- with(motorsim, matrix(times, ncol = 1))
y <- motorsim$accel
mu <- gp_post_mean_factory(X, y, kern, nugget)
gp_dist <- function(x1, x2) gp_arclen(mu, x1, x2, h = 1e-6)

###Plot the predictions, with some specified cutoff points
cutoffs <- c(0.18, 0.7)
cols <- c('red', 'black', 'blue')
zones <- with(motorsim, (times > cutoffs[1]) + (times > cutoffs[2]) + 1) 

with(motorsim, plot(times, accel, col = cols[zones]))
nn <- 1000
XX <- with(motorsim, seq(min(times), max(times), length.out = nn))
points(XX, sapply(as.numeric(XX), mu), type= 'l')
abline(v=cutoffs[1], col = 'red')
abline(v=cutoffs[2], col = 'blue')

### Do MDS usng Arclength
low_arclen <- smacof_forward_mds(high_d = X, weights = rep(1, 1), 
                                 dist.func = gp_dist, n.inits = 1)
low_euc <- smacof_forward_mds(high_d = X, weights = rep(1, 1), 
                                 dist.func = euclidean.dist, n.inits = 1)
