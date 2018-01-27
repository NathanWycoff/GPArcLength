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

## To load previously run results
load('motorcycle_save.RData')

#TODO: REMOVE
#Take away most of the observations
motorsim <- motorsim[seq(1,nrow(motorsim), by = 3),]

## Make the predictor data lie in [0,1]
motorsim$times <- (motorsim$times - min(motorsim$times)) / max(motorsim$times)

#Fit a GP to the data and get predictions
nn <- 10
XX <- with(motorsim, seq(min(times), max(times), length.out = nn))
fit <- with(motorsim, bgp(times, accel, XX, meanfn = 'constant', corr = 'exp', 
          trace = T, m0r1 = FALSE)) 
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

with(motorsim, plot(times, accel, lwd = 0))
with(motorsim, text(times, accel, 1:length(times), col = cols[zones]))
nn <- 1000
XX <- with(motorsim, seq(min(times), max(times), length.out = nn))
points(XX, sapply(as.numeric(XX), mu), type= 'l')
abline(v=cutoffs[1], col = 'red')
abline(v=cutoffs[2], col = 'blue')

### Do MDS usng Arclength
Rprof()
low_arclen <- smacof_forward_mds(high_d = X, weights = rep(1, 1), 
                                 dist.func = gp_dist, n.inits = 1)
Rprof(NULL)

euclidean.dist <- function(x1, x2) norm(x1-x2, '2')
low_euc <- smacof_forward_mds(high_d = X, weights = rep(1, 1), 
                                 dist.func = euclidean.dist, n.inits = 1)

# Compare the output
low_d1 <- low_arclen$par
low_d2 <- low_euc$par

cols <- c('red', 'black', 'blue')

par(mfrow=c(1,2))
plot(low_d1, lwd = 0, main = paste('GPArclength MDS', 'stress:',low_arclen$value), 
     lty=0)
text(low_d1[,1], low_d1[,2], rank(X), col = cols[zones])

plot(low_d2, lwd = 0, main = paste('Euclidean MDS', 'stress:',low_euc$value), lty=0)
text(low_d2[,1], low_d2[,2], rank(X), col = cols[zones])

save.image('motorcycle_save.RData')
