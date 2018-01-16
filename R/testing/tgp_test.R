#!/usr/bin/Rscript
#  tgp_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.16.2018

## Test that we can recover model parameters with tgp
require(tgp)
source('../lib/some_gp_funcs.R')

#Generate from a GP with known hyperparams
seed <- sample(1:1000,1)
set.seed(seed)
n <- 200
p <- 1
X <- matrix(runif(n*p), ncol = p)

kern <- kernel_factory(lengthscale=0.1)
nugget <- 0.01
y <- gen_gp(X, kern, nugget)
#Make mean 0 range 1
#y <- y - mean(y)

#Fit a GP to the data and recover params
nn <- 1
XX <- seq(0, 1, length.out = nn)
fit <- bgp(X, y, XX, meanfn = 'constant', corr = 'exp', 
          trace = T, m0r1 = FALSE, zcov = FALSE)
res <- fit$trace$XX[[1]]

#Extract posterior quantities
lengthscale <- mean(res$d)
nugget <- mean(res$nug)
covar <- mean(res$s2)

print(lengthscale)
print(nugget)
print(covar)

#Check if, visually, the predictions are reasonable.
if (p==1) {
    plot(X, y)
    inf_kern <- kernel_factory(lengthscale=lengthscale, covariance = covar)
    inf_mu <- gp_post_mean_factory(X, y, inf_kern, nugget)
    true_mu <- gp_post_mean_factory(X, y, kern, nugget)
    Xs <- sort(X)
    points(Xs, sapply(Xs, inf_mu), col = 'red', type = 'l')
    points(Xs, sapply(Xs, true_mu), col = 'blue', type = 'l')
}

#Looks good!
