#!/usr/bin/Rscript
#  test_mat_func.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.10.2018

## Give our matrix functions a test drive

source('../lib/some_matrix_funcs.R')
source('../lib/some_gp_funcs.R')
require(Matrix)

n <- 100
p <- 1
nugget <- 0.01

set.seed(123)
#X <- matrix(runif(n*p), ncol = p)
X <- matrix(seq(0,1, length.out=n), ncol = p)

cp <- 0.5

kern1 <- kernel_factory(lengthscale=0.1, covariance = 1)
kern2 <- kernel_factory(lengthscale=0.01, covariance = 1)
#par(mfrow=c(1,2))

kerns <- list(kern1, kern2)
kern_to_use <- (X > cp) + 1

#Create the covariance matrix
n <- nrow(X)
SIGMA <- matrix(NA, nrow = n, ncol = n)
for (i in 1:n) {
    for (j in 1:n) {
        SIGMA[i, j] <- sqrt(kerns[[kern_to_use[i]]](X[i,], X[j,]) * 
            kerns[[kern_to_use[j]]](X[i,], X[j,]))
        SIGMA[i, j] <- SIGMA[i, j] + (i==j) * nugget
    }
}

AM_APPROX <- absmax_psd(SIGMA)
FROB_APPROX <- as.matrix(nearPD(SIGMA)$mat)

par(mfrow=c(1,2))
image(abs(AM_APPROX - SIGMA), main = 'mine')
image(abs(FROB_APPROX - SIGMA), main = 'theirs')
