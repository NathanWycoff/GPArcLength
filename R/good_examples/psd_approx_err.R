#!/usr/bin/Rscript
#  ../dev_examples/psd_approx_err.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.18.2018

source('../lib/some_gp_funcs.R')

## Generate some data
n <- 100
p <- 1
X <- matrix(seq(0,1,length.out=n), ncol = p)

cp <- 0.5

## Some params to set
kern1 <- kernel_factory(lengthscale=0.1)
kern2 <- kernel_factory(lengthscale=0.01)
nugget <- 0.01

#Determine the base kernel for each point
kerns <- list(kern1, kern2)
kern_to_use <- (X > cp) + 1

#Create the covariance matrix
n <- nrow(X)
SIGMA <- matrix(NA, nrow = n, ncol = n)
for (i in 1:n) {
    for (j in 1:n) {
        SIGMA[i, j] <- (kerns[[kern_to_use[i]]](X[i,], X[j,]) * 
                        kerns[[kern_to_use[j]]](X[i,], X[j,]))
        SIGMA[i, j] <- SIGMA[i, j] + (i==j) * nugget
        #Get a little bit bigger of a nugget if we're on the boundry
        if ( (X[i,] - cp) * (X[j,] - cp) < 0)  {
            SIGMA[i, j] <- SIGMA[i, j] + (i==j) * nugget
        }
    }
}

#Sometimes this bad boi won't be PD
PSD_APPROX <- as.matrix(Matrix::nearPD(SIGMA)$mat)
GAMMA <- chol(PSD_APPROX)

#Plot the result sometimes
image(SIGMA[nrow(SIGMA):1,])
image(as.matrix(PSD_APPROX[nrow(PSD_APPROX):1,]))
image(as.matrix((PSD_APPROX - SIGMA)[nrow(PSD_APPROX):1,]))

print(norm(SIGMA - PSD_APPROX))

## The covariance on the border is much less than desired
f <- ceiling(n / 2)
d <- ceiling(n / 2) + 1
SIGMA[f,d]
PSD_APPROX[f,d]

# Generate some data
#Sometimes it's pretty OK
set.seed(1237)
z <- rnorm(n)
y <- t(GAMMA) %*% z

plot(X, y)
abline(v=cp, col = 'red')

#Sometimes it's less OK
set.seed(1236)
z <- rnorm(n)
y <- t(GAMMA) %*% z

plot(X, y)
abline(v=cp, col = 'red')
