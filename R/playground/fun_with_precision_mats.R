#!/usr/bin/Rscript
#  fun_with_precision_mats.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.09.2018

##What happens if we specify precision matrices instead of covariance matrices in GP's?

source('../lib/some_gp_funcs.R')

###### See if old sam knows what he's talking about regarding the "screen effect"
n <- 10
p <- 1
x <- matrix(seq(0,1,length.out=n), ncol = p)

kern <- kernel_factory()
nugget <- 1

y <- gen_gp(x, kern, nugget)

plot(x, y)


n <- nrow(x)
#Create a default kernel
if (missing(kern)) {
    kern <- kernel_factory()
}
#Create the covariance matrix
SIGMA <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
    for (j in 1:n) {
        SIGMA[i,j] <- kern(x[i,],x[j,]) + (i==j) * nugget
    }
}

print(solve(SIGMA))
plot(x, y)
#Doesn't look right to me

##### Try to specify a precision function
n <- 100
p <- 1
x <- matrix(runif(n*p), ncol = p)

kern <- function(a, b) 2*(-exp(-norm(a-b, '2') / 2) + 0.5)
nugget <- 100

#Create the covariance matrix
OMEGA <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
    for (j in 1:n) {
        OMEGA[i,j] <- kern(x[i,],x[j,]) + (i==j) * nugget
    }
}

SIGMA <- solve(OMEGA)

#Create a MVNormal draw from the implied covariance matrix
GAMMA <- chol(SIGMA)
z <- rnorm(n)
y <- t(GAMMA) %*% z

plot(x, y)
