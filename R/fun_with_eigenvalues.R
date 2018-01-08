#!/usr/bin/Rscript
#  fun_with_eigenvalues.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.08.2018

source('some_funcs.R')

#Mess around with making arbitrary symmetric matrices PSD matrices
n <- 300
p <- 2
nugget <- 1
X <- matrix(runif(n*p), ncol = p)

kern1 <- kernel_factory(lengthscale=0.1, covariance = 1)
kern2 <- kernel_factory(lengthscale=0.01, covariance = 1)
kern12 <- kernel_factory(lengthscale=0.1, covariance = 1)

A1 <- c(0,0.4)
A2 <- c(0.6, 1.0)

#Create the covariance matrix
#What kernel we use depends on where the points are in space
SIGMA <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
    for (j in 1:n) {
        if (A1[1] <= X[i,]  && X[i,] <= A1[2] && A1[1] <= X[j,] && X[j,] <= A1[2]) {
            #If both points are in A1
            SIGMA[i,j] <- kern1(X[i,],X[j,]) + (i==j) * nugget
        } else if (A2[1] <= X[i,]  && X[i,] <= A2[2] && A2[1] <= X[j,] && X[j,] <= A2[2]) {
            #If both points are in A2
            SIGMA[i,j] <- kern2(X[i,],X[j,]) + (i==j) * nugget
        } else {
            #Otherwise, use the global link
            SIGMA[i,j] <- kern12(X[i,],X[j,]) + (i==j) * nugget
        }
    }
}

#Make a PD projection of that matrix
ed <- eigen(SIGMA)
Q <- ed$vectors
ev <- ed$values
ev[ev < 0] <- 0
L <- diag(ev)

PSD_approx <- Q %*% L %*% solve(Q)
eigen(PSD_approx)$values
