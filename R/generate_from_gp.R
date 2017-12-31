#!/usr/bin/Rscript
#  generate_from_gp.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.31.2017

#' Create a kernel function to specs
#' @param type What kind of kernel? Currently only 'squared exponential' is allowed (and is implied to be seperable)
#' @param lengthscale For the squared exponential kernel, the separable lengthscale parameter
#' @param covariance For the SE kernel, covariance parameter
kernel_factory <- function(type = 'squared exponential', lengthscale = 1, 
                           covariance = 1) {
    if (type == 'squared exponential')
        return(function(x1, x2) covariance * exp(-sum((x1 - x2)^2) / lengthscale))
    error("kernel type not implemented")
}

#' Generate data from a GP model
#' @param X The locations where data were observed, an n by p matrix for n observations in p dimensional spcae.
#' @param kern The covariance kernel function which induces the covariance matrix
#' @param nugget the scalar nugget effect 
#' @return A vector of responses of length nrow(X)
gen_GP <- function(X, kern, nugget) {
    #Create a defalut kernel
    if (missing(kern)) {
        kern <- kernel_factory
    }
    #Create the covariance matrix
    SIGMA <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:n) {
            SIGMA[i,j] <- kern(X[i,],X[j,]) + nugget
        }
    }

    #Create a MVNormal draw from the implied covariance matrix
    GAMMA <- solve(chol(SIGMA))
    z <- rnorm(n)
    y <- GAMMA %*% z

    return(y)
}

#' Get the parametric function for a line passing through 2 points in arbitrary dimensional space. 
#' @param x1 The "from" point
#' @param x2 The "to" point
#' @return A parametric function returning coordinates along the line. For instance, ret(1) gives the coordinates 1 unit in the direction from x1 to x2.
#' @example
#' #Get the line segment from 1,1 to 2,2
#' f <- conn_line_seg(c(1,1), c(2,2))
#'
#' #Evaluating it at zero should give the first point:
#' f(0)
#' 
#' #But we can evaluate it arbitrarily along that interval
#' f(0.34)
#' 
#' #And even beyond, so be careful:
#' f(20)
conn_line_seg <- function(x1, x2) {
    #Get the normalized direction of the line
    d <- x2 - x1
    d <- d / norm(d,'2')

    return(function(t) x1 + t * d)
}

#' Numerically Calculate arc length of a function along a line - useful for double checking analytic answers. Simply breaks the function into many small lines and sums their lengths.
#' @param target The target function to evaluate, should accept a single, vector-valued argument.
#' @param a The "from" vector 
#' @param b the "to vector 
#' @param n The number of lines to use in the approximation
#' @return The scalar arc length
#' @example 
#' 
#' #Calculate the arc length of (log(sec(x)), log(sec(y))) from 0,0 to pi/4, pi/4
#' #Wolfram Alpha tells me this is 1.3612
#' target <- function(x) -log(cos(x[1])) -log(cos(x[2]))
#' num_arclength(target, c(0,0), c(pi/4, pi/4), 1000)
#' 
#' #Now try the same function from 0,0 to pi/8, pi/4, answer should be 1.0046
#' target <- function(x) -log(cos(x[1])) -log(cos(x[2]))
#' num_arclength(target, c(0,0), c(pi/8, pi/4), 1000)

num_arclength <- function(target, a, b, n) {
    #TODO: Gain speed by reimplementing this first part for higher D
    lens <- rep(0,n)
    h <- norm(b-a, '2') / n
    h2 <- h^2

    #Get the line connecting the points
    l <- conn_line_seg(a, b)

    #Get the lengths of a bunch of little lines.
    for (i in 1:n) {
        t1 <- l(a + (i-1) * h)
        t2 <- l(a + i * h)
        y1 <- target(t1)
        y2 <- target(t2)
        lens[i] <- norm(c(y2 - y1, t2 - t1), '2')
    }
    return(sum(lens))
}
