#!/usr/bin/Rscript
#  some_geom_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.09.2018

#Some functions related to geometry in R^k with p-norms

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
arclength_limdef <- function(target, a, b, n) {
    lens <- rep(0,n)
    h <- norm(b-a, '2') / n
    h2 <- h^2

    #Get the line connecting the points
    l <- conn_line_seg(a, b)

    #Get the lengths of a bunch of little lines.
    for (i in 1:n) {
        #t1 <- l(a + (i-1) * h)
        #t2 <- l(a + i * h)
        t1 <- l((i-1) * h)
        t2 <- l(i * h)
        y1 <- target(t1)
        y2 <- target(t2)
        #lens[i] <- norm(c(y2 - y1, t2 - t1), '2')
        lens[i] <- sqrt(h2 + (y2 - y1)^2)
    }
    return(sum(lens))
}

#' Get the arclength of line segments connecting a point cloud
#' @param X An n by p matrix representing n many points in p dimensional space. These should be ordered according to the path of which arc length calculation is desired.
#' @return Scalar sum of the lengths (Eucliden 2 norm) of the line segments connecting each of the points in X, ordered by the order of X's rows.
lineseg_arclength <- function(X) {
    ds <- diff(X)
    lens <- apply(ds, 1, function(i) norm(i, '2'))
    return(sum(lens))
}

#' Calculate the a 1D derivative use (forward) finite differences
#' @param f The function to differentiate
#' @param x The point at which to evaluate the derivative (a scalar)
#' @param h How far to move forards (a scalar)
#' @return The approximate value of the derivative
scalar_fd <- function(f, x, h = 1e-6) {
    num <- f(x + h) - f(x)
    return(num / h)
}

#' Calculate the directional gradient using (forward) finite differences
#' @param f The function to differentiate
#' @param x The point at which to evaluate the derivative (a vector)
#' @param d The direction in which to evaluate the derivative (a vector)
#' @param h How far to move in the d direction (a scalar)
#' @return The approximate value of the derivative
direct_fd <- function(f, x, d, h = 1e-6) {
    #TODO: Require unit vector to avoid constant rescaling
    d <- d / norm(d, '2')
    num <- f(x + h * d) - f(x)
    return(num / h)
}

#' Cubic Function Interpolation 
#'
#' A cubic interpolant to interpolate a function and its derivative at 2 points in 1D.
#' @param gamma1 The first input value
#' @param gamma2 The second input value
#' @param mu1 The first output value
#' @param mu2 The second output value
#' @param mu1p The first derivative value
#' @param mu2p The second output value
#' @return A scalar function giving the interpolating values
#' @example
#' #Points at which to interpolate
#' gamma1 <- 0.4
#' gamma2 <- 0.6
#' 
#' #Function values at those points, as well as derivative values
#' mu1 <- 1
#' mu2 <- 2
#' mu1p <- -1
#' mu2p <- 1
#' 
#' xs <- seq(0.39, 0.61, length.out = 100)
#' polynom <- get_interp()
#' plot(xs, sapply(xs, polynom))
#' abline(v=c(gamma1, gamma2), col = 'red')
#' points(c(gamma1, gamma2), c(mu1, mu2), col='red')
get_interp <- function(gamma1, gamma2, mu1, mu2, mu1p, mu2p) {
    #Create the linear system
    A <- matrix(c(gamma1^3, gamma2^3, 3 * gamma1^2, 3 * gamma2^2,
                  gamma1^2, gamma2^2, 2 * gamma1, 2 * gamma2, 
                  gamma1, gamma2, 1, 1, 
                  1, 1, 0, 0), ncol = 4)
    b <- c(mu1, mu2, mu1p, mu2)

    #Create the interpolating cubic
    coefs <- solve(A, b)
    polynom <- function(x) coefs %*% c(x^3, x^2, x, 1)

    return(polynom)
}
