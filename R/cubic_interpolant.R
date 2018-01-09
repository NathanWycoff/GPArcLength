#!/usr/bin/Rscript
#  cubic_interpolant.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.08.2018

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
