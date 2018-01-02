#!/usr/bin/Rscript
#  generate_from_gp.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.31.2017

#' Create a kernel function to specs
#' @param type What kind of kernel? Currently only 'squared exponential' is allowed (and is implied to be seperable)
#' @param lengthscale For the squared exponential kernel, the separable lengthscale parameter
#' @param covariance For the SE kernel, covariance parameter
#' @example
#' #Generate some data from a very standard kernel
#' set.seed(123)
#' n <- 20
#' p <- 10
#' x <- matrix(rnorm(n*p), ncol = p)
#' 
#' kern <- kernel_factory()
#' nugget <- 1
#' y <- gen_gp(x, kern, nugget)
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

#' Create a function which will give the posterior GP mean surface
#' @param X The locations where data were observed, an n by p matrix for n observations in p dimensional spcae.
#' @param y The reponse corresponding to the X matrix locations.
#' @param kern The covariance kernel function which induces the covariance matrix
#' @param nugget the scalar nugget effect 
#' @return A function of 1 parameter which will in turn return the predictive mean at that location
gp_post_mean_factory <- function(X, y, kern, nugget) {
    #Create the covariance matrix
    SIGMA <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:n) {
            SIGMA[i,j] <- kern(X[i,],X[j,]) + nugget
        }
    }
    
    #Prepare all constant terms (wrt x, the predictive location)
    alpha <- solve(SIGMA, y)
    
    #Prepare function to return
    f <- function(x) sapply(1:nrow(X), function(i) kern(x, X[i,])) %*% alpha
    return(f)
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

#' Calculate the Arc Length along a line for a GP mean surface using finite differenfces
#' to evaluate the directional derivative of the gradient.
#' @param post_mean A function of 1 parameter which evaluates the posterior mean surface at the location specified by its parameter, that is, which returns a scalar.
#' @param a The point "from" which to evaluate the arclen (since arclength is symmetric, interchanging a and b is not of any consequence).
#' @param b The point "to" which to evaluate the arclen (since arclength is symmetric, interchanging a and b is not of any consequence).
#' @param h The spacing used for finite difference evaluation of the post_mean's directional derivative.
#' @return The scalar arclength from a to b along the GP posterior mean surface
gp_arclen <- function(post_mean, a, b, h = 1e-6) {
    #Set up a function to get our desired directional derivative
    d <- b - a
    dkern_dir <- function(t) direct_fd(post_mean, a + t * d, d, h)

    #Set up the integrand for arc length
    f <- function(ins) sapply(ins, function(t) norm(d,'2') * sqrt(1 + dkern_dir(t)^2))

    return(integrate(f, 0, 1))
}

#' Calculate Arc Length using a simplified limit
#gp_lim_simplified <- function(post_mean, a, b, n = 1e5) {
#    d <- b - a
#    h <- norm(d, '2')
#    s <- 0
#    for (i in 1:n) {
#        tim1 <- a + (b - a) / n * i
#        f1 <- post_mean(tim1)
#        f2 <- post_mean(tim1 + 1/n * d) 
#        s <- s + sqrt(1 +  ((f2 - f1) / (h/n))^2)
#    }
#
#    return(h/n * s)
#}
#
##' Calculate Arc Length using a simplified limit and conn_line_segs
#gp_lim_simplified_conn <- function(post_mean, a, b, n = 1e5) {
#    d <- b - a
#    h <- norm(d, '2')
#    l <- conn_line_seg(a, b)
#    f1 <- post_mean(a)
#    s <- 0
#    for (i in 1:n) {
#        f2 <- post_mean(l(h*i/n))
#        s <- s + sqrt(1 +  ((f2 - f1) / (h/n))^2)
#        f1 <- f2
#    }
#
#    return(h/n * s)
#}
#
##' Calculate Arc length using 
