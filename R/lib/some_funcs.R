#!/usr/bin/Rscript
#  some_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.31.2017

require(Matrix)

source('cubic_interpolant.R')

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
gen_gp <- function(X, kern, nugget) {
    n <- nrow(X)
    #Create a defalut kernel
    if (missing(kern)) {
        kern <- kernel_factory()
    }
    #Create the covariance matrix
    SIGMA <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:n) {
            SIGMA[i,j] <- kern(X[i,],X[j,]) + (i==j) * nugget
        }
    }

    #Create a MVNormal draw from the implied covariance matrix
    GAMMA <- chol(SIGMA)
    z <- rnorm(n)
    y <- t(GAMMA) %*% z

    return(y)
}

#TODO: Not just in 1D, use Convex hulls
#' Generate Data from a multiple GP model
#' @param X The locations where data were observed, an n by p matrix for n observations in p dimensional spcae.
#' @param a1 The end of the first set
#' @param a2 The beginning of the second set
#' @param kern1
#' @param kern2
#' @return A vector of responses of length nrow(X)
# Verdict: issues with larger overwhelming smaller
gen_mult_gp_smooth <- function(X, a1, a2, kern1, kern2) {
    n <- nrow(X)
    #Create the covariance matrix
    #What kernel we use depends on where the points are in space
    SIGMA <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:n) {
            #Get the first point's kernel:
            if (X[i,] <= a1) {
                kl <- kern1(X[i,], X[j,])
            } else if (X[i,] >= a2) {
                kl <- kern2(X[i,], X[j,])
            } else {
                #Prepare interpolant
                cp1 <- a1
                cp2 <- a2
                k1 <- kern1(X[i,], X[j,])
                k2 <- kern2(X[i,], X[j,])
                k1p <- scalar_fd(function(a) kern1(a, X[j,]), cp1)
                k2p <- scalar_fd(function(a) kern2(a, X[j,]), cp2)
                interp <- get_interp(cp1, cp2, k1, k2, k1p, k2p)

                kl <- interp(X[i,])
            }
            #Get the second point's kernel:
            if (X[j,] <= a1) {
                kr <- kern1(X[i,], X[j,])
            } else if (X[j,] >= a2) {
                kr <- kern2(X[i,], X[j,])
            } else {
                #Prepare interpolant
                cp1 <- a1
                cp2 <- a2
                k1 <- kern1(X[i,], X[j,])
                k2 <- kern2(X[i,], X[j,])
                k1p <- scalar_fd(function(a) kern1(X[i,], a), cp1)
                k2p <- scalar_fd(function(a) kern2(X[i,], a), cp2)
                interp <- get_interp(cp1, cp2, k1, k2, k1p, k2p)

                kr <- interp(X[j,])
            }

            #Finally, assign a value to the covariance matrix
            SIGMA[i, j] <- max(kl, kr) + (i==j) * nugget
            #SIGMA[i, j] <- 2/(1/kl + 1/kr) + (i==j) * nugget
        }
    }

    #This SIGMA may not be PSD, so we find the nearest PSD matrix
    #This is done by simply setting any negative eigenvalues to 0
    PSD_approx <- nearPD(SIGMA)$mat

    #Create a MVNormal draw from the implied covariance matrix
    GAMMA <- chol(PSD_approx)
    z <- rnorm(n)
    y <- t(GAMMA) %*% z

    return(y)
}

#TODO: Not just in 1D, use Convex hulls
#' Generate Data from a multiple GP model
#' @param X The locations where data were observed, an n by p matrix for n observations in p dimensional spcae.
#' @param A1 A tuple c(x1, x2), the location of the first set, with x1 < x2.
#' @param A2 A tuple C(x1, x2), the location of the second set, with x1 < x2. Should not overlap with A1
#' @param kern1
#' @param kern2
#' @param kern12
#' @return A vector of responses of length nrow(X)
gen_mult_gp <- function(X, A1, A2, kern1, kern2, kern12) {
    n <- nrow(X)
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

    #This SIGMA may not be PSD, so we find the nearest PSD matrix
    #This is done by simply setting any negative eigenvalues to 0
    PSD_approx <- nearPD(SIGMA)$mat


    #Create a MVNormal draw from the implied covariance matrix
    GAMMA <- chol(PSD_approx)
    z <- rnorm(n)
    y <- t(GAMMA) %*% z

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

#' Get the arclength of line segments connecting a point cloud
#' @param X An n by p matrix representing n many points in p dimensional space. These should be ordered according to the path of which arc length calculation is desired.
#' @return Scalar sum of the lengths (Eucliden 2 norm) of the line segments connecting each of the points in X, ordered by the order of X's rows.
lineseg_arclength <- function(X) {
    ds <- diff(X)
    lens <- apply(ds, 1, function(i) norm(i, '2'))
    return(sum(lens))
}

#' Create a function which will give the posterior GP mean surface
#' @param X The locations where data were observed, an n by p matrix for n observations in p dimensional spcae.
#' @param y The reponse corresponding to the X matrix locations.
#' @param kern The covariance kernel function which induces the covariance matrix
#' @param nugget the scalar nugget effect 
#' @return A function of 1 parameter, representing a single predictive location, which will in turn return the predictive mean at that location
gp_post_mean_factory <- function(X, y, kern, nugget) {
    #Create the covariance matrix
    SIGMA <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        for (j in 1:n) {
            SIGMA[i,j] <- kern(X[i,],X[j,]) + (i==j) * nugget
        }
    }
    
    #Prepare all constant terms (wrt x, the predictive location)
    alpha <- solve(SIGMA, y)
    
    #Prepare function to return
    f <- function(x) sapply(1:nrow(X), function(i) kern(x, X[i,])) %*% alpha
    return(f)
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

    #Check to see if the points are right next to each other
    if (norm(d, '2') < 1e-5) {
        return(0)
    }
    dkern_dir <- function(t) direct_fd(post_mean, a + t * d, d, h)

    #Set up the integrand for arc length
    f <- function(ins) sapply(ins, function(t) norm(d,'2') * sqrt(1 + dkern_dir(t)^2))

    return(integrate(f, 0, 1)$value)
}
