#!/usr/bin/Rscript
#  some_gp_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.09.2018

## Some functions related to Gaussian Processes

source('../lib/some_geom_funcs.R')

cool_image <- function(X) {
    image(X[nrow(X):1,])
}

#' Create a kernel function to specs
#' @param type What kind of kernel? Currently only 'squared exponential' is allowed (and is implied to be seperable)
#' @param lengthscale For the squared exponential kernel, the separable lengthscale parameter
#' @param covariance For the SE kernel, covariance parameter
#' @param alpha For the RQ kernel, the relative large scale weighting
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
                           covariance = 1, alpha = 1) {
    #Verify params in proper region
    if (lengthscale <= 0 || covariance <= 0 || alpha <= 0) {
        stop("Lengthscale, covariance, and alpha ought to be positive numbers")
    }
    if (type == 'squared exponential')
        return(function(x1, x2) covariance * exp(-norm(x1 - x2, '2')^2 / lengthscale))
    else if (type == 'rational quadratic') {
        return(function(x1, x2)  {
               dist_part <- norm(x1 - x2, '2')^2 / (2 * alpha * lengthscale^2)
               return(covariance * (1 + dist_part)^(-alpha))
           })
    }
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
#' Generate Data from a multiple GP model with space between the local zones
#' i.e., zones do not form a partition of the space
#' @param X The locations where data were observed, an n by p matrix for n observations in p dimensional spcae.
#' @param a1 The end of the first set
#' @param a2 The beginning of the second set
#' @param kern1
#' @param kern2
#' @return A vector of responses of length nrow(X)
# Verdict: issues with larger overwhelming smaller
gen_nomansland <- function(X, a1, a2, kern1, kern2) {
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
            #SIGMA[i, j] <- max(kl, kr) + (i==j) * nugget
            SIGMA[i, j] <- (kr * kl) + (i==j) * nugget
        }
    }

    #This SIGMA may not be PSD, so we find the nearest PSD matrix
    #This is done by simply setting any negative eigenvalues to 0
    PSD_APPROX <- Matrix::nearPD(SIGMA)$mat
    #PSD_APPROX <- SIGMA

    #Plot the result sometimes
    image(SIGMA[nrow(SIGMA):1,])
    image(as.matrix(PSD_APPROX[nrow(PSD_APPROX):1,]))

    #Create a MVNormal draw from the implied covariance matrix
    GAMMA <- chol(PSD_APPROX)
    z <- rnorm(n)
    y <- t(GAMMA) %*% z

    return(y)
}

#TODO: Not just in 1D, not just 2 sections
#' Generate data from a multiple GP model
#'
#' Generate data from multiple GP model, the space being partitioned into sections.
gen_partition_gp <- function(X, cp, kern1, kern2) {
    #Determine the base kernel for each point
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

    print(norm(SIGMA - PSD_APPROX))

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
    PSD_approx <- Matrix::nearPD(SIGMA)$mat


    #Create a MVNormal draw from the implied covariance matrix
    GAMMA <- chol(PSD_approx)
    z <- rnorm(n)
    y <- t(GAMMA) %*% z

    return(y)
}

#' Create a function which will give the posterior GP mean surface
#' @param X The locations where data were observed, an n by p matrix for n observations in p dimensional spcae.
#' @param y The reponse corresponding to the X matrix locations.
#' @param kern The covariance kernel function which induces the covariance matrix
#' @param nugget the scalar nugget effect 
#' @return A function of 1 parameter, representing a single predictive location, which will in turn return the predictive mean at that location
gp_post_mean_factory <- function(X, y, kern, nugget) {
    n <- nrow(X)
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

#' Calculate the Arc Length along a line for a GP mean surface using finite differenfces
#' to evaluate the directional derivative of the gradient.
#' @param post_mean A function of 1 parameter which evaluates the posterior mean surface at the location specified by its parameter, that is, which returns a scalar.
#' @param a The point "from" which to evaluate the arclen (since arclength is symmetric, interchanging a and b is not of any consequence).
#' @param b The point "to" which to evaluate the arclen (since arclength is symmetric, interchanging a and b is not of any consequence).
#' @param h The spacing used for finite difference evaluation of the post_mean's directional derivative.
#' @param ... Additional parameters to be passed to the 'integrate' function to calculate arc length, e.g. rel.tol. 
#' @return The scalar arclength from a to b along the GP posterior mean surface
gp_arclen <- function(post_mean, a, b, h = 1e-6, ...) {
    #Set up a function to get our desired directional derivative
    d <- b - a

    #Check to see if the points are right next to each other
    if (norm(d, '2') < 1e-5) {
        return(0)
    }
    dkern_dir <- function(t) direct_fd(post_mean, a + t * d, d, h)

    #Set up the integrand for arc length
    f <- function(ins) sapply(ins, function(t) norm(d,'2') * sqrt(1 + dkern_dir(t)^2))

    return(integrate(f, 0, 1, ...)$value)
}
