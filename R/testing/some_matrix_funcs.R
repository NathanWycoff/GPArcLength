#!/usr/bin/Rscript
#  some_matrix_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.10.2018

## Some functions related to matrix theory

#' Yaboi's PSD Approx
#'
#' Get yaboi's Positive Semidefinite Approximation to an aribtrary matrix by setting negative eigenvalues to 0
#' Might minimize the frobenius norm or something like that
yabois_psd <- function(X) {
    er <- eigen(X)
    vals <- er$values
    vals[vals < 1e-5] <- 1e-5
    return(er$vectors %*% diag(vals) %*% t(er$vectors))
}

#' Returns the average squared difference along the middle area
max_abs_diff <- function(SIGMA, approx) {
    h <- nrow(SIGMA) / 2
    r <- c(h-1, h+1)
    mean((SIGMA[r, r] - approx[r, r])^2)
}

#' Returns a matrix formed by using eigenvalus "vals" with eigenvectors "vecs"
mat_given_eigenvals <- function(vecs, vals) {
    vecs %*% diag(vals) %*% t(vecs)
}


#' Yaboi's Other PSD Approx
#'
#' Get yaboi's other PSD approx to an arbitrary symmetric matrix by minimizing the difference of the max element of each matrix, keeping the eigenvalues the same
absmax_psd <- function(SIGMA) {
    ed <- eigen(SIGMA)
    init_vals <- ed$values
    init_vals[init_vals<0] <- 1e-5

    my_cost <- function(vals) max_abs_diff(SIGMA, mat_given_eigenvals(ed$vectors, vals))

    res <- optim(init_vals, my_cost, lower = 1e-5)

    absmax_approx <- ed$vectors %*% diag(res$par) %*% t(ed$vectors)
    return(absmax_approx)
}
