#!/usr/bin/Rscript
#  predict_2_kernels.r Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.08.2018

require(mds.methods)
source('../../lib/some_gp_funcs.R')

## Seeing what predictions looked like on the nonsmooth multiscale GP

######### Generate some data with specified weirdness
set.seed(1234)
n <- 50
p <- 1
X <- matrix(runif(n*p), ncol = p)
#X <- matrix(seq(0,1,length.out=n), ncol = p)

## Denote the top s which have the largest inner product with the p-vector of 1's as being in the different kernel space
sp <- 0.25#proportion of things in different kernel space
s <- ceiling(sp * n)
r <- rank(-X %*% rep(1,p))

##Xn -- x normal, obeys the kernel for most of the space
##Xd -- x different, obeys a kernel with a tenth of the lengthscale
Xn <- as.matrix(X[r > s,], ncol = 1)
Xd <- as.matrix(X[r <= s,], ncol = 1)

#Create the kernels and the response
kernn <- kernel_factory(lengthscale=0.1)
kernd <- kernel_factory(lengthscale=0.01, covariance = 10)
nugget <- 0.01
yn <- gen_gp(Xn, kernn, nugget)
yd <- gen_gp(Xd, kernd, nugget)

#Store them in one vector
y <- rep(NA, n)
y[r > s] <- yn
y[r <= s] <- yd

#Only works in 1D so far.
funky_kern <- function(x, y) {
    cp1 <- max(Xn)
    cp2 <- min(Xd)
    if (x <= cp1) {
        return(kernn(x, y))
    } else if (x >= cp2) {
        return(kernd(x, y))
    } else {
        l <- conn_line_seg(c(cp1, kernn(x, y)), c(cp2, kernd(x, y)))
        return(l(x - cp1)[2])
    }
}

##For 1D only, plot the points as well as the normal GP fit.
if (p == 1) {
    quartz()
    cols <- c('red', 'blue')
    plot(X, y, lwd=0)
    text(X, y, 1:n, col = cols[(r > s) + 1])
    XX <- as.matrix(seq(0,1,length.out=200), ncol = p)
    mu <- gp_post_mean_factory(X, y, funky_kern, nugget)
    yy <- sapply(1:nrow(XX), function(xx) mu(XX[xx,]))
    points(XX, yy, col = 'red', type = 'l')
}

