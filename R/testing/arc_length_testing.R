#!/usr/bin/Rscript
#  arc_length_testing.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.01.2018

## This file exists to compare various methods of arc length calculation for GP's

require(microbenchmark)
source('some_funcs.R')


######### Generate some example data
set.seed(123)
n <- 20
p <- 10
X <- matrix(rnorm(n*p), ncol = p)

kern <- kernel_factory()
nugget <- 1
y <- gen_gp(X, kern, nugget)

######### Now let's compute the arc length between points 1 and 2.
### Along the mean predictive surface:
## Using limit definition:
post_mean <- gp_post_mean_factory(X, y, kern, nugget)

lim_def_len <- arclength_limdef(post_mean, X[1,], X[2,], 1e4)
print(lim_def_len)

## Using integral of directional deriv with finite diffs for directional deriv calc
post_mean <- gp_post_mean_factory(X, y, kern, nugget)

int_def_len <- gp_arclen(post_mean, X[1,], X[2,])
print(int_def_len)

## Use Microbenchmark to time the two
microbenchmark(arclength_limdef(post_mean, X[2,], X[2,], 1e4), gp_arclen(post_mean, X[1,], X[2,]))
##Varying evaluations will vary precisions, but it seems like the integral definition is much faster

### Now, calculate the expected arc length
#Some params
a <- X[1,]
b <- X[2,]
d <- b - a
nn <- 1e2#How many points to use along arc

