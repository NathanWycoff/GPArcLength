#!/usr/bin/Rscript
#  motorcycle_download.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.13.2018

#Download the motorcycle data from 
###        Silverman,  B.W.  (1985)  Some  aspects  of   the   spline
###        smoothing   approach   to  non-parametric  curve  fitting.
###        Journal of the Royal Statistical Society, B, 47, 1-52.
#Thanks to some friends at CMU

url <- 'http://www.stat.cmu.edu/~larry/all-of-statistics/=data/motor.dat'
data <- read.table(url, header = TRUE)
motorsim <- data[,c('times', 'accel')]
save(motorsim, file='motorscyle.RData')
