#!/usr/bin/Rscript
#  frontend_fun.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.08.2018

require(scatterD3)
require(shiny)

#Play around with potential frontends

### ScatterD3 Standalone
n <- 20
X <- matrix(rnorm(n*2), ncol = 2)


scatterD3(x = X[,1], y = X[,2], lasso = TRUE)

### Shiny + ScatterD3

#Run the shiny app
runApp('./GPArcLen/')
