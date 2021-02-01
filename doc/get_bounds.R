## -----------------------------------------------------------------------------
set.seed(1406)
(thetas <- runif(3))
(gammas <- runif(3))

## -----------------------------------------------------------------------------
library(ACEBounds)

get_bounds(thetas = thetas,
           gammas = gammas)

## -----------------------------------------------------------------------------
probs <- c(.83, .05, .11, .01, 
           .88, .06, .05, .01, 
           .72, .05, .10, .13,
           .85, .02, .07, .06)

zetas <- array(probs, dim = c(2, 2, 4),
               dimnames = list(x = c(0, 1),
                               y = c(0, 1),
                               z = c(0, 1, 2, 3)))
zetas

get_bounds(zetas = zetas)

## -----------------------------------------------------------------------------
invalid_probs <- c(.83, .05, .11, .01, 
                   .88, .06, .05, .01, 
                   .72, .05, .20, .03)
invalid_zetas <- array(invalid_probs, dim = c(2, 2, 3),
                       dimnames = list(x = c(0, 1),
                                       y = c(0, 1),
                                       z = c(0, 1, 2)))
invalid_zetas

## ----eval = FALSE-------------------------------------------------------------
#  get_bounds(zetas = invalid_zetas)

## ----echo = FALSE, results="hide"---------------------------------------------
get_bounds(zetas = invalid_zetas, stop = FALSE)

## -----------------------------------------------------------------------------
get_bounds(zetas = invalid_zetas, stop = FALSE)

## -----------------------------------------------------------------------------
get_bounds(zetas = invalid_zetas, stop = FALSE, warning = FALSE)

