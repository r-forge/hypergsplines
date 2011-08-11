#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[hyp2f1.R] by DSB Mon 09/05/2011 18:07 (CEST)>
##
## Description:
## Interface to the hypergeometric function in C code.
##
## History:
## 13/09/2010    file creation
#####################################################################################


##' Gauss hypergeometric function
##'
##' This is an R interface to the C function taken from the Cephes library.
##' Note that all parameters must be numeric and must not be complex. Only
##' \code{x} may be a vector.
##' 
##' @param a first function parameter
##' @param b second function parameter
##' @param c third function parameter
##' @param x vector of abscissae with absolute value not greater than 1
##' @return the vector of function values
##'
##' @example hypergsplines/examples/hyp2f1.R
##' 
##' @export 
##' @keywords math
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
hyp2f1 <- function(a,
                   b,
                   c,
                   x)
{ 
    ## checks:
    stopifnot(is.numeric(a),
              identical(length(a), 1L),
              is.numeric(b),
              identical(length(b), 1L),
              is.numeric(c),
              identical(length(c), 1L),
              is.numeric(x),
              all(abs(x) <= 1))

    ## call the C++ function
    res <- .External(cpp_hyp2f1,
                     a, b, c, x)

    ## and return the result
    return(res)
}


##' Laplace approximation to the Gauss hypergeometric function
##'
##' This is an R interface to the C function implementing the Laplace
##' approximation to the Gauss hypergeometric function. Note that all parameters
##' must be numeric and must not be complex. Only \code{x} may be a vector.
##' 
##' @param a first function parameter
##' @param b second function parameter
##' @param c third function parameter
##' @param x vector of abscissae with absolute value not greater than 1
##' @return the logarithm of Laplace approximated values
##'
##' @export 
##' @keywords math
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
logHyp2f1Laplace <- function(a,
                             b,
                             c,
                             x)
{ 
    ## checks:
    stopifnot(is.numeric(a),
              identical(length(a), 1L),
              is.numeric(b),
              identical(length(b), 1L),
              is.numeric(c),
              identical(length(c), 1L),
              is.numeric(x),
              all(abs(x) <= 1))

    ## call the C++ function
    res <- .External(cpp_log_hyp2f1_laplace,
                     a, b, c, x)

    ## and return the result
    return(res)
}
