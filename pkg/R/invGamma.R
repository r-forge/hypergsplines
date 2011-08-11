#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
## 
## Time-stamp: <[invGamma.R] by DSB Die 14/09/2010 14:28 (CEST)> 
##
## Description:
## Functions for the inverse gamma distribution.
##
## History:
## 14/09/2010   copy from dsb package, and add Roxygen chunks.
#####################################################################################

##' Density function for the inverse gamma distribution
##'
##' @param x vector of values in the positive support
##' @param a shape parameter
##' @param b rate parameter
##' @param log logical; if \code{TRUE}, the log density is returned
##' @param normalize normalize the density function? (default)
##' @return the density values at \code{x}
##' 
##' @export 
##' @keywords distribution
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
dinvGamma <- function (x,
                       a,
                       b,
                       log = FALSE,
                       normalize = TRUE)
{
    ret <- - (a + 1) * log (x) - b / x
    if (normalize)
        ret <- ret + a * log (b) - lgamma (a)
    if (log)
        return (ret)
    else
        return (exp (ret))
}


##' Cumulative distribution function for the inverse gamma distribution
##'
##' @param q vector of quantiles in the positive support
##' @param a shape parameter
##' @param b rate parameter
##' @param lower.tail logical; if \code{TRUE} (default), probabilities
##' are P[X <= x], otherwise, P[X > x]
##' @param log.p logical; if \code{TRUE}, the log probabilities are returned 
##' @return the cdf values at \code{q}
##' 
##' @export 
##' @keywords distribution
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
pinvGamma <- function (q,
                       a,
                       b,
                       lower.tail = TRUE,
                       log.p = FALSE)
{
    pgamma (q = 1 / q,
            shape = a,
            rate = b,
            lower.tail = ! lower.tail,
            log.p = log.p) 
}


##' Quantile function for the inverse gamma distribution
##'
##' @param p vector of (log) probabilities
##' @param a shape parameter
##' @param b rate parameter
##' @param lower.tail logical; if \code{TRUE} (default), probabilities
##' \code{p} are P[X <= x], otherwise, P[X > x]
##' @param log.p logical; if \code{TRUE}, \code{p} are the log probabilities 
##' @return the quantiles values
##' 
##' @export 
##' @keywords distribution
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
qinvGamma <- function (p,
                       a,
                       b,
                       lower.tail = TRUE,
                       log.p = FALSE)
{
    1 / qgamma (p = p,
                shape = a,
                rate = b,
                lower.tail = ! lower.tail,
                log.p = log.p)
}


##' Generate random numbers from the inverse gamma distribution
##'
##' @param n number of observations. If \code{length(n) > 1},
##' the length is taken to be the number required
##' @param a shape parameter
##' @param b rate parameter
##' @return the random variates
##' 
##' @export 
##' @keywords distribution
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
rinvGamma <- function (n,
                       a,
                       b)
{
    1 / rgamma (n,
                shape = a,
                rate = b)
}

