#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[GPrior-classes.R] by DSB Don 20/09/2012 16:58 (CEST)>
##
## Description:
## Hopefully a clean class system for the different priors on g.
##
## History:
## 20/09/2012   - adapt from R-package "glmBfp",
##              - add HypergnPrior
#####################################################################################

## ----------------------------------------------------------------------------

##' The virtual g-prior class
##'
##' This is the virtual g-prior class from which other g-prior classes derive.
##' The slots are:
##' \describe{
##' \item{logDens}{the prior log density}
##' }
##'
##' @seealso \code{\linkS4class{HypergPrior}}, \code{\linkS4class{InvGammaGPrior}},
##' \code{\linkS4class{CustomGPrior}}
##' 
##' @name GPrior-class
##' @keywords classes internal
setClass(Class="GPrior",
         representation=
         representation(logDens="function"),
         contains=list("VIRTUAL"),
         validity=
         function(object){
             ## check that the exp of the log density function
             ## is a valid density function 

             XMIN <- .Machine$double.xmin
             EPS <- sqrt(.Machine$double.eps)
        
             integrand <- function(g) exp(object@logDens(g))
             integral <- integrate(f=integrand,
                                   lower=XMIN,
                                   upper=Inf)
        
             if(integral$message != "OK")
             {
                 return(integral$message)
             }
             else 
             {
                 if(abs(integral$value - 1) > EPS)
                 {
                     warning("density must be proper and normalized: (numerical) integral is ",
                             integral$value)
                 }

                 return(TRUE)
             }})

## ----------------------------------------------------------------------------

##' The hyper-g prior class
##'
##' The slots are:
##' \describe{
##' \item{a}{the hyperparameter}
##' }
##'
##' @seealso the constructor \code{\link{HypergPrior}}
##'
##' @name HypergPrior-class
##' @keywords classes
##' @export
setClass(Class="HypergPrior",
         representation=
         representation(a="numeric"),
         contains=list("GPrior"),
         validity=           
         function(object){
             if(object@a <= 3)
             {
                 return("the parameter a must be larger than 3 for proper posteriors")
             }
             else
             {
                 return(TRUE)
             }})


##' Initialization method for the "HypergPrior" class
##'
##' @usage \S4method{initialize}{HypergPrior}(.Object, a, \dots)
##' @param .Object the \code{\linkS4class{HypergPrior}} we want to initialize
##' @param a the hyperparameter value
##' @param \dots unused
##' @return the initialized object
##'
##' @name HypergPrior-initialize
##' @aliases HypergPrior-initialize initialize,HypergPrior-method
##' 
##' @keywords methods internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
setMethod("initialize",
    signature(.Object = "HypergPrior"),
    function (.Object, a, ...) 
    {
       .Object@logDens <- function(g)
       {
           return(log(a - 2) - log(2) - a/2 * log1p(g))
       }
       callNextMethod(.Object, a=a, ...)
    })


##' Constructor for the hyper-g prior class
##'
##' @param a the hyperparameter which must be larger than 3, and should not be larger than 4
##' in order not to favour too much shrinkage a priori (default: 4)
##' @return a new \code{\linkS4class{HypergPrior}} object
##'
##' @keywords classes
##' @export
HypergPrior <- function(a=4)
{
    return(new("HypergPrior",
               a=a))
}

## ----------------------------------------------------------------------------

##' The hyper-g/n prior class
##'
##' The slots are:
##' \describe{
##' \item{a}{the hyperparameter}
##' \item{n}{the sample size}
##' }
##'
##' @seealso the constructor \code{\link{HypergnPrior}}
##'
##' @name HypergnPrior-class
##' @keywords classes
##' @export
setClass(Class="HypergnPrior",
         representation=
         representation(a="numeric",
                        n="integer"),
         contains=list("GPrior"),
         validity=           
         function(object){
             if(object@a <= 3)
             {
                 return("the parameter a must be larger than 3 for proper posteriors")
             }
             else if(object@n <= 0L)
             {
                 return("the parameter n must be a positive integer")
             }
             else
             {
                 return(TRUE)
             }})


##' Initialization method for the "HypergnPrior" class
##'
##' @usage \S4method{initialize}{HypergnPrior}(.Object, a, n, \dots)
##' @param .Object the \code{\linkS4class{HypergnPrior}} we want to initialize
##' @param a the hyperparameter value
##' @param n the sample size
##' @param \dots unused
##' @return the initialized object
##'
##' @name HypergnPrior-initialize
##' @aliases HypergnPrior-initialize initialize,HypergnPrior-method
##' 
##' @keywords methods internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
setMethod("initialize",
    signature(.Object = "HypergnPrior"),
    function (.Object, a, n, ...) 
    {
       .Object@logDens <- function(g)
       {
           return(log(a - 2) - log(2) - log(n) - a/2 * log1p(g/n))
       }
       callNextMethod(.Object, a=a, n=n, ...)
    })


##' Constructor for the "HypergnPrior" class
##'
##' @param a the hyperparameter which must be larger than 3, and should not be
##' larger than 4 in order not to favour too much shrinkage a priori (default:
##' 4) 
##' @param n the sample size (positive integer) 
##' @return a new \code{\linkS4class{HypergnPrior}} object
##'
##' @keywords classes
##' @export
HypergnPrior <- function(a=4, n)
{
    return(new("HypergnPrior",
               a=a,
               n=n))
}



## ----------------------------------------------------------------------------

##' The inverse gamma g-prior class
##'
##' The slots are:
##' \describe{
##' \item{a}{the first hyperparameter}
##' \item{b}{the second hyperparameter}
##' }
##'
##' @seealso the constructor \code{\link{InvGammaGPrior}}
##' 
##' @name InvGammaGPrior-class
##' @keywords classes
##' @export
setClass(Class="InvGammaGPrior",
         representation=
         representation(a="numeric",
                        b="numeric"),
         contains=list("GPrior"),
         validity=           
         function(object){
             if((object@a <= 0) || (object@b <= 0))
             {
                 return("the parameters a and b must be positive")
             }
             else
             {
                 return(TRUE)
             }})

##' Initialization method for the "InvGammaGPrior" class
##'
##' @usage \S4method{initialize}{InvGammaGPrior}(.Object, a, b, \dots)
##' @param .Object the \code{\linkS4class{InvGammaGPrior}} we want to initialize
##' @param a the first hyperparameter value
##' @param b the second hyperparameter value
##' @param \dots unused
##' @return the initialized object
##'
##' @name InvGammaGPrior-initialize
##' @aliases InvGammaGPrior-initialize initialize,InvGammaGPrior-method
##' 
##' @keywords methods internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
setMethod("initialize",
    signature(.Object = "InvGammaGPrior"),
    function (.Object, a, b, ...) 
    {
       .Object@logDens <- function(g)
       {
           return(-(a + 1) * log(g) - b/g + a * log(b) - lgamma(a))
       }
       callNextMethod(.Object, a=a, b=b, ...)
    })

##' Constructor for the inverse gamma g-prior class
##'
##' @param a the first positive hyperparameter (default: 0.001)
##' @param b the second positive hyperparameter (default: 0.001)
##' @return a new \code{\linkS4class{InvGammaGPrior}} object
##'
##' @keywords classes
##' @export
InvGammaGPrior <- function(a=0.001, b=0.001)
{
    return(new("InvGammaGPrior",
               a=a,
               b=b))
}


## ----------------------------------------------------------------------------

##' The custom g-prior class
##'
##' This class wraps around a custom log prior density for the covariance factor g. 
##'
##' @seealso the constructor \code{\link{CustomGPrior}}
##' 
##' @name CustomGPrior-class
##' @keywords classes
##' @export
setClass(Class="CustomGPrior",
         contains=list("GPrior"))


##' Constructor for the custom g-prior class
##'
##' @param logDens the log prior density function for g
##' @return a new \code{\linkS4class{CustomGPrior}} object
##'
##' @keywords classes
##' @export
CustomGPrior <- function(logDens)
{
    return(new("CustomGPrior",
               logDens=logDens))
}

## ----------------------------------------------------------------------------
