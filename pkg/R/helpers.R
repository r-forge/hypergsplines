#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[helpers.R] by DSB Don 10/03/2011 11:40 (CET)>
##
## Description:
## Helper functions
##
## History:
## 10/03/2011   file creation: copy from glmBfp and add stuff
#####################################################################################

##' Predicate checking for a scalar
##'
##' @param x the object being checked
##' @return Returns \code{TRUE} if \code{x} is a length one vector, i.e. a
##' scalar.
##'
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
is.scalar <- function(x)
{
    return(identical(length(x), 1L))
}

##' Predicate checking for a boolean option
##'
##' @param x the object being checked
##' @return Returns \code{TRUE} if \code{x} is a logical scalar
##' 
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
is.bool <- function(x)
{
    return(is.scalar(x) &&
           is.logical(x))
}

##' Predicate checking for a scalar positive integer
##'
##' @param x the object being checked
##' @return Returns \code{TRUE} if \code{x} is a scalar positive integer.
##'
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
is.posInt <- function(x)
{
    return(is.scalar(x) &&
           is.integer(x) &&
           x > 0)
}
