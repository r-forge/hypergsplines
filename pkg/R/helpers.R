#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[helpers.R] by DSB Fre 11/05/2012 12:22 (CEST)>
##
## Description:
## Helper functions
##
## History:
## 10/03/2011   file creation: copy from glmBfp and add stuff
## 10/05/2012   add "checkModelConfigs"
## 11/05/2012   add "expWithConst"
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

##' Check coherency of model configurations with data
##'
##' @param modelConfigs the matrix of model configurations
##' @param modelData the data necessary for model estimation, which is the 
##' result from \code{\link{modelData}} or \code{\link{glmModelData}}
##' @return returns \code{TRUE} if model configurations are coherent with
##' \code{modelData}, otherwise throws error messages
##'
##' @keywords internal
##' @author Daniel Sabanes Bove
checkModelConfigs <- function(modelConfigs,
                              modelData)
{    
    stopifnot(is.matrix(modelConfigs),
              (nModels <- nrow(modelConfigs)) > 0L,
              (nCovs <- ncol(modelConfigs)) > 0L)
    mode(modelConfigs) <- "integer"
    
    ## check for coherency with modelData
    stopifnot(identical(nCovs, modelData$nCovs),
              identical(colnames(modelConfigs), colnames(modelData$X)))

    for(j in seq_len(modelData$nCovs))
    {
        if(modelData$continuous[j])
            stopifnot(all(modelConfigs[, j] %in% modelData$degrees))
        else
            stopifnot(all(modelConfigs[, j] %in% 0:1))
    }

    return(TRUE)
}

##' Exponentiate a vector of values plus a constant
##'
##' This functions tries to choose the constant such that the
##' smallest result is larger than zero and the largest result
##' is smaller than infinity, in the numerical representation.
##' If this is not possible, results numerically identical to
##' zero are tolerated. Only values which are finite are taken
##' into account for this considerations, because infinite values
##' will be either zero or infinity anyway, regardless of the
##' constant.
##' 
##' @param values argument values, of which at least one
##' must be finite
##' @return the exp(values + constant) results, see details
##'
##' @keywords internal
##' @author Daniel Sabanes Bove
expWithConst <- function(values)
{
    ## consider only finite values
    stopifnot(any(is.finite(values)))
    finValues <- values[is.finite(values)]    
    
    ## determine bounds for constant
    minConst <- max(finValues) - log(.Machine$double.xmax)
    maxConst <- min(finValues) - log(.Machine$double.xmin)

    ## determine constant
    if(minConst < maxConst)
    {
        constant <- mean(c(minConst, maxConst))
    } else {
        constant <- minConst
        warning("results numerically identical to zero tolerated") 
    }

    ## results
    results <- exp(values - constant)
    return(results)
}
