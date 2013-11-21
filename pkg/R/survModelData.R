#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[survModelData.R] by DSB Mon 26/08/2013 16:46 (CEST)>
##
## Description:
## Prepare model data for proportional hazards regression.
##
## History:
## 27/08/2012   file creation, adaptation from glmModelData
## 28/08/2012   keep attributes in X and Z[[j]]'s when appending the pseudo zero
##              obs
## 20/09/2012   move to class-based handling of hyperprior on g
## 06/12/2012   correct Poisson model approximation
#####################################################################################

##' @include helpers.R
##' @include getQmatrix.R
##' @include GPrior-classes.R
##' @include glmModelData.R
{}


##' Process the data needed for survival models
##'
##' @param times the numeric vector of survival times
##' @param X the numeric matrix of covariates (not including time)
##' @param observed the logical vector of observation indicators, \code{TRUE}
##' entries represent truly observed survival times, \code{FALSE} entries
##' represent censored survival times.
##' @param continuous see \code{\link{glmModelData}} for details
##' @param nKnots see \code{\link{glmModelData}} for details
##' @param splineType see \code{\link{glmModelData}} for details
##' @param gPrior see \code{\link{glmModelData}} for details. Defaults to
##' the hyper-g/n prior where n is chosen as the number of events, instead
##' of the number of observations. 
##' @return a list with the internally needed results.
##'
##' @example examples/survModelData.R
##'
##' @export
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
survModelData <- function(times,
                          X,
                          observed,
                          continuous=rep.int(TRUE, ncol(X)),
                          nKnots=4L,
                          splineType="linear",
                          gPrior=HypergnPrior(a=4, n=sum(observed)))
{
    ## checks and extracts:
    stopifnot(is.vector(times),
              is.numeric(times),
              is.matrix(X),
              is.numeric(X),
              is.vector(observed),
              is.logical(observed),
              all(! is.na(times)),
              all(! is.na(X)),
              all(! is.na(observed)),
              all(times > 0),
              is.logical(continuous),
              is.posInt(nKnots))

    if(is(gPrior, "character"))
    {
        gPrior <- match.arg(gPrior,
                            choices=c("hyper-g/n", "hyper-g"))
        gPrior <-
            if(gPrior == "hyper-g/n")
                HypergnPrior(a=4, n=sum(observed))
            else
                HypergPrior(a=4)
    } else {
        stopifnot(is(gPrior, "GPrior"))
    }

    ## ensure that X has column names
    if(is.null(colnames(X)))
        colnames(X) <- paste("V", seq_len(ncol(X)),
                             sep="")

    ## more checks
    stopifnot(identical(length(times),
                        nrow(X)),
              identical(length(times),
                        length(observed)),
              identical(length(continuous),
                        ncol(X)),
              nKnots > 1L && nKnots < length(times) - 2L)
    
    ## order everything after the survival times:
    timesOrder <- order(times)

    if(! all(timesOrder == seq_along(times)))
    {
        warning("Input data were reordered so that the survival times are sorted") 
    }
    
    times <- times[timesOrder]
    X <- X[timesOrder, ]
    observed <- observed[timesOrder]

    ## Get the matrix Q
    Q <- getQmatrix(times)

    ## Define the matrix to select observations with positive Q entry
    select <- (Q > 0)

    ## Define the pseudo response
    z <- diag(c(0L, as.integer(observed)))[select]

    ## Define the time design column for estimating the baseline hazard
    xTimes <- c(0, times)[col(Q)[select]]

    ## Define the covariates design matrix
    xCovariates <- rbind(0, X)[row(Q)[select], ]

    ## Define the offsets
    offsets <- log(Q[select])

    ## append time column to covariate matrix 
    xWhole <- cbind(time=xTimes,
                    xCovariates)
    continuous <- c(TRUE,
                    continuous)
    
    return(glmModelData(y=z,
                        X=xWhole,
                        continuous=continuous,
                        nKnots=nKnots,
                        splineType=splineType,
                        gPrior=gPrior,
                        offsets=offsets,
                        family=poisson(link="log")))
}
