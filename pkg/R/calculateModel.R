#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[calculateModel.R] by DSB Die 26/07/2011 11:26 (CEST)>
##
## Description:
## For one specific model configuration, do the marginalization over the
## spline coefficients, calculate the resulting coefficient of determination
## etc. We have a separate function for this because both model exploration
## and posterior parameter sampling need this functionality.
##
## History:
## 15/09/2010   file creation
## 27/09/2010   keep the column names of X.lin also when "rescaling"
##              with backsolve is necessary
## 13/10/2010   done: transfer this whole function to C++ via RcppArmadillo
##              and check if the speed is increased well enough
## 14/10/2010   use Z.tcrossprod.list of modelData
## 04/04/2011   Cleanup:
##              Include here only the C++ version. The old R version can be
##              found in the archive packages.
## 26/07/2011   allow for non-integer vector config
#####################################################################################

##' @include modelData.R
{}


##' Calculate intermediate information for a specific model
##'
##' @param config the model configuration vector
##' @param modelData the result from \code{\link{modelData}}
##' @return A list with necessary intermediate information
##'
##' @example hypergsplines/examples/calculateModel.R
##' 
##' @export 
##' @keywords regression internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
calculateModel <- function(config,
                           modelData)
{
    ## coerce config to double vector
    config <- as.double(config)
    
    ## checks:
    stopifnot(identical(length(config), modelData$nCovs))

    ## then directly go C++
    ret <- .Call("cpp_calculateModel",
                 config,
                 modelData)

    ## attach names
    colnames(ret$X.lin) <- names(ret$betaOLS) <-
        colnames(modelData$X)[ret$whichLinear]

    ## and return
    return(ret)
}
