#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[glmGetSamples.R] by DSB Don 28/07/2011 15:26 (CEST)>
##
## Description:
## Get posterior samples for a specific GLM configuration.
##
## History:
## 17/03/2011   modify normal model function
## 05/04/2011   We now already calculate the acceptance ratio and log marg lik 
##              estimate in C++.
## 28/07/2011   allow for non-integer vector config
#####################################################################################


##' @include glmModelData.R
##' @include options.R
{}

##' Get posterior samples for a specific model configuration for
##' generalised response
##'
##' @param config the model configuration vector
##' @param modelData the data necessary for model estimation, which
##' is the result from \code{\link{glmModelData}}
##' @param mcmc MCMC options, result from \code{\link{getMcmc}}
##' @param computation computation options produced by
##' \code{\link{getComputation}}
##' 
##' @return A list with samples from the shrinkage hyperparameter
##' t = g / (g + 1) and the (linear and spline)
##' coefficients, and further information.
##'
##' @example examples/glmGetSamples.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
glmGetSamples <- function(config,
                          modelData,
                          mcmc=getMcmc(),
                          computation=getComputation())
{
    ## checks and extracts:
    stopifnot(all(config >= 0),
              identical(length(config), modelData$nCovs))
    
    ## go C++
    ret <-
        .Call("cpp_glmGetSamples",
              as.double(config),
              modelData,
              mcmc,
              computation)
  
    ## return the list
    return(ret)
}
