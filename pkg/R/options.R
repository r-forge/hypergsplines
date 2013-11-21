#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[options.R] by DSB Die 28/08/2012 09:54 (CEST)>
##
## Description:
## Helper functions to collect options.
##
## History:
## 10/03/2011   file creation
## 14/03/2011   add gaussHermite to return list
## 17/03/2011   add getMcmc
## 08/04/2011   add option "nIwlsIterations" to getMcmc
## 19/05/2011   generalise "binaryLogisticCorrection" to "higherOrderCorrection"
#####################################################################################

##' @include helpers.R
{}


##' Collect the computation options
##' 
##' @param verbose should information on computation progress be given?
##' (default)
##' @param debug print debugging information? (not default)
##' @param nGaussHermite number of quantiles used in Gauss Hermite quadrature
##' for marginal likelihood approximation (and later in the MCMC sampler for the
##' approximation of the marginal covariance factor density).
##' @param useOpenMP shall OpenMP be used to accelerate the computations?
##' (default)
##' @param higherOrderCorrection should a higher-order correction of the
##' Laplace approximation be used for models with canonical links? (default)
##'
##' @export
##' @importFrom statmod gauss.quad
##' @keywords utilities
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getComputation <- function(verbose=TRUE,
                           debug=FALSE,
                           nGaussHermite=20L,
                           useOpenMP=TRUE,
                           higherOrderCorrection=TRUE)
{
    stopifnot(is.bool(verbose),
              is.bool(debug),
              is.posInt(nGaussHermite),
              is.bool(useOpenMP),
              is.bool(higherOrderCorrection))

    gaussHermite <- statmod::gauss.quad(n=nGaussHermite,
                                        kind="hermite")
    
    return(list(verbose=verbose,
                debug=debug,
                nGaussHermite=nGaussHermite,
                gaussHermite=gaussHermite,
                useOpenMP=useOpenMP,
                higherOrderCorrection=higherOrderCorrection))  
}



##' Collect the search options, todo: really use this in the package!
##' 
##' @param nModels how many best models should be saved? (default: 1\% of the
##' total number of (cached) models). Must not be larger than \code{nCache}.
##' @param nCache maximum number of best models to be cached at the same time
##' during the model sampling
##' @param chainlength length of the model sampling chain
##'
##' @export
##' @keywords utilities
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getSearch <- function(nModels=as.integer(max(1L, floor(nCache / 100))),
                      nCache=1e9L,
                      chainlength=1e4L)
{
    stopifnot(is.posInt(nModels),
              is.posInt(nCache),
              is.posInt(chainlength),
              nModels <= nCache)
    
    return(list(nModels=nModels,
                nCache=nCache,
                chainlength=chainlength))
}

##' Collect the MCMC options
##'
##' @param samples number of resulting samples (default: \code{10,000})
##' @param burnin number of burn-in iterations which are not saved (default:
##' \code{10,000}) 
##' @param step only every step-th iteration is saved after the burn-in
##' (default: \code{2})
##' @param nIwlsIterations maximum number of IWLS iterations in the proposal
##' step (default: \code{1}). Set that to a higher value to achieve higher
##' acceptance rates.
##' 
##' @export
##' @keywords utilities
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getMcmc <- function(samples=1e4L,
                    burnin=1e4L,
                    step=2L,
                    nIwlsIterations=1L)
{
    stopifnot(is.posInt(samples),
              is.posInt(burnin),
              is.posInt(step),
              is.posInt(nIwlsIterations))
    
    iterations <- as.integer(burnin + (step * samples))
    
    return(list(samples=samples,
                burnin=burnin,
                step=step,
                iterations=iterations,
                nIwlsIterations=nIwlsIterations))
}
