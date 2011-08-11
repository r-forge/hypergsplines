#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[stochSearch.R] by DSB Fre 05/08/2011 13:44 (CEST)>
##
## Description:
## Does a stochastic search for models with high posterior probability.
##
## History:
## 20/10/2010   file creation
## 22/10/2010   do post-processing of the C++ results to get something very
##              similar to the "exhaustive" output
## 25/01/2011   Add new prior types (as in getLogModelPrior before)
## 04/04/2011   Do not return R2 values for each model, because that does not
##              generalise to the GLM case, which we now include here (similarly
##              as in "exhaustive")
## 05/08/2011   now the configuration columns in the model data frame are of
##              mode "integer".
#####################################################################################

##' @include modelData.R
##' @include getLogModelPrior.R
{}

##' Stochastic search for models with high posterior probability
##' 
##' @param modelData the data necessary for model estimation, which is the
##' result from \code{\link{modelData}} or \code{\link{glmModelData}} 
##' @param modelPrior either \dQuote{flat} (default), \dQuote{exponential},
##' \dQuote{independent} or \dQuote{dependent}, see
##' \code{\link{getLogModelPrior}} for details. 
##' @param chainlength length of the model sampling chain (default: 100,000)
##' @param nCache maximum number of best models to be cached at the same time
##' during the model sampling (by default equal to \code{chainlength})
##' @param nModels how many best models should be saved? (default: 1\% of the
##' total number of \code{nCache}). Must not be larger than \code{nCache}.
##' @param computation computation options produced by
##' \code{\link{getComputation}}, only matters for generalised response models. 
##' @return a list with the data frame \dQuote{models} comprising the model
##' configurations, log marginal likelihoods / priors / posteriors and
##' hits in the MCMC run, the inclusion probabilities matrix
##' \dQuote{inclusionProbs}, the number of total visited models
##' \dQuote{numVisited} and the log normalization constant
##' \dQuote{logNormConst}.
##'
##' @example hypergsplines/examples/stochSearch.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
stochSearch <- function(modelData,
                        modelPrior=
                        c("flat",
                          "exponential",
                          "independent",
                          "dependent"),
                        chainlength=100000L,
                        nCache=chainlength,
                        nModels=as.integer(max(nCache / 100, 1)),
                        computation=getComputation())
{
    modelPrior <- match.arg(modelPrior)
    
    ## checks
    stopifnot(chainlength > 1,
              nModels >= 1,
              nModels <= nCache)

    ## bundle the search settings
    searchSettings <- list(chainlength=as.integer(chainlength),
                           nCache=as.integer(nCache),
                           nModels=as.integer(nModels))

    ## decide if this is a normal model
    isNormalModel <- is.null(modelData$family)
    
    ## then go C++
    ret <-
        if(isNormalModel)
            .Call("cpp_stochSearch",
                  modelData,
                  modelPrior,
                  searchSettings)
        else
            .Call("cpp_glmStochSearch",
                  modelData,
                  modelPrior,
                  searchSettings,
                  computation)

    ## postprocessing:
    configMatrix <- t(sapply(ret, "[[", "configuration"))
    colnames(configMatrix) <- colnames(modelData$X)
    mode(configMatrix) <- "integer"
    
    infoMatrix <- t(sapply(ret, function(x) unlist(x$information)))
    
    inclusionProbs <- attr(ret, "inclusionProbs")
    inclusionMatrix <- matrix(nrow=3L,
                              ncol=modelData$nCovs,
                              dimnames=
                              list(c("not included", "linear", "non-linear"),
                                   colnames(modelData$X)))
    inclusionMatrix["not included", ] <- 1 - inclusionProbs$included
    inclusionMatrix["non-linear", ] <- inclusionProbs$smooth
    inclusionMatrix["linear", ] <-
        inclusionProbs$included - inclusionProbs$smooth

    retList <- list(models=cbind(as.data.frame(configMatrix), infoMatrix),
                    inclusionProbs=inclusionMatrix,
                    numVisited=attr(ret, "numVisited"),
                    logNormConst=attr(ret, "logNormConst"))
    
    ## return the results
    return(retList)
}

