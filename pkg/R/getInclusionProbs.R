#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getInclusionProbs.R] by DSB Mon 07/02/2011 16:50 (CET)>
##
## Description:
## Extract posterior inclusion probabilities from the models table.
##
## History:
## 13/10/2010   file creation
## 07/02/2011   take the minimum only within the finite values of logPost
#####################################################################################


##' Extract posterior inclusion probabilities from the models table
##'
##' @param models the models table (a data frame), the result from
##' \code{\link{exhaustive}}
##' @param modelData data used for model estimation
##' @param logMargLiks vector of log marginal likelihoods (defaults
##' to the \code{logMargLik} column of \code{models})
##' @param logPriors vector of log prior model probabilities (defaults
##' to the \code{logPrior} column of \code{models})
##' @return a nice matrix with the probabilities for exclusion, linear 
##' and non-linear inclusion of the covariates
##'
##' @example examples/getInclusionProbs.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getInclusionProbs <- function(models,
                              modelData,
                              logMargLiks=models$logMargLik,
                              logPriors=models$logPrior)
{
    ## checks:
    stopifnot(identical(nrow(models),
                        length(logMargLiks)),
              identical(nrow(models),
                        length(logPriors)))

    ## compute posterior model probabilities
    logPost <- logMargLiks + logPriors
    post <- exp(logPost - min(logPost[is.finite(logPost)]))
    post <- post / sum(post)

    ## extract the covariate names
    covNames <- colnames(modelData$X)
    
    ## construct the empty inclusion probs matrix
    inclusion <- matrix(data=NA,
                        nrow=3L,
                        ncol=modelData$nCovs,
                        dimnames=
                        list(c("not included", "linear", "non-linear"),
                             covNames))

    ## fill it
    for(cov in covNames)
    {
        inclusion["not included", cov] <-
            sum(post[models[, cov] == 0])
        
        inclusion["linear", cov] <-
            sum(post[models[, cov] == 1])
        
        inclusion["non-linear", cov] <-
            sum(post[models[, cov] > 1])
    }

    ## return it
    return(inclusion)
}
