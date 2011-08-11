#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getLogModelPrior.R] by DSB Mit 29/06/2011 16:42 (CEST)>
##
## Description:
## Get the log prior probability for a model.
##
## History:
## 13/10/2010   file creation
## 18/10/2010   also allow a matrix argument "config", for better speed
## 07/01/2011   add new prior types: "independent" and "dependent"
## 10/01/2011   - fix copy/paste bug: variable "d" -> "config"
##              - allow additional elements on the right of the vector/matrix
## 12/01/2011   - correct nCovs extraction for the "dependent" prior
##              - be careful with linear covariates which cannot be
##              spline transformed
## 29/06/2011   revise model probabilities
#####################################################################################

##' Get the log prior probability for a model
##'
##' The \dQuote{exponential} prior assumes \eqn{f(\boldsymbol{d}) \propto
##' \exp(-2 \sum_{j=1}^{p} d_{j})}{f(d) ~ exp(-2 * sum(d_j))} which penalizes 
##' models with many degrees of freedom regardless of linear (\eqn{d_j = 1}) 
##' or smooth (\eqn{d_j > 1}) inclusion of covariates.
##'
##' The \dQuote{independent} prior (as the exponential prior) that the
##' different covariates are independent,
##' \eqn{f(\boldsymbol{d}) = \prod_{j=1}^{p} f(d_{j})}{f(d) = prod(f(d_j))}, 
##' where \eqn{f(d_{j})}{f(d_j)} is built hierarchically:
##' The prior probability for exclusion is 0.5, thus f(0)=0.5. Given inclusion,
##' the other half is split up uniformly between all possible degrees of freedom
##' larger than 0. 
##'
##' The \dQuote{dependent} prior imposes a uniform distribution for the number
##' of included variables. The \dQuote{position} of the  included variables,
##' i.e. which of the  variables are included, is uniformly distributed on all
##' possible configurations. The degrees of freedoms are then conditionally
##' independent and uniformly distributed on the respective possible values.
##'
##' @param config the model configuration vector, or a matrix of
##' model configurations in rows. There may be additional elements
##' to the right, e.g. the log marginal likelihood column/value,
##' if \code{modelData} is provided. These are then ignored internally.
##' @param type either \dQuote{flat} (default), \dQuote{exponential},
##' \dQuote{independent} or \dQuote{dependent}. See the details.
##' @param modelData the model setup from \code{\link{modelData}}
##' @return the (unnormalized) log prior probability
##'
##' @example hypergsplines/examples/getLogModelPrior.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getLogModelPrior <- function(config,
                             type=
                             c("flat",
                               "exponential",
                               "independent",
                               "dependent"),
                             modelData)
{
    ## determine the prior type
    type <- match.arg(type)

    ## determine arg type
    isMatrix <- identical(length(dim(config)),
                          2L)

    ## coerce config to correct form if we have modelData
    if(! missing(modelData))
    {
        config <- 
            if(isMatrix)
            {
                as.matrix(config[, 1:modelData$nCovs])
            }
            else
            {
                config[1:modelData$nCovs]
            }
    }
    
    ## then compute and return the corresponding log prior probability/ies
    logProb <-
        if(type == "flat")
        {
            if(isMatrix)
            {
                rep.int(0, nrow(config))
            }
            else
            {
                0
            }                
        }
        else if(type == "exponential")
        {
            if(isMatrix)
            {
                -2 * rowSums(config)
            }
            else
            {
                -2 * sum(config)
            }
        }
        else if(type == "dependent")
        {
            nCovs <- ncol(modelData$X)
            nInclDegrees <- length(setdiff(modelData$degrees, 0L))

            included <- config > 0
            
            if(isMatrix)
            {                
                nIncluded <- rowSums(config > 0)
                nSplinePossible <- rowSums(included &
                                           rep(modelData$continuous,
                                               each=nrow(included)))
            }
            else
            {
                nIncluded <- sum(included)
                nSplinePossible <- sum(included & modelData$continuous)
            }
            
            - lchoose(nCovs, nIncluded) - log1p(nCovs) -
                    nSplinePossible * log(nInclDegrees)
        }
        else ## if(type=="independent")
        {
            nInclDegrees <- length(setdiff(modelData$degrees, 0L))

            exclLogProb <- log(1/2)
            inclDegreeLogProb <- log(1/2 / nInclDegrees)

            ## log prior probabilities for continuous covs:
            smoothLogProbs <-
                c(exclLogProb,                  
                  rep.int(inclDegreeLogProb, nInclDegrees))

            ## and for linear covs:
            linearLogProbs <-
                c(exclLogProb,
                  exclLogProb)          # only because we have log(1/2)
            
            if(isMatrix)
            {
                smoothConfig <- config[, modelData$continuous, drop=FALSE]
                linearConfig <- config[, ! modelData$continuous, drop=FALSE]

                smoothTmp <- smoothLogProbs[smoothConfig + 1]
                dim(smoothTmp) <- dim(smoothConfig)

                linearTmp <- linearLogProbs[linearConfig + 1]
                dim(linearTmp) <- dim(linearConfig)                

                rowSums(smoothTmp) + rowSums(linearTmp)
            }
            else
            {            
                smoothConfig <- config[modelData$continuous]
                linearConfig <- config[! modelData$continuous]

                sum(smoothLogProbs[smoothConfig + 1]) +
                    sum(linearLogProbs[linearConfig + 1])
            }
        } 

    return(logProb)            
}

