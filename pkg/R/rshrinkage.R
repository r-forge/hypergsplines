#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
## 
## Time-stamp: <[rshrinkage.R] by DSB Mon 18/06/2012 11:14 (CEST)>
##
## Description:
## Sample from the model-specific posterior of the shrinkage factor
## t = g / (1 + g). 
##
## History:
## 04/09/2008   file creation
## 05/09/2008   correct innerFraction computation
## 05/10/2009   some beautifications
## 14/09/2010   copy from bfp package
## 29/03/2012   - different interface
##              - include sampling for the hyper-g/n prior
#####################################################################################

##' Sample from the model-specific posterior of the shrinkage factor 
##'
##' Sample from the model-specific posterior of the shrinkage factor
##' t = g / (1 + g), using either inverse sampling (for the hyper-g
##' prior) or numerical approximations (for the hyper-g/n prior).
##' 
##' @param n number of samples
##' @param marginal the marginal model, which is the result from
##' \code{\link{calculateModel}}
##' @param modelData the data necessary for model estimation, which
##' is the result from \code{\link{modelData}}
##' @return \code{n} posterior shrinkage factor samples
##'
##' @importFrom Runuran pinv.new ur
##' @export 
##' @keywords distribution
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
rshrinkage <- function(n, 
                       marginal,
                       modelData)
{
    ## extract contents
    R2 <- min(marginal$coefR2, 1)
    nObs <- modelData$nObs
    p <- marginal$dim.lin
    alpha <- 4

    ## decide computation depending on hyperprior on g
    if(identical(modelData$gPrior,
                 "hyper-g"))
    {
        ## parameters for the beta functions
        shape1 <- (nObs - p - alpha + 1) / 2
        shape2 <- (p + alpha - 2) / 2

        ## uniforms needed for inverse sampling
        uniforms <- runif(n=n)
        
        ## evaluate the quantile function at these values to obtain the samples
        oneMinusR2 <- 1 - R2

        innerFraction <-
            oneMinusR2 /
                qbeta(uniforms +
                      (1 - uniforms) * pbeta(oneMinusR2, shape1, shape2),
                      shape1, shape2)
        
        ret <- (1 - innerFraction) / R2
        
    } else { ## identical(modelData$gPrior, "hyper-g/n")
          
        ## unnormalised density function for z = log(g / n):
        unDensz <- function(z, log=FALSE, logOffset=0)
        {
            logDens <- dlogis(z, log=TRUE) +
                (nObs - 1 - p) / 2 * log1p(nObs * exp(z)) -
                    (nObs - 1) / 2 * log1p(nObs * exp(z) * (1 - R2)) +
                        logOffset

            if(log){
                return(logDens)
            } else {
                return(exp(logDens))
            }
        }

        ## maximize it
        optimRes <- optimize(unDensz,
                             lower=-100, upper=200,
                             maximum=TRUE,
                             log=TRUE)
        zMax <- optimRes$maximum

        ## now sample
        unuranObject <- Runuran::pinv.new(pdf=unDensz,
                                          lb=-Inf,
                                          ub=Inf,
                                          islog=TRUE,
                                          center=zMax,
                                          log=TRUE,
                                          logOffset=-optimRes$objective)
        ret <- Runuran::ur(unuranObject,
                           n=n)

        ## and finally transform z = log(g/n) to t = g / (1 + g)
        ret <- exp(ret) * n
        ret <- ret / (1 + ret)
    }

    ## return the shrinkage samples
    return(ret)
}
