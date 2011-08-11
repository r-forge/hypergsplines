#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
## 
## Time-stamp: <[rshrinkage.R] by DSB Mit 15/09/2010 14:52 (CEST)>
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
#####################################################################################

##' Sample from the model-specific posterior of the shrinkage factor 
##'
##' Sample from the model-specific posterior of the shrinkage factor
##' t = g / (1 + g), using inverse sampling.
##' 
##' @param n number of samples
##' @param R2 coefficient of determination in the model
##' @param nObs number of observations used to fit the model
##' @param p number of design matrix columns without counting
##' the intercept
##' @param alpha used hyperparameter for hyper-g prior
##' @return \code{n} posterior shrinkage factor samples
##' 
##' @export 
##' @keywords distribution
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
rshrinkage <- function(n, 
                       R2,
                       nObs,
                       p,   
                       alpha)
{
    ## todo: add checks of input parameters?
    
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

    ## return the shrinkage samples
    return(ret)    
}
