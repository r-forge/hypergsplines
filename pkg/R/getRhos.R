#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getRhos.R] by DSB Don 28/07/2011 11:38 (CEST)>
##
## Description:
## Compute the penalty parameters corresponding to degrees of freedom
## for a specific spline basis.
##
## History:
## 10/09/2010   file creation
## 13/09/2010   fix brace bug, include example and note.
## 09/03/2011   add weights argument needed for the GLM case.
## 22/07/2011   - use SVD of sqrt(weights) * Z instead of
##              eigendecomposition of crossprod(Z, weights * Z) to get the
##              eigenvalues, seems to produce better results!
##              - and use a wider range for finding the root
##              - exploit that the function is strictly increasing
## 28/07/2011   different interface: now lambdas are expected as input,
##              in order to not compute them too often without need
#####################################################################################


##' Compute the penalty parameters corresponding to degrees of freedom
##'
##' @param lambdas the eigenvalues of Z^T * W * Z, where Z is the
##' spline basis matrix
##' @param degrees the vector of degrees of freedom
##' @return the vector of rho parameters corresponding to the \code{degrees} 
##'
##' @example examples/getRhos.R
##' @note All \code{degrees} must be positive, and smaller than the dimension
##' of the spline basis (which is the number of columns of \code{Z}).
##' 
##' @export 
##' @keywords programming
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getRhos <- function(lambdas,
                    degrees)
{
    ## checks and extracts
    dimBasis <- length(lambdas)
    nDegrees <- length(degrees)

    stopifnot(max(degrees) < dimBasis,
              all(degrees > 0))
   
    ## so the degree function is:
    getDegree <- function(rho)
    {
        sum(lambdas / (lambdas + 1 / rho))
    }

    ## and the inverse is:
    getRho <- function(d,
                       lower=.Machine$double.xmin,
                       upper=.Machine$double.xmax)
    {
        uniroot(f=function(rho) getDegree(rho) - d,
                lower=lower,
                upper=upper,
                maxiter=1e4L)$root
    }

    ## compute and return the rhos:

    ## very naive would be:
    ## ret <- sapply(degrees, getRho)

    ## but we can exploit that the function is
    ## strictly increasing:

    ## order degrees
    degreesOrder <- order(degrees)
    degrees <- degrees[degreesOrder]

    ## calculate rhos
    ret <- numeric(length(degrees))
    ret[1L] <- getRho(degrees[1L])
    for(i in seq_along(degrees)[-1L])
    {
        ret[i] <- getRho(degrees[i],
                         lower=ret[i-1])
    }

    ## and order back the result
    ret <- ret[order(degreesOrder)]

    ## then return it
    return(ret)
}
