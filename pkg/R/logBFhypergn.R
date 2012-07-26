#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[logBFhypergn.R] by DSB Fre 15/06/2012 17:07 (CEST)>
##
## Description:
## Internal function which computes the log Bayes Factor against
## the null model if the hyper-g/n prior is used, via the Appell F1 function from
## the "appell" package. This is then called from C++ via Rcpp.
##
## History:
## 30/03/2012   file creation
## 15/06/2012   check size of model
#####################################################################################

##' Compute the log Bayes Factor against the null model under hyper-g/n prior.
##'
##' @param n number of observations
##' @param p number of covariates (excluding the intercept)
##' @param R2 coefficient of determination
##' @return the (double) value of the log BF, NA if the computation failed. 
##'
##' @importFrom appell appellf1
##' @keywords internal math
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
logBFhypergn <- function(n, p, R2)
{
    ## checks
    if((p + 1 >= n) | (R2 >= 1))
    {
        ## then p + 1 >= n or R2 >= 1, so the model is too large
        return(NA_real_)
    }
    
    ## convenient:
    oneMinusR2inv <- 1 / (1 - R2)

    ## call the Appell F1 function
    appellRes <- try(appell::appellf1(1 + p / 2,
                                      (1 + p - n) / 2,
                                      (n - 1) / 2,
                                      2 + p / 2,
                                      (n - 1) / n,
                                      (n - oneMinusR2inv) / n,
                                      debug=FALSE,
                                      userflag=1L))

    ## check if OK
    if(inherits(appellRes, "try-error") || is.na(appellRes$val) || (Re(appellRes$val) < 0))
    {
        return(NA_real_)
    } else {
        ## if it is OK, then add the other log terms
        logRet <- log(abs(appellRes$val)) - log(p + 2) + log(2) -
            p / 2 * log(n) + (n - 1) / 2 * log(oneMinusR2inv)
            
        return(logRet)
    }
}
