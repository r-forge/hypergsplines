#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[makeBasis.R] by DSB Son 21/11/2010 14:41 (CET)>
##
## Description:
## Construct the quantile-based linear spline basis for a numeric covariate
## vector. 
##
## History:
## 10/09/2010   file creation
## 13/09/2010   include example
## 30/09/2010   alter the function to (alternatively) directly accept knot
##              locations, this is needed for prediction e.g.!
## 26/10/2010   calculate the quantiles on the *unique* values of x, to avoid
##              problems with not really continuous covariates
## 19/11/2010   now optionally use cubic O-splines!
##              This re-uses code from supplementary material (Wand and Ormerod,
##              2010).
## 21/11/2010   redesign the function to be more flexible with predefined settings
#####################################################################################


##' Construction of quantile-based linear spline basis
##'
##' @param x the covariate vector
##' @param type the type of splines to be used, \dQuote{linear} splines
##' modelled by a truncated power series basis, or \dQuote{cubic} splines
##' modelled by O'Sullivan splines.
##' @param nKnots number of quantile-based knots to be used
##' @param settings list of settings. By default, this is \code{NULL}
##' and a whole new basis depending on \code{type} and \code{nKnots}
##' is built. Otherwise this is expected to be \code{attributes(basis)}
##' and the settings of the old \code{basis} are used to compute the
##' spline basis vectors at new grid points \code{x}.
##' @return the spline basis in a matrix with \code{length(x)} rows. As an
##' attribute, a list with the settings is provided.
##'
##' @example examples/makeBasis.R
##'
##' @importFrom splines bs spline.des
##' @export 
##' @keywords programming
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
makeBasis <- function(x,
                      type=c("linear", "cubic"),
                      nKnots,
                      settings=NULL)
{
    ## if a list of settings is provided, overwrite the other
    ## arguments accordingly
    if(! is.null(settings))
    {
        type <- settings$type
        knots <- settings$knots
        nKnots <- length(knots)
    }
    else ## otherwise determine type and knots
    {
        type <- match.arg(type)
        knots <- quantile(unique(x),
                          probs=
                          seq(from=0,
                              to=1,
                              length=nKnots + 2L))[2L:(nKnots + 1L)]
    }
    
    ## then we can construct the spline basis matrix, depending on type:  
    if(type=="linear")
    {
        basis <- t(sapply(x,
                          FUN="-",
                          knots))
        basis <- pmax(basis, 0)
        
        ## return with additional attributes
        return(structure(basis,
                         type=type,
                         knots=knots))
    }
    else ## type == "cubic"
    {
        ## use these bounds for function estimation
        bounds <- 
            if(is.null(settings))
                c(min(x) - 0.5,
                  max(x) + 0.5)
            else
                settings$bounds
        
        ## form the basis matrix
        basis <- splines::bs(x,
                             knots=knots,
                             degree=3L,
                             Boundary.knots=bounds,
                             intercept=TRUE)
        
        ## obtain the scaling matrix:
        LZ <-
            if(is.null(settings))
            {                
                ## construct the penalty matrix Omega
                allKnots <- c(rep(bounds[1L], 4L),
                              knots,
                              rep(bounds[2L], 4L))
                
                L <- 3L * (nKnots + 8L)
                
                xtilde <- (rep(allKnots, each=3L)[- c(1L, (L - 1L), L)] + 
                           rep(allKnots, each=3L)[- c(1L, 2L, L)]) / 2
                
                wts <- rep(diff(allKnots), each=3L) * rep(c(1, 4, 1) / 6, nKnots + 7L) 
                
                Bdd <- splines::spline.des(allKnots,
                                           xtilde,
                                           derivs=rep(2L, length(xtilde)),
                                           outer.ok=TRUE)$design  
                Omega <- t(Bdd * wts) %*% Bdd

                ## Obtain the spectral decomposition of Omega:
                eigOmega <- eigen(Omega)
                
                ## Obtain the matrix for linear transformation of B to Z:
                indsZ <- seq_len(nKnots + 2L)
                UZ <- eigOmega$vectors[, indsZ]
                LZ <- t( t(UZ) / sqrt(eigOmega$values[indsZ]) )
                
                ## Perform stability check:   
                indsX <- (nKnots + 3L):(nKnots + 4L)
                UX <- eigOmega$vectors[, indsX]   
                L <- cbind(UX, LZ)
                stabCheck <- t(crossprod(L, t(crossprod(L, Omega))))          
                if (sum(stabCheck^2) > 1.0001 * (nKnots + 2))
                {
                    warning("numerical instability arising from spectral decomposition")
                }

                LZ
            }
            else
            {
                settings$LZ
            }
        
        ## scale the basis matrix with LZ:
        basis <- basis %*% LZ
        
        ## return with additional attributes
        return(structure(basis,
                         type=type,
                         knots=knots,
                         bounds=bounds,
                         LZ=LZ))
    }
}    

