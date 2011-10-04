#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getFunctionSamples.R] by DSB Mon 04/04/2011 15:42 (CEST)>
##
## Description:
## Extract function values samples from posterior coefficients samples.
##
## History:
## 27/09/2010   file creation
## 30/09/2010   take the same knots as used in the model fitting, not new quantile-
##              based ones!
## 21/11/2010   adapt to new signature of makeBasis: we pass now the whole list
##              of attributes instead of only the knot locations
## 11/01/2011   adapt to change of samples$linearCoefs
## 04/04/2011   Note that we can also use this function for GLMs!
#####################################################################################

##' @include makeBasis.R
{}


##' Extract function values samples from posterior coefficients samples
##'
##' @param x the vector of abscissa values (on the \emph{original} scale!)
##' @param covName string with the name of the covariate
##' @param samples the samples object (either from \code{\link{getSamples}} or
##' the \code{samples} element from \code{\link{glmGetSamples}})
##' @param modelData the corresponding model data object
##' @return the function values samples as a matrix.
##'
##' @example examples/getFunctionSamples.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getFunctionSamples <- function(x,
                               covName,
                               samples,
                               modelData)
{
    ## checks and extracts:
    stopifnot(is.numeric(x),
              is.character(covName),
              identical(length(covName), 1L),
              covName %in% names(samples$linearCoefs))

    ## scale the x vector
    x.scaled <- (x - attr(modelData$X, "scaled:center")[covName]) /
        attr(modelData$X, "scaled:scale")[covName]

    ## center the x vector
    x.scaled.centered <- x.scaled -
        attr(modelData$X, "scaled:colMeans")[covName] 

    ## compute the resulting function value samples:
    funSamples <- tcrossprod(x.scaled.centered,
                             samples$linearCoefs[[covName]])

    ## if the covariate is also included nonlinearly, add that part:
    if(covName %in% names(samples$splineCoefs))
    {
        zAttrs <- attributes(modelData$Z.list[[covName]])

        ## make and center the spline basis 
        z <- makeBasis(x=x.scaled,
                       settings=zAttrs)

        z <- z - tcrossprod(rep(1,
                                length(x.scaled)),
                            zAttrs[["scaled:colMeans"]])

        z <- z - tcrossprod(x.scaled.centered,
                            zAttrs[["scaled:crossX"]])

        ## resulting nonlinear part samples
        z.samples <- z %*% samples$splineCoefs[[covName]]

        ## add to total
        funSamples <- funSamples + z.samples
    }

    ## return the samples
    return(funSamples)
}
