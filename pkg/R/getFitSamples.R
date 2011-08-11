#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getFitSamples.R] by DSB Don 19/05/2011 10:37 (CEST)>
##
## Description:
## Extract fit samples on linear predictor scale
## from posterior coefficients samples.
##
## History:
## 30/09/2010   file creation
## 11/01/2011   adapt to change of samples$linearCoefs
## 04/04/2011   - rename to "getFitSamples" because prediction samples would
##              comprise the additional observation error.
##              - Note that we can also use this function for GLMs!
#####################################################################################

##' @include modelData.R
##' @include getFunctionSamples.R
{}


##' Extract fit samples on linear predictor scale
##' from posterior coefficients samples
##'
##' @param X the numeric matrix with new covariate values on the
##' \emph{original} scale, with the same column layout as originally
##' provided to \code{\link{modelData}} or \code{\link{glmModelData}}. 
##' @param samples the samples object (either from \code{\link{getSamples}} or 
##' the \code{samples} element from \code{\link{glmGetSamples}})
##' @param modelData the corresponding model data object
##' @return the fit samples as a matrix.
##'
##' @example hypergsplines/examples/getFitSamples.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getFitSamples <- function(X,
                          samples,
                          modelData)
{
    ## checks and extracts:
    stopifnot(is.matrix(X),
              is.numeric(X),
              all(! is.na(X)),
              identical(ncol(X), modelData$nCovs))

    nNew <- nrow(X)
    covNames <- colnames(modelData$X)

    ## start with the intercept samples
    fitSamples <- tcrossprod(rep(1L, nNew),
                              samples$intercept)
    ## so the layout will be nNew x nSamples

    ## then sequentially add the function samples, from all covariates
    ## included in the model which has been sampled:
    for(thisCov in names(samples$linearCoefs))
    {
        ## note that we do only rely on the fact that X has the
        ## same column order as modelData$X, but not the same
        ## column names possibly.
        fitSamples <- fitSamples +
            getFunctionSamples(x=X[, match(thisCov, covNames)],
                               covName=thisCov,
                               samples=samples,
                               modelData=modelData)
    }
    
    ## return the fit samples
    return(fitSamples)
}
