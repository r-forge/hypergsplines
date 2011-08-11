#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getSamples.R] by DSB Don 28/07/2011 11:43 (CEST)>
##
## Description:
## Get posterior samples for a specific model configuration.
##
## History:
## 14/09/2010   file creation
## 15/09/2010   first complete version
## 27/09/2010   sort the spline coefficients samples according to the
##              covariates they belong to
## 28/09/2010   get the samples.linearCoefs proper variable (row)names
## 29/09/2010   correct the precision matrix calculation: *inverse* of diagRho
##              is added!
## 21/11/2010   adapt to the generalization dimSplineBasis of nKnots:
##              with cubic O-Splines, nKnots is smaller than the dimension of the
##              spline basis matrix.
## 11/01/2011   change format of linearCoefs to list with covariate names
##              instead of matrix with row names (necessary for compatibility
##              of output with that from getBmaSamples).
## 26/07/2011   start revision to allow for non-integer vector config
## 28/07/2011   different interface for getRhos
#####################################################################################


##' @include rshrinkage.R
##' @include invGamma.R
##' @include modelData.R
##' @include calculateModel.R
{}

##' Get posterior samples for a specific model configuration
##'
##' @param config the model configuration vector
##' @param nSamples number of samples to simulate
##' @param modelData the data necessary for model estimation, which
##' is the result from \code{\link{modelData}}
##' @return A list with samples from the shrinkage hyperparameter
##' t = g / (g + 1), the regression variance, and the (linear and spline)
##' coefficients.
##'
##' @example hypergsplines/examples/getSamples.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getSamples <- function(config,
                       nSamples,
                       modelData)
{
    ## checks and extracts:
    stopifnot(is.integer(nSamples),
              nSamples > 0L,
              identical(length(nSamples),
                        1L),
              identical(length(config), modelData$nCovs))
    config <- as.numeric(config)
    
    ## calculate the (marginal) model
    marginal <- calculateModel(config=config,
                               modelData=modelData)

    ## sample the shrinkage hyperparameter t,
    ## yet only if this is not the null model:
    samples.t <-
        if(marginal$isNullModel)
            numeric(0)
        else
            rshrinkage(n=nSamples,
                       R2=marginal$coefR2,
                       nObs=modelData$nObs,
                       p=marginal$dim.lin,
                       alpha=modelData$a)

    ## compute corresponding parameters of the inverse gamma posterior
    aStar <- (modelData$nObs - 1) / 2

    sst <- sum((marginal$y - marginal$meanY)^2)
    bStar <- (1 - samples.t * marginal$coefR2) * sst / 2

    ## now sample the regression variance
    samples.sigma2 <- rinvGamma(n=nSamples,
                                a=aStar,
                                b=bStar)

    ## the intercept
    samples.intercept <- rnorm(n=nSamples,
                               mean=mean(modelData$y), 
                               sd=sqrt(samples.sigma2 / modelData$nObs))

    ## the coefficients of linear effects
    samples.linearCoefs <-
        if(marginal$isNullModel)
            list()
        else
        {
            ## cholesky decomposition of XtX (without intercept)
            XtX <- crossprod(marginal$X.lin)
            XtXroot <- chol(XtX)
            
            ## then sample the coefficients:

            ## first N(0, I_{dim.lin})
            simCoefs <- matrix(data=rnorm(n=marginal$dim.lin * nSamples),
                               nrow=marginal$dim.lin,
                               ncol=nSamples)

            ## then N(0, t * sigma2 * I_{dim.lin})
            simCoefs <- sweep(simCoefs,
                              MARGIN=2,
                              STATS=sqrt(samples.t * samples.sigma2),
                              FUN="*")
            
            ## then N(0, t * sigma2 * (XtX)^{-1}) 
            simCoefs <- backsolve(r=XtXroot,
                                  x=simCoefs,
                                  k=marginal$dim.lin)

            ## and finally N(t * betaOLS, t * sigma2 * (XtX)^{-1})
            simCoefs + sapply(samples.t,
                              FUN="*",
                              marginal$betaOLS)
        }

    ## finally the spline coefficients samples
    samples.splineCoefs <-
        if(marginal$hasOnlyLinear)
            list()
        else
        {
            ## number of columns of the spline design matrix:
            dim.spline <- modelData$dimSplineBasis * length(marginal$whichSpline)
                
            ## diagonal rho matrix entries will come in here:
            diagRho <- numeric(dim.spline)

            ## the whole spline design matrix will come in here:
            splineDesign <- matrix(nrow=modelData$nObs,
                                   ncol=dim.spline)

            ## start filling
            nSplineTerms <- length(marginal$whichSpline)
            splineNames <- character(nSplineTerms)
            for(i in seq_len(nSplineTerms))
            {
                ## which covariate?
                s <- marginal$whichSpline[i]

                ## save the name for the list construction below
                splineNames[i] <- names(modelData$Z.list[s])
                
                ## get the correct rho
                d <- config[s]

                thisRho <- 
                    if(identical(d, floor(d)))
                    {
                        dNumber <- match(d - 1, modelData$splineDegrees)
                        modelData$rho.list[[s]][dNumber]
                    }
                    else
                    {
                        getRhos(lambdas=modelData$lambdas.list[[s]],
                                degrees=d - 1)
                    }
                
                ## positions of rho and spline basis columns
                positions <- (i - 1) * modelData$dimSplineBasis + (1:modelData$dimSplineBasis)
                
                ## put rho into the diagonal (dimSplineBasis times)
                diagRho[positions] <- thisRho 
                
                ## put the spline basis into the design matrix
                splineDesign[, positions] <- modelData$Z.list[[s]] 
            }

            ## so we can construct the precision matrix:
            precMatrix <- crossprod(splineDesign)
            diag(precMatrix) <- diag(precMatrix) + 1 / diagRho

            ## make a Cholesky decomposition of it
            precRoot <- chol(precMatrix)

            ## then normal variate generation:

            ## first N(0, I_{dim.spline})
            simCoefs <- matrix(data=rnorm(n=dim.spline * nSamples),
                               nrow=dim.spline,
                               ncol=nSamples)

            ## then N(0, sigma2 * I_{dim.spline})
            simCoefs <- sweep(simCoefs,
                              MARGIN=2,
                              STATS=sqrt(samples.sigma2),
                              FUN="*")
            
            ## then N(0, sigma2 * precMatrix^{-1}) 
            simCoefs <- backsolve(r=precRoot, x=simCoefs, k=dim.spline)

            ## compute the mean vector samples:
            w <- modelData$y -
                tcrossprod(rep(1, modelData$nObs), samples.intercept) -
                    modelData$X[, marginal$whichLinear, drop=FALSE] %*% 
                        samples.linearCoefs
            
            means <- forwardsolve(l=precRoot,
                                  x=
                                  crossprod(splineDesign,
                                            w), 
                                  upper.tri=TRUE,
                                  transpose=TRUE) 
            means <- backsolve(r=precRoot, x=means)
            
            ## and finally we have N(means, sigma2 * precMatrix^{-1})
            ## random vectors:
            simCoefs <- simCoefs + means

            ## now sort these into a list corresponding to the
            ## covariates they belong to:
            splineCoefs <- list()
            for(i in seq_len(nSplineTerms))
            {
                ## positions of spline coefficients
                positions <- (i - 1) * modelData$dimSplineBasis + (1:modelData$dimSplineBasis)

                ## put the spline coefficients into the list
                splineCoefs[[splineNames[i]]] <-
                    simCoefs[positions, , drop=FALSE] 
            }

            ## "return" the list
            splineCoefs
        }

    if(! marginal$isNullModel)
    {
        ## convert to list with covariates as names
        samples.linearCoefs <-
            split(samples.linearCoefs,
                  colnames(marginal$X.lin)[row(samples.linearCoefs)])
    }
    
    ## return the list with all samples
    return(list(t=samples.t,
                sigma2=samples.sigma2,
                intercept=samples.intercept,
                linearCoefs=samples.linearCoefs,
                splineCoefs=samples.splineCoefs))
}
