#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getBmaSamples.R] by DSB Die 04/10/2011 16:20 (CEST)>
##
## Description:
## Get posterior samples from the Bayesian Model Average (BMA).
##
## History:
## 11/01/2011   file creation
## 12/01/2011   - change to *log*PostProbs argument
##              - correct modelFreqs computation and bind it together with the
##              postProbs to the configuration matrix and return that in the list.
## 04/04/2011   extended to GLM case
## 19/05/2011   now "higherOrderCorrection" instead of
##              "binaryLogisticCorrection" is passed to "getComputation"
## 04/10/2011   include zeroes for models where the covariate is not included.
##              Note that this is necessary for correct computation of fit
##              samples based on a model average. However, this changes the
##              interpretation of the samples, and therefore curve estimates
##              based on these samples: it is no longer conditional on inclusion
##              of the covariate, but marginally over all models, also those
##              not including the covariate.
#####################################################################################


##' @include getSamples.R
{}

##' Get posterior samples from the Bayesian Model Average (BMA)
##'
##' @param config the data frame/matrix with model specifications, e.g.
##' the result from \code{\link{exhaustive}}
##' @param logPostProbs vector of log posterior probabilites (will be
##' exponentiated and normalized within the function) for the weighting of
##' the models in \code{config} 
##' @param nSamples number of samples to simulate
##' @param modelData the data necessary for model estimation, which
##' is the result from \code{\link{modelData}} or \code{\link{glmModelData}}
##' @param mcmc MCMC options produced by
##' \code{\link{getMcmc}}, only matters for generalised response models. Then,
##' the burn-in and thinning parameters will be applied for each sampled model. 
##' @param computation computation options produced by
##' \code{\link{getComputation}}, only matters for generalised response models. 
##' @return A list with samples from the shrinkage hyperparameter, the
##' regression variance, and the (linear and spline) coefficients, analogous to
##' the return value from \code{\link{getSamples}} or
##' \code{\link{glmGetSamples}}. The only difference is that
##' \dQuote{linearCoefs} and \dQuote{splineCoefs} contain zeroes
##' for samples where the model did not contain that covariate linearly or
##' smoothly. This is necessary to ensure compatibility with
##' \code{\link{getFunctionSamples}} and \code{\link{getFitSamples}}. Moreover,
##' the model specifications matrix is appended with columns \dQuote{postProb}
##' and \dQuote{sampleFreq}, containing the posterior probability and the
##' sampling frequency, respectively.
##'
##' @example examples/getBmaSamples.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getBmaSamples <- function(config,
                          logPostProbs,
                          nSamples,
                          modelData,
                          mcmc=getMcmc(),
                          computation=getComputation())
{    
    ## checks and extracts:
    config <- as.matrix(config[, 1:modelData$nCovs])
    nModels <- nrow(config)

    postProbs <- exp(logPostProbs - min(logPostProbs))
    
    stopifnot(all(postProbs >= 0),
              any(postProbs > 0),
              all(is.finite(postProbs)),
              nModels > 1L,
              identical(length(postProbs),
                        nModels),
              is.integer(nSamples),
              nSamples > 0L,
              identical(length(nSamples),
                        1L),
              all(config %in% modelData$degrees),
              identical(ncol(config), modelData$nCovs))

    covnames <- colnames(modelData$X)

    ## decide if this is a normal model
    isNormalModel <- is.null(modelData$family)
        
    ## setup samples containers:

    ## shrinkage hyperparameter
    samples.t <- numeric(0L)

    ## the regression variance
    samples.sigma2 <-
        if(isNormalModel)
            numeric(0L)
        else
            NULL
    
    ## the intercept
    samples.intercept <- numeric(0L)
    ## the coefficients of linear effects
    samples.linearCoefs <- list()
    ## finally the spline coefficients samples
    samples.splineCoefs <- list()
    
    ## Distribute samples to models:
    postProbs <- postProbs / sum(postProbs)

    modelDraws <- sample(seq_len(nModels),
                         size=nSamples,
                         replace=TRUE,
                         prob=postProbs)
    
    modelFreqs <- tabulate(modelDraws,
                           nbins=nModels)

    ## progress bar setup:
    nBars <- 20L
    steps <- floor(seq(from=0, to=100, length=nBars + 1L))[- 1L]

    lastWasTwoDigits <- FALSE
    for(k in seq_len(nBars))
    {
        if(k %% 5L == 0L)
        {
            cat(steps[k])
            
            if(steps[k] >= 10)
                lastWasTwoDigits <- TRUE
        }
        else
        {
            if(lastWasTwoDigits)
                lastWasTwoDigits <- FALSE
            else
                cat("-")
        }
    }
    cat("\n")
    
    for(k in seq_len(nBars))
    {
        if(k %% 5L == 0L)
            cat("|")
        else
            cat("-")
    }
    cat("\n")

    barInterval <- floor(nSamples / nBars)

    
    ## initialize samples counters:
    count <- 0L
    oldIntervalCount <- 0L
    
    ## start sampling:
    for(i in which(modelFreqs > 0L))
    {
        ## sample this model
        thisModel <- config[i, ]
        thisSamples <-
            if(isNormalModel)
                getSamples(config=thisModel,
                           nSamples=modelFreqs[i],
                           modelData=modelData)
            else
                glmGetSamples(config=thisModel,
                              modelData=modelData,
                              mcmc=
                              getMcmc(samples=modelFreqs[i],
                                      burnin=mcmc$burnin,
                                      step=mcmc$step),
                              computation=
                              getComputation(verbose=FALSE,
                                             debug=
                                             computation$debug,
                                             nGaussHermite=
                                             computation$nGaussHermite,
                                             useOpenMP=
                                             computation$useOpenMP,
                                             higherOrderCorrection=
                                             computation$higherOrderCorrection))$samples  
        
        count <- count + modelFreqs[i]

        ## and put results into containers
        samples.t <- c(samples.t,
                       thisSamples$t)
        if(isNormalModel)
        {
            samples.sigma2 <- c(samples.sigma2,
                                thisSamples$sigma2)
        }
        samples.intercept <- c(samples.intercept,
                               thisSamples$intercept)

        ## for all included covariates:
        for(cov in covnames)
        {
            ## if there are linear coefs
            linearCoefs <-
                if(! is.null(thisSamples$linearCoefs[[cov]]))
                    ## also include these
                    thisSamples$linearCoefs[[cov]]
                else
                    ## otherwise include zeroes
                    rep.int(0L, times=modelFreqs[i])
            
            samples.linearCoefs[[cov]] <- c(samples.linearCoefs[[cov]],
                                            linearCoefs)
          
            ## and if there are spline coefs
            splineCoefs <- 
                if(! is.null(thisSamples$splineCoefs[[cov]]))
                    ## also include these
                    thisSamples$splineCoefs[[cov]]
                else
                    ## otherwise include zeroes
                    matrix(data=0,
                           nrow=modelData$dimSplineBasis,
                           ncol=modelFreqs[i])
            
            samples.splineCoefs[[cov]] <- cbind(samples.splineCoefs[[cov]],
                                                splineCoefs)
        }

        ## echo progress of sampling
        newIntervalCount <- floor(count / barInterval)
        if((countDiff <- newIntervalCount - oldIntervalCount) > 0L)
        {
            cat(rep.int("=", countDiff),
                sep="")
        }
        oldIntervalCount <- newIntervalCount
    }
    cat("\n")

    ## bind to configuration matrix
    config <- cbind(config,
                    postProb=postProbs,
                    sampleFreq=modelFreqs)
    
    ## return the list with all samples
    return(list(t=samples.t,
                sigma2=samples.sigma2,
                intercept=samples.intercept,
                linearCoefs=samples.linearCoefs,
                splineCoefs=samples.splineCoefs,
                config=config))
}
