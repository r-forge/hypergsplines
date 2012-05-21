#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[exhaustive.R] by DSB Don 10/05/2012 18:32 (CEST)>
##
## Description:
## Implement the first function in this package. This evaluates all possible
## models for given hyperparameters, and returns a data frame / list with
## the corresponding marginal likelihood values.
##
## History:
## 10/09/2010   file creation
## 13/09/2010   first complete version
## 15/09/2010   next very compact version which avoids code duplication through 
##              intermediate functions.
## 30/09/2010   add column names to modelConfigs
## 13/10/2010   include progress bar
## 22/10/2010   adapt for the possibility of binary non-continuous covariates
## 04/04/2011   Cleanup:
##              Include here only the general C++ functions. The old function
##              can be found in the archive packages.
## 10/06/2011   Add optional argument "modelConfigs"
## 29/06/2011   Catch errors in expand.grid, which stem from too large model
##              spaces.
## 10/05/2012   use new helper function "checkModelConfigs"
#####################################################################################

##' @include modelData.R
##' @include options.R
##' @include helpers.R
{}

##' Evaluate all possible models
##' 
##' @param modelData the data necessary for model estimation, which is the
##' result from \code{\link{modelData}} or \code{\link{glmModelData}}
##' @param modelConfigs optional matrix of model configurations, which are then
##' evaluated instead of all possible configurations. It is check for coherency
##' with \code{modelData}.
##' @param algorithm either \dQuote{2} (default, fast for small dimension of
##' spline coefficients) or \dQuote{1} (fast for small number of observations),
##' specifying the algorithm version. Only matters for normal models, for GLMs
##' always a type 2 algorithm is used.
##' @param computation computation options produced by
##' \code{\link{getComputation}}, only matters for generalised response models. 
##' @return a data frame with model specifications and marginal likelihood
##' values
##'
##' @example examples/exhaustive.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
exhaustive <- function(modelData,
                       modelConfigs=NULL,
                       algorithm=c("2", "1"),
                       computation=getComputation())
{
    ## check if modelConfigs is provided
    if(is.null(modelConfigs))
    {
        ## if not, generate the data frame of all model configurations:
        linearDegrees <- 0:1

        possibilities <- list()
        for(j in seq_len(modelData$nCovs))
        {
            possibilities[[j]] <-
                if(modelData$continuous[j])
                    modelData$degrees
                else
                    linearDegrees
        }
        
        modelConfigs <-
            tryCatch(as.matrix(expand.grid(possibilities)),
                     error=
                     function(e){
                         print(e)
                         nDegs <- length(modelData$degrees)
                         nCont <- sum(modelData$continuous)
                         nBin <- sum(!modelData$continuous)
                         nModels <- (nDegs^nCont) * (2^nBin) 
                         myError <-
                             simpleError(paste("Too many possible models (",
                                               nModels, "), please supply ",
                                               "an argument modelConfigs with ",
                                               "a smaller subset of models.",
                                               sep=""))
                         stop(myError)
                     })

        colnames(modelConfigs) <- colnames(modelData$X)    
    }
    else
    {
        ## if it is provided, check for consistency
        checkModelConfigs(modelConfigs=modelConfigs,
                          modelData=modelData)
    }

    
    ## decide if this is a normal model
    isNormalModel <- is.null(modelData$family)
    
    ## go C++
    ret <-
        if(isNormalModel)
        {
            algorithm <- match.arg(algorithm)
            
            if(algorithm == 1L)
                .Call("cpp_exhaustive", ## slow for large n
                      modelData,
                      modelConfigs)
            else ## if(algorithm == 2L)
                .Call("cpp_exhaustive2", ## slow for large p * dimSplineBasis
                      modelData,
                      modelConfigs)    
        }
        else
        {
            .Call("cpp_glmExhaustive",
                  modelData,
                  modelConfigs,
                  computation)
        }
        
    ## return the model configurations along with the results
    return(cbind(modelConfigs,
                 ret))
}

