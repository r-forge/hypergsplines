#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[postOptimize.R] by DSB Don 28/07/2011 15:20 (CEST)>
##
## Description:
## After finding a good model, do an "optim" run on the
## smoothing parameters of the continuous variables, as a form of
## postprocessing.
##
## History:
## 21/06/2011   file creation
## 22/07/2011   use svd() instead of eigen(crossprod())
## 26/07/2011   larger range for rho optimisation
## 28/07/2011   now also works for GAMs
#####################################################################################

##' Post-optimizing a good model
##'
##' After finding a good model, do an "optim" of the marginal likelihood with
##' respect the smoothing parameters of the continuous variables, as a form of
##' postprocessing.
##'
##' @param modelData the data necessary for model estimation, which is the
##' result from \code{\link{modelData}} or \code{\link{glmModelData}}
##' @param modelConfig the model configuration
##' @param computation computation options produced by
##' \code{\link{getComputation}}, only matters for generalised response models.  
##' @return the optimized model configuration, which is non-integer for included
##' continuous variables
##'
##' @example examples/postOptimize.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
postOptimize <- function(modelData,
                         modelConfig,
                         computation=getComputation())
{
    ## which parameters do we need to optimize?
    whichOptimize <- which((modelConfig > 0) & (modelData$continuous))
    nParameters <- length(whichOptimize)
    if(! (nParameters > 0))
    {
        warning("no parameters to be optimized")
        return(modelConfig)
    }
    
    ## translate the start degrees of freedom into rhos
    startRhos <- numeric(nParameters)
    for(i in seq_len(nParameters))
    {
        startRhos[i] <-
            getRhos(lambdas=modelData$lambdas.list[[whichOptimize[i]]],
                    degrees=
                    max(modelConfig[whichOptimize[i]] - 1L,
                        0.01))
        
    }

    ## function to translate back into degrees of freedom:
    getDegrees <- function(rhos)
    {
        degrees <- modelConfig

        for(i in seq_len(nParameters))
        {
            lambdas <- modelData$lambdas.list[[whichOptimize[i]]]
            degrees[whichOptimize[i]] <-
                sum(lambdas / (lambdas + 1 / rhos[i])) + 1L 
        }

        return(degrees)
    }
    
    
    ## define the function to be optimized
    ## depending on if this is a normal model
    target <- 
        if(is.null(modelData$family))
            ## the normal response case.
            function(rhos)
            {
                config <- as.double(getDegrees(rhos))
                .Call("cpp_logMargLik",
                      config,
                      modelData)
            }
        else
            ## the generalised response case.
            function(rhos)
            {
                config <- as.double(getDegrees(rhos))
                .Call("cpp_glmLogMargLik",
                      config,
                      modelData,                      
                      computation)                
            }

    ## now do the optim run
    optimResult <- optim(par=startRhos,
                         fn=target,
                         method="L-BFGS-B",
                         lower=
                         rep(.Machine$double.xmin,
                             nParameters),
                         upper=
                         rep(.Machine$double.xmax,
                             nParameters),
                         control=list(fnscale=-1))

    ## check convergence
    if(optimResult$convergence != 0L)
    {
        stop(optimResult$message)
    }

    ## translate result into degrees of freedom
    finalRhos <- optimResult$par
    finalDegrees <- getDegrees(finalRhos)
    
    ## return the vector of optimized degrees of freedom
    return(finalDegrees)     
}
