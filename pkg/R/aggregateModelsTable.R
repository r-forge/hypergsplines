#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[aggregateModelsTable.R] by DSB Don 04/08/2011 17:28 (CEST)>
##
## Description:
## Aggregate the models found by exhaustive or stochastic search.
##
## History:
## 04/08/2011   file creation
#####################################################################################

##' Aggregate the models found by exhaustive or stochastic search
##'
##' This function reduces the model configurations to meta-models which
##' distinguish only 0, 1, 2, ..., \code{cut} degrees of freedom, that is,
##' degrees of freedom greater or equal than \code{cut} are seen as identical.
##' The function returns the meta-model strings for each single model,
##' and the posterior probabilities of these meta-models, ordered from
##' top to bottom.
##' 
##' @param modelsTable the model configurations in a data frame
##' @param posterior the posterior probabilities of the models
##' @param cut the (integer) cutpoint (see details, default is 1)
##' @return a list with elements \code{metaConfig} and \code{metaProb}
##' containing the meta-model strings for each single model and the
##' posterior probabilities of the meta-models, respectively.
##'
##' @example hypergsplines/examples/aggregateModelsTable.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
aggregateModelsTable <- function(modelsTable,
                                 posterior,
                                 cut=1L)
{
    ## checks etc.
    stopifnot(is.data.frame(modelsTable),
              all(sapply(modelsTable, is.integer)),
              is.numeric(posterior),
              identical(nrow(modelsTable),
                        length(posterior)),
              cut > 0L)

    ## go C++
    ret <- .Call("cpp_aggregateModelsTable",
                 modelsTable,
                 posterior,
                 as.integer(cut))

    ## sort post probs of meta-models
    ret$metaProb <- sort.int(x=ret$metaProb,
                             decreasing=TRUE,
                             method="quick")
    
    ## and return
    return(ret)
}
