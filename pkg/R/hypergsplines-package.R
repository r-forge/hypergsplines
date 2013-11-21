#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[hypergsplines-package.R] by DSB Don 06/12/2012 12:14 (CET)>
##
## Description:
## Package description.
##
## History:
## 10/09/2010   file creation
## 13/09/2010   use dynamic library and cpp_hyp2f1 inside it
## 13/10/2010   add cpp_calculateModel
## 15/10/2010   add cpp_exhaustive and cpp_exhaustive2
## 20/10/2010   add cpp_stochSearch
## 10/03/2011   add cpp_glmExhaustive and "generalised linear models" to
##              description
## 18/03/2011   add cpp_glmGetSamples and cpp_linApproxDens
## 04/04/2011   add cpp_glmStochSearch
## 09/05/2011   add cpp_log_hyp2f1_laplace
## 22/06/2011   add cpp_logMargLik, cpp_glmLogMargLik
## 04/08/2011   add cpp_aggregateModelsTable
## 28/03/2012   add cpp_logMargLikHypergn
## 29/03/2012   hyper-g priors, not only one and the hyper-g prior
## 20/09/2012   add imports from methods package
## 06/12/2012   also import "is" from methods package
#####################################################################################

##' Bayesian model selection with penalized splines and hyper-g priors
##'
##' This R-package implements Bayesian model selection with penalized
##' splines, using hyper-g priors for the linear effects coefficients and the
##' regression variance in the normal linear and generalised linear models.
##' 
##' @name hypergsplines-package
##' @aliases hypergsplines
##' @docType package
##'
##' @useDynLib hypergsplines
##' cpp_hyp2f1 cpp_calculateModel cpp_exhaustive cpp_exhaustive2 cpp_stochSearch
##' cpp_glmExhaustive cpp_glmGetSamples cpp_linApproxDens cpp_glmStochSearch
##' cpp_log_hyp2f1_laplace cpp_logMargLik cpp_glmLogMargLik
##' cpp_aggregateModelsTable cpp_logMargLikHypergn
##' @importFrom methods setClass setOldClass setGeneric setMethod representation
##' signature prototype initialize new is
##' 
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
##' @keywords package regression
{}

