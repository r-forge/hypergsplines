#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Hyper-g Splines
##        
## Time-stamp: <[linApproxDens.R] by DSB Fre 18/03/2011 17:02 (CET)>
##
## Description:
## Test the C++ object class LinApproxDens
##
## History:
## 18/03/2011   file creation
#####################################################################################

##' Test the C++ object class LinApproxDens
##'
##' @param args the x values
##' @param logDens the unnormalised log density values
##' @param grid the grid on which to evaluate the linear spline approximation
##' @param nSamples number of samples to produce
##' 
##' @return list with approximate density values and samples
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
linApproxDens <- function(args,
                          logDens,
                          grid,
                          nSamples)
{
    .Call("cpp_linApproxDens",
          args,
          logDens,
          grid,
          nSamples)
}
