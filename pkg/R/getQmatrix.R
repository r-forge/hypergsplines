#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[getQmatrix.R] by DSB Don 06/12/2012 09:30 (CET)>
##
## Description:
## Compute the Q matrix necessary for the Poisson likelihood offsets for the
## survival modelling. 
##
## History:
## 27/08/2012   file creation
## 28/08/2012   do not require ordered survival times
## 06/12/2012   rename function and compute only the Q matrix here
#####################################################################################

##' Compute the Q matrix necessary for the offsets for survival regression 
##'
##' This function uses Algorithm 1 from our Reader Reaction article
##' to compute the Q matrix.
##' 
##' @param times vector of \eqn{n} ordered survival times, not including the 
##' pseudo-zero! 
##' @return the \eqn{(n+1) \times (n+1)} matrix Q
##'
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getQmatrix <- function(times)
{
    ## checks
    stopifnot(! is.unsorted(times),
              times[1] > 0)

    ## how many times?
    nObs <- length(times)
    
    ## add 0 in front for times
    times <- c(0, times)
    
    ## allocate vectors
    mainDiag <- belowDiag <- offsets <- numeric(nObs + 1)

    ## initialize mainDiag
    mainDiag[1] <- 0

    ## the index to which the last difference was positive
    ## (start with index 1, that is the zero in front)
    lastDiffIndex <- 1

    ## loop
    for(i in 2:(nObs + 1))
    {
        timeDiff <- times[i] - times[i - 1]

        if(timeDiff > 0)
        {
            lastDiffIndex <- i - 1
            
            mainDiag[i] <- 1/2 * timeDiff
            belowDiag[i - 1] <- mainDiag[i - 1] + mainDiag[i]
            
        } else {
            mainDiag[i] <- 1/2 * (times[i] - times[lastDiffIndex])
            belowDiag[i - 1] <- 0
        }        
    }

    ## insert into matrix
    Q <- matrix(data=0,
                nrow=nObs + 1,
                ncol=nObs + 1)

    diag(Q) <- mainDiag
    for(i in seq_len(nObs))
    {
        Q[(i + 1):(nObs + 1), i] <- belowDiag[i]
    }

    return(Q)
}

