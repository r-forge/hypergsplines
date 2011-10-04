#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[modelData.R] by DSB Don 28/07/2011 11:43 (CEST)>
##
## Description:
## Calculate the information from the input data which is used for both
## model search and posterior parameter sampling.
##
## History:
## 15/09/2010   file creation
## 27/09/2010   - save column means of X and Z's and also the crossproducts of x_j with
##              Z_j because they are needed in the construction of the design matrices
##              for new data.
##              - ensure that X has column names (give it some if it has not)
## 30/09/2010   also save the original, unscaled X matrix
## 14/10/2010   precompute here the tcrossprod's of the Z's. This needs more memory,
##              but avoids re-computing them many times in the exhaustive model
##              space exploration.
## 20/10/2010   include information whether a covariate is only to be included
##              linearly at most, in the logical vector "continuous".
## 22/10/2010   adapt spline construction to this possibility
## 26/10/2010   catch error from getRhos and translate into correct error message:
##              there are too many knots for a covariate (so the covariate has
##              fewer unique values than needed)
## 21/11/2010   add option for the spline type
## 21/01/2011   try splineDegrees up to dimSplineBasis - 1 instead of nKnots - 1,
##              this is a difference for O'Sullivan splines.
## 28/07/2011   save lambdas.list, different interface for getRhos
#####################################################################################

##' @include makeBasis.R
##' @include getRhos.R
{}

##' Process the data needed for modelling
##'
##' @param y the numeric response vector
##' @param X the numeric matrix of covariates
##' @param continuous logical vector specifying which covariates
##' really are continous and can be included nonlinearly in a model
##' (default: all covariates are continuous)
##' @param nKnots number of (quantile-based) spline knots (default: 4)
##' @param splineType type of splines to be used (default: \dQuote{linear}),
##' see \code{\link{makeBasis}} for possible types.
##' @param a hyperparameter for the hyper-g prior (default: 4)
##' @return a list with the internally needed results.
##'
##' @example examples/modelData.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
modelData <- function(y,
                      X,
                      continuous=rep.int(TRUE, ncol(X)),
                      nKnots=4L,
                      splineType="linear",
                      a=4)
{ 
    ## checks and extracts:
    stopifnot(is.vector(y),
              is.numeric(y),
              is.matrix(X),
              is.numeric(X),
              all(! is.na(y)),
              all(! is.na(X)),
              is.logical(continuous),
              is.integer(nKnots),
              is.numeric(a))
   
    nKnots <- nKnots[1L]
    a <- a[1L]
    
    nObs <- length(y)
    nCovs <- ncol(X)
    
    stopifnot(identical(nObs,
                        nrow(X)),
              identical(length(continuous),
                        nCovs),
              nKnots > 1L && nKnots < nObs - 2L,
              a > 3 && a <= 4)

    ## ensure that X has column names
    if(is.null(colnames(X)))
        colnames(X) <- paste("V", seq_len(nCovs),
                             sep="")

    ## save the original X (with column names)
    origX <- X

    ## rescale the covariates to [0, 1]:
    X <- scale(X,
               center=apply(X, 2L, min),
               scale=apply(X, 2L, function(x) diff(range(x))))

    ## what is now the vector of column means of X?
    attr(X, "scaled:colMeans") <-
        xColMeans <-
            colMeans(X)
   
    ## construct the spline basis for each continuous covariate, then
    ## center covariates and spline bases, and make them orthogonal:
    Z.list <- list()
    for(j in seq_len(nCovs))
    {
        ## construct spline basis
        Z.list[[j]] <- 
            if(continuous[j])
                makeBasis(X[, j],
                          type=splineType,
                          nKnots=nKnots)
            else ## input a dummy matrix
                matrix(nrow=0, ncol=0)
        
        ## center X[, j]:
        X[, j] <- X[, j] - xColMeans[j]

        if(continuous[j])
        {
            ## for Z_j we must yet compute and save the column means
            ## and cross-products with X[, j]
            attr(Z.list[[j]], "scaled:colMeans") <-
                zColMeans <-
                    colMeans(Z.list[[j]])

            attr(Z.list[[j]], "scaled:crossX") <-
                zCrossX <-
                    (crossprod(X[, j], Z.list[[j]]) /
                     crossprod(X[, j])[1])[1, ]

            ## then use this for "centering"
            Z.list[[j]] <- Z.list[[j]] -
                tcrossprod(rep(1, nObs), zColMeans) -
                    tcrossprod(X[, j], zCrossX)
        }
    }    
    names(Z.list) <- colnames(X)

    ## get the spline basis dimension (if there is any continuous covariate)
    dimSplineBasis <-
        if(any(continuous))
            ncol(Z.list[[min(which(continuous))]])
        else
            nKnots
    
    ## precompute the tcrossprod of the spline basis for each cov: 
    Z.tcrossprod.list <- lapply(Z.list, tcrossprod)
    
    ## what are the possible degrees of freedom for each (continuous) covariate?
    splineDegrees <- 1L:(dimSplineBasis - 1L)
    
    degrees <- c(0L, 1L, splineDegrees + 1L)

    ## for each continuous covariate, we must compute the 
    ## penalty parameters (rho_j) corresponding to the spline degrees:
    lambdas.list <- list()
    rho.list <- list()
    for(j in seq_len(nCovs))
    {
        lambdas.list[[j]] <-
            if(continuous[j])
                svd(Z.list[[j]],
                    nu=0L,
                    nv=0L)$d^2
            else
                numeric(0L)
            
        rho.list[[j]] <-
            if(continuous[j])
                tryCatch(getRhos(lambdas=lambdas.list[[j]],
                                 degrees=splineDegrees),
                         error=
                         function(e){
                             print(e)
                             myError <- 
                                 simpleError(paste("Too many knots", 
                                                   "for covariate",
                                                   names(Z.list)[j])) 
                             stop(myError)
                         })
            else ## dummy vector
                numeric(0L)
    }
    names(rho.list) <- names(Z.list)

    ## return the list with the results
    return(list(nKnots=nKnots,
                dimSplineBasis=dimSplineBasis,
                a=a,
                nObs=nObs,
                nCovs=nCovs,
                y=y,
                X=X,
                continuous=continuous,
                origX=origX,
                Z.list=Z.list,
                Z.tcrossprod.list=Z.tcrossprod.list,
                degrees=degrees,
                splineDegrees=splineDegrees,
                rho.list=rho.list,
                lambdas.list=lambdas.list))
}

