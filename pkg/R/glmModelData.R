#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[glmModelData.R] by DSB Don 28/07/2011 11:42 (CEST)>
##
## Description:
## Calculate the information from the input data which is used for both
## model search and posterior parameter sampling for GLMs.
##
## History:
## 09/03/2011   start by modifying the normal linear model version.
## 15/03/2011   in the zColMeans construction we also need to coerce it to a
##              vector
## 19/05/2011   add "offsets" argument, which is checked and returned inside the 
##              family list
## 22/07/2011   - use offsets argument to compute intercept estimate in the null
##              model
##              - also use the offsets in constructing the fixed weight matrix
##              - remove dependency on glmNullModelInfo.R which is no longer
##              needed
## 28/07/2011   save lambdas.list, different interface for getRhos
#####################################################################################

##' @include getFamily.R
##' @include makeBasis.R
##' @include getRhos.R
##' @include helpers.R
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
##' @param weights optionally a vector of positive weights (if not provided, a 
##' vector of ones)
##' @param offsets this can be used to specify an _a priori_ known component to
##' be included in the linear predictor during fitting. This must be a numeric
##' vector of length equal to the number of cases (if not provided, a vector of
##' zeroes)
##' @param family distribution and link (as in the glm function)
##' @param phi value of the dispersion parameter (defaults to 1)
##' @return a list with the internally needed results.
##'
##' @example hypergsplines/examples/glmModelData.R
##' 
##' @export 
##' @keywords regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
glmModelData <- function(y,
                         X,
                         continuous=rep.int(TRUE, ncol(X)),
                         nKnots=4L,
                         splineType="linear",
                         a=4,
                         weights=rep.int(1L, length(y)),
                         offsets=rep.int(0L, length(y)),
                         family=gaussian,       
                         phi=1)
{ 
    ## checks and extracts:
    stopifnot(is.vector(y),
              is.numeric(y),
              is.matrix(X),
              is.numeric(X),
              all(! is.na(y)),
              all(! is.na(X)),
              is.logical(continuous),
              is.posInt(nKnots),
              is.numeric(a),
              is.numeric(weights),
              all(weights >= 0),
              is.numeric(offsets))
  
    nKnots <- nKnots[1L]
    a <- a[1L]
    
    nObs <- length(y)
    nCovs <- ncol(X)
    
    stopifnot(identical(nObs,
                        nrow(X)),
              identical(length(continuous),
                        nCovs),
              nKnots > 1L && nKnots < nObs - 2L,
              a > 3 && a <= 4,
              identical(nObs,
                        length(offsets)))

    ## family stuff, similar as in glmBayesMfp:
    family <- getFamily(family, phi)    
    init <- family$init(y=y, weights=weights)
    Y <- init$y
    family$weights <- as.double(init$weights)
    family$offsets <- as.double(offsets)
    family$dispersions <- as.double(family$phi / init$weights) # and zero
                                        # weights ?! 
    family$linPredStart <- as.double(init$linPredStart)
    family$loglik <- function(mu)
    {
        return(- 0.5 * sum(family$dev.resids(y=Y, mu=mu, wt=weights)) / phi)
    }

    ## fit the null model:
    nullModelFit <- glm(Y ~ 1,
                        weights=family$weights,
                        family=family$family,
                        offset=family$offsets)

    intercept <- coef(nullModelFit)["(Intercept)"]

    ## compute fixed weight matrix from the null model:
    weightMatrixEntries <-
        family$mu.eta(intercept + family$offsets)^2 /
        family$variance(family$linkinv(intercept + family$offsets)) /
            family$dispersions 

    ## compute log marginal likelihood in the null model
    ## via the Laplace approximation:
    nullModelLogMargLik <-
        - nullModelFit$deviance / 2 + 0.5 * log(2 * pi * vcov(nullModelFit))
    
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

    ## build the vector of column "means":
    attr(X, "scaled:colMeans") <- 
        xColMeans <-
            (crossprod(weightMatrixEntries, X) /
             sum(weightMatrixEntries))[1L, ]
    
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
        
        ## "center" X[, j]:       
        X[, j] <- X[, j] - xColMeans[j]

        if(continuous[j])
        {
            ## for Z_j we must yet compute and save the column "means"
            ## and cross-products with X[, j]
            attr(Z.list[[j]], "scaled:colMeans") <-
                zColMeans <-
                    (crossprod(weightMatrixEntries, Z.list[[j]]) /
                        sum(weightMatrixEntries))[1L, ]

            attr(Z.list[[j]], "scaled:crossX") <-
                zCrossX <-
                    (crossprod(X[, j] * weightMatrixEntries, Z.list[[j]]) /
                     sum(X[, j]^2 * weightMatrixEntries))[1L, ]

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
    
    ## what are the possible degrees of freedom for each (continuous) covariate?
    splineDegrees <- 1:(dimSplineBasis - 1L)
    
    degrees <- c(0L, 1L, splineDegrees + 1L)

    ## for each continuous covariate, we must compute the 
    ## penalty parameters (rho_j) corresponding to the spline degrees:
    lambdas.list <- list()
    rho.list <- list()
    for(j in seq_len(nCovs))
    {
        lambdas.list[[j]] <-
            if(continuous[j])
                svd(sqrt(weightMatrixEntries) * Z.list[[j]],
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
                numeric(0)
    }
    names(rho.list) <- names(Z.list)

    ## return the list with the results
    return(list(nKnots=nKnots,
                dimSplineBasis=dimSplineBasis,
                a=a,
                nObs=nObs,
                nCovs=nCovs,
                Y=Y,
                X=X,
                family=family,
                weightMatrixEntries=weightMatrixEntries,
                nullModelLogMargLik=nullModelLogMargLik,
                continuous=continuous,
                origX=origX,
                Z.list=Z.list,
                degrees=degrees,
                splineDegrees=splineDegrees,
                rho.list=rho.list,
                lambdas.list=lambdas.list))
}

