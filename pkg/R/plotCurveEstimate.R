#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: hypergsplines
##        
## Time-stamp: <[plotCurveEstimate.R] by DSB Mon 04/04/2011 16:02 (CEST)>
##
## Description:
## Plot function value estimates.
##
## History:
## 29/09/2010   file creation
## 30/09/2010   add partial residual option (which is on by default)
## 13/10/2010   add better default ylim, so that all partial residuals are shown
## 29/10/2010   update to better match bfp plots
## 11/01/2011   adapt to change of samples$linearCoefs
## 21/01/2011   add new options "hpd" and "estimate" to also allow quantile-based
##              (and only quantile-based) summaries of the samples
## 04/04/2011   Note that we can also use this function for GLMs! Only
##              difference is that we cannot plot partial residuals.
#####################################################################################

##' @include getFunctionSamples.R
##' @include getFitSamples.R
##' @include hpds.R
{}


##' Plot function value estimates
##'
##' @param covName string with the name of the covariate
##' @param samples the samples object (either from \code{\link{getSamples}} or 
##' the \code{samples} element from \code{\link{glmGetSamples}})
##' @param modelData the corresponding model data object
##' @param nGrid number of abscissa values for the grid (default: 200)
##' @param plevel credible level for the pointwise credible intervals (default:
##' 0.95, and \code{NULL} suppresses it)
##' @param slevel credible level for simultaneous credible band (defaults to
##' \code{plevel}, and \code{NULL} suppresses it)
##' @param partialResids add partial residuals to the plot? (default,
##' is only possible for normal models)
##' @param hpd use HPD intervals / bands? (default) Otherwise equi-tailed
##' intervals / bands are computed and plotted.
##' @param estimate type of the estimated curve
##' @param ylim y axis limits (has a sensible default to include all points in
##' the plot) 
##' @param lty line type for (1) mean curve and (2, 3) credible interval bounds
##' (default: \code{c(1, 2, 2, 2, 2)})
##' @param col line color(s) (default: black, blue, blue, green, green)
##' @param xlab x axis label (default: \code{covName})
##' @param ylab y axis label (default: f(\code{covName}) )
##' @param \dots \dots additional plotting parameters
##' @return currently nothing.
##'
##' @example hypergsplines/examples/plotCurveEstimate.R
##' 
##' @export 
##' @keywords hplot
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
plotCurveEstimate <- function(covName,
                              samples,
                              modelData,
                              nGrid=200L,
                              plevel=0.95,
                              slevel=plevel,
                              partialResids=TRUE,
                              hpd=TRUE,
                              estimate=c("mean", "median"),
                              ...,
                              ylim=NULL,
                              lty=c(1, 2, 2, 2, 2),
                              col=
                              c("black",
                                "blue",
                                "blue",
                                "green",
                                "green"),
                              xlab=covName,
                              ylab=
                              paste("f(", covName, ")",
                                    sep=""))
{
    ## checks and extracts:
    estimate <- match.arg(estimate)
    stopifnot(is.character(covName),
              identical(length(covName), 1L),
              covName %in% names(samples$linearCoefs),
              nGrid > 1L,
              is.logical(hpd) && identical(length(hpd), 1L))

    ## get a grid of x values, on the original scale
    x <- modelData$origX[, covName]
    
    x.grid <- seq(from=min(x),
                  to=max(x),
                  length=nGrid)

    ## get the samples
    curveSamples <- getFunctionSamples(x=x.grid,
                                       covName=covName,
                                       samples=samples,
                                       modelData=modelData)

    ## which lines are we going to plot?
    lineData <-
        if(estimate == "mean")
            list(mean=rowMeans(curveSamples))
        else ## estimate == "median"
            list(median=
                 apply(curveSamples,
                       1L,
                       median))
    
    ## get the (optional) pointwise and simultaneous bounds:
    if(! is.null(plevel))
    {
        
        stopifnot(identical(length(plevel), 1L),
                  plevel > 0,
                  plevel < 1)
        
        tmp <-
            if(hpd)
                apply(curveSamples,
                      1L,
                      empiricalHpd,
                      level=plevel)
            else
                apply(curveSamples,
                      1L,
                      quantile,
                      p=
                      c((1 - plevel) / 2,
                        (1 + plevel) / 2))

        lineData$plower <- tmp[1L, ]
        lineData$pupper <- tmp[2L, ]
    }

    if(! is.null(slevel))
    {
        
        stopifnot(identical(length(slevel), 1L),
                  slevel > 0,
                  slevel < 1)
        
        tmp <-
            if(hpd)
                scrHpd(curveSamples,
                       level=slevel,
                       mode=lineData$mean)
            else
                scrBesag(curveSamples,
                         level=slevel)
        
        lineData$slower <- tmp[, "lower"]
        lineData$supper <- tmp[, "upper"]
    }

    lineData <- as.data.frame(lineData)

    ## partial residuals computation?
    isNormalModel <- is.null(modelData$family)
    if(isNormalModel && partialResids)
    {
        ## get the residuals:
        resids <- modelData$y -
            rowMeans(getFitSamples(X=modelData$origX,
                                   samples=samples,
                                   modelData=modelData))

        ## get estimated function values at the original x's:
        funEstimates <- rowMeans(getFunctionSamples(x=x,
                                                    covName=covName,
                                                    samples=samples,
                                                    modelData=modelData))

        ## so what are the partial residuals?
        pResids <- funEstimates + resids

        ## get a sensible default ylim
        if(is.null(ylim))
        {
            ylim <- range(c(lineData, pResids))
        }
    }
    
    ## setup the plot
    matplot(x=x.grid,
            y=lineData,
            type="n",                   # empty!
            ...,
            ylim=ylim,
            lty=lty,
            col=col,
            xlab=xlab,
            ylab=ylab)

    ## optionally add the partial residuals
    if(isNormalModel && partialResids)
    {
        points(x,
               pResids,
               col="gray",
               cex=0.5)        
    }

    ## plot the lines now to avoid overplotting by the points
    matplot(x=x.grid,
            y=lineData,
            type="l",                   # draw the lines here
            ...,
            ylim=ylim,
            lty=lty,
            col=col,
            xlab=xlab,
            ylab=ylab,
            add=TRUE)
}
