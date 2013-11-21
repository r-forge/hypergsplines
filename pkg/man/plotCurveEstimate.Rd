\name{plotCurveEstimate}
\alias{plotCurveEstimate}
\title{Plot function value estimates}
\usage{
  plotCurveEstimate(covName, samples, modelData,
    nGrid = 200L, plevel = 0.95, slevel = plevel,
    partialResids = TRUE, hpd = TRUE,
    estimate = c("mean", "median"), plot = TRUE, ...,
    ylim = NULL, lty = c(1, 2, 2, 2, 2),
    col = c("black", "blue", "blue", "green", "green"),
    xlab = covName,
    ylab = paste("f(", covName, ")", sep = ""))
}
\arguments{
  \item{covName}{string with the name of the covariate}

  \item{samples}{the samples object (either from
  \code{\link{getSamples}} or the \code{samples} element
  from \code{\link{glmGetSamples}})}

  \item{modelData}{the corresponding model data object}

  \item{nGrid}{number of abscissa values for the grid
  (default: 200)}

  \item{plevel}{credible level for the pointwise credible
  intervals (default: 0.95, and \code{NULL} suppresses it)}

  \item{slevel}{credible level for simultaneous credible
  band (defaults to \code{plevel}, and \code{NULL}
  suppresses it)}

  \item{partialResids}{add partial residuals to the plot?
  (default, is only possible for normal models)}

  \item{hpd}{use HPD intervals / bands? (default) Otherwise
  equi-tailed intervals / bands are computed and plotted.}

  \item{plot}{produce a plot? (default: yes) Otherwise you
  could use the returned information to produce the plot
  yourself}

  \item{estimate}{type of the estimated curve}

  \item{ylim}{y axis limits (has a sensible default to
  include all points in the plot)}

  \item{lty}{line type for (1) mean curve and (2, 3)
  credible interval bounds (default: \code{c(1, 2, 2, 2,
  2)})}

  \item{col}{line color(s) (default: black, blue, blue,
  green, green)}

  \item{xlab}{x axis label (default: \code{covName})}

  \item{ylab}{y axis label (default: f(\code{covName}) )}

  \item{\dots}{\dots additional plotting parameters}
}
\value{
  A list with the elements \code{x.grid} and
  \code{lineData}, containing the x axis grid and the
  different y values for plotting, respectively.
}
\description{
  Plot function value estimates
}
\examples{
## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                gPrior="hyper-g")

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)
str(samples)
summary(samples$t)

## and plot resulting curve estimates:
par(mfrow=c(1, 2))

plotCurveEstimate(covName="GNP",
                  samples=samples,
                  modelData=md)

plotCurveEstimate(covName="Armed.Forces",
                  samples=samples,
                  modelData=md)

plotCurveEstimate(covName="Armed.Forces",
                  samples=samples,
                  modelData=md,
                  hpd=FALSE,
                  estimate="median")
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{hplot}

