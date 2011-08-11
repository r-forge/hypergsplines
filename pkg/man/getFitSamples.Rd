\name{getFitSamples}
\alias{getFitSamples}
\title{Extract fit samples on linear predictor scale...}
\usage{getFitSamples(X, samples, modelData)}
\description{Extract fit samples on linear predictor scale
from posterior coefficients samples}
\value{the fit samples as a matrix.}
\keyword{regression}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{X}{the numeric matrix with new covariate values on the
\emph{original} scale, with the same column layout as originally
provided to \code{\link{modelData}} or \code{\link{glmModelData}}.}
\item{samples}{the samples object (either from \code{\link{getSamples}} or 
the \code{samples} element from \code{\link{glmGetSamples}})}
\item{modelData}{the corresponding model data object}}
\examples{## get some data
attach(longley)

## get model data
X <- cbind(GNP, Armed.Forces)
md <- modelData(y=Employed,
                X=X)

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)

## and then get fit samples at the original X:
res <- getFitSamples(X=X,
                     samples=samples,
                     modelData=md)
str(res)}
