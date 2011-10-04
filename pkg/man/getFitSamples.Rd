\name{getFitSamples}
\alias{getFitSamples}
\title{Extract fit samples on linear predictor scale
from posterior coefficients samples}
\usage{
  getFitSamples(X, samples, modelData)
}
\arguments{
  \item{X}{the numeric matrix with new covariate values on
  the \emph{original} scale, with the same column layout as
  originally provided to \code{\link{modelData}} or
  \code{\link{glmModelData}}.}

  \item{samples}{the samples object (either from
  \code{\link{getSamples}} or the \code{samples} element
  from \code{\link{glmGetSamples}})}

  \item{modelData}{the corresponding model data object}
}
\value{
  the fit samples as a matrix.
}
\description{
  Extract fit samples on linear predictor scale from
  posterior coefficients samples
}
\examples{
## get some data
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
str(res)

## this must also work for BMA samples:
tab <- exhaustive(modelData=md)
samples <- getBmaSamples(config=tab,
                         logPostProbs=tab$logMargLik,
                         nSamples=1000L,
                         modelData=md)
str(samples)
res <- getFitSamples(X=X,
                     samples=samples,
                     modelData=md)
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

