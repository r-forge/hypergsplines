\name{getFunctionSamples}
\alias{getFunctionSamples}
\title{Extract function values samples from posterior coefficients samples...}
\usage{getFunctionSamples(x, covName, samples, modelData)}
\description{Extract function values samples from posterior coefficients samples}
\value{the function values samples as a matrix.}
\keyword{regression}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{x}{the vector of abscissa values (on the \emph{original} scale!)}
\item{covName}{string with the name of the covariate}
\item{samples}{the samples object (either from \code{\link{getSamples}} or
the \code{samples} element from \code{\link{glmGetSamples}})}
\item{modelData}{the corresponding model data object}}
\examples{## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)

## and then get function samples:
res <- getFunctionSamples(x=
                          seq(from=min(GNP),
                              to=max(GNP),
                              length=100L),
                          covName="GNP",
                          samples=samples,
                          modelData=md)
str(res)


## get model data with cubic splines
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                nKnots=10L,
                splineType="cubic")

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)

## and then get function samples:
res <- getFunctionSamples(x=
                          seq(from=min(GNP),
                              to=max(GNP),
                              length=100L),
                          covName="GNP",
                          samples=samples,
                          modelData=md)
str(res)}
