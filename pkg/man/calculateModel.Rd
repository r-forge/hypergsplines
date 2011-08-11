\name{calculateModel}
\alias{calculateModel}
\title{Calculate intermediate information for a specific model...}
\usage{calculateModel(config, modelData)}
\description{Calculate intermediate information for a specific model}
\value{A list with necessary intermediate information}
\keyword{regression}
\keyword{internal}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{config}{the model configuration vector}
\item{modelData}{the result from \code{\link{modelData}}}}
\examples{## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## calculate a specific model
res <- calculateModel(config=c(2, 3),
                      modelData=md)

## look at the result
str(res)}
