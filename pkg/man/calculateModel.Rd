\name{calculateModel}
\alias{calculateModel}
\title{Calculate intermediate information for a specific model}
\usage{
  calculateModel(config, modelData)
}
\arguments{
  \item{config}{the model configuration vector}

  \item{modelData}{the result from \code{\link{modelData}}}
}
\value{
  A list with necessary intermediate information
}
\description{
  Calculate intermediate information for a specific model
}
\examples{
## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## calculate a specific model
res <- calculateModel(config=c(2, 3),
                      modelData=md)

## look at the result
str(res)
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{internal}
\keyword{regression}

