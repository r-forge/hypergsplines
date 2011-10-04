\name{getSamples}
\alias{getSamples}
\title{Get posterior samples for a specific model configuration}
\usage{
  getSamples(config, nSamples, modelData)
}
\arguments{
  \item{config}{the model configuration vector}

  \item{nSamples}{number of samples to simulate}

  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{modelData}}}
}
\value{
  A list with samples from the shrinkage hyperparameter t =
  g / (g + 1), the regression variance, and the (linear and
  spline) coefficients.
}
\description{
  Get posterior samples for a specific model configuration
}
\examples{
## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get posterior samples for a specific model configuration
res <- getSamples(config=c(2, 2),
                  nSamples=1000L,
                  modelData=md)
str(res)
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

