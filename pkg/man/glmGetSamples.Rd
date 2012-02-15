\name{glmGetSamples}
\alias{glmGetSamples}
\title{Get posterior samples for a specific model configuration for
generalised response}
\usage{
  glmGetSamples(config, modelData, mcmc = getMcmc(),
    computation = getComputation())
}
\arguments{
  \item{config}{the model configuration vector}

  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{glmModelData}}}

  \item{mcmc}{MCMC options, result from
  \code{\link{getMcmc}}}

  \item{computation}{computation options produced by
  \code{\link{getComputation}}}
}
\value{
  A list with samples from the shrinkage hyperparameter t =
  g / (g + 1) and the (linear and spline) coefficients, and
  further information.
}
\description{
  Get posterior samples for a specific model configuration
  for generalised response
}
\examples{
## get some data
attach(longley)

## get the model data
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=cbind(GNP, Armed.Forces),
                   family=binomial)

## and get samples for one specific model
samples <- glmGetSamples(config=c(2L, 3L),
                         modelData=md,
                         computation=
                         getComputation(higherOrderCorrection=FALSE))

str(samples)

}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

