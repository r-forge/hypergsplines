\name{getInclusionProbs}
\alias{getInclusionProbs}
\title{Extract posterior inclusion probabilities from the models table}
\usage{
  getInclusionProbs(models, modelData,
    logMargLiks = models$logMargLik,
    logPriors = models$logPrior)
}
\arguments{
  \item{models}{the models table (a data frame), the result
  from \code{\link{exhaustive}}}

  \item{modelData}{data used for model estimation}

  \item{logMargLiks}{vector of log marginal likelihoods
  (defaults to the \code{logMargLik} column of
  \code{models})}

  \item{logPriors}{vector of log prior model probabilities
  (defaults to the \code{logPrior} column of
  \code{models})}
}
\value{
  a nice matrix with the probabilities for exclusion,
  linear and non-linear inclusion of the covariates
}
\description{
  Extract posterior inclusion probabilities from the models
  table
}
\examples{
## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get a list of all possible models with this data
models <- exhaustive(md)

## attach log prior probabilities
models$logPrior <- apply(models[, 1:2],
                         1L,
                         getLogModelPrior,
                         type="exponential",
                         modelData=md)

## then we can compute the inclusion probabilities
getInclusionProbs(models=models,
                  modelData=md)
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

