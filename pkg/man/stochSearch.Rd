\name{stochSearch}
\alias{stochSearch}
\title{Stochastic search for models with high posterior probability}
\usage{
  stochSearch(modelData,
    modelPrior = c("flat", "exponential", "independent", "dependent", "dep.linear"),
    chainlength = 100000L, nCache = chainlength,
    nModels = as.integer(max(nCache/100, 1)),
    computation = getComputation())
}
\arguments{
  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{modelData}} or
  \code{\link{glmModelData}}}

  \item{modelPrior}{either \dQuote{flat} (default),
  \dQuote{exponential}, \dQuote{independent},
  \dQuote{dependent}, or \dQuote{dep.linear}, see
  \code{\link{getLogModelPrior}} for details.}

  \item{chainlength}{length of the model sampling chain
  (default: 100,000)}

  \item{nCache}{maximum number of best models to be cached
  at the same time during the model sampling (by default
  equal to \code{chainlength})}

  \item{nModels}{how many best models should be saved?
  (default: 1\% of the total number of \code{nCache}). Must
  not be larger than \code{nCache}.}

  \item{computation}{computation options produced by
  \code{\link{getComputation}}, only matters for
  generalised response models.}
}
\value{
  a list with the data frame \dQuote{models} comprising the
  model configurations, log marginal likelihoods / priors /
  posteriors and hits in the MCMC run, the inclusion
  probabilities matrix \dQuote{inclusionProbs}, the number
  of total visited models \dQuote{numVisited} and the log
  normalization constant \dQuote{logNormConst}.
}
\description{
  Stochastic search for models with high posterior
  probability
}
\examples{
## get some data
attach(longley)

## get large model data
md <- modelData(y=Employed,
                X=
                cbind(GNP,
                      Armed.Forces,
                      Population,
                      Year))

## do a stochastic search over the model space
res <- stochSearch(md)
res

## now the same, but with cubic splines:

## get large model data
md <- modelData(y=Employed,
                X=
                cbind(GNP,
                      Armed.Forces,
                      Population,
                      Year),
                splineType="cubic")

## do a stochastic search over the model space
res <- stochSearch(md)
res

## and now for generalised response:

## get the model data
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=
                   cbind(GNP,
                         Armed.Forces,
                         Population,
                         Year),
                   family=binomial)

## do a stochastic search over the model space
res <- stochSearch(md,
                   chainlength=1000L,
                   computation=
                   getComputation(higherOrderCorrection=FALSE))
res
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

