\name{getBmaSamples}
\alias{getBmaSamples}
\title{Get posterior samples from the Bayesian Model Average (BMA)}
\usage{
  getBmaSamples(config, logPostProbs, nSamples, modelData,
    mcmc = getMcmc(), computation = getComputation())
}
\arguments{
  \item{config}{the data frame/matrix with model
  specifications, e.g. the result from
  \code{\link{exhaustive}}}

  \item{logPostProbs}{vector of log posterior probabilites
  (will be exponentiated and normalized within the
  function) for the weighting of the models in
  \code{config}}

  \item{nSamples}{number of samples to simulate}

  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{modelData}} or
  \code{\link{glmModelData}}}

  \item{mcmc}{MCMC options produced by
  \code{\link{getMcmc}}, only matters for generalised
  response models. Then, the burn-in and thinning
  parameters will be applied for each sampled model.}

  \item{computation}{computation options produced by
  \code{\link{getComputation}}, only matters for
  generalised response models.}
}
\value{
  A list with samples from the shrinkage hyperparameter,
  the regression variance, and the (linear and spline)
  coefficients, analogous to the return value from
  \code{\link{getSamples}} or \code{\link{glmGetSamples}}.
  The only difference is that \dQuote{linearCoefs} and
  \dQuote{splineCoefs} contain zeroes for samples where the
  model did not contain that covariate linearly or
  smoothly. This is necessary to ensure compatibility with
  \code{\link{getFunctionSamples}} and
  \code{\link{getFitSamples}}. Moreover, the model
  specifications matrix is appended with columns
  \dQuote{postProb} and \dQuote{sampleFreq}, containing the
  posterior probability and the sampling frequency,
  respectively.
}
\description{
  Get posterior samples from the Bayesian Model Average
  (BMA)
}
\examples{
## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                gPrior="hyper-g/n")

## get models table
tab <- exhaustive(modelData=md)
tab

## get posterior samples from the BMA assuming
## a flat model prior
res <- getBmaSamples(config=tab,
                     logPostProbs=tab$logMargLik,
                     nSamples=1000L,
                     modelData=md)
str(res)

summary(res$t)
hist(res$t, nclass=100)

## now for generalised response:

## get the model data
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=cbind(GNP, Armed.Forces),
                   family=binomial)

## get models table
tab <- exhaustive(modelData=md,
                  computation=
                  getComputation(higherOrderCorrection=FALSE))

## get posterior samples from the BMA assuming
## a flat model prior
res <- getBmaSamples(config=tab,
                     logPostProbs=tab$logMargLik,
                     nSamples=1000L,
                     modelData=md,
                     mcmc=
                     getMcmc(burnin=10L,
                             step=1L),
                     computation=
                     getComputation(higherOrderCorrection=FALSE))
str(res)

hist(res$t, nclass=100)
res$config
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

