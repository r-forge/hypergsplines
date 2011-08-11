\name{glmGetSamples}
\alias{glmGetSamples}
\title{Get posterior samples for a specific model configuration for...}
\usage{glmGetSamples(config, modelData, mcmc=getMcmc(),
    computation=getComputation())}
\description{Get posterior samples for a specific model configuration for
generalised response}
\value{A list with samples from the shrinkage hyperparameter
t = g / (g + 1) and the (linear and spline)
coefficients, and further information.}
\keyword{regression}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{config}{the model configuration vector}
\item{modelData}{the data necessary for model estimation, which
is the result from \code{\link{glmModelData}}}
\item{mcmc}{MCMC options, result from \code{\link{getMcmc}}}
\item{computation}{computation options produced by
\code{\link{getComputation}}}}
\examples{## get some data
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

str(samples)}
