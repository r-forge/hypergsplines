\name{exhaustive}
\alias{exhaustive}
\title{Evaluate all possible models}
\usage{
  exhaustive(modelData,
    modelPrior = c("flat", "exponential", "independent", "dependent", "dep.linear"),
    modelConfigs = NULL, algorithm = c("2", "1"),
    computation = getComputation(), order = FALSE)
}
\arguments{
  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{modelData}} or
  \code{\link{glmModelData}}}

  \item{modelPrior}{either \dQuote{flat} (default),
  \dQuote{exponential}, \dQuote{independent},
  \dQuote{dependent}, or \dQuote{dep.linear}, see
  \code{\link{getLogModelPrior}} for details.}

  \item{modelConfigs}{optional matrix of model
  configurations, which are then evaluated instead of all
  possible configurations. It is check for coherency with
  \code{modelData}.}

  \item{algorithm}{either \dQuote{2} (default, fast for
  small dimension of spline coefficients) or \dQuote{1}
  (fast for small number of observations), specifying the
  algorithm version. Only matters for normal models, for
  GLMs always a type 2 algorithm is used.}

  \item{computation}{computation options produced by
  \code{\link{getComputation}}, only matters for
  generalised response models.}

  \item{order}{should the models be ordered after their
  posterior probability? (default: \code{FALSE} for
  backwards compatibility)}
}
\value{
  a list with the data frame \dQuote{models} comprising the
  model configurations, (R2 for normal models) / log
  marginal likelihoods / log priors / posteriors; and the
  inclusion probabilities matrix \dQuote{inclusionProbs}.
}
\description{
  Evaluate all possible models
}
\examples{
## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get a list of all possible models with this data
res <- exhaustive(md)

res

## now the same, but with cubic splines and algorithm 1:

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                splineType="cubic")

## get a list of all possible models with this data
res <- exhaustive(md,
                  algorithm="1")

res

## now only compute for two certain model configurations:
configs <- cbind(GNP=c(1L, 3L),
                 Armed.Forces=c(2L, 3L))
res <- exhaustive(md,
                  modelConfigs=configs)

## now for generalised response:

## get the model data
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=cbind(GNP, Armed.Forces),
                   family=binomial)

## and do the exhaustive search
res <- exhaustive(md,
                  modelPrior="dependent",
                  computation=
                  getComputation(higherOrderCorrection=FALSE,
                                 debug=FALSE))

res$models <- res$models[order(res$models$post, decreasing=TRUE), ]
res

res1 <- exhaustive(md,
                  computation=
                  getComputation(higherOrderCorrection=FALSE,
                                 debug=FALSE))
res2 <- exhaustive(md,
                  computation=
                  getComputation(higherOrderCorrection=FALSE,
                                 debug=FALSE))
res3 <- exhaustive(md,
                  computation=
                  getComputation(higherOrderCorrection=FALSE,
                                 debug=TRUE))

str(res1)
identical(res1, res2)
identical(res1, res3)


## now with offsets:
set.seed(93)
offsets <- rnorm(n=length(Employed))

md <- glmModelData(y=round(Employed / 10),
                   X=cbind(GNP, Armed.Forces),
                   family=poisson,
                   offsets=offsets)

res <- exhaustive(md,
                  computation=
                  getComputation(higherOrderCorrection=TRUE,
                                 debug=TRUE))
res




res <- exhaustive(md,
                  computation=
                  getComputation(higherOrderCorrection=TRUE,
                                 debug=TRUE))
res
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

