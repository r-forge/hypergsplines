\name{postOptimize}
\alias{postOptimize}
\title{Post-optimizing a good model}
\usage{
  postOptimize(modelData, modelConfig,
    computation = getComputation())
}
\arguments{
  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{modelData}} or
  \code{\link{glmModelData}}}

  \item{modelConfig}{the model configuration}

  \item{computation}{computation options produced by
  \code{\link{getComputation}}, only matters for
  generalised response models.}
}
\value{
  the optimized model configuration, which is non-integer
  for included continuous variables
}
\description{
  After finding a good model, do an "optim" of the marginal
  likelihood with respect the smoothing parameters of the
  continuous variables, as a form of postprocessing.
}
\examples{
## get some data
attach(longley)


## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                splineType="cubic")

## get a list of all possible models with this data
res <- exhaustive(md)

## now optimize the best model (best wrt to marg lik)
bestConfig <- res[which.max(res$logMargLik), 1:2]
bestConfig

optimConfig <- postOptimize(modelData=md,
                            modelConfig=bestConfig)
optimConfig


## now for binary response:

## get the model data
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=cbind(GNP, Armed.Forces),
                   family=binomial)

## and do the exhaustive search
res <- exhaustive(md,
                  computation=getComputation(higherOrderCorrection=FALSE))

## now optimize the best model (best wrt to marg lik)
bestConfig <- res[which.max(res$logMargLik), 1:2]
bestConfig

optimConfig <- postOptimize(modelData=md,
                            modelConfig=bestConfig,
                            computation=getComputation(higherOrderCorrection=FALSE))
optimConfig
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

