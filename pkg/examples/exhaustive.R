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
                  computation=
                  getComputation(higherOrderCorrection=FALSE,
                                 debug=FALSE))
res$logPrior <-
    getLogModelPrior(config=res[, 1:2],
                     type="dependent",
                     modelData=md)

resIncProbs <-
    getInclusionProbs(models=res,
                      modelData=md)

res$unLogPost <- res$logMargLik + res$logPrior
res$post <- exp(res$unLogPost - min(res$unLogPost[is.finite(res$unLogPost)]))
res$post <- res$post / sum(res$post)

res <- res[order(res$unLogPost, decreasing=TRUE), ]

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
