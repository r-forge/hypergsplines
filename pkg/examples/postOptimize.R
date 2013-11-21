## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                splineType="cubic",
                gPrior="hyper-g/n")

## get a list of all possible models with this data
res <- exhaustive(md)$models

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
                  computation=getComputation(higherOrderCorrection=FALSE))$models

## now optimize the best model (best wrt to marg lik)
bestConfig <- res[which.max(res$logMargLik), 1:2]
bestConfig

## todo!
## optimConfig <- postOptimize(modelData=md,
##                             modelConfig=bestConfig,
##                             computation=getComputation(higherOrderCorrection=FALSE))
## optimConfig
