## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get a list of all possible models with this data
models <- exhaustive(md,
                     modelPrior="exponential")$models

## attach log prior probabilities
models$logPrior2 <- apply(models[, 1:2],
                          1L,
                          getLogModelPrior,
                          type="exponential",
                          modelData=md)

stopifnot(all.equal(models$logPrior,
                    models$logPrior2))

## then we can compute the inclusion probabilities
getInclusionProbs(models=models,
                  modelData=md)
