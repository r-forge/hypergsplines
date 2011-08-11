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

