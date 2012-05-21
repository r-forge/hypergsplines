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
