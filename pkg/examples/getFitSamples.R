## get some data
attach(longley)

## get model data
X <- cbind(GNP, Armed.Forces)
md <- modelData(y=Employed,
                X=X)

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)

## and then get fit samples at the original X:
res <- getFitSamples(X=X,
                     samples=samples,
                     modelData=md)
str(res)

## this must also work for BMA samples:
tab <- exhaustive(modelData=md)
samples <- getBmaSamples(config=tab,
                         logPostProbs=tab$logMargLik,
                         nSamples=1000L,
                         modelData=md)
str(samples)
res <- getFitSamples(X=X,
                     samples=samples,
                     modelData=md)
