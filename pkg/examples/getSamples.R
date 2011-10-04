## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get posterior samples for a specific model configuration
res <- getSamples(config=c(2, 2),
                  nSamples=1000L,
                  modelData=md)
str(res)
