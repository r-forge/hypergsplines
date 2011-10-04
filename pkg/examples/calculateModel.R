## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## calculate a specific model
res <- calculateModel(config=c(2, 3),
                      modelData=md)

## look at the result
str(res)
