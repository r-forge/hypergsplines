## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)

## and then get function samples:
res <- getFunctionSamples(x=
                          seq(from=min(GNP),
                              to=max(GNP),
                              length=100L),
                          covName="GNP",
                          samples=samples,
                          modelData=md)
str(res)


## get model data with cubic splines
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                nKnots=10L,
                splineType="cubic")

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)

## and then get function samples:
res <- getFunctionSamples(x=
                          seq(from=min(GNP),
                              to=max(GNP),
                              length=100L),
                          covName="GNP",
                          samples=samples,
                          modelData=md)
str(res)
