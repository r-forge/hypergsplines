## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                gPrior="hyper-g")

## get posterior samples for a specific model configuration
samples <- getSamples(config=c(2, 1),
                      nSamples=1000L,
                      modelData=md)
str(samples)
summary(samples$t)

## and plot resulting curve estimates:
par(mfrow=c(1, 2))

plotCurveEstimate(covName="GNP",
                  samples=samples,
                  modelData=md)

plotCurveEstimate(covName="Armed.Forces",
                  samples=samples,
                  modelData=md)

plotCurveEstimate(covName="Armed.Forces",
                  samples=samples,
                  modelData=md,
                  hpd=FALSE,
                  estimate="median")
