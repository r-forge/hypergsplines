## get some data
attach(longley)

## get large model data
md <- modelData(y=Employed,
                X=
                cbind(GNP,
                      Armed.Forces,
                      Population,
                      Year))

## do a stochastic search over the model space
res <- stochSearch(md)
res

## now the same, but with cubic splines:

## get large model data
md <- modelData(y=Employed,
                X=
                cbind(GNP,
                      Armed.Forces,
                      Population,
                      Year),
                splineType="cubic")

## do a stochastic search over the model space
res <- stochSearch(md)
res

## and now for generalised response:

## get the model data
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=
                   cbind(GNP,
                         Armed.Forces,
                         Population,
                         Year),
                   family=binomial)

## do a stochastic search over the model space
res <- stochSearch(md,
                   chainlength=1000L,
                   computation=
                   getComputation(higherOrderCorrection=FALSE))
res
