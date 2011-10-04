## get some data
attach(longley)

## try the function
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=cbind(GNP, Armed.Forces),
                   family=binomial)

## look at the results
str(md)


## try again with cubic splines
md <- glmModelData(y=as.numeric(Employed > 64),
                   X=cbind(GNP, Armed.Forces),
                   nKnots=10L,
                   splineType="cubic",
                   family=binomial)

## look at the results
str(md)
