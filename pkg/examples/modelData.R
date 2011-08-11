## get some data
attach(longley)


## try the function
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## look at the results
str(md)


## try again with cubic splines
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces),
                nKnots=10L,
                splineType="cubic")

## look at the results
str(md)
