## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get a list of all possible models with this data
exRes <- exhaustive(md)
exRes

## add post prob
exRes$post <- exp(exRes$logMargLik) / sum(exp(exRes$logMargLik))

## get meta-model table
aggRes <- aggregateModelsTable(modelsTable=exRes[, 1:2],
                               posterior=exRes$post)

aggRes$metaProb

## the top meta model
topMeta <- names(aggRes$metaProb)[1L]
topMeta

## models corresponding to top meta-model
exRes[aggRes$metaConfig == topMeta, ]


