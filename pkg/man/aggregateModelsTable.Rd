\name{aggregateModelsTable}
\alias{aggregateModelsTable}
\title{Aggregate the models found by exhaustive or stochastic search}
\usage{
  aggregateModelsTable(modelsTable, posterior, cut = 1L)
}
\arguments{
  \item{modelsTable}{the model configurations in a data
  frame}

  \item{posterior}{the posterior probabilities of the
  models}

  \item{cut}{the (integer) cutpoint (see details, default
  is 1)}
}
\value{
  a list with elements \code{metaConfig} and
  \code{metaProb} containing the meta-model strings for
  each single model and the posterior probabilities of the
  meta-models, respectively.
}
\description{
  This function reduces the model configurations to
  meta-models which distinguish only 0, 1, 2, ...,
  \code{cut} degrees of freedom, that is, degrees of
  freedom greater or equal than \code{cut} are seen as
  identical. The function returns the meta-model strings
  for each single model, and the posterior probabilities of
  these meta-models, ordered from top to bottom.
}
\examples{
## get some data
attach(longley)

## get model data
md <- modelData(y=Employed,
                X=cbind(GNP, Armed.Forces))

## get a list of all possible models with this data
exRes <- exhaustive(md)$models
exRes

## get meta-model table
aggRes <- aggregateModelsTable(modelsTable=exRes[, 1:2],
                               posterior=exRes$post)

aggRes$metaProb

## the top meta model
topMeta <- names(aggRes$metaProb)[1L]
topMeta

## models corresponding to top meta-model
exRes[aggRes$metaConfig == topMeta, ]


}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

