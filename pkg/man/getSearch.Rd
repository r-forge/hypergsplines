\name{getSearch}
\alias{getSearch}
\title{Collect the search options...}
\usage{getSearch(nModels=as.integer(max(1L, floor(nCache/100))),
    nCache=1000000000, chainlength=10000)}
\description{Collect the search options}
\keyword{utilities}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{nModels}{how many best models should be saved? (default: 1\% of the
total number of (cached) models). Must not be larger than \code{nCache}.}
\item{nCache}{maximum number of best models to be cached at the same time
during the model sampling}
\item{chainlength}{length of the model sampling chain}}
