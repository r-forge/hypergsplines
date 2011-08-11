\name{modelData}
\alias{modelData}
\title{Process the data needed for modelling...}
\usage{modelData(y, X, continuous=rep.int(TRUE, ncol(X)), nKnots=4,
    splineType="linear", a=4)}
\description{Process the data needed for modelling}
\value{a list with the internally needed results.}
\keyword{regression}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{y}{the numeric response vector}
\item{X}{the numeric matrix of covariates}
\item{continuous}{logical vector specifying which covariates
really are continous and can be included nonlinearly in a model
(default: all covariates are continuous)}
\item{nKnots}{number of (quantile-based) spline knots (default: 4)}
\item{splineType}{type of splines to be used (default: \dQuote{linear}),
see \code{\link{makeBasis}} for possible types.}
\item{a}{hyperparameter for the hyper-g prior (default: 4)}}
\examples{## get some data
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
str(md)}
