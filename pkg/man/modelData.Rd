\name{modelData}
\alias{modelData}
\title{Process the data needed for modelling}
\usage{
  modelData(y, X, continuous = rep.int(TRUE, ncol(X)),
    nKnots = 4L, splineType = "linear",
    gPrior = c("hyper-g/n", "hyper-g"))
}
\arguments{
  \item{y}{the numeric response vector}

  \item{X}{the numeric matrix of covariates}

  \item{continuous}{logical vector specifying which
  covariates really are continous and can be included
  nonlinearly in a model (default: all covariates are
  continuous)}

  \item{nKnots}{number of (quantile-based) spline knots
  (default: 4)}

  \item{splineType}{type of splines to be used (default:
  \dQuote{linear}), see \code{\link{makeBasis}} for
  possible types.}

  \item{gPrior}{hyperprior for g, either \dQuote{hyper-g/n}
  (default) or \dQuote{hyper-g}.}
}
\value{
  a list with the internally needed results.
}
\description{
  Process the data needed for modelling
}
\examples{
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
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

