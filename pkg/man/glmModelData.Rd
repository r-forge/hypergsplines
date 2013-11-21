\name{glmModelData}
\alias{glmModelData}
\title{Process the data needed for modelling}
\usage{
  glmModelData(y, X, continuous = rep.int(TRUE, ncol(X)),
    nKnots = 4L, splineType = "linear",
    gPrior = HypergnPrior(a = 4, n = length(y)),
    weights = rep.int(1L, length(y)),
    offsets = rep.int(0L, length(y)), family = gaussian,
    phi = 1)
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

  \item{gPrior}{A g-prior class object. Defaults to a
  hyper-g/n prior. See \code{\linkS4class{GPrior}} for more
  information. Deprecated but still possible for
  backwards-compatibility is the use of the strings
  \dQuote{hyper-g/n} or \dQuote{hyper-g}.}

  \item{weights}{optionally a vector of positive weights
  (if not provided, a vector of ones)}

  \item{offsets}{this can be used to specify an _a priori_
  known component to be included in the linear predictor
  during fitting. This must be a numeric vector of length
  equal to the number of cases (if not provided, a vector
  of zeroes)}

  \item{family}{distribution and link (as in the glm
  function)}

  \item{phi}{value of the dispersion parameter (defaults to
  1)}
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
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

