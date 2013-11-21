\name{survModelData}
\alias{survModelData}
\title{Process the data needed for survival models}
\usage{
  survModelData(times, X, observed,
    continuous = rep.int(TRUE, ncol(X)), nKnots = 4L,
    splineType = "linear",
    gPrior = HypergnPrior(a = 4, n = sum(observed)))
}
\arguments{
  \item{times}{the numeric vector of survival times}

  \item{X}{the numeric matrix of covariates (not including
  time)}

  \item{observed}{the logical vector of observation
  indicators, \code{TRUE} entries represent truly observed
  survival times, \code{FALSE} entries represent censored
  survival times.}

  \item{continuous}{see \code{\link{glmModelData}} for
  details}

  \item{nKnots}{see \code{\link{glmModelData}} for details}

  \item{splineType}{see \code{\link{glmModelData}} for
  details}

  \item{gPrior}{see \code{\link{glmModelData}} for details.
  Defaults to the hyper-g/n prior where n is chosen as the
  number of events, instead of the number of observations.}
}
\value{
  a list with the internally needed results.
}
\description{
  Process the data needed for survival models
}
\examples{
## get some data
library(survival)
pbc <- na.omit(pbc)
pbc$sex <- as.numeric(pbc$sex == "f")

## try the function
md <- survModelData(times=pbc$time,
                    X=
                    as.matrix(subset(pbc,
                                     select=
                                     c(trt,
                                       age,
                                       sex,
                                       ascites,
                                       hepato,
                                       spiders,
                                       edema,
                                       bili,
                                       chol,
                                       albumin,
                                       copper,
                                       alk.phos,
                                       ast,
                                       trig,
                                       platelet,
                                       protime,
                                       stage))),
                    observed=
                    pbc$status == 2,
                    continuous=
                    c(FALSE,
                      TRUE,
                      FALSE,
                      FALSE,
                      FALSE,
                      FALSE,
                      FALSE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      TRUE,
                      FALSE),
                    nKnots=6L,
                    splineType="cubic",
                    gPrior="hyper-g/n")

## look at the results
str(md)
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

