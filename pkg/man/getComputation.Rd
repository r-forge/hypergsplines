\name{getComputation}
\alias{getComputation}
\title{Collect the computation options...}
\usage{getComputation(verbose=TRUE, debug=FALSE, nGaussHermite=20,
    useOpenMP=TRUE, higherOrderCorrection=TRUE)}
\description{Collect the computation options}
\keyword{utilities}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{verbose}{should information on computation progress be given?
(default)}
\item{debug}{print debugging information? (not default)}
\item{nGaussHermite}{number of quantiles used in Gauss Hermite quadrature
for marginal likelihood approximation (and later in the MCMC sampler for the
approximation of the marginal covariance factor density).}
\item{useOpenMP}{shall OpenMP be used to accelerate the computations?
(default)}
\item{higherOrderCorrection}{should a higher-order correction of the
Laplace approximation be used for models with canonical links? (default)}}
