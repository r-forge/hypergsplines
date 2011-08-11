\name{rshrinkage}
\alias{rshrinkage}
\title{Sample from the model-specific posterior of the shrinkage factor...}
\usage{rshrinkage(n, R2, nObs, p, alpha)}
\description{Sample from the model-specific posterior of the shrinkage factor}
\details{Sample from the model-specific posterior of the shrinkage factor
t = g / (1 + g), using inverse sampling.}
\value{\code{n} posterior shrinkage factor samples}
\keyword{distribution}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{n}{number of samples}
\item{R2}{coefficient of determination in the model}
\item{nObs}{number of observations used to fit the model}
\item{p}{number of design matrix columns without counting
the intercept}
\item{alpha}{used hyperparameter for hyper-g prior}}
