\name{logBFhypergn}
\alias{logBFhypergn}
\title{Compute the log Bayes Factor against the null model under hyper-g/n prior.}
\usage{
  logBFhypergn(n, p, R2)
}
\arguments{
  \item{n}{number of observations}

  \item{p}{number of covariates (excluding the intercept)}

  \item{R2}{coefficient of determination}
}
\value{
  the (double) value of the log BF, NA if the computation
  failed.
}
\description{
  Compute the log Bayes Factor against the null model under
  hyper-g/n prior.
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{internal}
\keyword{math}

