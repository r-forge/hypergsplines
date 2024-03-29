\name{qinvGamma}
\alias{qinvGamma}
\title{Quantile function for the inverse gamma distribution}
\usage{
  qinvGamma(p, a, b, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{p}{vector of (log) probabilities}

  \item{a}{shape parameter}

  \item{b}{rate parameter}

  \item{lower.tail}{logical; if \code{TRUE} (default),
  probabilities \code{p} are P[X <= x], otherwise, P[X >
  x]}

  \item{log.p}{logical; if \code{TRUE}, \code{p} are the
  log probabilities}
}
\value{
  the quantiles values
}
\description{
  Quantile function for the inverse gamma distribution
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{distribution}

