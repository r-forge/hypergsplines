\name{pinvGamma}
\alias{pinvGamma}
\title{Cumulative distribution function for the inverse gamma distribution}
\usage{
  pinvGamma(q, a, b, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{q}{vector of quantiles in the positive support}

  \item{a}{shape parameter}

  \item{b}{rate parameter}

  \item{lower.tail}{logical; if \code{TRUE} (default),
  probabilities are P[X <= x], otherwise, P[X > x]}

  \item{log.p}{logical; if \code{TRUE}, the log
  probabilities are returned}
}
\value{
  the cdf values at \code{q}
}
\description{
  Cumulative distribution function for the inverse gamma
  distribution
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{distribution}

