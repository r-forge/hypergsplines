\name{logHyp2f1Laplace}
\alias{logHyp2f1Laplace}
\title{Laplace approximation to the Gauss hypergeometric function}
\usage{
  logHyp2f1Laplace(a, b, c, x)
}
\arguments{
  \item{a}{first function parameter}

  \item{b}{second function parameter}

  \item{c}{third function parameter}

  \item{x}{vector of abscissae with absolute value not
  greater than 1}
}
\value{
  the logarithm of Laplace approximated values
}
\description{
  This is an R interface to the C function implementing the
  Laplace approximation to the Gauss hypergeometric
  function. Note that all parameters must be numeric and
  must not be complex. Only \code{x} may be a vector.
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{math}

