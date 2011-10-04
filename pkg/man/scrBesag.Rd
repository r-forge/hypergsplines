\name{scrBesag}
\alias{scrBesag}
\title{Simultaneous credible band computation}
\usage{
  scrBesag(samples, level = 0.95)
}
\arguments{
  \item{samples}{m by n matrix where m is the number of
  parameters, n is the number of samples and hence each
  (multivariate) sample is a column in the matrix
  \code{samples}}

  \item{level}{the credible level (default: 0.95)}
}
\value{
  matrix with \sQuote{lower} and \sQuote{upper} columns
}
\description{
  Calculates a series of simultaneous credible bounds for
  one parameter type.
}
\author{
  Thomas Kneib
}
\references{
  J. Besag, P. Green, D. Higdon, K. Mengersen (1995):
  Bayesian Methods for Nonlinear Classification and
  Regression, \emph{Statistical Science} \bold{10}/1,
  3--66,
  \url{http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.27.8854}
}
\keyword{htest}

