\name{dinvGamma}
\alias{dinvGamma}
\title{Density function for the inverse gamma distribution}
\usage{
  dinvGamma(x, a, b, log = FALSE, normalize = TRUE)
}
\arguments{
  \item{x}{vector of values in the positive support}

  \item{a}{shape parameter}

  \item{b}{rate parameter}

  \item{log}{logical; if \code{TRUE}, the log density is
  returned}

  \item{normalize}{normalize the density function?
  (default)}
}
\value{
  the density values at \code{x}
}
\description{
  Density function for the inverse gamma distribution
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{distribution}

