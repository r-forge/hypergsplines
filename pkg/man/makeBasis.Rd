\name{makeBasis}
\alias{makeBasis}
\title{Construction of quantile-based linear spline basis}
\usage{
  makeBasis(x, type = c("linear", "cubic"), nKnots,
  settings = NULL)
}
\arguments{
  \item{x}{the covariate vector}

  \item{type}{the type of splines to be used,
  \dQuote{linear} splines modelled by a truncated power
  series basis, or \dQuote{cubic} splines modelled by
  O'Sullivan splines.}

  \item{nKnots}{number of quantile-based knots to be used}

  \item{settings}{list of settings. By default, this is
  \code{NULL} and a whole new basis depending on
  \code{type} and \code{nKnots} is built. Otherwise this is
  expected to be \code{attributes(basis)} and the settings
  of the old \code{basis} are used to compute the spline
  basis vectors at new grid points \code{x}.}
}
\value{
  the spline basis in a matrix with \code{length(x)} rows.
  As an attribute, a list with the settings is provided.
}
\description{
  Construction of quantile-based linear spline basis
}
\examples{
x <- seq(from=-5,
         to=5,
         length=1000)

Z <- makeBasis(x=x,
               nKnots=4L,
               type="linear")
str(Z)

plot(x,
     Z \%*\% c(1, -3, 2, 1),
     type="l")
abline(v=attr(Z, "knots"),
       col="gray",
       lty=2)

Z <- makeBasis(x=x,
               nKnots=10L,
               type="cubic")
dimZ <- ncol(Z)

coefs <- rnorm(dimZ)

plot(x,
     Z \%*\% coefs,
     type="l")
abline(v=attr(Z, "knots"),
       col="gray",
       lty=2)
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{programming}

