\name{getRhos}
\alias{getRhos}
\title{Compute the penalty parameters corresponding to degrees of freedom...}
\usage{getRhos(lambdas, degrees)}
\description{Compute the penalty parameters corresponding to degrees of freedom}
\value{the vector of rho parameters corresponding to the \code{degrees}}
\note{All \code{degrees} must be positive, and smaller than the dimension
of the spline basis (which is the number of columns of \code{Z}).}
\keyword{programming}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{lambdas}{the eigenvalues of Z^T * W * Z, where Z is the
spline basis matrix}
\item{degrees}{the vector of degrees of freedom}}
\examples{x <- seq(from=-5,
         to=5,
         length=100)
Z <- makeBasis(x=x,
               nKnots=4L)
lambdas <- svd(Z, nu=0L, nv=0L)$d^2
getRhos(lambdas=lambdas,
        degrees=c(1, 2, 3, 3.9))}
