\name{hyp2f1}
\alias{hyp2f1}
\title{Gauss hypergeometric function...}
\usage{hyp2f1(a, b, c, x)}
\description{Gauss hypergeometric function}
\details{This is an R interface to the C function taken from the Cephes library.
Note that all parameters must be numeric and must not be complex. Only
\code{x} may be a vector.}
\value{the vector of function values}
\keyword{math}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{a}{first function parameter}
\item{b}{second function parameter}
\item{c}{third function parameter}
\item{x}{vector of abscissae with absolute value not greater than 1}}
\examples{## examples from the Digital Library of Mathematical Functions
## (http://dlmf.nist.gov/15.3)

x <- seq(from=-0.023,
         to=1,
         length=1000)
plot(x,
     hyp2f1(a=5,
            b=-10,
            c=1,
            x=x),
     type="l")


x <- seq(from=-1,
         to=0.022,
         length=1000)
plot(x,
     hyp2f1(a=5,
            b=10,
            c=1,
            x=x),
     type="l")}
