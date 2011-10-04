\name{linApproxDens}
\alias{linApproxDens}
\title{Test the C++ object class LinApproxDens}
\usage{
  linApproxDens(args, logDens, grid, nSamples)
}
\arguments{
  \item{args}{the x values}

  \item{logDens}{the unnormalised log density values}

  \item{grid}{the grid on which to evaluate the linear
  spline approximation}

  \item{nSamples}{number of samples to produce}
}
\value{
  list with approximate density values and samples
}
\description{
  Test the C++ object class LinApproxDens
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{internal}

