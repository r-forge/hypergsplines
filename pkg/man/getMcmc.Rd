\name{getMcmc}
\alias{getMcmc}
\title{Collect the MCMC options}
\usage{
  getMcmc(samples = 10000L, burnin = 10000L, step = 2L,
    nIwlsIterations = 1L)
}
\arguments{
  \item{samples}{number of resulting samples (default:
  \code{10,000})}

  \item{burnin}{number of burn-in iterations which are not
  saved (default: \code{10,000})}

  \item{step}{only every step-th iteration is saved after
  the burn-in (default: \code{2})}

  \item{nIwlsIterations}{maximum number of IWLS iterations
  in the proposal step (default: \code{1}). Set that to a
  higher value to achieve higher acceptance rates.}
}
\description{
  Collect the MCMC options
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{utilities}

