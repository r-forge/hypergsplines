\name{getMcmc}
\alias{getMcmc}
\title{Collect the MCMC options...}
\usage{getMcmc(samples=10000, burnin=10000, step=2, nIwlsIterations=1)}
\description{Collect the MCMC options}
\keyword{utilities}
\author{Daniel Sabanes Bove \email{daniel.sabanesbove@ifspm.uzh.ch}}
\arguments{\item{samples}{number of resulting samples (default: \code{10,000})}
\item{burnin}{number of burn-in iterations which are not saved (default:
\code{10,000})}
\item{step}{only every step-th iteration is saved after the burn-in
(default: \code{2})}
\item{nIwlsIterations}{maximum number of IWLS iterations in the proposal
step (default: \code{1}). Set that to a higher value to achieve higher
acceptance rates.}}
