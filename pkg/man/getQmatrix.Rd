\name{getQmatrix}
\alias{getQmatrix}
\title{Compute the Q matrix necessary for the offsets for survival regression}
\usage{
  getQmatrix(times)
}
\arguments{
  \item{times}{vector of \eqn{n} ordered survival times,
  not including the pseudo-zero!}
}
\value{
  the \eqn{(n+1) \times (n+1)} matrix Q
}
\description{
  This function uses Algorithm 1 from our Reader Reaction
  article to compute the Q matrix.
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{internal}

