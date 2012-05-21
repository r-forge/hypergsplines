\name{rshrinkage}
\alias{rshrinkage}
\title{Sample from the model-specific posterior of the shrinkage factor}
\usage{
  rshrinkage(n, marginal, modelData)
}
\arguments{
  \item{n}{number of samples}

  \item{marginal}{the marginal model, which is the result
  from \code{\link{calculateModel}}}

  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{modelData}}}
}
\value{
  \code{n} posterior shrinkage factor samples
}
\description{
  Sample from the model-specific posterior of the shrinkage
  factor t = g / (1 + g), using either inverse sampling
  (for the hyper-g prior) or numerical approximations (for
  the hyper-g/n prior).
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{distribution}

