\name{HypergnPrior}
\alias{HypergnPrior}
\title{Constructor for the "HypergnPrior" class}
\usage{
  HypergnPrior(a = 4, n)
}
\arguments{
  \item{a}{the hyperparameter which must be larger than 3,
  and should not be larger than 4 in order not to favour
  too much shrinkage a priori (default: 4)}

  \item{n}{the sample size (positive integer)}
}
\value{
  a new \code{\linkS4class{HypergnPrior}} object
}
\description{
  Constructor for the "HypergnPrior" class
}
\keyword{classes}

