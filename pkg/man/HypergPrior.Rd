\name{HypergPrior}
\alias{HypergPrior}
\title{Constructor for the hyper-g prior class}
\usage{
  HypergPrior(a = 4)
}
\arguments{
  \item{a}{the hyperparameter which must be larger than 3,
  and should not be larger than 4 in order not to favour
  too much shrinkage a priori (default: 4)}
}
\value{
  a new \code{\linkS4class{HypergPrior}} object
}
\description{
  Constructor for the hyper-g prior class
}
\keyword{classes}

