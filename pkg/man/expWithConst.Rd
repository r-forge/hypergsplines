\name{expWithConst}
\alias{expWithConst}
\title{Exponentiate a vector of values plus a constant}
\usage{
  expWithConst(values)
}
\arguments{
  \item{values}{argument values, of which at least one must
  be finite}
}
\value{
  the exp(values + constant) results, see details
}
\description{
  This functions tries to choose the constant such that the
  smallest result is larger than zero and the largest
  result is smaller than infinity, in the numerical
  representation. If this is not possible, results
  numerically identical to zero are tolerated. Only values
  which are finite are taken into account for this
  considerations, because infinite values will be either
  zero or infinity anyway, regardless of the constant.
}
\author{
  Daniel Sabanes Bove
}
\keyword{internal}

