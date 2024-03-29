\name{empiricalHpd}
\alias{empiricalHpd}
\title{Construct an empirical HPD interval from samples}
\usage{
  empiricalHpd(theta, level)
}
\arguments{
  \item{theta}{the vector of samples}

  \item{level}{the credible level}
}
\value{
  A vector with the estimated lower and upper bounds of the
  HPD interval.
}
\description{
  Construct an empirical highest posterior density (HPD)
  interval from samples which have been drawn from the
  distribution of a quantity of interest.
}
\seealso{
  \code{\link{scrHpd}}
}
\keyword{htest}

