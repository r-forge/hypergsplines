\name{getFamily}
\alias{getFamily}
\title{Helper function for glmModelData: Extracts an S3 family object}
\usage{
  getFamily(family, phi)
}
\arguments{
  \item{family}{the family argument passed to
  \code{\link{glmModelData}}}

  \item{phi}{the dispersion argument passed to
  \code{\link{glmModelData}}}
}
\value{
  The returned family object also includes a custom
  \sQuote{init} function, which takes the response vector
  (or matrix) (\sQuote{y}) and the corresponding weight
  vector (\sQuote{weights}), processes them to response
  vector \sQuote{y} and possibly altered weights
  \sQuote{weights}, and includes starting values
  \sQuote{mustart} for the IWLS algorithm. For example,
  here the binomial special case of two-column-response
  matrix is treated exactly in the same way as in
  \code{\link{glm}}.
}
\description{
  Extracts an S3 family object, which (at least) the usual
  elements "family", "link", "linkfun", "linkinv",
  "variance", "mu.eta", "dev.resids", plus the additional
  elements: "phi" which includes just the dispersion
  parameter to be used, and "simulate" which can generate
  random variates for a given linear predictor and weight
  vector (this function of course also uses the "phi"
  value)
}
\keyword{internal}

