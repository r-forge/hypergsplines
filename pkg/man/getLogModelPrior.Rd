\name{getLogModelPrior}
\alias{getLogModelPrior}
\title{Get the log prior probability for a model}
\usage{
  getLogModelPrior(config, type = c("flat", "exponential",
  "independent", "dependent", "dep.linear"), modelData)
}
\arguments{
  \item{config}{the model configuration vector, or a matrix
  of model configurations in rows. There may be additional
  elements to the right, e.g. the log marginal likelihood
  column/value, if \code{modelData} is provided. These are
  then ignored internally.}

  \item{type}{either \dQuote{flat} (default),
  \dQuote{exponential}, \dQuote{independent} or
  \dQuote{dependent}. See the details.}

  \item{modelData}{the model setup from
  \code{\link{modelData}}}
}
\value{
  the (unnormalized) log prior probability
}
\description{
  The \dQuote{exponential} prior assumes
  \eqn{f(\boldsymbol{d}) \propto \exp(-2 \sum_{j=1}^{p}
  d_{j})}{f(d) ~ exp(-2 * sum(d_j))} which penalizes models
  with many degrees of freedom regardless of linear
  (\eqn{d_j = 1}) or smooth (\eqn{d_j > 1}) inclusion of
  covariates.
}
\details{
  The \dQuote{independent} prior (as the exponential prior)
  that the different covariates are independent,
  \eqn{f(\boldsymbol{d}) = \prod_{j=1}^{p} f(d_{j})}{f(d) =
  prod(f(d_j))}, where \eqn{f(d_{j})}{f(d_j)} is built
  hierarchically: The prior probability for exclusion is
  0.5, thus f(0)=0.5. Given inclusion, the other half is
  split up uniformly between all possible degrees of
  freedom larger than 0.

  The \dQuote{dependent} prior imposes a uniform
  distribution for the number of included variables. The
  \dQuote{position} of the included variables, i.e. which
  of the variables are included, is uniformly distributed
  on all possible configurations. The degrees of freedoms
  are then conditionally independent and uniformly
  distributed on the respective possible values.

  The \dQuote{dep.linear} prior is similar, but imposes a
  uniform distribution for the number of smoothly included
  variables conditional on inclusion. The \dQuote{position}
  of the (smoothly) included variables, i.e. which of the
  (linearly included) variables are included (as splines),
  is uniformly distributed on the possible configurations
  at both steps. The spline variance parameters are then
  conditionally independent and uniformly distributed on
  the respective possible values. This leads to a fixed
  prior probability of linear inclusion of 1/4.
}
\examples{
## artificial configuration:
config <- c(a=0L, b=1L, c=3L, d=0L)

## get log prior probs:
getLogModelPrior(config=config,
                 type="flat")

getLogModelPrior(config=config,
                 type="exponential")
}
\author{
  Daniel Sabanes Bove
  \email{daniel.sabanesbove@ifspm.uzh.ch}
}
\keyword{regression}

