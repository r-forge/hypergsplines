\name{checkModelConfigs}
\alias{checkModelConfigs}
\title{Check coherency of model configurations with data}
\usage{
  checkModelConfigs(modelConfigs, modelData)
}
\arguments{
  \item{modelConfigs}{the matrix of model configurations}

  \item{modelData}{the data necessary for model estimation,
  which is the result from \code{\link{modelData}} or
  \code{\link{glmModelData}}}
}
\value{
  returns \code{TRUE} if model configurations are coherent
  with \code{modelData}, otherwise throws error messages
}
\description{
  Check coherency of model configurations with data
}
\author{
  Daniel Sabanes Bove
}
\keyword{internal}

