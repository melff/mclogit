\name{simulate.mclogit}
\alias{simulate.mclogit}
\alias{simulate.mblogit}
\alias{simulate.mmclogit}
\alias{simulate.mmblogit}
\title{
  Simulating responses from baseline-category and conditional logit models.
}
\description{
  The \code{simulate()} methods allow to simulate responses from models
  fitted with \code{mclogit()} and \code{mblogit()}. Currently only
  models \emph{without} random effects are supported for this.
}
\usage{
\method{simulate}{mblogit}(object, nsim = 1, seed = NULL, ...)
\method{simulate}{mclogit}(object, nsim = 1, seed = NULL, ...)

# These methods are currently just 'stubs', causing an error
# message stating that simulation from models with random
# effects are not supported yet
\method{simulate}{mmblogit}(object, nsim = 1, seed = NULL, ...)
\method{simulate}{mmclogit}(object, nsim = 1, seed = NULL, ...)
}
\arguments{
  \item{object}{an object from the relevant class}
  \item{nsim}{a number, specifying the number of simulated responses
    for each observation.}
  \item{seed}{an object specifying if and how the random number
    generator should be initialized (‘seeded’). The interpetation of
  this argument follows the default method, see \code{link[stats]{simulate}}
  }
}
\value{
  The result of the \code{simulate()} method
  for objects created by \code{mclogit()} is a matrix matrices
  with the number of columns given by the \code{nsim=} argument and a
  number of rows corresponding to the observations used to fit the model
  to. The contents of the columns are counts (or zero-one values), with
  group-wise multinomial distribution (within choice sets) just like it
  is assumed for the original response.

  The shape of the result of the \code{simulate()} method
  for objects created by \code{mclogit()} depends on the kind of
  response to which the model was fitted. If this response is a matrix,
  the result of \code{simulate()} is an array with the first two dimensions
  equal to the dimensions of the matrix and the third dimension equal to
  the \code{nsim=} argument. If this response is a factor and the model
  was fitted with \code{from.table=FALSE}, then the result is a matrix
  with as many rows as observations and as many columns as given by
  \code{nsim=}, where the columns contain the numeric codes of the
  factor levels of the original response. If the response is a factor
  and the model was fitted with \code{from.table=TRUE}, the result is
  matrix with counts that have a multinomial distribution within groups
  of rows that correspond to covariate classes.
}
