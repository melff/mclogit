\name{mclogit.fit}
\alias{mclogit.fit}
\alias{mmclogit.fitPQL}
\title{
  Internal functions used for model fit.
}
\description{
  These functions are exported and documented for use by other packages. They are not
  intended for end users.
}
\usage{
mclogit.fit(y, s, w, X, start = NULL, offset = NULL, control = mclogit.control())

mmclogit.fitPQL(y, s, w, X, Z, G, groups, start,
                  offset = NULL, control = mclogit.control())
}
\arguments{
  \item{y}{a response vector. Should be binary.}
  \item{s}{a vector identifying individuals or covariate strata}
  \item{w}{a vector with observation weights.}
  \item{X}{a model matrix; required.}
  \item{Z}{the random effects design matrix.}
  \item{G}{a list of design matrices for the (co-)variance parameters.}
  \item{groups}{a list of grouping factors.}
  \item{start}{an optional numerical vector of starting values
    for the coefficients.
  }

  \item{offset}{an optional model offset. Currently only supported
    for models without random effects.}
  \item{control}{a list of parameters for the fitting process.
    See \code{\link{mclogit.control}} }

}
\value{
  A list with components describing the fitted model.
}
