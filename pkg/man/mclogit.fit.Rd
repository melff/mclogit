\name{mclogit.fit}
\alias{mclogit.fit}
\alias{mmclogit.fitPQLMQL}
\title{
  Internal functions used for model fit.
}
\description{
  These functions are exported and documented for use by other packages. They are not
  intended for end users.
}
\usage{
mclogit.fit(y, s, w, X,
            overdispersion=FALSE,
            start = NULL, offset = NULL,
            control = mclogit.control())

mmclogit.fitPQLMQL(y, s, w, X, Z, groups, start,
                   offset = NULL, method=c("PQL","MQL"),
                   control = mmclogit.control())
}
\arguments{
  \item{y}{a response vector. Should be binary.}
  \item{s}{a vector identifying individuals or covariate strata}
  \item{w}{a vector with observation weights.}
  \item{X}{a model matrix; required.}
  \item{overdispersion}{a logical value or a character string; whether and how
     an overdispersion parameter should be estimated. For details see \code{\link{overdispersion}}.}
  \item{Z}{the random effects design matrix.}
  \item{groups}{a list of grouping factors.}
  \item{start}{an optional numerical vector of starting values
    for the coefficients.
  }
  \item{offset}{an optional model offset. Currently only supported
    for models without random effects.}
  \item{method}{a character string, either "PQL" or "MQL", specifies
    the type of the quasilikelihood approximation.}
  \item{control}{a list of parameters for the fitting process.
    See \code{\link{mclogit.control}} }

}
\value{
  A list with components describing the fitted model.
}
