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
            dispersion=FALSE,
            start = NULL, offset = NULL,
            control = mclogit.control(),
            groups=NULL)

mmclogit.fitPQLMQL(y, s, w, X, Z, d, 
                   start = NULL,
                   start.Phi = NULL,
                   start.b = NULL,
                   offset = NULL, method=c("PQL","MQL"),
                   estimator = c("ML","REML"),
                   control = mmclogit.control())
}
\arguments{
  \item{y}{a response vector. Should be binary.}
  \item{s}{a vector identifying individuals or covariate strata}
  \item{w}{a vector with observation weights.}
  \item{X}{a model matrix; required.}
  \item{dispersion}{a logical value or a character string; whether and how
     a dispersion parameter should be estimated. For details see \code{\link{dispersion}}.}
  \item{groups}{a vector that identifies the groups relevant 
     for the estimation of the dispersion parameter.}
  \item{Z}{the random effects design matrix.}
  \item{d}{dimension of random effects. Typically $d=1$ for random intercepts
  only, $d>1$ for models with random intercepts.}
  \item{start}{an optional numerical vector of starting values
    for the coefficients.
  }
  \item{offset}{an optional model offset. Currently only supported
    for models without random effects.}
  \item{start.Phi}{an optional matrix of strarting values for the
    (co-)variance parameters.}
  \item{start.b}{an optional list of vectors with starting values
    for the random effects.}
  \item{method}{a character string, either "PQL" or "MQL", specifies
    the type of the quasilikelihood approximation.}
  \item{estimator}{a character string; either "ML" or "REML",
    specifies which estimator is to be used/approximated.}
  \item{control}{a list of parameters for the fitting process.
    See \code{\link{mclogit.control}} }

}
\value{
  A list with components describing the fitted model.
}
