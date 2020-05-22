\name{overdispersion}
\alias{overdispersion.mclogit}
\alias{update.mclogit}
\title{Overdispersion in Multinomial Logit Models}
\description{
  The function 
}
\usage{
overdispersion(object,method, \dots)
\method{overdispersion}{mclogit}(object,method=NULL, \dots)
\method{update}{mclogit}(object, formula., overdispersion, \dots)
}
\arguments{
  \item{object}{an object that inherits class \code{"mclogit"}.
    When passed to \code{overdispersion()}, it 
    should be the result of a call of \code{mclogit()} of
    \code{mblogit()}, \emph{without} random effects.
    The same restriction applies to a call of \code{update()} if an
    \code{overdispersion=} argument is given that is different from
    \code{NULL}.
  }
  \item{method}{a character string, either \code{"Afroz"},
    \code{"Fletcher"}, \code{"Pearson"}, or \code{"Deviance"}, that
    specifies the estimator of the overdispersion; or
    \code{NULL}, in which case the default estimator, \code{"Afroz"} 
    is used. The estimators are discussed in Afroz et al. (2019).
  }
  \item{formula.}{a changes to the model formula,
    see \code{\link[stats]{update.default}} and
    \code{\link[stats]{update.formula}}.}
  \item{overdispersion}{a real number used as dispersion parameter;
    a character vector that specifies the method to compute the overdispersion;
    a logical value -- if \code{TRUE} the default method
    (\code{"Afroz"}) is used, if \code{FALSE}, the dispersion parameter
    is set to 1, that is, no overdispersion.}
  \item{\dots}
}
\references{
  Afroz, Farzana, Matt Parry, and David Fletcher. (2019).
  "Estimating Overdispersion in Sparse Multinomial Data."
  \emph{Biometrics}, early view. \url{https://doi.org/10.1111/biom.13194}.
}
