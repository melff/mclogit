\name{dispersion}
\alias{dispersion}
\alias{dispersion.mclogit}
\alias{update.mclogit}
\title{Overdispersion in Multinomial Logit Models}
\description{
  The function 
}
\usage{
dispersion(object,method, \dots)
\method{dispersion}{mclogit}(object,method=NULL, \dots)
\method{update}{mclogit}(object, formula., dispersion, \dots)
}
\arguments{
  \item{object}{an object that inherits class \code{"mclogit"}.
    When passed to \code{dispersion()}, it 
    should be the result of a call of \code{mclogit()} of
    \code{mblogit()}, \emph{without} random effects.
    The same restriction applies to a call of \code{update()} if an
    \code{dispersion=} argument is given that is different from
    \code{NULL}.
  }
  \item{method}{a character string, either \code{"Afroz"},
    \code{"Fletcher"}, \code{"Pearson"}, or \code{"Deviance"}, that
    specifies the estimator of the dispersion; or
    \code{NULL}, in which case the default estimator, \code{"Afroz"} 
    is used. The estimators are discussed in Afroz et al. (2019).
  }
  \item{formula.}{a changes to the model formula,
    see \code{\link[stats]{update.default}} and
    \code{\link[stats]{update.formula}}.}
  \item{dispersion}{a real number used as dispersion parameter;
    a character vector that specifies the method to compute the dispersion;
    a logical value -- if \code{TRUE} the default method
    (\code{"Afroz"}) is used, if \code{FALSE}, the dispersion parameter
    is set to 1, that is, no dispersion.}
  \item{\dots}{other arguments, ignored or passed to other methods.}
}
\references{
  Afroz, Farzana, Matt Parry, and David Fletcher. (2019).
  "Estimating Overdispersion in Sparse Multinomial Data."
  \emph{Biometrics}, early view. \url{https://doi.org/10.1111/biom.13194}.
}
\examples{
library(MASS) # For 'housing' data

house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, weights = Freq,
                         data = housing)

dispersion(house.mblogit,method="Afroz")
dispersion(house.mblogit,method="Deviance")

summary(house.mblogit)

phi.Afroz <- dispersion(house.mblogit,method="Afroz")
summary(house.mblogit, dispersion=phi.Afroz)

summary(update(house.mblogit, dispersion="Afroz"))

}
