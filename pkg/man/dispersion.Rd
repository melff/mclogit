\name{dispersion}
\alias{dispersion}
\alias{dispersion.mclogit}
\title{Overdispersion in Multinomial Logit Models}
\description{
  The function \code{dispersion()} extracts the dispersion parameter
  from a multinomial logit model or computes a dispersion parameter
  estimate based on a given method. This dispersion parameter can
  be attached to a model using \code{update()}. It can also given as an
  argument to \code{summary()}.
}
\usage{
dispersion(object,method, \dots)
\method{dispersion}{mclogit}(object,method=NULL, \dots)
}
\arguments{
  \item{object}{an object that inherits class \code{"mclogit"}.
    When passed to \code{dispersion()}, it 
    should be the result of a call of \code{mclogit()} of
    \code{mblogit()}, \emph{without} random effects.
  }
  \item{method}{a character string, either \code{"Afroz"},
    \code{"Fletcher"}, \code{"Pearson"}, or \code{"Deviance"}, that
    specifies the estimator of the dispersion; or
    \code{NULL}, in which case the default estimator, \code{"Afroz"} 
    is used. The estimators are discussed in Afroz et al. (2019).
  }
  \item{\dots}{other arguments, ignored or passed to other methods.}
}
\references{
  Afroz, Farzana, Matt Parry, and David Fletcher. (2020).
  "Estimating Overdispersion in Sparse Multinomial Data."
  \emph{Biometrics} 76(3): 834-842. \doi{10.1111/biom.13194}.
}
\examples{
library(MASS) # For 'housing' data

# Note that with a factor response and frequency weighted data,
# Overdispersion will be overestimated:
house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, weights = Freq,
                         data = housing)

dispersion(house.mblogit,method="Afroz")
dispersion(house.mblogit,method="Deviance")

summary(house.mblogit)

phi.Afroz <- dispersion(house.mblogit,method="Afroz")
summary(house.mblogit, dispersion=phi.Afroz)

summary(update(house.mblogit, dispersion="Afroz"))

# In order to be able to estimate overdispersion accurately,
# data like the above (which usually comes from applying
# 'as.data.frame' to a contingency table) the model has to be
# fitted with the optional argument 'from.table=TRUE':
house.mblogit.corrected <- mblogit(Sat ~ Infl + Type + Cont, weights = Freq,
                                   data = housing, from.table=TRUE,
                                   dispersion="Afroz")
# Now the estimated dispersion parameter is no longer larger than 20,
# but just bit over 1.0.
summary(house.mblogit.corrected)
}
