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
dispersion(object, method, \dots)
\method{dispersion}{mclogit}(object, method=NULL, \dots)
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
house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, 
                         weights = Freq,
                         data = housing)
dispersion(house.mblogit, method = "Afroz")
dispersion(house.mblogit, method = "Deviance")

# In order to be able to estimate overdispersion accurately,
# data like the above (which usually comes from applying
# 'as.data.frame' to a contingency table) the model has to be
# fitted with the optional argument 'aggregate=TRUE' or 
# by requesting the dispersion in advance.
house.mblogit.agg <- mblogit(Sat ~ Infl + Type + Cont, 
                             weights = Freq,
                             data = housing, 
                             aggregate = TRUE)
# Now the estimated dispersion parameter is no longer larger than 20,
# but just bit over 1.0.
dispersion(house.mblogit.agg, method = "Afroz")
dispersion(house.mblogit.agg, method = "Deviance")

# It is possible to obtain the dispersion after estimating the coefficients:
phi.Afroz <- dispersion(house.mblogit.agg, method = "Afroz")
summary(house.mblogit.agg, dispersion = phi.Afroz)

summary(update(house.mblogit.agg, dispersion = "Afroz"))

# If an estimate of the (over-)dispersion is requested, 'aggregate' is set to
# TRUE by default:
house.mblogit.odsp <- mblogit(Sat ~ Infl + Type + Cont, 
                              weights = Freq,
                              data = housing, 
                              dispersion = "Afroz")
summary(house.mblogit.odsp)
dispersion(house.mblogit.odsp, method = "Deviance")

# Note that aggregation (either implicitly or explicitly required) affects
# the reported deviance in surprising ways:
house.mblogit.o.00 <- mblogit(Sat ~ Infl, 
                              weights = Freq,
                              data = housing, 
                              dispersion = TRUE)
deviance(house.mblogit.o.00)
dispersion(house.mblogit.o.00)
# The deviance is (almost) zero, because aggregation leads to a two-way
# table and a single-predictor model is already saturated.

# In order to make models comparable, one will need to set the groups:
house.mblogit.o.0 <- mblogit(Sat ~ Infl, 
                             weights = Freq,
                             data = housing, 
                             groups = ~ Infl + Type + Cont,
                             dispersion = TRUE)
deviance(house.mblogit.o.0)
dispersion(house.mblogit.o.0)

anova(house.mblogit.o.0,
      house.mblogit.odsp)

# These complications with the deviances do not arise if no aggregation is 
# requested:
house.mblogit.0 <- mblogit(Sat ~ Infl, 
                           weights = Freq,
                           data = housing)
anova(house.mblogit.0,
      house.mblogit)


# Using frequences on the left-hand side is perhaps the safest option:
housing.wide <- memisc::Aggregate(table(Sat) ~ Infl + Type + Cont,
                                  data = housing) # Note that 'Aggegate' uses
                                                # variable 'Freq' for weighting.
house.mblogit.wide <- mblogit(cbind(Low,Medium,High) ~ Infl + Type + Cont, 
                              data = housing.wide)
summary(house.mblogit.wide)
dispersion(house.mblogit.wide, method = "Afroz")

house.mblogit.wide.0 <- mblogit(cbind(Low,Medium,High) ~ Infl, 
                                data = housing.wide)
summary(house.mblogit.wide.0)
dispersion(house.mblogit.wide.0, method="Afroz")

anova(house.mblogit.wide.0,
      house.mblogit.wide)
}
