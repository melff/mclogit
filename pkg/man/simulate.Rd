\name{simulate.mclogit}
\alias{simulate.mclogit}
\alias{simulate.mblogit}
\alias{simulate.mmclogit}
\alias{simulate.mmblogit}
\title{
  Simulating responses from baseline-category and conditional logit models
}
\description{
  The \code{simulate()} methods allow to simulate responses from models
  fitted with \code{mclogit()} and \code{mblogit()}. Currently only
  models \emph{without} random effects are supported for this.
}
\usage{
\method{simulate}{mblogit}(object, nsim = 1, seed = NULL, \dots)
\method{simulate}{mclogit}(object, nsim = 1, seed = NULL, \dots)

# These methods are currently just 'stubs', causing an error
# message stating that simulation from models with random
# effects are not supported yet
\method{simulate}{mmblogit}(object, nsim = 1, seed = NULL, \dots)
\method{simulate}{mmclogit}(object, nsim = 1, seed = NULL, \dots)
}
\arguments{
  \item{object}{an object from the relevant class}
  \item{nsim}{a number, specifying the number of simulated responses
    for each observation.}
  \item{seed}{an object specifying if and how the random number
    generator should be initialized ('seeded'). The interpetation of
  this argument follows the default method, see \code{link[stats]{simulate}}
  }
  \item{\dots}{other arguments, ignored.}
}
\value{
  The result of the \code{\link[stats]{simulate}} method for objects
  created by \code{\link{mclogit}} is a data frame with one variable for
  each requested simulation run (their number is given by the
  \code{nsim=} argument). The contents of the columns are counts (or
  zero-one values), with group-wise multinomial distribution (within
  choice sets) just like it is assumed for the original response.

  The shape of the result of the \code{\link[stats]{simulate}} method
  for objects created by \code{\link{mblogit}} is also a data frame.
  The variables within the data frame have a mode or shape that
  corresponds to the response to which the model was fitted. If the
  response is a matrix of counts, then the variables in the data frame
  are also matrices of counts.  If the response is a factor and
  \code{\link{mblogit}} was called with an argument
  \code{from.table=FALSE}, the variables in the data frame are factors
  with the same factor levels as the response to which the model was
  fitted. If instead the function was called with
  \code{from.table=TRUE}, the variables in the data frame are counts,
  which represent frequency weights that would result from applying
  \code{\link[base]{as.data.frame}} to a contingency table of simulated
  frequency counts.
}
\examples{
library(MASS)
(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, 
                          data = housing,
                          weights=Freq,
                          from.table=TRUE))
sm <- simulate(house.mblogit,nsim=7)

housing.long <- housing[rep(seq.int(nrow(housing)),housing$Freq),]
(housel.mblogit <- mblogit(Sat ~ Infl + Type + Cont,
                           data=housing.long))
sml <- simulate(housel.mblogit,nsim=7)

housing.table <- xtabs(Freq~.,data=housing)
housing.mat <- memisc::to.data.frame(housing.table)
head(housing.mat)

(housem.mblogit <- mblogit(cbind(Low,Medium,High) ~
                               Infl + Type + Cont,
                           data=housing.mat))
smm <- simulate(housem.mblogit,nsim=7)

str(sm)
str(sml)
str(smm)

head(smm[[1]])
}
