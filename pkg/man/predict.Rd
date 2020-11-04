\name{predict}
\alias{predict.mblogit}
\alias{predict.mmblogit}
\alias{predict.mclogit}
\alias{predict.mmclogit}

\title{Predicting responses or linear parts of the baseline-category and
  conditional logit models}

\description{
  The \code{predict()} methods allow to obtain within-sample and
  out-of-sample predictions from models
  fitted with \code{mclogit()} and \code{mblogit()}.

  For models with random effecs fitted using the PQL-method, it is
  possible to obtain responses that are conditional on the reconstructed
  random effects.
}

\usage{

\method{predict}{mblogit}(object, newdata=NULL,type=c("link","response"),se.fit=FALSE, \dots)
\method{predict}{mclogit}(object, newdata=NULL,type=c("link","response"),se.fit=FALSE, \dots)
\method{predict}{mmblogit}(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,
                             conditional=TRUE, \dots)
\method{predict}{mmclogit}(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,
                             conditional=TRUE, \dots)
}
\arguments{
 \item{object}{an object in class "mblogit", "mmblogit", "mclogit", or
   "mmclogit"}
 \item{newdata}{an optional data frame with new data}
 \item{type}{a character string specifying the kind of prediction}
 \item{se.fit}{a logical value; whether predictions should be
   accompanied with standard errors}
 \item{conditional}{a logical value; whether predictions should be made
   conditional on the random effects (or whether they are set to zero,
   i.e. their expectation). This argument is consequential only if
   the "mmblogit" or "mmclogit" object was created with \code{method="PQL"}.}
 \item{\dots}{other arguments, ignored.}
}
\value{
  The \code{predict} methods return either a matrix (unless called with
  \code{se.fit=TRUE}) or a list with two matrix-valued elements
  \code{"fit"} and \code{"se.fit"}.
}
\examples{
library(MASS)
(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, 
                          data = housing,
                          weights=Freq))

head(pred.house.mblogit <- predict(house.mblogit))
str(pred.house.mblogit <- predict(house.mblogit,se=TRUE))

head(pred.house.mblogit <- predict(house.mblogit,
                                   type="response"))
str(pred.house.mblogit <- predict(house.mblogit,se=TRUE,
                                  type="response"))
\donttest{ # This takes a bit longer.
data(electors)
(mcre <- mclogit(
    cbind(Freq,interaction(time,class))~econ.left/class+welfare/class+auth/class,
    random=~1|party.time,
    data=within(electors,party.time<-interaction(party,time))))

str(predict(mcre))
str(predict(mcre,type="response"))

str(predict(mcre,se.fit=TRUE))
str(predict(mcre,type="response",se.fit=TRUE))
}
}
