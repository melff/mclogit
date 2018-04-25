\name{mclogit}
\alias{mclogit}
\alias{predict.mclogit}
\alias{anova.mclogit}
\alias{print.mclogit}
\alias{vcov.mclogit}
\alias{deviance.mclogit}
\alias{logLik.mclogit}
\alias{summary.mclogit}
\alias{print.summary.mclogit}
\alias{fitted.mclogit}
\alias{predict.mclogit}
\alias{residuals.mclogit}
\alias{AIC.mclogit}
\alias{BIC.mclogit}
\alias{anova.mclogit}

\title{Conditional Logit Models and Mixed Conditional Logit Models}
\description{
  \code{mclogit} fits conditional logit models and mixed conditional
  logit models to count data and individual choice data,
  where the choice set may vary across choice occasions.

  Conditional logit models without random effects are fitted by
  Fisher-scoring/IWLS. Models with random effects 
  (mixed conditional logit models) are estimated via maximum likelihood
  with a simple Laplace aproximation (aka PQL).
}
\usage{

mclogit(formula, data=parent.frame(), random=NULL,
        subset, weights = NULL, offset=NULL, na.action = getOption("na.action"),
        model = TRUE, x = FALSE, y = TRUE, contrasts=NULL,
        start=NULL,
        control=mclogit.control(\dots), \dots)
}
\arguments{
  \item{formula}{a model formula: a symbolic description of the
    model to be fitted. The left-hand side contains is expected
    to be a two-column matrix. The first column contains
    the choice counts or choice indicators (alternative is
    chosen=1, is not chosen=0). The second column contains
    unique numbers for each choice set.

    If individual-level data is used, choice sets correspond
    to the individuals, if aggregated data with choice counts are used,
    choice sets may e.g. correspond to covariate classes within clusters.

    The right-hand of the formula contains choice predictors. It should be noted
    that constants are deleted from the formula as are predictors that do not vary
    within choice sets.
    }

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{glm} is called.}

  \item{random}{an optional formula that specifies the random-effects structure or
    NULL.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.  Should be \code{NULL} or a numeric vector.}

  \item{offset}{an optional model offset. Currently only supported
    for models without random effects.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{start}{an optional numerical vector of starting values
    for the conditional logit parameters.
  }

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{x, y}{
    logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.
  }

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{control}{a list of parameters for the fitting process.
    See \code{\link{mclogit.control}} }

  \item{\dots}{
    arguments to be passed to \code{mclogit.control}
    }

}
\value{
  \code{mclogit} returns an object of class "mclogit", which has almost the
  same structure as an object of class "\link[stats]{glm}".
}

\note{
  Covariates that are constant within choice sets are automatically
  dropped from the model formula specified by the \code{formula}
  argument of \code{mclogit}.

  If the model contains random effects, these should
  \itemize{
    \item either vary within choice sets (e.g. the levels of a factor
    that defines the choice sets should not be nested within the levels
    of factor)
    \item or be random coefficients of covariates that vary within
    choice sets. 
  }
  
  In earlier versions of the package (prior to 0.6) it will lead to a
  failure of the model fitting algorithm if these conditions are not
  satisfied. Since version 0.6 of the package, the function
  \code{mclogit} will complain about such model a misspecification
  explicitely.
}

\examples{
data(Transport)

summary(mclogit(
  cbind(resp,suburb)~distance+cost,
  data=Transport
  ))

data(electors)

summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left/class+welfare/class+auth/class,
  random=~1|party.time,
  data=within(electors,party.time<-interaction(party,time))))
}
\keyword{models}
\keyword{regression}