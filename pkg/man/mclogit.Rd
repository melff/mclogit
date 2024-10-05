\name{mclogit}
\alias{mclogit}
\alias{anova.mclogit}
\alias{print.mclogit}
\alias{vcov.mclogit}
\alias{deviance.mclogit}
\alias{logLik.mclogit}
\alias{summary.mclogit}
\alias{print.summary.mclogit}
\alias{fitted.mclogit}
\alias{residuals.mclogit}
\alias{weights.mclogit}
\alias{AIC.mclogit}
\alias{BIC.mclogit}
\alias{update.mclogit}
\alias{anova.mclogit}
\alias{summary.mmclogit}
\alias{print.summary.mmclogit}
\alias{ranef.mmclogit}

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
        method = NULL, estimator=c("ML","REML"),
        dispersion = FALSE,
        start=NULL,
        groups = NULL,
        control=if(length(random))
                    mmclogit.control(\dots)
                else mclogit.control(\dots), \dots)

\method{update}{mclogit}(object, formula., dispersion, \dots)

\method{summary}{mclogit}(object, dispersion = NULL, correlation = FALSE,
        symbolic.cor = FALSE,  \dots)    
}
\arguments{
  \item{formula}{a model formula: a symbolic description of the
    model to be fitted. The left-hand side should result in
    a two-column matrix. The first column contains
    the choice counts or choice indicators (alternative is
    chosen=1, is not chosen=0). The second column contains
    unique numbers for each choice set.

    The left-hand side can either take the form \code{cbind(choice,set)}
    or (from version 0.9.1) \code{choice|set}

    If individual-level data is used, choice sets correspond
    to individuals, if aggregated data with choice counts are used,
    choice sets usually correspond to covariate classes.

    The right-hand of the formula contains choice predictors. It should be noted
    that constants are deleted from the formula as are predictors that do not vary
    within choice sets.
    }

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{glm} is called.}

  \item{random}{an optional formula or list of formulas that specify the
    random-effects structure or NULL.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.  Should be \code{NULL} or a numeric vector.}

  \item{offset}{an optional model offset.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{start}{an optional numerical vector of starting values for the
    conditional logit parameters.  If the model has random effects, the
    vector should have a "VarCov" attribute wtih starting values for the
    random effects (co-)variances. If the random effects model is
    estimated with the "PQL" method, the starting values matrix should
    also have a "random.effects" attribute, which should have the same
    structure as the "random.effects" component of an object returned by
    \code{mblogit()}.  }

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{x, y}{
    logical values indicating whether the response vector and model
    matrix used in the fitting process should be returned as components
    of the returned value.
  }

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}

  \item{method}{\code{NULL} or a character string, either "PQL" or "MQL", specifies
      the type of the quasilikelihood approximation to be used if
      a random-effects model is to be estimated.}

  \item{estimator}{a character string; either "ML" or "REML",
    specifies which estimator is to be used/approximated.}
    
  \item{dispersion}{a real number used as dispersion parameter;
    a character vector that specifies the method to compute the dispersion;
    a logical value -- if \code{TRUE} the default method
    (\code{"Afroz"}) is used, if \code{FALSE}, the dispersion parameter
    is set to 1, that is, no dispersion. For details see \code{\link{dispersion}}.}

  \item{groups}{an optional formula that specifies groups of observations
    relevant for the estimation of overdispersion. Covariates should be
    constant within groups, otherwise a warning is generated
    since the overdispersion estimate may
    be imprecise.
  }
  
  \item{control}{a list of parameters for the fitting process.
    See \code{\link{mclogit.control}} }

  \item{\dots}{
    arguments to be passed to \code{mclogit.control} or \code{mmclogit.control}
    }

  \item{object}{an object that inherits class \code{"mclogit"}.
    When passed to \code{dispersion()}, it 
    should be the result of a call of \code{mclogit()} of
    \code{mblogit()}, \emph{without} random effects.
  }
  \item{formula.}{a changes to the model formula,
    see \code{\link[stats:update]{update.default}} and
    \code{\link[stats]{update.formula}}.}
  \item{correlation}{logical; see \code{\link[stats]{summary.lm}}.}
  \item{symbolic.cor}{logical; see \code{\link[stats]{summary.lm}}.}
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

  From version 0.9.7 it is possible to choose the optimization
  technique used for the inner iterations of the PQL/MQL: either
  \code{\link[stats]{nlminb}} (the default), \code{\link[stats]{nlm}},
  or any of the algorithms (other than "Brent" supported by
  \code{\link[stats]{optim}}).  To choose the optimizer, use the
  appropriate argument for \code{\link{mmclogit.control}} .
}

\references{
    Agresti, Alan (2002).
    \emph{Categorical Data Analysis.} 2nd ed, Hoboken, NJ: Wiley.
    \doi{10.1002/0471249688}

    Breslow, N.E. and D.G. Clayton (1993).
    "Approximate Inference in Generalized Linear Mixed Models".
    \emph{Journal of the American Statistical Association} 88 (421): 9-25.
    \doi{10.1080/01621459.1993.10594284}
    
    Elff, Martin (2009).
    "Social Divisions, Party Positions, and Electoral Behaviour".
    \emph{Electoral Studies} 28(2): 297-308.
    \doi{10.1016/j.electstud.2009.02.002}
    
    McFadden, D. (1973).
    "Conditionial Logit Analysis of Qualitative Choice Behavior".
    Pp. 105-135 in P. Zarembka (ed.).
    \emph{Frontiers in Econometrics}.
    New York: Wiley.
    \url{https://eml.berkeley.edu/reprints/mcfadden/zarembka.pdf}
}

\examples{
data(Transport)

summary(mclogit(
  cbind(resp,suburb)~distance+cost,
  data=Transport
  ))
# New syntactic sugar:
summary(mclogit(
  resp|suburb~distance+cost,
  data=Transport
  ))


\dontrun{ # This takes a bit longer.
data(electors)

electors <- within(electors,{
    party.time <-interaction(party,time)
    time.class <- interaction(time,class)
})

# Time points nested within parties
summary(mclogit(
  Freq|time.class~econ.left/class+welfare/class+auth/class,
  random=~1|party/time,
  data=electors))

# Party-level random intercepts and random slopes varying over time points
summary(mclogit(
  Freq|time.class~econ.left/class+welfare/class+auth/class,
  random=list(~1|party,~econ.left+0|time),
  data=electors))
}
}
\keyword{models}
\keyword{regression}

\seealso{
Conditional logit models are also supported by \pkg{gmnl}, \pkg{mlogit}, and \pkg{survival}.
\pkg{survival} supports conditional logit models for binary panel data and case-control studies.     
\pkg{mlogit} and \pkg{gmnl} treat conditional logit models from an econometric perspective.
Unlike the present package, they focus on the random utility interpretation of discrete choice models
and support generalisations of conditional logit models, such as nested logit models, that are intended
to overcome the IIA (indipendence from irrelevant alterantives) assumption. Mixed multinomial models are
also supported and estimated using simulation-based techniques. Unlike the present package,
mixed or random-effects extensions are mainly intended to fit repeated choices of the same individuals and not
aggregated choices of many individuals facing identical alternatives.
}
