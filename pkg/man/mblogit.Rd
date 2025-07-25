% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mblogit.R
\name{mblogit}
\alias{mblogit}
\alias{print.mblogit}
\alias{summary.mblogit}
\alias{print.summary.mblogit}
\alias{fitted.mblogit}
\alias{weights.mblogit}
\alias{print.mmblogit}
\alias{summary.mmblogit}
\alias{print.summary.mmblogit}
\title{Baseline-Category Logit Models for Categorical and Multinomial Responses}
\usage{
mblogit(
  formula,
  data = parent.frame(),
  random = NULL,
  catCov = c("free", "diagonal", "single"),
  subset,
  weights = NULL,
  offset = NULL,
  na.action = getOption("na.action"),
  model = TRUE,
  x = FALSE,
  y = TRUE,
  contrasts = NULL,
  method = NULL,
  estimator = c("ML", "REML"),
  dispersion = FALSE,
  start = NULL,
  aggregate = !isFALSE(dispersion),
  groups = NULL,
  from.table = FALSE,
  control = if (length(random)) mmclogit.control(...) else mclogit.control(...),
  ...
)
}
\arguments{
\item{formula}{the model formula. The response must be a factor or a matrix
of counts.}

\item{data}{an optional data frame, list or environment (or object coercible
by \code{\link{as.data.frame}} to a data frame) containing the variables
in the model.  If not found in \code{data}, the variables are taken from
\code{environment(formula)}, typically the environment from which
\code{glm} is called.}

\item{random}{an optional formula or list of formulas that specify the
random-effects structure or NULL.}

\item{catCov}{a character string that specifies optional restrictions
on the covariances of random effects between the logit equations.
"free" means no restrictions, "diagonal" means that random effects
pertinent to different categories are uncorrelated, while "single" means
that the random effect variances pertinent to all categories are identical.}

\item{subset}{an optional vector specifying a subset of observations to be
used in the fitting process.}

\item{weights}{an optional vector of weights to be used in the fitting
process.  Should be \code{NULL} or a numeric vector.}

\item{offset}{an optional model offset. If not NULL, must be a matrix
if as many columns as the response has categories or one less.}

\item{na.action}{a function which indicates what should happen when the data
contain \code{NA}s.  The default is set by the \code{na.action} setting
of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.
The \sQuote{factory-fresh} default is \code{\link{na.omit}}.  Another
possible value is \code{NULL}, no action.  Value \code{\link{na.exclude}}
can be useful.}

\item{model}{a logical value indicating whether \emph{model frame} should be
included as a component of the returned value.}

\item{x, y}{logical values indicating whether the response vector and model
matrix used in the fitting process should be returned as components of
the returned value.}

\item{contrasts}{an optional list. See the \code{contrasts.arg} of
\code{model.matrix.default}.}

\item{method}{\code{NULL} or a character string, either "PQL" or "MQL",
specifies the type of the quasilikelihood approximation to be used if a
random-effects model is to be estimated.}

\item{estimator}{a character string; either "ML" or "REML", specifies which
estimator is to be used/approximated.}

\item{dispersion}{a logical value or a character string; whether and how a
dispersion parameter should be estimated. For details see
\code{\link{dispersion}}.}

\item{start}{an optional matrix of starting values (with as many rows
as logit equations). If the model has random effects, the matrix
should have a "VarCov" attribute wtih starting values for
the random effects (co-)variances. If the random effects model
is estimated with the "PQL" method, the starting values matrix
should also have a "random.effects" attribute, which should have
the same structure as the "random.effects" component of an object
returned by \code{mblogit()}.}

\item{aggregate}{a logical value; whether to aggregate responses by
covariate classes and groups before estimating the model
if the response variable is a factor. 

This will not affect the estimates, but the dispersion and the
residual degrees of freedom. If \code{aggregate=TRUE}, the 
dispersion will be relative to a saturated model; it will be much 
smaller than with \code{aggregate=TRUE}. In particular, with only
a single covariate and no grouping, the deviance will be close to
zero. If \code{dispersion} is not \code{FALSE}, then the
default value of \code{aggregate} will be \code{TRUE}. For details see
\code{\link{dispersion}}.

This argument has consequences only if the response in \code{formula}
is a factor.}

\item{groups}{an optional formula that specifies groups of observations
relevant for the estimation of overdispersion. For details see
\code{\link{dispersion}}.}

\item{from.table}{a logical value; should be FALSE. This argument
only exists for the sake of compatibility and will be removed
in the next relase.}

\item{control}{a list of parameters for the fitting process.  See
\code{\link{mclogit.control}}}

\item{\dots}{arguments to be passed to \code{mclogit.control} or
\code{mmclogit.control}}
}
\value{
\code{mblogit} returns an object of class "mblogit", which has almost
    the same structure as an object of class "\link[stats]{glm}". The
    difference are the components \code{coefficients}, \code{residuals},
    \code{fitted.values}, \code{linear.predictors}, and \code{y}, which are
    matrices with number of columns equal to the number of response
    categories minus one.
}
\description{
The function \code{mblogit} fits baseline-category logit models for categorical
and multinomial count responses with fixed alternatives.
}
\details{
The function \code{mblogit} internally rearranges the data into a
    'long' format and uses \code{\link{mclogit.fit}} to compute
    estimates. Nevertheless, the 'user data' are unaffected.
}
\examples{
library(MASS) # For 'housing' data
library(nnet)
library(memisc)

(house.mult<- multinom(Sat ~ Infl + Type + Cont, weights = Freq,
                       data = housing))


(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, weights = Freq,
                         data = housing))

summary(house.mult)

summary(house.mblogit)

mtable(house.mblogit)
}
\references{
Agresti, Alan. 2002.
   \emph{Categorical Data Analysis.} 2nd ed, Hoboken, NJ: Wiley.
   \doi{10.1002/0471249688}

   Breslow, N.E. and D.G. Clayton. 1993.
   "Approximate Inference in Generalized Linear Mixed Models".
   \emph{Journal of the American Statistical Association} 88 (421): 9-25.
   \doi{10.1080/01621459.1993.10594284}
}
\seealso{
The function \code{\link[nnet]{multinom}} in package \pkg{nnet} also
    fits multinomial baseline-category logit models, but has a slightly less
    convenient output and does not support overdispersion or random
    effects. However, it provides some other options. Baseline-category logit
    models are also supported by the package \pkg{VGAM}, as well as some
    reduced-rank and (semi-parametric) additive generalisations.  The package
    \pkg{mnlogit} estimates logit models in a way optimized for large numbers
    of alternatives.
}
