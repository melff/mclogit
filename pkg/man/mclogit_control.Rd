\name{mclogit.control}
\alias{mclogit.control}
\alias{mmclogit.control}

\title{Control Parameters for the Fitting Process}
\description{
  \code{mclogit.control} returns a list of default parameters
  that control the fitting process of \code{mclogit}.
}
\usage{
mclogit.control(epsilon = 1e-08,
                maxit = 25, trace=TRUE)
mmclogit.control(epsilon = 1e-08,
                 maxit = 25, trace=TRUE,
                 trace.inner=FALSE,
                 avoid.increase = FALSE,
                 break.on.increase = FALSE,
                 break.on.infinite = FALSE,
                 break.on.negative = FALSE,
                 inner.optimizer = "nlminb",
                 maxit.inner = switch(inner.optimizer,
                                      SANN          = 10000,
                                      `Nelder-Mead` = 500,
                                      100),
                 CG.type = 1,
                 NM.alpha = 1,
                 NM.beta = 0.5,
                 NM.gamma = 2.0,
                 SANN.temp = 10,
                 SANN.tmax = 10,
                 grtol = 1e-6,
                 xtol = 1e-8,
                 maxeval = 100,
                 gradstep = c(1e-6, 1e-8),
                 use.gradient = c("analytic","numeric"))
}
\arguments{
   \item{epsilon}{positive convergence tolerance \eqn{\epsilon};
    the iterations converge when
    \eqn{|dev - dev_{old}|/(|dev| + 0.1) < \epsilon}{|dev - devold|/(|dev| + 0.1) < \epsilon}.}
   \item{maxit}{integer giving the maximal number of IWLS or PQL iterations.}
   \item{trace}{logical indicating if output should be produced for each
    iteration.}
   \item{trace.inner}{logical; indicating if output should be produced for each
     inner iteration of the PQL method.}
   \item{avoid.increase}{logical; should an increase of the deviance
     be avoided by step truncation?}
   \item{break.on.increase}{logical; should an increase of the deviance
     be avoided by stopping the algorithm?}
   \item{break.on.infinite}{logical; should an infinite deviance
     stop the algorithm instead of leading to step truncation?}
   \item{break.on.negative}{logical; should a negative deviance
     stop the algorithm?}
   \item{inner.optimizer}{a character string, one of
     "nlminb", "nlm", "ucminf", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN".
     See \code{\link[stats]{nlminb}}, \code{\link[stats]{nlm}}, 
     \code{\link[ucminf]{ucminf}}, or
     \code{\link[stats]{optim}}.
   }
   \item{maxit.inner}{integer; the maximum number of inner iterations}
   \item{CG.type}{integer; the \code{type} argument passed to
     \code{\link{optim}} if
     "CG" is selected as inner optimizer.}
   \item{NM.alpha}{integer; the \code{alpha} argument passed to
     \code{\link{optim}} if
     "Nelder-Mead" is selected as inner optimizer.}
   \item{NM.beta}{integer; the \code{beta} argument passed to
     \code{\link{optim}} if
     "Nelder-Mead" is selected as inner optimizer.}
   \item{NM.gamma}{integer; the \code{gamma} argument passed to
     \code{\link{optim}} if
     "Nelder-Mead" is selected as inner optimizer.}
   \item{SANN.temp}{integer; the \code{temp} argument passed to
     \code{\link{optim}} if
     "SANN" is selected as inner optimizer.}
   \item{SANN.tmax}{integer; the \code{tmax} argument passed to
     \code{\link{optim}} if
     "SANN" is selected as inner optimizer.}
   \item{grtol}{numeric; the \code{grtol} control parameter
     for \code{ucminf} if "ucminf" is selected as inner optimizer.}
   \item{xtol}{numeric; the \code{xtol} control parameter
     for \code{ucminf} if "ucminf" is selected as inner optimizer.}
   \item{maxeval}{integer; the \code{maxeval} control parameter
     for \code{ucminf} if "ucminf" is selected as inner optimizer.}
   \item{gradstep}{a numeric vector of length; the \code{gradstep} control parameter
     for \code{ucminf} if "ucminf" is selected as inner optimizer.}
   \item{use.gradient}{a character string; whether the gradient should
     be computed analytically or whether a finite-difference approximation
     should be used.}
}
\value{
  A list.
}
