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
                 break.on.negative = FALSE)
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
}
\value{
  A list.
}
