\name{mclogit.control}
\alias{mclogit.control}

\title{Control Parameters for the Fitting Process}
\description{
  \code{mclogit.control} returns a list of default parameters
  that control the fitting process of \code{mclogit}.
}
\usage{
mclogit.control(epsilon = 1e-08,
                maxit = 25, trace=TRUE)
}
\arguments{
   \item{epsilon}{positive convergence tolerance \eqn{\epsilon};
    the iterations converge when
    \eqn{|dev - dev_{old}|/(|dev| + 0.1) < \epsilon}{|dev - devold|/(|dev| + 0.1) < \epsilon}.}
   \item{maxit}{integer giving the maximal number of IWLS or PQL iterations.}
   \item{trace}{logical indicating if output should be produced for each
    iteration.}
}
\value{
  A list.
}
