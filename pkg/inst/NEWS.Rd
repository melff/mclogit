\name{NEWS}
\title{\emph{memisc} News}
\encoding{UTF-8}

\section{Version 0.4}{
  \subsection{NEW FEATURES}{
  \itemize{
    \item New \code{nobs} and \code{extractAIC} methods for \code{mclogit} objects, so that \code{drop1.default} should work with these.
    \item New function \code{mblogit} to fit multinomial baseline logit models.
    \item \code{mclogit} \code{mclogit.fit} now allow user-provided starting values.
    }
  }
  \subsection{BUGFIXES}{
  \itemize{
    \item \code{getSummary} methods now return "contrasts" and "xlevels" components.
    \item Fixed prediction method for \code{mclogit} results.
    \item Corrected handling of weights and standard errors of prediction.
    \item Matrices returned by the \code{mclogit} method of \code{vcov()} have
      row and column names.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
   \itemize{
     \item \code{mclogit.fit} and \code{mclogit.fit.rePQL} are exported to enable their use by other packages.
   }
  }
}

