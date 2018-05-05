\name{NEWS}
\title{\emph{mclogit} News}
\encoding{UTF-8}

\section{Version 0.6}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{mclogit} now complains about (i.e. throws an error
      exception) when the random effects structure cannot be estimated,
      e.g. because random effects are constant within choice sets
      and therefore drop out by the formation of conditional logits.
    }
  }
}

\section{Version 0.5}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{mclogit} now supports conditional logit models with
      random slopes.
      \item \code{mblogit} now supports multinomial baseline models with
      random intercept and random slopes.
    }
  }
  \subsection{BUGFIXES}{
  \itemize{
    \item \code{predict} methods of objects created by \code{mclogit}
    and \code{mblogit} are better in handling missing data.
    }
  }
}

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
    \item The number of observations is now displayed where it was not before.
    \item \code{nobs} is defined in \code{mclogit.fit.rePQL}.
    }
  }
  \subsection{USER-VISIBLE CHANGES}{
   \itemize{
     \item \code{mclogit.fit} and \code{mclogit.fit.rePQL} are exported to enable their use by other packages.
   }
  }
}

