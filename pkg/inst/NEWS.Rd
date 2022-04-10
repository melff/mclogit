\name{NEWS}
\title{\emph{mclogit} News}
\encoding{UTF-8}

\section{Version 0.9}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item It is now possible to estimate models with non-nested
      (e.g. crossed) random effects. Such models can be specified
      by providing a list of formulas as \code{random=} argument to the
      \code{mclogit()} or \code{mblogit()} function.
      \item The left-hand side of conditional logit models can now
      more conveniently specified using the vertical-bar (\code{|})
      operator.
    }
  }
  \subsection{BUGFIXES}{
    \itemize{
      \item Singular initial covariance matrices no longer cause errors.
      \item A warning about unreliable results is issued if
      \code{anova()} is applied to models with random effects.
    }
  }
}

\section{Version 0.8}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item It is now possible to use the MQL estimation technique as an
      alternative to PQL.
      \item As an alternative to extending a logit model with random
      effects, it is now possible to add an over-dispersion parameter to
      the model.
      \item In addition to approximate the ML estimator, MQL and PQL
      have a variant that approximates the REML estimator.
      \item There is now a \code{simulate()} method for objects returned
      by \code{mblogit()} or \code{mclogit()} (but only for those
      without random effects).
      \item Predictions from random-effects models estimated using the
      PQL technique now are now conditional on the random effects
      (unless requested otherwise).
    }
  }
  \subsection{BUGFIXES}{
    \itemize{
      \item \code{mclogit()} now handles empty responses (i.e. counts that
      sum to zero) correclty.
      \item \code{mclogit()} now flags non-numeric response vectors as an
      error.
      \item \code{predict()} now handles scaled independent variables correcty.
    }
  }
  \subsection{IMPROVEMENTS}{
    \itemize{
      \item \code{summary()} shows the number of groups per random
      effects level (if present).
      \item \code{mclogit()} and \code{mblogit()} with random effects now work with
      \code{formula=}-argumements passed in variables.
    }
  }
}

\section{Version 0.7}{
  \subsection{IMPROVEMENTS}{
    \itemize{
      \item The algorithm for fitting random-effects models tended to
      stop prematurely returning the starting values obtained using a
      methods of moments. It has been completely refactored and proceeds
      similar to the PQL algorithm in Professor Brian Ripley's MASS
      package:
      
      Now an inner step, in which a linear mixed model is fitted to a
      working dependent variable is nested into outer step iterations,
      in which the the working dependent variable is
      updated.

      \item Also, the PQL algorithm no longer builds on starting values from a
      no-random-effects model, because surprisingly this makes the
      algorithm more stable and not less.
      
      As a consequence, the algorithm does a much better job at avoiding
      divergence or running into numerical difficulties.

      \item The PQL estimator for random-effects model uses a (inverse)
      Cholesky factor parametrisation, which makes sure that
      random-effects (Co-)Variance matrices are always positive
      (semi-)definite.
    }
  }
}

\section{Version 0.6}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{mclogit} now complains about (i.e. throws an error
      exception) when the random effects structure cannot be estimated,
      e.g. because random effects are constant within choice sets
      and therefore drop out by the formation of conditional logits.
    }
  }
  \subsection{BUGFIXES}{
  \itemize{
    \item \code{mblogit} now handles responses with only two columns.
    \item \code{mblogit} now can deal with matrix responses that have
    rows that sum to zero.
    \item \code{mclogit} and \code{mblogit} now return a component
    named "df.residual" instead of "residual.df".
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

