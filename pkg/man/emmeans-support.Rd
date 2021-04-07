\name{emmeans-support}
\alias{emmeans-support}
\title{emmeans support for \code{mblogit} models}
\description{
  Support for \code{mblogit} and \code{mmblogit} models is provided
  for interfacing with the \bold{emmeans} package, if it is installed. An
  optional \code{mode} argument allows for choosing two modes of estimation:
  \code{"prob"} (the default), which provides estimated multinomial
  probabilities, and \code{"latent"}, which provides the linear predictor, on
  the log scale. For either mode, the variables in the reference grid comprise the
  predictors in the right-hand side of the model, plus the response variable
  (whose levels delineate the multinomial responses).
  
  With \code{mode = "prob"}, the user may optionally specify an iteger value for
  \code{N.sim}, in which case responses are simulated using the model equation.
  The \code{emmeans} results are then summarized using the simulated values
  as if it were a model fitted by MCMC methods.
}
\examples{
if (requireNamespace("MASS") && requireNamespace("emmeans")) {
  house.mbl <- mblogit(Sat ~ Type + Infl, weights = Freq,
      data = MASS::housing)
  print(emmeans::emmeans(house.mbl, ~ Sat | Type)) # defaults to mode = "prob"
  
  set.seed(123)
  emmeans::emmeans(house.mbl, ~ Sat | Type, N.sim = 1000)
}
}




