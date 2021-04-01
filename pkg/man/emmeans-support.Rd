\name{emmeans-support}
\alias{emmeans-support}
\title{emmeans support for \code{mblogit} models}
\description{
  Support for \code{mblogit} and \code{mmblogit} models is provided
  for interfacing with the \bold{emmeans} package, if it is installed. An
  optional \code{mode} argument allows for choosing two modes of estimation:
  \code{"prob"} (the default), which provides estimated multinomial
  probabilities, and \code{"latent"}, which provides the linear predictor, on
  the log scale. These operate the same way as for \code{nnet::multinom}
  models.
}
\examples{
if (requireNamespace("MASS") && requireNamespace("emmeans")) {
  house.mbl <- mblogit(Sat ~ Type + Infl, weights = Freq,
      data = MASS::housing)
  emmeans::emmeans(house.mbl, ~ Sat | Type, mode = "prob")
}
}




