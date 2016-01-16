
# mclogit: Mixed conditional logit models in R [![Build Status](https://travis-ci.org/melff/mclogit.svg?branch=master)](https://travis-ci.org/melff/mclogit) [![CRAN](http://www.r-pkg.org/badges/version/mclogit)](http://cran.rstudio.com/package=mclogit)


This packages provides allows to estimate conditional logit models of binary responses and multinomial counts, with or without alternative-specific random effects (random intercepts only, no random slopes yet). The current implementation of the estimator for random effects variants of the model uses a Laplace approximation (or PQL) approach and thus should be used only if groups sizes are large.
