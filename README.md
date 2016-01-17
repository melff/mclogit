
# mclogit: Mixed conditional logit models in R 

[![CRAN](http://www.r-pkg.org/badges/version/mclogit)](http://cran.rstudio.com/package=mclogit)
[![Build Status](https://travis-ci.org/melff/mclogit.svg?branch=master)](https://travis-ci.org/melff/mclogit) 
[![Build status](https://ci.appveyor.com/api/projects/status/289k656f3jsbotd2?svg=true)](https://ci.appveyor.com/project/melff/mclogit)

This packages provides allows to estimate conditional logit models of binary responses and multinomial counts, with or without alternative-specific random effects (random intercepts only, no random slopes yet). The current implementation of the estimator for random effects variants of the model uses a Laplace approximation (or PQL) approach and thus should be used only if groups sizes are large.

Releases will be published on [CRAN](http://cran.r-project.org/package=mclogit). Pre-releases will be available here in source and binary form. To install from the sources on GitHub you can use `devtools::install_github("melff/mclogit",subdir="pkg")`.

More information about the package can be found at [my site](http://www.elff.eu/software/mclogit/)
