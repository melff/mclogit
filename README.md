
# mclogit: Mixed conditional logit models in R 

[![CRAN](http://www.r-pkg.org/badges/version/mclogit)](http://cran.rstudio.com/package=mclogit)
[![Total downloads from RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/grand-total/mclogit)](http://cran.r-project.org/web/packages/mclogit/index.html)
[![Current release on GitHub](http://img.shields.io/github/release/melff/mclogit.svg)](http://github.com/melff/mclogit/releases/)
[![Build Status](https://travis-ci.org/melff/mclogit.svg?branch=master)](https://travis-ci.org/melff/mclogit) 

<!-- [![Build status](https://ci.appveyor.com/api/projects/status/289k656f3jsbotd2?svg=true)](https://ci.appveyor.com/project/melff/mclogit) one CI service is enough -->

This packages provides estimators for multinomial logit models in their conditional logit and baseline logit variants, with or without random effects, with or without overdispersion. Random effects models are estimated using the PQL technique (based on a Laplace approximation) or the MQL technique (based on a Solomon-Cox approximation). Estimates should be treated with caution if the group sizes are small.

Releases will be published on [CRAN](http://cran.r-project.org/package=mclogit). Pre-releases will be available here in source and binary form. To install from the sources on GitHub you can use `devtools::install_github("melff/mclogit",subdir="pkg")`.

More information about the package can be found at [my site](http://www.elff.eu/software/mclogit/)
