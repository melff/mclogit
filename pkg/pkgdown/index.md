
#  mclogit: Multinomial Logit Models for Categorical Responses and Discrete Choices

[![CRAN](http://www.r-pkg.org/badges/version/mclogit)](http://cran.rstudio.com/package=mclogit)
[![Total downloads from RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/grand-total/mclogit)](http://cran.r-project.org/web/packages/mclogit/index.html)
[![Current release on GitHub](http://img.shields.io/github/release/melff/mclogit.svg)](http://github.com/melff/mclogit/releases/)
[![Total downloads from RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/mclogit)](http://cran.r-project.org/web/packages/mclogit/index.html)

This packages provides estimators for multinomial logit models in their
conditional logit (for discrete choices) and baseline logit variants (for
categorical responses), optionally with overdispersion or random effects.
Random effects models are estimated using the PQL technique (based on a Laplace
approximation) or the MQL technique (based on a Solomon-Cox approximation).
Estimates should be treated with caution if the group sizes are small.

### The statistical models

- [Conditional logit models](articles/conditional-logit.html)
- [Baseline-category logit models](articles/baseline-logit.html)
- [The relation between baseline logit and conditional logit models](articles/baseline-and-conditional-logit.html)
- [Random effects in baseline logit models and conditional logit models](articles/random-effects.html)

### Technical aspects of model fitting

- [The IWLS algorithm used to fit conditional logit models](articles/fitting-mclogit.html)
- [Approximate Inference for Multinomial Logit Models with Random Effects](articles/approximations.html)
- [Bias reduction using Firth's penalized likelihood technique](articles/Firth-bias-reduction.html)
