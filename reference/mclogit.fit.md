# Internal functions used for model fit.

These functions are exported and documented for use by other packages.
They are not intended for end users.

## Usage

``` r
mclogit.fit(y, s, w, X,
            dispersion=FALSE,
            start = NULL, offset = NULL,
            control = mclogit.control(),
            Firth=FALSE)

mmclogit.fitPQLMQL(y, s, w, X, Z, d, 
                   start = NULL,
                   start.Phi = NULL,
                   start.b = NULL,
                   offset = NULL, method=c("PQL","MQL"),
                   estimator = c("ML","REML"),
                   control = mmclogit.control())
```

## Arguments

- y:

  a response vector. Should be binary.

- s:

  a vector identifying individuals or covariate strata

- w:

  a vector with observation weights.

- X:

  a model matrix; required.

- dispersion:

  a logical value or a character string; whether and how a dispersion
  parameter should be estimated. For details see
  [`dispersion`](https://melff.github.io/mclogit/reference/dispersion.md).

- Z:

  the random effects design matrix.

- d:

  dimension of random effects. Typically \$d=1\$ for random intercepts
  only, \$d\>1\$ for models with random intercepts.

- start:

  an optional numerical vector of starting values for the coefficients.

- offset:

  an optional model offset. Currently only supported for models without
  random effects.

- start.Phi:

  an optional matrix of strarting values for the (co-)variance
  parameters.

- start.b:

  an optional list of vectors with starting values for the random
  effects.

- method:

  a character string, either "PQL" or "MQL", specifies the type of the
  quasilikelihood approximation.

- estimator:

  a character string; either "ML" or "REML", specifies which estimator
  is to be used/approximated.

- control:

  a list of parameters for the fitting process. See
  [`mclogit.control`](https://melff.github.io/mclogit/reference/mclogit_control.md)

- Firth:

  a logical value; whether to use Firth's bias correction

## Value

A list with components describing the fitted model.
