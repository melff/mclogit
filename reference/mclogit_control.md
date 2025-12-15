# Control Parameters for the Fitting Process

`mclogit.control` returns a list of default parameters that control the
fitting process of `mclogit`.

## Usage

``` r
mclogit.control(epsilon = 1e-08,
                maxit = 25, trace=TRUE)
mmclogit.control(epsilon = 1e-08,
                 maxit = 25, trace=TRUE,
                 trace.inner=FALSE,
                 avoid.increase = FALSE,
                 break.on.increase = FALSE,
                 break.on.infinite = FALSE,
                 break.on.negative = FALSE,
                 inner.optimizer = "nlminb",
                 maxit.inner = switch(inner.optimizer,
                                      SANN          = 10000,
                                      `Nelder-Mead` = 500,
                                      100),
                 CG.type = 1,
                 NM.alpha = 1,
                 NM.beta = 0.5,
                 NM.gamma = 2.0,
                 SANN.temp = 10,
                 SANN.tmax = 10,
                 grtol = 1e-6,
                 xtol = 1e-8,
                 maxeval = 100,
                 gradstep = c(1e-6, 1e-8),
                 use.gradient = c("analytic","numeric"))
```

## Arguments

- epsilon:

  positive convergence tolerance \\\epsilon\\; the iterations converge
  when \\\|dev - dev\_{old}\|/(\|dev\| + 0.1) \< \epsilon\\.

- maxit:

  integer giving the maximal number of IWLS or PQL iterations.

- trace:

  logical indicating if output should be produced for each iteration.

- trace.inner:

  logical; indicating if output should be produced for each inner
  iteration of the PQL method.

- avoid.increase:

  logical; should an increase of the deviance be avoided by step
  truncation?

- break.on.increase:

  logical; should an increase of the deviance be avoided by stopping the
  algorithm?

- break.on.infinite:

  logical; should an infinite deviance stop the algorithm instead of
  leading to step truncation?

- break.on.negative:

  logical; should a negative deviance stop the algorithm?

- inner.optimizer:

  a character string, one of "nlminb", "nlm", "ucminf", "Nelder-Mead",
  "BFGS", "CG", "L-BFGS-B", "SANN". See
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html),
  [`nlm`](https://rdrr.io/r/stats/nlm.html),
  [`ucminf`](https://rdrr.io/pkg/ucminf/man/ucminf.html), or
  [`optim`](https://rdrr.io/r/stats/optim.html).

- maxit.inner:

  integer; the maximum number of inner iterations

- CG.type:

  integer; the `type` argument passed to
  [`optim`](https://rdrr.io/r/stats/optim.html) if "CG" is selected as
  inner optimizer.

- NM.alpha:

  integer; the `alpha` argument passed to
  [`optim`](https://rdrr.io/r/stats/optim.html) if "Nelder-Mead" is
  selected as inner optimizer.

- NM.beta:

  integer; the `beta` argument passed to
  [`optim`](https://rdrr.io/r/stats/optim.html) if "Nelder-Mead" is
  selected as inner optimizer.

- NM.gamma:

  integer; the `gamma` argument passed to
  [`optim`](https://rdrr.io/r/stats/optim.html) if "Nelder-Mead" is
  selected as inner optimizer.

- SANN.temp:

  integer; the `temp` argument passed to
  [`optim`](https://rdrr.io/r/stats/optim.html) if "SANN" is selected as
  inner optimizer.

- SANN.tmax:

  integer; the `tmax` argument passed to
  [`optim`](https://rdrr.io/r/stats/optim.html) if "SANN" is selected as
  inner optimizer.

- grtol:

  numeric; the `grtol` control parameter for `ucminf` if "ucminf" is
  selected as inner optimizer.

- xtol:

  numeric; the `xtol` control parameter for `ucminf` if "ucminf" is
  selected as inner optimizer.

- maxeval:

  integer; the `maxeval` control parameter for `ucminf` if "ucminf" is
  selected as inner optimizer.

- gradstep:

  a numeric vector of length; the `gradstep` control parameter for
  `ucminf` if "ucminf" is selected as inner optimizer.

- use.gradient:

  a character string; whether the gradient should be computed
  analytically or whether a finite-difference approximation should be
  used.

## Value

A list.
