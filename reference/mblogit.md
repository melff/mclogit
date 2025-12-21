# Baseline-Category Logit Models for Categorical and Multinomial Responses

The function `mblogit` fits baseline-category logit models for
categorical and multinomial count responses with fixed alternatives.

## Usage

``` r
mblogit(
  formula,
  data = parent.frame(),
  random = NULL,
  catCov = c("free", "diagonal", "single"),
  subset,
  weights = NULL,
  offset = NULL,
  na.action = getOption("na.action"),
  model = TRUE,
  x = FALSE,
  y = TRUE,
  contrasts = NULL,
  method = NULL,
  estimator = c("ML", "REML"),
  dispersion = FALSE,
  start = NULL,
  aggregate = FALSE,
  groups = NULL,
  from.table = FALSE,
  Firth = FALSE,
  control = if (length(random)) mmclogit.control(...) else mclogit.control(...),
  ...
)
```

## Arguments

- formula:

  the model formula. The response must be a factor or a matrix of
  counts.

- data:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in `data`,
  the variables are taken from `environment(formula)`, typically the
  environment from which `glm` is called.

- random:

  an optional formula or list of formulas that specify the
  random-effects structure or NULL.

- catCov:

  a character string that specifies optional restrictions on the
  covariances of random effects between the logit equations. "free"
  means no restrictions, "diagonal" means that random effects pertinent
  to different categories are uncorrelated, while "single" means that
  the random effect variances pertinent to all categories are identical.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector.

- offset:

  an optional model offset. If not NULL, must be a matrix if as many
  columns as the response has categories or one less.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. The default is set by the `na.action` setting of
  [`options`](https://rdrr.io/r/base/options.html), and is
  [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.
  The ‘factory-fresh’ default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html). Another possible
  value is `NULL`, no action. Value
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) can be useful.

- model:

  a logical value indicating whether *model frame* should be included as
  a component of the returned value.

- x, y:

  logical values indicating whether the response vector and model matrix
  used in the fitting process should be returned as components of the
  returned value.

- contrasts:

  an optional list. See the `contrasts.arg` of `model.matrix.default`.

- method:

  `NULL` or a character string, either "PQL" or "MQL", specifies the
  type of the quasilikelihood approximation to be used if a
  random-effects model is to be estimated.

- estimator:

  a character string; either "ML" or "REML", specifies which estimator
  is to be used/approximated.

- dispersion:

  a logical value or a character string; whether and how a dispersion
  parameter should be estimated. For details see
  [`dispersion`](https://melff.github.io/mclogit/reference/dispersion.md).

- start:

  an optional matrix of starting values (with as many rows as logit
  equations). If the model has random effects, the matrix should have a
  "VarCov" attribute wtih starting values for the random effects
  (co-)variances. If the random effects model is estimated with the
  "PQL" method, the starting values matrix should also have a
  "random.effects" attribute, which should have the same structure as
  the "random.effects" component of an object returned by `mblogit()`.

- aggregate:

  a logical value; whether to aggregate responses by covariate classes
  and groups before estimating the model if the response variable is a
  factor.

  This will not affect the estimates, but the dispersion and the
  residual degrees of freedom. If `aggregate=TRUE`, the dispersion will
  be relative to a saturated model; it will be much smaller than with
  `aggregate=TRUE`. In particular, with only a single covariate and no
  grouping, the deviance will be close to zero. If `dispersion` is not
  `FALSE`, then the default value of `aggregate` will be `TRUE`. For
  details see
  [`dispersion`](https://melff.github.io/mclogit/reference/dispersion.md).

  This argument has consequences only if the response in `formula` is a
  factor.

- groups:

  an optional formula that specifies groups of observations relevant for
  the estimation of overdispersion. For details see
  [`dispersion`](https://melff.github.io/mclogit/reference/dispersion.md).

- from.table:

  a logical value; should be FALSE. This argument only exists for the
  sake of compatibility and will be removed in the next relase.

- Firth:

  a logical value; whether to use Firth's bias correction.

- control:

  a list of parameters for the fitting process. See
  [`mclogit.control`](https://melff.github.io/mclogit/reference/mclogit_control.md)

- ...:

  arguments to be passed to `mclogit.control` or `mmclogit.control`

## Value

`mblogit` returns an object of class "mblogit", which has almost the
same structure as an object of class
"[glm](https://rdrr.io/r/stats/glm.html)". The difference are the
components `coefficients`, `residuals`, `fitted.values`,
`linear.predictors`, and `y`, which are matrices with number of columns
equal to the number of response categories minus one.

## Details

The function `mblogit` internally rearranges the data into a 'long'
format and uses
[`mclogit.fit`](https://melff.github.io/mclogit/reference/mclogit.fit.md)
to compute estimates. Nevertheless, the 'user data' are unaffected.

## References

Agresti, Alan. 2002. *Categorical Data Analysis.* 2nd ed, Hoboken, NJ:
Wiley. [doi:10.1002/0471249688](https://doi.org/10.1002/0471249688)

Breslow, N.E. and D.G. Clayton. 1993. "Approximate Inference in
Generalized Linear Mixed Models". *Journal of the American Statistical
Association* 88 (421): 9-25.
[doi:10.1080/01621459.1993.10594284](https://doi.org/10.1080/01621459.1993.10594284)

Firth, David. 1993. "Bias Reduction of Maximum Likelihood Estimates".
*Biometrika* 80 (1): 27–38.
[doi:10.1093/biomet/80.1.27](https://doi.org/10.1093/biomet/80.1.27)

## See also

The function [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) in
package nnet also fits multinomial baseline-category logit models, but
has a slightly less convenient output and does not support
overdispersion or random effects. However, it provides some other
options. Baseline-category logit models are also supported by the
package VGAM, as well as some reduced-rank and (semi-parametric)
additive generalisations. The package mnlogit estimates logit models in
a way optimized for large numbers of alternatives.

## Examples

``` r
library(MASS) # For 'housing' data
library(nnet)
library(memisc)
#> Loading required package: lattice
#> 
#> Attaching package: ‘memisc’
#> The following object is masked from ‘package:Matrix’:
#> 
#>     as.array
#> The following objects are masked from ‘package:stats’:
#> 
#>     contr.sum, contr.treatment, contrasts
#> The following object is masked from ‘package:base’:
#> 
#>     as.array

(house.mult<- multinom(Sat ~ Infl + Type + Cont, weights = Freq,
                       data = housing))
#> # weights:  24 (14 variable)
#> initial  value 1846.767257 
#> iter  10 value 1747.045232
#> final  value 1735.041933 
#> converged
#> Call:
#> multinom(formula = Sat ~ Infl + Type + Cont, data = housing, 
#>     weights = Freq)
#> 
#> Coefficients:
#>        (Intercept) InflMedium  InflHigh TypeApartment TypeAtrium TypeTerrace
#> Medium  -0.4192316  0.4464003 0.6649367    -0.4356851  0.1313663  -0.6665728
#> High    -0.1387453  0.7348626 1.6126294    -0.7356261 -0.4079808  -1.4123333
#>         ContHigh
#> Medium 0.3608513
#> High   0.4818236
#> 
#> Residual Deviance: 3470.084 
#> AIC: 3498.084 


(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, weights = Freq,
                         data = housing))
#> 
#> Iteration 1 - deviance = 3493.764 - criterion = 0.9614469
#> Iteration 2 - deviance = 3470.111 - criterion = 0.00681597
#> Iteration 3 - deviance = 3470.084 - criterion = 7.82437e-06
#> Iteration 4 - deviance = 3470.084 - criterion = 7.469596e-11
#> converged
#> 
#> Call: mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq)
#> 
#> Coefficients:
#>             Predictors
#> Logit eqn.    (Intercept)  InflMedium  InflHigh  TypeApartment  TypeAtrium
#>   Medium/Low  -0.4192       0.4464      0.6649   -0.4357         0.1314   
#>   High/Low    -0.1387       0.7349      1.6126   -0.7356        -0.4080   
#>             Predictors
#> Logit eqn.    TypeTerrace  ContHigh
#>   Medium/Low  -0.6666       0.3609 
#>   High/Low    -1.4123       0.4818 
#> 
#> Null Deviance:     3694 
#> Residual Deviance: 3470

summary(house.mult)
#> Call:
#> multinom(formula = Sat ~ Infl + Type + Cont, data = housing, 
#>     weights = Freq)
#> 
#> Coefficients:
#>        (Intercept) InflMedium  InflHigh TypeApartment TypeAtrium TypeTerrace
#> Medium  -0.4192316  0.4464003 0.6649367    -0.4356851  0.1313663  -0.6665728
#> High    -0.1387453  0.7348626 1.6126294    -0.7356261 -0.4079808  -1.4123333
#>         ContHigh
#> Medium 0.3608513
#> High   0.4818236
#> 
#> Std. Errors:
#>        (Intercept) InflMedium  InflHigh TypeApartment TypeAtrium TypeTerrace
#> Medium   0.1729344  0.1415572 0.1863374     0.1725327  0.2231065   0.2062532
#> High     0.1592295  0.1369380 0.1671316     0.1552714  0.2114965   0.2001496
#>         ContHigh
#> Medium 0.1323975
#> High   0.1241371
#> 
#> Residual Deviance: 3470.084 
#> AIC: 3498.084 

summary(house.mblogit)
#> 
#> Call:
#> mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq)
#> 
#> Equation for Medium vs Low:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    -0.4192     0.1729  -2.424 0.015342 *  
#> InflMedium      0.4464     0.1416   3.153 0.001613 ** 
#> InflHigh        0.6649     0.1863   3.568 0.000359 ***
#> TypeApartment  -0.4357     0.1725  -2.525 0.011562 *  
#> TypeAtrium      0.1314     0.2231   0.589 0.555980    
#> TypeTerrace    -0.6666     0.2063  -3.232 0.001230 ** 
#> ContHigh        0.3609     0.1324   2.726 0.006420 ** 
#> 
#> Equation for High vs Low:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    -0.1387     0.1592  -0.871 0.383570    
#> InflMedium      0.7349     0.1369   5.366 8.03e-08 ***
#> InflHigh        1.6126     0.1671   9.649  < 2e-16 ***
#> TypeApartment  -0.7356     0.1553  -4.738 2.16e-06 ***
#> TypeAtrium     -0.4080     0.2115  -1.929 0.053730 .  
#> TypeTerrace    -1.4123     0.2001  -7.056 1.71e-12 ***
#> ContHigh        0.4818     0.1241   3.881 0.000104 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Approximate residual Deviance: 3470 
#> Number of Fisher scoring iterations:  4 
#> Number of observations:  1681 
#> 

mtable(house.mblogit)
#> 
#> Calls:
#> house.mblogit: mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq)
#> 
#> ================================================
#>                          Medium/Low  High/Low   
#> ------------------------------------------------
#>   (Intercept)              -0.419*   -0.139     
#>                            (0.173)   (0.159)    
#>   Infl: Medium/Low          0.446**   0.735***  
#>                            (0.142)   (0.137)    
#>   Infl: High/Low            0.665***  1.613***  
#>                            (0.186)   (0.167)    
#>   Type: Apartment/Tower    -0.436*   -0.736***  
#>                            (0.173)   (0.155)    
#>   Type: Atrium/Tower        0.131    -0.408     
#>                            (0.223)   (0.211)    
#>   Type: Terrace/Tower      -0.667**  -1.412***  
#>                            (0.206)   (0.200)    
#>   Cont: High/Low            0.361**   0.482***  
#>                            (0.132)   (0.124)    
#> ------------------------------------------------
#>   Deviance               3470.1                 
#>   N                      1681                   
#> ================================================
#>   Significance: *** = p < 0.001;   
#>                 ** = p < 0.01; * = p < 0.05  
```
