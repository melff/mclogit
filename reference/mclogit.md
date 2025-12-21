# Conditional Logit Models and Mixed Conditional Logit Models

`mclogit` fits conditional logit models and mixed conditional logit
models to count data and individual choice data, where the choice set
may vary across choice occasions.

Conditional logit models without random effects are fitted by
Fisher-scoring/IWLS. Models with random effects (mixed conditional logit
models) are estimated via maximum likelihood with a simple Laplace
aproximation (aka PQL).

## Usage

``` r
mclogit(formula, data=parent.frame(), random=NULL,
        subset, weights = NULL, offset=NULL, na.action = getOption("na.action"),
        model = TRUE, x = FALSE, y = TRUE, contrasts=NULL,
        method = NULL, estimator=c("ML","REML"),
        dispersion = FALSE,
        start=NULL,
        groups = NULL,
        Firth = FALSE,
        control=if(length(random))
                    mmclogit.control(...)
                else mclogit.control(...), ...)

# S3 method for class 'mclogit'
update(object, formula., dispersion, ...)

# S3 method for class 'mclogit'
summary(object, dispersion = NULL, correlation = FALSE,
        symbolic.cor = FALSE,  ...)
```

## Arguments

- formula:

  a model formula: a symbolic description of the model to be fitted. The
  left-hand side should result in a two-column matrix. The first column
  contains the choice counts or choice indicators (alternative is
  chosen=1, is not chosen=0). The second column contains unique numbers
  for each choice set.

  The left-hand side can either take the form `cbind(choice,set)` or
  (from version 0.9.1) `choice|set`

  If individual-level data is used, choice sets correspond to
  individuals, if aggregated data with choice counts are used, choice
  sets usually correspond to covariate classes.

  The right-hand of the formula contains choice predictors. It should be
  noted that constants are deleted from the formula as are predictors
  that do not vary within choice sets.

- data:

  an optional data frame, list or environment (or object coercible by
  [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a data
  frame) containing the variables in the model. If not found in `data`,
  the variables are taken from `environment(formula)`, typically the
  environment from which `glm` is called.

- random:

  an optional formula or list of formulas that specify the
  random-effects structure or NULL.

- subset:

  an optional vector specifying a subset of observations to be used in
  the fitting process.

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector.

- offset:

  an optional model offset.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. The default is set by the `na.action` setting of
  [`options`](https://rdrr.io/r/base/options.html), and is
  [`na.fail`](https://rdrr.io/r/stats/na.fail.html) if that is unset.
  The ‘factory-fresh’ default is
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html). Another possible
  value is `NULL`, no action. Value
  [`na.exclude`](https://rdrr.io/r/stats/na.fail.html) can be useful.

- start:

  an optional numerical vector of starting values for the conditional
  logit parameters. If the model has random effects, the vector should
  have a "VarCov" attribute wtih starting values for the random effects
  (co-)variances. If the random effects model is estimated with the
  "PQL" method, the starting values matrix should also have a
  "random.effects" attribute, which should have the same structure as
  the "random.effects" component of an object returned by
  [`mblogit()`](https://melff.github.io/mclogit/reference/mblogit.md).

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

  a real number used as dispersion parameter; a character vector that
  specifies the method to compute the dispersion; a logical value – if
  `TRUE` the default method (`"Afroz"`) is used, if `FALSE`, the
  dispersion parameter is set to 1, that is, no dispersion. For details
  see
  [`dispersion`](https://melff.github.io/mclogit/reference/dispersion.md).

- groups:

  an optional formula that specifies groups of observations relevant for
  the estimation of overdispersion. Covariates should be constant within
  groups, otherwise a warning is generated since the overdispersion
  estimate may be imprecise.

- control:

  a list of parameters for the fitting process. See
  [`mclogit.control`](https://melff.github.io/mclogit/reference/mclogit_control.md)

- Firth:

  a logical value; whether to use Firth's bias correction

- ...:

  arguments to be passed to `mclogit.control` or `mmclogit.control`

- object:

  an object that inherits class `"mclogit"`. When passed to
  [`dispersion()`](https://melff.github.io/mclogit/reference/dispersion.md),
  it should be the result of a call of `mclogit()` of
  [`mblogit()`](https://melff.github.io/mclogit/reference/mblogit.md),
  *without* random effects.

- formula.:

  a changes to the model formula, see
  [`update.default`](https://rdrr.io/r/stats/update.html) and
  [`update.formula`](https://rdrr.io/r/stats/update.formula.html).

- correlation:

  logical; see [`summary.lm`](https://rdrr.io/r/stats/summary.lm.html).

- symbolic.cor:

  logical; see [`summary.lm`](https://rdrr.io/r/stats/summary.lm.html).

## Value

`mclogit` returns an object of class "mclogit", which has almost the
same structure as an object of class
"[glm](https://rdrr.io/r/stats/glm.html)".

## Note

Covariates that are constant within choice sets are automatically
dropped from the model formula specified by the `formula` argument of
`mclogit`.

If the model contains random effects, these should

- either vary within choice sets (e.g. the levels of a factor that
  defines the choice sets should not be nested within the levels of
  factor)

- or be random coefficients of covariates that vary within choice sets.

In earlier versions of the package (prior to 0.6) it will lead to a
failure of the model fitting algorithm if these conditions are not
satisfied. Since version 0.6 of the package, the function `mclogit` will
complain about such model a misspecification explicitely.

From version 0.9.7 it is possible to choose the optimization technique
used for the inner iterations of the PQL/MQL: either
[`nlminb`](https://rdrr.io/r/stats/nlminb.html) (the default),
[`nlm`](https://rdrr.io/r/stats/nlm.html), or any of the algorithms
(other than "Brent" supported by
[`optim`](https://rdrr.io/r/stats/optim.html)). To choose the optimizer,
use the appropriate argument for
[`mmclogit.control`](https://melff.github.io/mclogit/reference/mclogit_control.md)
.

## References

Agresti, Alan (2002). *Categorical Data Analysis.* 2nd ed, Hoboken, NJ:
Wiley. [doi:10.1002/0471249688](https://doi.org/10.1002/0471249688)

Breslow, N.E. and D.G. Clayton (1993). "Approximate Inference in
Generalized Linear Mixed Models". *Journal of the American Statistical
Association* 88 (421): 9-25.
[doi:10.1080/01621459.1993.10594284](https://doi.org/10.1080/01621459.1993.10594284)

Elff, Martin (2009). "Social Divisions, Party Positions, and Electoral
Behaviour". *Electoral Studies* 28(2): 297-308.
[doi:10.1016/j.electstud.2009.02.002](https://doi.org/10.1016/j.electstud.2009.02.002)

Firth, David. 1993. "Bias Reduction of Maximum Likelihood Estimates".
*Biometrika* 80 (1): 27–38.
[doi:10.1093/biomet/80.1.27](https://doi.org/10.1093/biomet/80.1.27)

McFadden, D. (1973). "Conditionial Logit Analysis of Qualitative Choice
Behavior". Pp. 105-135 in P. Zarembka (ed.). *Frontiers in
Econometrics*. New York: Wiley.
<https://eml.berkeley.edu/reprints/mcfadden/zarembka.pdf>

## See also

Conditional logit models are also supported by gmnl, mlogit, and
survival. survival supports conditional logit models for binary panel
data and case-control studies. mlogit and gmnl treat conditional logit
models from an econometric perspective. Unlike the present package, they
focus on the random utility interpretation of discrete choice models and
support generalisations of conditional logit models, such as nested
logit models, that are intended to overcome the IIA (indipendence from
irrelevant alterantives) assumption. Mixed multinomial models are also
supported and estimated using simulation-based techniques. Unlike the
present package, mixed or random-effects extensions are mainly intended
to fit repeated choices of the same individuals and not aggregated
choices of many individuals facing identical alternatives.

## Examples

``` r
data(Transport)

summary(mclogit(
  cbind(resp,suburb)~distance+cost,
  data=Transport
  ))
#> 
#> Iteration 1 - deviance = 39.74973 - criterion = 0.8590917
#> Iteration 2 - deviance = 10.50328 - criterion = 2.758244
#> Iteration 3 - deviance = 9.231325 - criterion = 0.1363107
#> Iteration 4 - deviance = 9.227742 - criterion = 0.0003840654
#> Iteration 5 - deviance = 9.227742 - criterion = 3.446459e-09
#> converged
#> 
#> Call:
#> mclogit(formula = cbind(resp, suburb) ~ distance + cost, data = Transport)
#> 
#>          Estimate Std. Error z value Pr(>|z|)    
#> distance -1.43940    0.05318  -27.07   <2e-16 ***
#> cost     -0.97753    0.03987  -24.52   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Null Deviance:     2734 
#> Residual Deviance: 9.228 
#> Number of Fisher Scoring iterations:  5 
#> Number of observations:  1994 
#> 
#> 
# New syntactic sugar:
summary(mclogit(
  resp|suburb~distance+cost,
  data=Transport
  ))
#> 
#> Iteration 1 - deviance = 39.74973 - criterion = 0.8590917
#> Iteration 2 - deviance = 10.50328 - criterion = 2.758244
#> Iteration 3 - deviance = 9.231325 - criterion = 0.1363107
#> Iteration 4 - deviance = 9.227742 - criterion = 0.0003840654
#> Iteration 5 - deviance = 9.227742 - criterion = 3.446459e-09
#> converged
#> 
#> Call:
#> mclogit(formula = resp | suburb ~ distance + cost, data = Transport)
#> 
#>          Estimate Std. Error z value Pr(>|z|)    
#> distance -1.43940    0.05318  -27.07   <2e-16 ***
#> cost     -0.97753    0.03987  -24.52   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Null Deviance:     2734 
#> Residual Deviance: 9.228 
#> Number of Fisher Scoring iterations:  5 
#> Number of observations:  1994 
#> 
#> 


if (FALSE)  # This takes a bit longer.
data(electors)

electors <- within(electors,{
    party.time <-interaction(party,time)
    time.class <- interaction(time,class)
})

# Time points nested within parties
summary(mclogit(
  Freq|time.class~econ.left/class+welfare/class+auth/class,
  random=~1|party/time,
  data=electors))
#> 
#> Warning: Inner iterations did not coverge - nlminb message: singular convergence (7)
#> 
#> Iteration 1 - deviance = 495.6469 - criterion = 0.1640698
#> Iteration 2 - deviance = 379.0387 - criterion = 0.02944294
#> Iteration 3 - deviance = 363.3644 - criterion = 0.006445485
#> Iteration 4 - deviance = 362.7738 - criterion = 0.0003382211
#> Iteration 5 - deviance = 362.7685 - criterion = 6.930491e-07
#> Iteration 6 - deviance = 362.7684 - criterion = 2.672469e-12
#> converged
#> 
#> Call:
#> mclogit(formula = Freq | time.class ~ econ.left/class + welfare/class + 
#>     auth/class, data = electors, random = ~1 | party/time)
#> 
#> Coefficents:
#>                           Estimate Std. Error z value Pr(>|z|)    
#> econ.left                 -0.04035    0.73004  -0.055   0.9559    
#> welfare                    1.95936    1.16585   1.681   0.0928 .  
#> auth                       0.16811    0.62796   0.268   0.7889    
#> econ.left:classnew.middle -2.04700    0.11816 -17.324   <2e-16 ***
#> econ.left:classold.middle -3.39457    0.17416 -19.491   <2e-16 ***
#> classnew.middle:welfare   -0.75181    0.07501 -10.023   <2e-16 ***
#> classold.middle:welfare   -1.27119    0.14517  -8.757   <2e-16 ***
#> classnew.middle:auth      -1.49715    0.05188 -28.855   <2e-16 ***
#> classold.middle:auth       1.40982    0.05997  23.510   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Co-)Variances:
#> Grouping level: party 
#>          Estimate   Std.Err.
#>          (Const.)   (Const.)
#> (Const.)  2.002      2.143  
#> 
#> Grouping level: party:time 
#>          Estimate    Std.Err.
#>          (Const.)    (Const.)
#> (Const.) 2.372e-08   1.09e-24
#> 
#> Approximate residual deviance: 362.8 
#> Number of Fisher scoring iterations:  6
#> Number of observations
#>   Groups by party: 6
#>   Groups by party:time: 150
#>   Individual observations:  37500
#> 

# Party-level random intercepts and random slopes varying over time points
summary(mclogit(
  Freq|time.class~econ.left/class+welfare/class+auth/class,
  random=list(~1|party,~econ.left+0|time),
  data=electors))
#> 
#> Iteration 1 - deviance = 495.6162 - criterion = 0.1640703
#> Iteration 2 - deviance = 378.9601 - criterion = 0.02945638
#> Iteration 3 - deviance = 363.1957 - criterion = 0.006453599
#> Iteration 4 - deviance = 362.5831 - criterion = 0.0003390485
#> Iteration 5 - deviance = 362.5743 - criterion = 6.960937e-07
#> Iteration 6 - deviance = 362.574 - criterion = 2.692628e-12
#> converged
#> 
#> Call:
#> mclogit(formula = Freq | time.class ~ econ.left/class + welfare/class + 
#>     auth/class, data = electors, random = list(~1 | party, ~econ.left + 
#>     0 | time))
#> 
#> Coefficents:
#>                           Estimate Std. Error z value Pr(>|z|)    
#> econ.left                 -0.04034    0.73006  -0.055   0.9559    
#> welfare                    1.95937    1.16587   1.681   0.0928 .  
#> auth                       0.16811    0.62797   0.268   0.7889    
#> econ.left:classnew.middle -2.04703    0.11816 -17.324   <2e-16 ***
#> econ.left:classold.middle -3.39462    0.17416 -19.492   <2e-16 ***
#> classnew.middle:welfare   -0.75181    0.07501 -10.023   <2e-16 ***
#> classold.middle:welfare   -1.27119    0.14517  -8.757   <2e-16 ***
#> classnew.middle:auth      -1.49715    0.05188 -28.855   <2e-16 ***
#> classold.middle:auth       1.40982    0.05997  23.510   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Co-)Variances:
#> Grouping level: party 
#>          Estimate   Std.Err.
#>          (Const.)   (Const.)
#> (Const.)  2.002      2.143  
#> 
#> Grouping level: time 
#>           Estimate    Std.Err. 
#>           econ.left   econ.left
#> econ.left 0.0002345   2.579e-12
#> 
#> Approximate residual deviance: 362.6 
#> Number of Fisher scoring iterations:  6
#> Number of observations
#>   Groups by party: 6
#>   Groups by time: 25
#>   Individual observations:  37500
#> 
 # \dontrun{}
```
