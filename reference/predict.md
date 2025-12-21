# Predicting responses or linear parts of the baseline-category and conditional logit models

The [`predict()`](https://rdrr.io/r/stats/predict.html) methods allow to
obtain within-sample and out-of-sample predictions from models fitted
with [`mclogit()`](https://melff.github.io/mclogit/reference/mclogit.md)
and [`mblogit()`](https://melff.github.io/mclogit/reference/mblogit.md).

For models with random effecs fitted using the PQL-method, it is
possible to obtain responses that are conditional on the reconstructed
random effects.

## Usage

``` r
# S3 method for class 'mblogit'
predict(object, newdata=NULL,type=c("link","response"),se.fit=FALSE, ...)
# S3 method for class 'mclogit'
predict(object, newdata=NULL,type=c("link","response"),se.fit=FALSE, ...)
# S3 method for class 'mmblogit'
predict(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,
                             conditional=TRUE, ...)
# S3 method for class 'mmclogit'
predict(object, newdata=NULL,type=c("link","response"),se.fit=FALSE,
                             conditional=TRUE, ...)
```

## Arguments

- object:

  an object in class "mblogit", "mmblogit", "mclogit", or "mmclogit"

- newdata:

  an optional data frame with new data

- type:

  a character string specifying the kind of prediction

- se.fit:

  a logical value; whether predictions should be accompanied with
  standard errors

- conditional:

  a logical value; whether predictions should be made conditional on the
  random effects (or whether they are set to zero, i.e. their
  expectation). This argument is consequential only if the "mmblogit" or
  "mmclogit" object was created with `method="PQL"`.

- ...:

  other arguments, ignored.

## Value

The `predict` methods return either a matrix (unless called with
`se.fit=TRUE`) or a list with two matrix-valued elements `"fit"` and
`"se.fit"`.

## Examples

``` r
library(MASS)
(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, 
                          data = housing,
                          weights=Freq))
#> 
#> Iteration 1 - deviance = 3493.764 - criterion = 0.9614469
#> Iteration 2 - deviance = 3470.111 - criterion = 0.00681597
#> Iteration 3 - deviance = 3470.084 - criterion = 7.82437e-06
#> Iteration 4 - deviance = 3470.084 - criterion = 7.46957e-11
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

head(pred.house.mblogit <- predict(house.mblogit))
#>        Medium       High
#> 1 -0.41922874 -0.1387428
#> 2 -0.41922874 -0.1387428
#> 3 -0.41922874 -0.1387428
#> 4  0.02716715  0.5961205
#> 5  0.02716715  0.5961205
#> 6  0.02716715  0.5961205
str(pred.house.mblogit <- predict(house.mblogit,se=TRUE))
#> List of 2
#>  $ fit   : num [1:72, 1:2] -0.4192 -0.4192 -0.4192 0.0272 0.0272 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:72] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:2] "Medium" "High"
#>  $ se.fit: num [1:72, 1:2] 0.173 0.173 0.173 0.17 0.17 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "Medium" "High"

head(pred.house.mblogit <- predict(house.mblogit,
                                   type="response"))
#>         Low    Medium      High
#> 1 0.3955687 0.2601077 0.3443236
#> 2 0.3955687 0.2601077 0.3443236
#> 3 0.3955687 0.2601077 0.3443236
#> 4 0.2602403 0.2674072 0.4723526
#> 5 0.2602403 0.2674072 0.4723526
#> 6 0.2602403 0.2674072 0.4723526
str(pred.house.mblogit <- predict(house.mblogit,se=TRUE,
                                  type="response"))
#> List of 2
#>  $ fit   : num [1:72, 1:3] 0.396 0.396 0.396 0.26 0.26 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:72] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  $ se.fit: num [1:72, 1:3] 0.0343 0.0343 0.0343 0.0273 0.0273 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:72] "1" "2" "3" "4" ...
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
 # This takes a bit longer.
data(electors)
(mcre <- mclogit(
    cbind(Freq,interaction(time,class))~econ.left/class+welfare/class+auth/class,
    random=~1|party.time,
    data=within(electors,party.time<-interaction(party,time))))
#> 
#> Iteration 1 - deviance = 1070.463 - criterion = 0.1596265
#> Iteration 2 - deviance = 965.7808 - criterion = 0.0253941
#> Iteration 3 - deviance = 949.163 - criterion = 0.005112972
#> Iteration 4 - deviance = 947.9638 - criterion = 0.0002016053
#> Iteration 5 - deviance = 947.8468 - criterion = 2.470044e-07
#> Iteration 6 - deviance = 947.8431 - criterion = 4.405134e-13
#> converged
#> mclogit(formula = cbind(Freq, interaction(time, class)) ~ econ.left/class + 
#>     welfare/class + auth/class, data = within(electors, party.time <- interaction(party, 
#>     time)), random = ~1 | party.time)
#> 
#> Coefficients:
#>                 econ.left                    welfare  
#>                  -0.17380                    2.05525  
#>                      auth  econ.left:classnew.middle  
#>                   0.08059                   -1.66428  
#> econ.left:classold.middle    classnew.middle:welfare  
#>                  -2.96667                   -0.99252  
#>   classold.middle:welfare       classnew.middle:auth  
#>                  -1.62032                   -1.39064  
#>      classold.middle:auth  
#>                   1.45728  
#> 
#> (Co-)Variances:
#> Grouping level: party.time 
#>           (Const.)
#> (Const.)  1.604   
#> 
#> Approximate residual deviance: 947.8

str(predict(mcre))
#>  num [1:450] 1.168 4.367 0.148 -1.285 -1.378 ...
str(predict(mcre,type="response"))
#>  num [1:450] 0.03789 0.92862 0.01366 0.00326 0.00297 ...

str(predict(mcre,se.fit=TRUE))
#> List of 2
#>  $ fit   : num [1:450] 1.168 4.367 0.148 -1.285 -1.378 ...
#>  $ se.fit: Named num [1:450] 0.533 0.531 0.609 0.536 0.539 ...
#>   ..- attr(*, "names")= chr [1:450] "1" "2" "3" "4" ...
str(predict(mcre,type="response",se.fit=TRUE))
#> List of 2
#>  $ fit   : num [1:450] 0.03789 0.92862 0.01366 0.00326 0.00297 ...
#>  $ se.fit: num [1:450] 0.004138 0.007417 0.004969 0.000562 0.000509 ...
```
