# Simulating responses from baseline-category and conditional logit models

The [`simulate()`](https://rdrr.io/r/stats/simulate.html) methods allow
to simulate responses from models fitted with
[`mclogit()`](https://melff.github.io/mclogit/reference/mclogit.md) and
[`mblogit()`](https://melff.github.io/mclogit/reference/mblogit.md).
Currently only models *without* random effects are supported for this.

## Usage

``` r
# S3 method for class 'mblogit'
simulate(object, nsim = 1, seed = NULL, ...)
# S3 method for class 'mclogit'
simulate(object, nsim = 1, seed = NULL, ...)

# These methods are currently just 'stubs', causing an error
# message stating that simulation from models with random
# effects are not supported yet
# S3 method for class 'mmblogit'
simulate(object, nsim = 1, seed = NULL, ...)
# S3 method for class 'mmclogit'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  an object from the relevant class

- nsim:

  a number, specifying the number of simulated responses for each
  observation.

- seed:

  an object specifying if and how the random number generator should be
  initialized ('seeded'). The interpetation of this argument follows the
  default method, see `link[stats]{simulate}`

- ...:

  other arguments, ignored.

## Value

The result of the [`simulate`](https://rdrr.io/r/stats/simulate.html)
method for objects created by
[`mclogit`](https://melff.github.io/mclogit/reference/mclogit.md) is a
data frame with one variable for each requested simulation run (their
number is given by the `nsim=` argument). The contents of the columns
are counts (or zero-one values), with group-wise multinomial
distribution (within choice sets) just like it is assumed for the
original response.

The shape of the result of the
[`simulate`](https://rdrr.io/r/stats/simulate.html) method for objects
created by
[`mblogit`](https://melff.github.io/mclogit/reference/mblogit.md) is
also a data frame. The variables within the data frame have a mode or
shape that corresponds to the response to which the model was fitted. If
the response is a matrix of counts, then the variables in the data frame
are also matrices of counts. If the response is a factor and
[`mblogit`](https://melff.github.io/mclogit/reference/mblogit.md) was
called with an argument `from.table=FALSE`, the variables in the data
frame are factors with the same factor levels as the response to which
the model was fitted. If instead the function was called with
`from.table=TRUE`, the variables in the data frame are counts, which
represent frequency weights that would result from applying
[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html) to a
contingency table of simulated frequency counts.

## Examples

``` r
library(MASS)
(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, 
                          data = housing,
                          weights=Freq,
                          aggregate=TRUE))
#> 
#> Iteration 1 - deviance = 38.84842 - criterion = 0.992521
#> Iteration 2 - deviance = 38.66222 - criterion = 0.004803721
#> Iteration 3 - deviance = 38.6622 - criterion = 3.782555e-07
#> Iteration 4 - deviance = 38.6622 - criterion = 2.749622e-15
#> converged
#> 
#> Call: mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq, 
#>     aggregate = TRUE)
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
#> Null Deviance:     262.1 
#> Residual Deviance: 38.66
sm <- simulate(house.mblogit,nsim=7)

housing.long <- housing[rep(seq.int(nrow(housing)),housing$Freq),]
(housel.mblogit <- mblogit(Sat ~ Infl + Type + Cont,
                           data=housing.long))
#> 
#> Iteration 1 - deviance = 3474.691 - criterion = 0.5057269
#> Iteration 2 - deviance = 3470.086 - criterion = 0.001326883
#> Iteration 3 - deviance = 3470.084 - criterion = 7.526912e-07
#> Iteration 4 - deviance = 3470.084 - criterion = 1.308345e-12
#> converged
#> 
#> Call: mblogit(formula = Sat ~ Infl + Type + Cont, data = housing.long)
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
sml <- simulate(housel.mblogit,nsim=7)

housing.table <- xtabs(Freq~.,data=housing)
housing.mat <- memisc::to.data.frame(housing.table)
head(housing.mat)
#>     Infl      Type Cont Low Medium High
#> 1    Low     Tower  Low  21     21   28
#> 2 Medium     Tower  Low  34     22   36
#> 3   High     Tower  Low  10     11   36
#> 4    Low Apartment  Low  61     23   17
#> 5 Medium Apartment  Low  43     35   40
#> 6   High Apartment  Low  26     18   54

(housem.mblogit <- mblogit(cbind(Low,Medium,High) ~
                               Infl + Type + Cont,
                           data=housing.mat))
#> 
#> Iteration 1 - deviance = 38.84842 - criterion = 0.992521
#> Iteration 2 - deviance = 38.66222 - criterion = 0.004803721
#> Iteration 3 - deviance = 38.6622 - criterion = 3.782555e-07
#> Iteration 4 - deviance = 38.6622 - criterion = 2.749622e-15
#> converged
#> 
#> Call: mblogit(formula = cbind(Low, Medium, High) ~ Infl + Type + Cont, 
#>     data = housing.mat)
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
#> Null Deviance:     262.1 
#> Residual Deviance: 38.66
smm <- simulate(housem.mblogit,nsim=7)

str(sm)
#> 'data.frame':    0 obs. of  7 variables:
#>  $ sim_1: int 
#>  $ sim_2: int 
#>  $ sim_3: int 
#>  $ sim_4: int 
#>  $ sim_5: int 
#>  $ sim_6: int 
#>  $ sim_7: int 
#>  - attr(*, "seed")= int [1:626] 10403 207 1980538363 -125047968 -820381145 -1073960099 -33499050 1164085917 1218464490 -1045685882 ...
str(sml)
#> 'data.frame':    1681 obs. of  7 variables:
#>  $ sim_1: Factor w/ 3 levels "Low","Medium",..: 2 3 3 3 3 1 3 1 1 1 ...
#>  $ sim_2: Factor w/ 3 levels "Low","Medium",..: 1 3 1 1 1 2 1 1 1 2 ...
#>  $ sim_3: Factor w/ 3 levels "Low","Medium",..: 1 1 2 1 3 2 3 2 3 1 ...
#>  $ sim_4: Factor w/ 3 levels "Low","Medium",..: 1 2 2 2 3 2 3 3 1 3 ...
#>  $ sim_5: Factor w/ 3 levels "Low","Medium",..: 1 3 1 2 1 1 3 3 1 2 ...
#>  $ sim_6: Factor w/ 3 levels "Low","Medium",..: 2 2 1 3 3 3 3 1 1 3 ...
#>  $ sim_7: Factor w/ 3 levels "Low","Medium",..: 2 1 3 1 1 3 2 2 2 3 ...
#>  - attr(*, "seed")= int [1:626] 10403 19 1492790388 -1629574902 -1817687082 -2080485428 -622102478 -1390801052 -1280540579 249411485 ...
str(smm)
#> 'data.frame':    24 obs. of  7 variables:
#>  $ sim_1: int [1:24, 1:3] 26 23 2 53 45 20 15 5 3 22 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  $ sim_2: int [1:24, 1:3] 30 24 10 55 41 25 19 9 3 17 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  $ sim_3: int [1:24, 1:3] 23 25 10 54 46 32 13 8 7 22 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  $ sim_4: int [1:24, 1:3] 22 30 9 52 46 32 13 9 3 15 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  $ sim_5: int [1:24, 1:3] 22 23 3 57 45 27 8 11 4 21 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  $ sim_6: int [1:24, 1:3] 24 22 8 49 36 29 14 10 5 19 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  $ sim_7: int [1:24, 1:3] 35 25 16 53 43 23 18 9 5 20 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "Low" "Medium" "High"
#>  - attr(*, "seed")= int [1:626] 10403 554 1019603289 -1305277535 -470137223 -57077173 1775921955 -181041711 861366705 -1087987206 ...

head(smm[[1]])
#>      Low Medium High
#> [1,]  26     17   27
#> [2,]  23     28   41
#> [3,]   2      7   48
#> [4,]  53     20   28
#> [5,]  45     36   37
#> [6,]  20     22   56
```
