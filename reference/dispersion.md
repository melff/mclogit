# Overdispersion in Multinomial Logit Models

The function `dispersion()` extracts the dispersion parameter from a
multinomial logit model or computes a dispersion parameter estimate
based on a given method. This dispersion parameter can be attached to a
model using [`update()`](https://rdrr.io/r/stats/update.html). It can
also given as an argument to
[`summary()`](https://rdrr.io/r/base/summary.html).

## Usage

``` r
dispersion(object, method, ...)
# S3 method for class 'mclogit'
dispersion(object, method=NULL, ...)
```

## Arguments

- object:

  an object that inherits class `"mclogit"`. When passed to
  `dispersion()`, it should be the result of a call of
  [`mclogit()`](https://melff.github.io/mclogit/reference/mclogit.md) of
  [`mblogit()`](https://melff.github.io/mclogit/reference/mblogit.md),
  *without* random effects.

- method:

  a character string, either `"Afroz"`, `"Fletcher"`, `"Pearson"`, or
  `"Deviance"`, that specifies the estimator of the dispersion; or
  `NULL`, in which case the default estimator, `"Afroz"` is used. The
  estimators are discussed in Afroz et al. (2019).

- ...:

  other arguments, ignored or passed to other methods.

## References

Afroz, Farzana, Matt Parry, and David Fletcher. (2020). "Estimating
Overdispersion in Sparse Multinomial Data." *Biometrics* 76(3): 834-842.
[doi:10.1111/biom.13194](https://doi.org/10.1111/biom.13194) .

## Examples

``` r
library(MASS) # For 'housing' data

# Note that with a factor response and frequency weighted data,
# Overdispersion will be overestimated:
house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, 
                         weights = Freq,
                         data = housing)
#> 
#> Iteration 1 - deviance = 3493.764 - criterion = 0.9614469
#> Iteration 2 - deviance = 3470.111 - criterion = 0.00681597
#> Iteration 3 - deviance = 3470.084 - criterion = 7.82437e-06
#> Iteration 4 - deviance = 3470.084 - criterion = 7.469596e-11
#> converged
dispersion(house.mblogit, method = "Afroz")
#> [1] 20.45129
dispersion(house.mblogit, method = "Deviance")
#> [1] 26.69295

# In order to be able to estimate overdispersion accurately,
# data like the above (which usually comes from applying
# 'as.data.frame' to a contingency table) the model has to be
# fitted with the optional argument 'aggregate=TRUE' or 
# by requesting the dispersion in advance.
house.mblogit.agg <- mblogit(Sat ~ Infl + Type + Cont, 
                             weights = Freq,
                             data = housing, 
                             aggregate = TRUE)
#> 
#> Iteration 1 - deviance = 38.84842 - criterion = 0.992521
#> Iteration 2 - deviance = 38.66222 - criterion = 0.004803721
#> Iteration 3 - deviance = 38.6622 - criterion = 3.782555e-07
#> Iteration 4 - deviance = 38.6622 - criterion = 3.666163e-15
#> converged
# Now the estimated dispersion parameter is no longer larger than 20,
# but just bit over 1.0.
dispersion(house.mblogit.agg, method = "Afroz")
#> [1] 1.121765
dispersion(house.mblogit.agg, method = "Deviance")
#> [1] 1.137124

# It is possible to obtain the dispersion after estimating the coefficients:
phi.Afroz <- dispersion(house.mblogit.agg, method = "Afroz")
summary(house.mblogit.agg, dispersion = phi.Afroz)
#> 
#> Call:
#> mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq, 
#>     aggregate = TRUE)
#> 
#> Equation for Medium vs Low:
#>               Estimate Std. Error t value Pr(>|t|)   
#> (Intercept)    -0.4192     0.1729  -2.424  0.02081 * 
#> InflMedium      0.4464     0.1416   3.153  0.00336 **
#> InflHigh        0.6649     0.1863   3.568  0.00109 **
#> TypeApartment  -0.4357     0.1725  -2.525  0.01639 * 
#> TypeAtrium      0.1314     0.2231   0.589  0.55987   
#> TypeTerrace    -0.6666     0.2063  -3.232  0.00273 **
#> ContHigh        0.3609     0.1324   2.726  0.01007 * 
#> 
#> Equation for High vs Low:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)    -0.1387     0.1592  -0.871 0.389681    
#> InflMedium      0.7349     0.1369   5.366 5.74e-06 ***
#> InflHigh        1.6126     0.1671   9.649 2.89e-11 ***
#> TypeApartment  -0.7356     0.1553  -4.738 3.75e-05 ***
#> TypeAtrium     -0.4080     0.2115  -1.929 0.062112 .  
#> TypeTerrace    -1.4123     0.2001  -7.056 3.79e-08 ***
#> ContHigh        0.4818     0.1241   3.881 0.000454 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Dispersion:  1.121765  on  34  degrees of freedom
#> Approximate residual Deviance: 38.66 
#> Number of Fisher scoring iterations:  4 
#> Number of observations:  1681 
#> 

summary(update(house.mblogit.agg, dispersion = "Afroz"))
#> 
#> Call:
#> mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq, 
#>     aggregate = TRUE)
#> 
#> Equation for Medium vs Low:
#>               Estimate Std. Error t value Pr(>|t|)   
#> (Intercept)    -0.4192     0.1832  -2.289  0.02842 * 
#> InflMedium      0.4464     0.1499   2.977  0.00533 **
#> InflHigh        0.6649     0.1974   3.369  0.00189 **
#> TypeApartment  -0.4357     0.1827  -2.384  0.02284 * 
#> TypeAtrium      0.1314     0.2363   0.556  0.58189   
#> TypeTerrace    -0.6666     0.2184  -3.051  0.00440 **
#> ContHigh        0.3609     0.1402   2.573  0.01461 * 
#> 
#> Equation for High vs Low:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)    -0.1387     0.1686  -0.823 0.416418    
#> InflMedium      0.7349     0.1450   5.067 1.41e-05 ***
#> InflHigh        1.6126     0.1770   9.110 1.20e-10 ***
#> TypeApartment  -0.7356     0.1645  -4.473 8.19e-05 ***
#> TypeAtrium     -0.4080     0.2240  -1.821 0.077370 .  
#> TypeTerrace    -1.4123     0.2120  -6.662 1.20e-07 ***
#> ContHigh        0.4818     0.1315   3.665 0.000837 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Dispersion:  1.121765  on  34  degrees of freedom
#> Approximate residual Deviance: 38.66 
#> Number of Fisher scoring iterations:  4 
#> Number of observations:  1681 
#> 

# If an estimate of the (over-)dispersion is requested, 'aggregate' is set to
# TRUE by default:
house.mblogit.odsp <- mblogit(Sat ~ Infl + Type + Cont, 
                              weights = Freq,
                              data = housing, 
                              dispersion = "Afroz")
#> 
#> Iteration 1 - deviance = 38.84842 - criterion = 0.992521
#> Iteration 2 - deviance = 38.66222 - criterion = 0.004803721
#> Iteration 3 - deviance = 38.6622 - criterion = 3.782555e-07
#> Iteration 4 - deviance = 38.6622 - criterion = 3.666163e-15
#> converged
summary(house.mblogit.odsp)
#> 
#> Call:
#> mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq, 
#>     dispersion = "Afroz")
#> 
#> Equation for Medium vs Low:
#>               Estimate Std. Error t value Pr(>|t|)   
#> (Intercept)    -0.4192     0.1832  -2.289  0.02842 * 
#> InflMedium      0.4464     0.1499   2.977  0.00533 **
#> InflHigh        0.6649     0.1974   3.369  0.00189 **
#> TypeApartment  -0.4357     0.1827  -2.384  0.02284 * 
#> TypeAtrium      0.1314     0.2363   0.556  0.58189   
#> TypeTerrace    -0.6666     0.2184  -3.051  0.00440 **
#> ContHigh        0.3609     0.1402   2.573  0.01461 * 
#> 
#> Equation for High vs Low:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)    -0.1387     0.1686  -0.823 0.416418    
#> InflMedium      0.7349     0.1450   5.067 1.41e-05 ***
#> InflHigh        1.6126     0.1770   9.110 1.20e-10 ***
#> TypeApartment  -0.7356     0.1645  -4.473 8.19e-05 ***
#> TypeAtrium     -0.4080     0.2240  -1.821 0.077370 .  
#> TypeTerrace    -1.4123     0.2120  -6.662 1.20e-07 ***
#> ContHigh        0.4818     0.1315   3.665 0.000837 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Dispersion:  1.121765  on  34  degrees of freedom
#> Approximate residual Deviance: 38.66 
#> Number of Fisher scoring iterations:  4 
#> Number of observations:  1681 
#> 
dispersion(house.mblogit.odsp, method = "Deviance")
#> [1] 1.137124

# Note that aggregation (either implicitly or explicitly required) affects
# the reported deviance in surprising ways:
house.mblogit.o.00 <- mblogit(Sat ~ Infl, 
                              weights = Freq,
                              data = housing, 
                              dispersion = TRUE)
#> 
#> Iteration 1 - deviance = 2.084495e-10 - criterion = 0.01687143
#> Iteration 2 - deviance = 6.794565e-14 - criterion = 2.083815e-09
#> converged
deviance(house.mblogit.o.00)
#> [1] 6.794565e-14
dispersion(house.mblogit.o.00)
#> [1] Inf
# The deviance is (almost) zero, because aggregation leads to a two-way
# table and a single-predictor model is already saturated.

# In order to make models comparable, one will need to set the groups:
house.mblogit.o.0 <- mblogit(Sat ~ Infl, 
                             weights = Freq,
                             data = housing, 
                             groups = ~ Infl + Type + Cont,
                             dispersion = TRUE)
#> Warning: Argument 'groups' is inconsequential unless aggregate=TRUE
#> 
#> Iteration 1 - deviance = 111.578 - criterion = 0.9973916
#> Iteration 2 - deviance = 111.0847 - criterion = 0.004437062
#> Iteration 3 - deviance = 111.0846 - criterion = 5.432493e-07
#> Iteration 4 - deviance = 111.0846 - criterion = 9.969427e-15
#> converged
deviance(house.mblogit.o.0)
#> [1] 111.0846
dispersion(house.mblogit.o.0)
#> [1] 2.567212

anova(house.mblogit.o.0,
      house.mblogit.odsp)
#> Analysis of Deviance Table
#> 
#> Model 1: Sat ~ Infl
#> Model 2: Sat ~ Infl + Type + Cont
#>   Resid. Df Resid. Dev Df Deviance
#> 1        42    111.085            
#> 2        34     38.662  8   72.422

# These complications with the deviances do not arise if no aggregation is 
# requested:
house.mblogit.0 <- mblogit(Sat ~ Infl, 
                           weights = Freq,
                           data = housing)
#> 
#> Iteration 1 - deviance = 3557.711 - criterion = 0.9621399
#> Iteration 2 - deviance = 3542.511 - criterion = 0.004290705
#> Iteration 3 - deviance = 3542.506 - criterion = 1.326243e-06
#> Iteration 4 - deviance = 3542.506 - criterion = 6.94199e-13
#> converged
anova(house.mblogit.0,
      house.mblogit)
#> Analysis of Deviance Table
#> 
#> Model 1: Sat ~ Infl
#> Model 2: Sat ~ Infl + Type + Cont
#>   Resid. Df Resid. Dev Df Deviance
#> 1       138     3542.5            
#> 2       130     3470.1  8   72.422


# Using frequences on the left-hand side is perhaps the safest option:
housing.wide <- memisc::Aggregate(table(Sat) ~ Infl + Type + Cont,
                                  data = housing) # Note that 'Aggegate' uses
                                                # variable 'Freq' for weighting.
house.mblogit.wide <- mblogit(cbind(Low,Medium,High) ~ Infl + Type + Cont, 
                              data = housing.wide)
#> 
#> Iteration 1 - deviance = 38.84842 - criterion = 0.992521
#> Iteration 2 - deviance = 38.66222 - criterion = 0.004803721
#> Iteration 3 - deviance = 38.6622 - criterion = 3.782555e-07
#> Iteration 4 - deviance = 38.6622 - criterion = 3.666163e-15
#> converged
summary(house.mblogit.wide)
#> 
#> Call:
#> mblogit(formula = cbind(Low, Medium, High) ~ Infl + Type + Cont, 
#>     data = housing.wide)
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
#> Approximate residual Deviance: 38.66 
#> Number of Fisher scoring iterations:  4 
#> Number of observations:  24 
#> 
dispersion(house.mblogit.wide, method = "Afroz")
#> [1] 1.121765

house.mblogit.wide.0 <- mblogit(cbind(Low,Medium,High) ~ Infl, 
                                data = housing.wide)
#> 
#> Iteration 1 - deviance = 111.578 - criterion = 0.9973916
#> Iteration 2 - deviance = 111.0847 - criterion = 0.004437062
#> Iteration 3 - deviance = 111.0846 - criterion = 5.432493e-07
#> Iteration 4 - deviance = 111.0846 - criterion = 9.969427e-15
#> converged
summary(house.mblogit.wide.0)
#> 
#> Call:
#> mblogit(formula = cbind(Low, Medium, High) ~ Infl, data = housing.wide)
#> 
#> Equation for Medium vs Low:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.5061     0.0971  -5.212 1.87e-07 ***
#> InflMedium    0.4200     0.1399   3.002  0.00268 ** 
#> InflHigh      0.6026     0.1832   3.288  0.00101 ** 
#> 
#> Equation for High vs Low:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.47712    0.09623  -4.958 7.12e-07 ***
#> InflMedium   0.72519    0.13380   5.420 5.96e-08 ***
#> InflHigh     1.54140    0.16213   9.507  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Approximate residual Deviance: 111.1 
#> Number of Fisher scoring iterations:  4 
#> Number of observations:  24 
#> 
dispersion(house.mblogit.wide.0, method="Afroz")
#> [1] 2.567212

anova(house.mblogit.wide.0,
      house.mblogit.wide)
#> Analysis of Deviance Table
#> 
#> Model 1: cbind(Low, Medium, High) ~ Infl
#> Model 2: cbind(Low, Medium, High) ~ Infl + Type + Cont
#>   Resid. Df Resid. Dev Df Deviance
#> 1        42    111.085            
#> 2        34     38.662  8   72.422
```
