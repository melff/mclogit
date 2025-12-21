# Class, Party Position, and Electoral Choice

This is an artificial data set on electoral choice as influenced by
class and party positions.

## Usage

``` r
data(electors)
```

## Format

A data frame containing the following variables:

- class:

  class position of voters

- party:

  party that runs for election

- Freq:

  freqency by which each party list is chosen by members of each class

- time:

  time variable, runs from zero to one

- econ.left:

  economic-policy "leftness" of each party

- welfare:

  emphasis of welfare expansion of each party

- auth:

  position on authoritarian issues

## Examples

``` r
data(electors)

summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left+welfare+auth,
  data=electors))
#> 
#> Iteration 1 - deviance = 85051.49 - criterion = 0.9989204
#> Iteration 2 - deviance = 76759.94 - criterion = 0.108019
#> Iteration 3 - deviance = 74896.56 - criterion = 0.02487934
#> Iteration 4 - deviance = 74890.9 - criterion = 7.559543e-05
#> Iteration 5 - deviance = 74890.9 - criterion = 1.726814e-09
#> converged
#> 
#> Call:
#> mclogit(formula = cbind(Freq, interaction(time, class)) ~ econ.left + 
#>     welfare + auth, data = electors)
#> 
#>            Estimate Std. Error z value Pr(>|z|)    
#> econ.left -0.507265   0.007495 -67.679  < 2e-16 ***
#> welfare    0.564650   0.010700  52.769  < 2e-16 ***
#> auth       0.030305   0.005749   5.271 1.36e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Null Deviance:     80580 
#> Residual Deviance: 74890 
#> Number of Fisher Scoring iterations:  5 
#> Number of observations:  37500 
#> 
#> 

summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left/class+welfare/class+auth/class,
  data=electors))
#> 
#> Iteration 1 - deviance = 7377.939 - criterion = 0.9875551
#> Iteration 2 - deviance = 4589.544 - criterion = 0.6075407
#> Iteration 3 - deviance = 4293.485 - criterion = 0.06895374
#> Iteration 4 - deviance = 4277.887 - criterion = 0.00364612
#> Iteration 5 - deviance = 4277.808 - criterion = 1.852771e-05
#> Iteration 6 - deviance = 4277.808 - criterion = 5.890774e-10
#> converged
#> 
#> Call:
#> mclogit(formula = cbind(Freq, interaction(time, class)) ~ econ.left/class + 
#>     welfare/class + auth/class, data = electors)
#> 
#>                           Estimate Std. Error z value Pr(>|z|)    
#> econ.left                 -0.77851    0.02312 -33.671  < 2e-16 ***
#> welfare                    3.43776    0.03170 108.431  < 2e-16 ***
#> auth                      -0.13740    0.03608  -3.808  0.00014 ***
#> econ.left:classnew.middle  0.44546    0.02588  17.212  < 2e-16 ***
#> econ.left:classold.middle -0.44082    0.10387  -4.244  2.2e-05 ***
#> classnew.middle:welfare   -3.12917    0.03696 -84.659  < 2e-16 ***
#> classold.middle:welfare   -5.27438    0.07286 -72.393  < 2e-16 ***
#> classnew.middle:auth      -0.86676    0.03947 -21.957  < 2e-16 ***
#> classold.middle:auth       1.39435    0.05615  24.831  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Null Deviance:     80580 
#> Residual Deviance: 4278 
#> Number of Fisher Scoring iterations:  6 
#> Number of observations:  37500 
#> 
#> 

if (FALSE) # This takes a bit longer.
summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left/class+welfare/class+auth/class,
  random=~1|party.time,
  data=within(electors,party.time<-interaction(party,time))))

summary(mclogit(
  cbind(Freq,interaction(time,class))~econ.left/(class*time)+welfare/class+auth/class,
  random=~1|party.time,
  data=within(electors,{
        party.time <-interaction(party,time)
        econ.left.sq <- (econ.left-mean(econ.left))^2
        })))
#> 
#> Iteration 1 - deviance = 1071.031 - criterion = 0.1597241
#> Iteration 2 - deviance = 965.6196 - criterion = 0.02540274
#> Iteration 3 - deviance = 948.8356 - criterion = 0.005154656
#> Iteration 4 - deviance = 947.6262 - criterion = 0.000205486
#> Iteration 5 - deviance = 947.5081 - criterion = 2.557134e-07
#> Iteration 6 - deviance = 947.5042 - criterion = 5.256386e-13
#> converged
#> 
#> Call:
#> mclogit(formula = cbind(Freq, interaction(time, class)) ~ econ.left/(class * 
#>     time) + welfare/class + auth/class, data = within(electors, 
#>     {
#>         party.time <- interaction(party, time)
#>         econ.left.sq <- (econ.left - mean(econ.left))^2
#>     }), random = ~1 | party.time)
#> 
#> Coefficents:
#>                                Estimate Std. Error z value Pr(>|z|)    
#> econ.left                      -0.13335    0.20837  -0.640    0.522    
#> welfare                         2.05552    0.21245   9.675   <2e-16 ***
#> auth                            0.08071    0.11717   0.689    0.491    
#> econ.left:classnew.middle      -1.69581    0.11631 -14.580   <2e-16 ***
#> econ.left:classold.middle      -3.04338    0.20351 -14.954   <2e-16 ***
#> econ.left:time                 -0.07782    0.30228  -0.257    0.797    
#> classnew.middle:welfare        -0.99267    0.06073 -16.346   <2e-16 ***
#> classold.middle:welfare        -1.62088    0.12850 -12.614   <2e-16 ***
#> classnew.middle:auth           -1.39056    0.04673 -29.754   <2e-16 ***
#> classold.middle:auth            1.45722    0.05814  25.063   <2e-16 ***
#> econ.left:classnew.middle:time  0.06049    0.14449   0.419    0.675    
#> econ.left:classold.middle:time  0.14723    0.26232   0.561    0.575    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Co-)Variances:
#> Grouping level: party.time 
#>          Estimate   Std.Err.
#>          (Const.)   (Const.)
#> (Const.)  1.604      0.3098 
#> 
#> Approximate residual deviance: 947.5 
#> Number of Fisher scoring iterations:  6
#> Number of observations
#>   Groups by party.time: 150
#>   Individual observations:  37500
#> 
 # \dontrun{}
```
