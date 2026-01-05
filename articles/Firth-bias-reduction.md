# Bias reduction using Firth's penalized likelihood technique

Firth’s penalized likelihood technique (Firth 1993) has two advantages
in comparison to the conventional maximum likelihood estimator:

1.  It has a smaller asymptotic bias than the MLE, i.e. it converges
    “faster” to its limit as the sample size approaches infinity: In
    somewhat imprecise terms, its bias is $O\left( n^{- 2} \right)$
    instead of $O\left( n^{- 1} \right)$, where $n$ is the sample size.

2.  In case of logistic regression and other models for categorical
    responses, tends to be less afflicted by the problem of complete
    separation or quasi-complete separation. That is, Firth’s penalized
    likelihood technique gives finite coefficient estimates in cases
    where a (finite) MLE for the coefficients does not exist.,

In case of the models supported by the *mclogit* package, the difference
between Firth’s PML technique and conventional ML estimation is that it
solves for each coefficient $\alpha_{q}$ the modified gradient equation
$$0 = U^{*}\left( \alpha_{q} \right) = \frac{\partial\ell}{\partial\alpha_{q}} + \frac{1}{2}{tr}\left( (\mathbf{X}\prime\mathbf{W}\mathbf{X})^{- 1}\mathbf{X}\prime\frac{\partial\mathbf{W}}{\partial\alpha_{q}}\mathbf{X} \right)$$
instead of
$$0 = U\left( \alpha_{q} \right) = \frac{\partial\ell}{\partial\alpha_{q}}$$
for the MLE, where $\mathbf{X}\prime\mathbf{W}\mathbf{X}$ is the
information matrix for the full coefficient vector. This technique is
equivalent to maximizing the penalized (log-)likelihood
$$\mathcal{L}^{*}({\mathbf{α}}) = \mathcal{L}({\mathbf{α}})\frac{1}{2}\det(\mathbf{X}\prime\mathbf{W}\mathbf{X})\quad\text{or}\quad\ell^{*}({\mathbf{α}}) = \ell({\mathbf{α}}) + \frac{1}{2}\ln\det(\mathbf{X}\prime\mathbf{W}\mathbf{X}).$$

Like the conventional MLE, the *mclogit* package uses an IWLS-algorithm,
however, using a modified working response vector with elements
$$y_{ij}^{*} = \mathbf{x}_{ij}\prime{\mathbf{α}} + \frac{y_{ij} - \pi_{ij}}{\pi_{ij}} + \frac{1}{2}\sum\limits_{r}\sum\limits_{s}\zeta_{ij,rs}I^{(r,s)}$$
where $I^{(r,s)}$ is the $r,s$ element of
$(\mathbf{X}\prime\mathbf{W}\mathbf{X})^{- 1}$ and $$\begin{aligned}
\zeta_{ij,rs} & {= \left( x_{ijr} - \sum\limits_{k}x_{ikr}\pi_{ik} \right)\left( x_{ijs} - \sum\limits_{l}\pi_{il}x_{ils} \right) - \sum\limits_{k}\pi_{ik}\left( x_{ikr} - \sum\limits_{l}x_{ilr}\pi_{il} \right)\left( x_{iks} - \sum\limits_{l}\pi_{il}x_{ils} \right)} \\
 & {= \Delta_{i,j,{\mathbf{π}}_{i}}\left( \Delta_{i,j,{\mathbf{π}}_{i}}x_{ijr}\Delta_{i,j,{\mathbf{π}}_{i}}x_{ijs} \right)}
\end{aligned}$$ where
$\Delta_{i,j,{\mathbf{π}}_{i}}a_{ij} = a_{ij} - \sum_{k}\pi_{ik}a_{ik}.$

### A comparison of results obtained with other packages

For a multinomial baseline logit model for a two-category response,
`mblogit(...,Firth=TRUE)` gives the identical result as function
[`brglm()`](https://rdrr.io/pkg/brglm/man/brglm.html) from Ioannis
Kosmidis’ package *brglm* (Kosmidis 2025):

``` r
if(require("brglm", quietly = TRUE)) {
    suppressMessages(library(mclogit))
    data(lizards)

  lizards.brglm <- brglm(cbind(grahami, opalinus) ~ height + diameter +
                         light + time, family = binomial(logit), data=lizards,
                         method = "brglm.fit")

    print(summary(lizards.brglm))

    lizards.mblogit_F <- mblogit(cbind(opalinus, grahami) ~ height + diameter + light + time,
                                 data = lizards, Firth = TRUE)
    print(summary(lizards.mblogit_F))
}
```

    ## 'brglm' will gradually be superseded by the 'brglm2' R package (https://cran.r-project.org/package=brglm2), which provides utilities for mean and median bias reduction for all GLMs.
    ##  Methods for the detection of separation and infinite estimates in binomial-response models are provided by the 'detectseparation' R package (https://cran.r-project.org/package=detectseparation).

    ## 
    ## Call:
    ## brglm(formula = cbind(grahami, opalinus) ~ height + diameter + 
    ##     light + time, family = binomial(logit), data = lizards, method = "brglm.fit")
    ## 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    1.9018     0.3374   5.638 1.73e-08 ***
    ## height>=5ft    1.1064     0.2544   4.349 1.37e-05 ***
    ## diameter>2in  -0.7536     0.2103  -3.584 0.000338 ***
    ## lightshady    -0.8177     0.3186  -2.566 0.010277 *  
    ## timemidday     0.2280     0.2488   0.916 0.359621    
    ## timelate      -0.7273     0.2975  -2.445 0.014482 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 68.132  on 22  degrees of freedom
    ## Residual deviance: 14.246  on 17  degrees of freedom
    ## Penalized deviance: -4.00651 
    ## AIC:  83.07 
    ## 
    ## 
    ## Iteration 1 - deviance = 14.43761 - criterion = 0.7002064
    ## Iteration 2 - deviance = 14.24815 - criterion = 0.01320475
    ## Iteration 3 - deviance = 14.24623 - criterion = 0.0001332939
    ## Iteration 4 - deviance = 14.24623 - criterion = 9.080537e-08
    ## Iteration 5 - deviance = 14.24623 - criterion = 2.278489e-10
    ## converged
    ## 
    ## Call:
    ## mblogit(formula = cbind(opalinus, grahami) ~ height + diameter + 
    ##     light + time, data = lizards, Firth = TRUE)
    ## 
    ## Equation for grahami vs opalinus:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    1.9018     0.3374   5.638 1.73e-08 ***
    ## height>=5ft    1.1064     0.2544   4.349 1.37e-05 ***
    ## diameter>2in  -0.7536     0.2103  -3.584 0.000338 ***
    ## lightshady    -0.8177     0.3186  -2.566 0.010277 *  
    ## timemidday     0.2280     0.2488   0.916 0.359621    
    ## timelate      -0.7273     0.2975  -2.445 0.014482 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate residual Deviance: 14.25 
    ## Number of Fisher scoring iterations:  5 
    ## Number of observations:  23

``` r
library(MASS) #For the housing data
print(house.mblogit <- mblogit(Sat ~ Infl + Type + Cont, weights = Freq,
                          data = housing,
                          Firth = TRUE))
```

    ## 
    ## Iteration 1 - deviance = 3494.982 - criterion = 0.9614604
    ## Iteration 2 - deviance = 3470.128 - criterion = 0.007162046
    ## Iteration 3 - deviance = 3470.092 - criterion = 1.040578e-05
    ## Iteration 4 - deviance = 3470.092 - criterion = 4.681856e-09
    ## converged
    ## 
    ## Call: mblogit(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq, 
    ##     Firth = TRUE)
    ## 
    ## Coefficients:
    ##             Predictors
    ## Logit eqn.    (Intercept)  InflMedium  InflHigh  TypeApartment  TypeAtrium
    ##   Medium/Low  -0.4169       0.4441      0.6615   -0.4338         0.1300   
    ##   High/Low    -0.1385       0.7314      1.6033   -0.7314        -0.4067   
    ##             Predictors
    ## Logit eqn.    TypeTerrace  ContHigh
    ##   Medium/Low  -0.6620       0.3587 
    ##   High/Low    -1.4036       0.4792 
    ## 
    ## Null Deviance:     3694 
    ## Residual Deviance: 3470

``` r
if(require("brglm2", quietly = TRUE)){
    print(brmultinom(Sat ~ Infl + Type + Cont, weights = Freq,
                   data = housing, ref = 1,
                   type = "AS_mean"))
}
```

    ## Call:
    ## brmultinom(formula = Sat ~ Infl + Type + Cont, data = housing, 
    ##     weights = Freq, ref = 1, type = "AS_mean")
    ## 
    ## Coefficients:
    ##         (Intercept)  InflMedium  InflHigh  TypeApartment  TypeAtrium
    ## Medium  -0.41687      0.44413     0.66149  -0.43385        0.13005  
    ## High    -0.13851      0.73136     1.60325  -0.73136       -0.40671  
    ##         TypeTerrace  ContHigh
    ## Medium  -0.66200      0.35873
    ## High    -1.40361      0.47923
    ## 
    ## Residual Deviance: 3470.092

## References

Firth, David. 1993. “Bias Reduction of Maximum Likelihood Estimates.”
*Biometrika* 80 (1): 27–38. <https://doi.org/10.1093/biomet/80.1.27>.

Kosmidis, Ioannis. 2025. *brglm: Bias Reduction in Binary-Response
Generalized Linear Models*. <https://cran.r-project.org/package=brglm>.
