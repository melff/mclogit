# Approximate Inference for Multinomial Logit Models with Random Effects

## The problem

A crucial problem for inference about non-linear models with random
effects is that the likelihood function for such models involves
integrals for which no analytical solution exists.

For given values $\mathbf{b}$ of the random effects the likelihood
function of a conditional logit model (and therefore also of a
baseline-logit model) can be written in the form

$$\mathcal{L}_{\text{cpl}}(\mathbf{y},\mathbf{b}) = \exp\left( \ell_{\text{cpl}}(\mathbf{y},\mathbf{b}) \right) = \exp\left( \ell\left( \mathbf{y}|\mathbf{b};{\mathbf{α}} \right) - \frac{1}{2}\ln\det(\mathbf{\Sigma}) - \frac{1}{2}\mathbf{b}\prime\mathbf{\Sigma}^{- 1}\mathbf{b} \right)$$

However, this “complete data” likelihood function cannot be used for
inference, because it depends on the unobserved random effects. To
arrive at a likelihood function that depends only on observed data, one
needs to used the following integrated likelihood function:

$$\mathcal{L}_{\text{obs}}(\mathbf{y}) = \int\exp\left( \ell_{\text{cpl}}(\mathbf{y},\mathbf{b}) \right)\partial\mathbf{b} = \int\exp\left( \ell\left( \mathbf{y}|\mathbf{b};{\mathbf{α}} \right) - \frac{1}{2}\ln\det(\mathbf{\Sigma}) - \frac{1}{2}\mathbf{b}\prime\mathbf{\Sigma}^{- 1}\mathbf{b} \right)\partial\mathbf{b}$$

In general, this integral cannot be “solved”, i.e. eliminated from the
formula by analytic means (it is “analytically untractable”). Instead,
one will compute it either using numeric techniques (e.g. using
numerical quadrature) or approximate it using analytical techniques.
Unless there is only a single level of random effects numerical
quadrature can become computationally be demanding, that is, the
computation of the (log-)likelihood function and its derivatives can
take a lot of time even on modern, state-of-the-art computer hardware.
Yet approximations based on analytical techniques hand may lead to
biased estimates in particular in samples where the number of
observations relative to the number of random offects is small, but at
least they are much easier to compute and sometimes making inference
possible after all.

The package “mclogit” supports to kinds of analytical approximations,
the Laplace approximation and what one may call the Solomon-Cox
appoximation. Both approximations are based on a quadratic expansion of
the integrand so that the thus modified integral does have a closed-form
solution, i.e. is analytically tractable.

## The Laplace approximation and PQL

### Laplace approximation

The (first-order) Laplace approximation is based on the quadratic
expansion the logarithm of the integrand, the complete-data
log-likelihood

$$\ell_{\text{cpl}}(\mathbf{y},\mathbf{b}) \approx \ell\left( \mathbf{y}|\widetilde{\mathbf{b}};{\mathbf{α}} \right) - \frac{1}{2}\left( \mathbf{b} - \widetilde{\mathbf{b}} \right)\prime\widetilde{\mathbf{H}}\left( \mathbf{b} - \widetilde{\mathbf{b}} \right) - \frac{1}{2}\ln\det(\mathbf{\Sigma}) - \frac{1}{2}\left( \mathbf{b} - \widetilde{\mathbf{b}} \right)\prime\mathbf{\Sigma}^{- 1}\left( \mathbf{b} - \widetilde{\mathbf{b}} \right)$$

where $\widetilde{\mathbf{b}}$ is the solution to

$$\frac{\partial\ell_{\text{cpl}}(\mathbf{y},\mathbf{b})}{\partial\mathbf{b}} = 0$$

and
$\widetilde{\mathbf{H}} = \mathbf{H}\left( \widetilde{\mathbf{b}} \right)$
is the value of the negative Hessian with respect to $\mathbf{b}$

$$\mathbf{H}(\mathbf{b}) = - \frac{\partial^{2}\ell\left( \mathbf{y}|\mathbf{b};{\mathbf{α}} \right)}{\partial\mathbf{b}\partial\mathbf{b}\prime}$$

for $\mathbf{b} = \widetilde{\mathbf{b}}$.

Since this quadratic expansion—let us call it
$\ell_{\text{Lapl}}^{*}(\mathbf{y},\mathbf{b})$—is a (multivariate)
quadratic function of $\mathbf{b}$, the integral of its exponential does
have a closed-form solution (the relevant formula can be found in
Harville (1997)).

For purposes of estimation, the resulting approximate log-likelihood is
more useful:

$$\ell_{\text{Lapl}}^{*} = \ln\int\exp\left( \ell_{\text{Lapl}}(\mathbf{y},\mathbf{b}) \right)\partial\mathbf{b} = \ell\left( \mathbf{y}|\widetilde{\mathbf{b}};{\mathbf{α}} \right) - \frac{1}{2}\widetilde{\mathbf{b}}\prime\mathbf{\Sigma}^{- 1}\widetilde{\mathbf{b}} - \frac{1}{2}\ln\det(\mathbf{\Sigma}) - \frac{1}{2}\ln\det\left( \widetilde{\mathbf{H}} + \mathbf{\Sigma}^{- 1} \right).$$

### Penalized quasi-likelihood (PQL)

If one disregards the dependence of $\widetilde{\mathbf{H}}$ on
$\mathbf{α}$ and $\mathbf{b}$, then $\widetilde{\mathbf{b}}$ maximizes
not only $\ell_{\text{cpl}}(\mathbf{y},\mathbf{b})$ but also
$\ell_{\text{Lapl}}^{*}$. This motivates the following IWLS/Fisher
scoring equations for $\widehat{\mathbf{α}}$ and
$\widetilde{\mathbf{b}}$ (see Breslow and Clayton (1993) and [this
page](https://melff.github.io/mclogit/articles/fitting-mclogit.md)):

$$\begin{array}{r}
{\begin{bmatrix}
{\mathbf{X}\prime\mathbf{W}\mathbf{X}} & {\mathbf{X}\prime\mathbf{W}\mathbf{Z}} \\
{\mathbf{Z}\prime\mathbf{W}\mathbf{X}} & {\mathbf{Z}\prime\mathbf{W}\mathbf{Z} + \mathbf{\Sigma}^{- 1}} \\
 & 
\end{bmatrix}\begin{bmatrix}
\widehat{\mathbf{α}} \\
\widetilde{\mathbf{b}} \\

\end{bmatrix} = \begin{bmatrix}
{\mathbf{X}\prime\mathbf{W}\mathbf{y}^{*}} \\
{\mathbf{Z}\prime\mathbf{W}\mathbf{y}^{*}}
\end{bmatrix}}
\end{array}$$

where

$$\mathbf{y}^{*} = \mathbf{X}{\mathbf{α}} + \mathbf{Z}\mathbf{b} + \mathbf{W}^{-}(\mathbf{y} - {\mathbf{π}})$$

is the IWLS “working dependend variable” with $\mathbf{α}$,
$\mathbf{b}$, $\mathbf{W}$, and $\mathbf{π}$ computed in an earlier
iteration.

Substitutions lead to the equations:

$$\left( \mathbf{X}\mathbf{V}^{-}\mathbf{X} \right)\widehat{\mathbf{α}} = \mathbf{X}\mathbf{V}^{-}\mathbf{y}^{*}$$

and

$$\left( \mathbf{Z}\prime\mathbf{W}\mathbf{Z} + \mathbf{\Sigma}^{- 1} \right)\mathbf{b} = \mathbf{Z}\prime\mathbf{W}\left( \mathbf{y}^{*} - \mathbf{X}{\mathbf{α}} \right)$$

which can be solved to compute $\widehat{\mathbf{α}}$ and
$\widetilde{\mathbf{b}}$ (for given $\mathbf{\Sigma}$)

Here

$$\mathbf{V} = \mathbf{W}^{-} + \mathbf{Z}\mathbf{\Sigma}\mathbf{Z}\prime$$

and

$$\mathbf{V}^{-} = \mathbf{W} - \mathbf{W}\mathbf{Z}\prime\left( \mathbf{Z}\prime\mathbf{W}\mathbf{Z} + \mathbf{\Sigma}^{- 1} \right)^{- 1}\mathbf{Z}\mathbf{W}$$

Following Breslow and Clayton (1993) the variance parameters in
$\mathbf{\Sigma}$ are estimated by minimizing

$$q_{1} = \det(\mathbf{V}) + \left( \mathbf{y}^{*} - \mathbf{X}{\mathbf{α}} \right)\mathbf{V}^{-}\left( \mathbf{y}^{*} - \mathbf{X}{\mathbf{α}} \right)$$

or the “REML” variant:

$$q_{2} = \det(\mathbf{V}) + \left( \mathbf{y}^{*} - \mathbf{X}{\mathbf{α}} \right)\mathbf{V}^{-}\left( \mathbf{y}^{*} - \mathbf{X}{\mathbf{α}} \right) + \det\left( \mathbf{X}\prime\mathbf{V}^{-}\mathbf{X} \right)$$

This motivates the following algorithm, which is strongly inspired by
the [`glmmPQL()`](https://rdrr.io/pkg/MASS/man/glmmPQL.html) function in
Brian Ripley’s *R* package
[MASS](https://cran.r-project.org/package=MASS) (Venables and Ripley
2002):

1.  Create some suitable starting values for $\mathbf{π}$, $\mathbf{W}$,
    and $\mathbf{y}^{*}$
2.  Construct the “working dependent variable” $\mathbf{y}^{*}$
3.  Minimize $q_{1}$ (quasi-ML) or $q_{2}$ (quasi-REML) iteratively
    (inner loop), to obtain an estimate of $\mathbf{\Sigma}$
4.  Obtain $hat{\mathbf{α}}$ and $\widetilde{\mathbf{b}}$ based on the
    current estimate of $\mathbf{\Sigma}$
5.  Compute updated
    ${\mathbf{η}} = \mathbf{X}{\mathbf{α}} + \mathbf{Z}\mathbf{b}$,
    $\mathbf{π}$, $\mathbf{W}$.
6.  If the change in $\mathbf{η}$ is smaller than a given tolerance
    criterion stop the algorighm and declare it as converged. Otherwise
    go back to step 2 with the updated values of $\widehat{\mathbf{α}}$
    and $\widetilde{\mathbf{b}}$.

This algorithm is a modification of the
[IWLS](https://melff.github.io/mclogit/articles/fitting-mclogit.md)
algorithm used to fit conditional logit models without random effects.
Instead of just solving a linear requatoin in step 3, it estimates a
weighted linear mixed-effects model. In contrast to
[`glmmPQL()`](https://rdrr.io/pkg/MASS/man/glmmPQL.html) it does not use
the `lme()` function from package
[nlme](https://cran.r-project.org/package=nlme) (Pinheiro and Bates
2000) for this, because the weighting matrix $\mathbf{W}$ is
non-diagonal. Instead, $q_{1}$ or $q_{2}$ are minimized using the
function `nlminb` from the standard *R* package “stats” or some other
optimizer chosen by the user.

## The Solomon-Cox approximation and MQL

### The Solomon-Cox approximation

The (first-order) Solomon approximation (Solomon and Cox 1992) is based
on the quadratic expansion the integrand

$$\ell_{\text{cpl}}(\mathbf{y},\mathbf{b}) \approx \ell\left( \mathbf{y}|\mathbf{0};{\mathbf{α}} \right) + \mathbf{g}_{0}\prime\mathbf{b} - \frac{1}{2}\mathbf{b}\prime\mathbf{H}_{0}\mathbf{b} - \frac{1}{2}\ln\det(\mathbf{\Sigma}) - \frac{1}{2}\mathbf{b}\prime\mathbf{\Sigma}^{- 1}\mathbf{b}$$

where $\mathbf{g}\_ 0 = \mathbf{g}(\mathbf{0})$ is the gradient of
$\ell(\mathbf{y} \parallel \mathbf{b};{\mathbf{α}})$

$$\mathbf{g}(\mathbf{b}) = - \frac{\partial\ell\left( \mathbf{y}|\mathbf{b};{\mathbf{α}} \right)}{\partial\mathbf{b}}$$

at $\mathbf{b} = \mathbf{0}$, while
$\mathbf{H}\_ 0 = \mathbf{H}(\mathbf{0})$ is the negative Hessian at
$\mathbf{b} = \mathbf{0}$.

Like before, the integral of the exponential this quadratic expansion
(which we refer to as $\ell_{\text{SC}}(\mathbf{y},\mathbf{b})$) has a
closed-form solution, as does its logarithm, which is:

$$\ln\int\exp\left( \ell_{\text{SC}}(\mathbf{y},\mathbf{b}) \right)\partial\mathbf{b} = \ell\left( \mathbf{y}|\mathbf{0};{\mathbf{α}} \right) - \frac{1}{2}\mathbf{g}_{0}\prime\left( \mathbf{H}_{0} + \mathbf{\Sigma}^{- 1} \right)^{- 1}\mathbf{g}_{0} - \frac{1}{2}\ln\det(\mathbf{\Sigma}) - \frac{1}{2}\ln\det\left( \mathbf{H}_{0} + \mathbf{\Sigma}^{- 1} \right).$$

### Marginal quasi-likelhood (MQL)

The resulting estimation technique is very similar to PQL (again, see
Breslow and Clayton 1993 for a discussion). The only difference is the
construction of the “working dependent” variable $\mathbf{y}^{*}$. With
PQL it is constructed as
$$\mathbf{y}^{*} = \mathbf{X}{\mathbf{α}} + \mathbf{Z}\mathbf{b} + \mathbf{W}^{-}(\mathbf{y} - {\mathbf{π}})$$
while the MQL working dependent variable is just

$$\mathbf{y}^{*} = \mathbf{X}{\mathbf{α}} + \mathbf{W}^{-}(\mathbf{y} - {\mathbf{π}})$$

so that the algorithm has the following steps:

1.  Create some suitable starting values for $\mathbf{π}$, $\mathbf{W}$,
    and $\mathbf{y}^{*}$
2.  Construct the “working dependent variable” $\mathbf{y}^{*}$
3.  Minimize $q_{1}$ (quasi-ML) or $q_{2}$ (quasi-REML) iteratively
    (inner loop), to obtain an estimate of $\mathbf{\Sigma}$
4.  Obtain $\widehat{\mathbf{α}}$ based on the current estimate of
    $\mathbf{\Sigma}$
5.  Compute updated ${\mathbf{η}} = \mathbf{X}{\mathbf{α}}$,
    $\mathbf{π}$, $\mathbf{W}$.
6.  If the change in $\mathbf{η}$ is smaller than a given tolerance
    criterion stop the algorighm and declare it as converged. Otherwise
    go back to step 2 with the updated values of $\widehat{\mathbf{α}}$.

## References

Breslow, Norman E., and David G. Clayton. 1993. “Approximate Inference
in Generalized Linear Mixed Models.” *Journal of the American
Statistical Association* 88 (421): 9–25.

Harville, David A. 1997. *Matrix Algebra from a Statistician’s
Perspective*. New York: Springer.

Pinheiro, José C., and Douglas M. Bates. 2000. *Mixed-Effects Models in
s and s-PLUS*. New York: Springer. <https://doi.org/10.1007/b98882>.

Solomon, P. J., and D. R. Cox. 1992. “Nonlinear Component of Variance
Models.” *Biometrika* 79 (1): 1–11.
<https://doi.org/10.1093/biomet/79.1.1>.

Venables, W. N., and B. D. Ripley. 2002. *Modern Applied Statistics with
s*. Fourth. New York: Springer. <https://www.stats.ox.ac.uk/pub/MASS4/>.
