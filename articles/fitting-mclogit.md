# The IWLS algorithm used to fit conditional logit models

The package “mclogit” fits conditional logit models using a maximum
likelihood estimator. It does this by maximizing the log-likelihood
function using an *iterative weighted least-squares* (IWLS) algorithm,
which follows the algorithm used by the
[`glm.fit()`](https://rdrr.io/r/stats/glm.html) function from the
“stats” package of *R* (Nelder and Wedderburn 1972; McCullagh and Nelder
1989; R Core Team 2023).

If $\pi_{ij}$ is the probability that individual $i$ chooses alternative
$j$ from his/her choice set $\mathcal{S}_{i}$, where

$$\pi_{ij} = \frac{\exp\left( \eta_{ij} \right)}{\sum\limits_{k \in \mathcal{S}_{i}}\exp\left( \eta_{ik} \right)}$$

and if $y_{ij}$ is the dummy variable with equals 1 if individual $i$
chooses alternative $j$ and equals 0 otherwise, the log-likelihood
function (given that the choices are identically independent distributed
given $\pi_{ij}$) can be written as

$$\ell = \sum\limits_{i,j}y_{ij}\ln\pi_{ij} = \sum\limits_{i,j}y_{ij}\eta_{ij} - \sum\limits_{i}\ln\left( \sum\limits_{j}\exp\left( \eta_{ij} \right) \right)$$

If the data are aggregated in the terms of counts such that $n_{ij}$ is
the number of individuals with the same choice set and the same choice
probabilities $\pi_{ij}$ that have chosen alternative $j$, the
log-likelihood is (given that the choices are identically independent
distributed given $\pi_{ij}$)

$$\ell = \sum\limits_{i,j}n_{ij}\ln\pi_{ij} = \sum\limits_{i,j}n_{ij}\eta_{ij} - \sum\limits_{i}n_{i +}\ln\left( \sum\limits_{j}\exp\left( \eta_{ij} \right) \right)$$

where $n_{i +} = \sum_{j \in \mathcal{S}_{i}}n_{ij}$.

If

$$\eta_{ij} = \alpha_{1}x_{1ij} + \cdots + \alpha_{r}x_{rij} = \mathbf{x}_{ij}\prime{\mathbf{α}}$$

then the gradient of the log-likelihood with respect to the coefficient
vector $\mathbf{α}$ is

$$\frac{\partial\ell}{\partial{\mathbf{α}}} = \sum\limits_{i,j}\frac{\partial\eta_{ij}}{\partial{\mathbf{α}}}\frac{\partial\ell}{\partial\eta_{ij}} = \sum\limits_{i,j}\mathbf{x}_{ij}\left( n_{ij} - n_{i +}\pi_{ij} \right) = \sum\limits_{i,j}\mathbf{x}_{ij}n_{i +}\left( y_{ij} - \pi_{ij} \right) = \mathbf{X}\prime\mathbf{N}(\mathbf{y} - {\mathbf{π}})$$

and the Hessian is

$$\frac{\partial^{2}\ell}{\partial{\mathbf{α}}\partial{\mathbf{α}}\prime} = \sum\limits_{i,j}\frac{\partial\eta_{ij}}{\partial{\mathbf{α}}}\frac{\partial\eta_{ij}}{\partial{\mathbf{α}}\prime}\frac{\partial\ell^{2}}{\partial\eta_{ij}^{2}} = - \sum\limits_{i,j,k}\mathbf{x}_{ij}n_{i +}\left( \delta_{jk} - \pi_{ij}\pi_{ik} \right)\mathbf{x}_{ij}\prime = - \mathbf{X}\prime\mathbf{W}\mathbf{X}$$

Here $y_{ij} = n_{ij}/n_{i +}$, while $\mathbf{N}$ is a diagonal matrix
with diagonal elements $n_{i +}$.

Newton-Raphson iterations then take the form

$${\mathbf{α}}^{(s + 1)} = {\mathbf{α}}^{(s)} - \left( \frac{\partial^{2}\ell}{\partial{\mathbf{α}}\partial{\mathbf{α}}\prime} \right)^{- 1}\frac{\partial\ell}{\partial{\mathbf{α}}} = {\mathbf{α}}^{(s)} + (\mathbf{X}\prime\mathbf{W}\mathbf{X})^{- 1}\mathbf{X}\prime\mathbf{N}(\mathbf{y} - {\mathbf{π}})$$

where $\mathbf{π}$ and $\mathbf{W}$ are evaluated at
${\mathbf{α}} = {\mathbf{α}}^{(s)}$.

Multiplying by $\mathbf{X}\prime\mathbf{W}\mathbf{X}$ gives

$$\mathbf{X}\prime\mathbf{W}\mathbf{X}{\mathbf{α}}^{(s + 1)} = \mathbf{X}\prime\mathbf{W}\mathbf{X}{\mathbf{α}}^{(s)} + \mathbf{X}\prime\mathbf{N}(\mathbf{y} - {\mathbf{π}}) = \mathbf{X}\prime\mathbf{W}\left( \mathbf{X}{\mathbf{α}}^{(s)} + \mathbf{W}^{-}\mathbf{N}(\mathbf{y} - {\mathbf{π}}) \right) = \mathbf{X}\prime\mathbf{W}\mathbf{y}^{*}$$

where $\mathbf{W}^{-}$ is a generalized inverse of $\mathbf{W}$ and
$\mathbf{y}^{*}$ is a “working response vector” with elements

$$y_{ij}^{*} = \mathbf{x}_{ij}\prime{\mathbf{α}}^{(s)} + \frac{y_{ij} - \pi_{ij}}{\pi_{ij}}$$

The IWLS algorithm thus involves the following steps:

1.  Create some suitable starting values for $\mathbf{π}$, $\mathbf{W}$,
    and $\mathbf{y}^{*}$

2.  Construct the “working dependent variable” $\mathbf{y}^{*}$

3.  Solve the equation

    $$\mathbf{X}\prime\mathbf{W}\mathbf{X}{\mathbf{α}} = \mathbf{X}\prime\mathbf{W}\mathbf{y}^{*}$$

    for $\mathbf{α}$.

4.  Compute updated $\mathbf{η}$, $\mathbf{π}$, $\mathbf{W}$, and
    $\mathbf{y}^{*}$.

5.  Compute the updated value for the log-likelihood or the deviance

    $$d = 2\sum\limits_{i,j}n_{ij}\ln\frac{y_{ij}}{\pi_{ij}}$$

6.  If the decrease of the deviance (or the increase of the
    log-likelihood) is smaller than a given tolerance criterian
    (typically $\Delta d \leq 10^{- 7}$) stop the algorighm and declare
    it as converged. Otherwise go back to step 2 with the updated value
    of $\mathbf{α}$.

The starting values for the algorithm used by the *mclogit* package are
constructe as follows:

1.  Set

    $$\eta_{ij}^{(0)} = \ln\left( n_{ij} + \frac{1}{2} \right) - \frac{1}{q_{i}}\sum\limits_{k \in \mathcal{S}_{i}}\ln\left( n_{ij} + \frac{1}{2} \right)$$

    (where $q_{i}$ is the size of the choice set $\mathcal{S}_{i}$)

2.  Compute the starting values of the choice probabilities
    $\pi_{ij}^{(0)}$ according to the equation at the beginning of the
    page

3.  Compute intial values of the working dependent variable according to

    $$y_{ij}^{*{(0)}} = \eta_{ij}^{(0)} + \frac{y_{ij} - \pi_{ij}^{(0)}}{\pi_{ij}^{(0)}}$$

## References

McCullagh, P., and J. A. Nelder. 1989. *Generalized Linear Models*.
Monographs on Statistics & Applied Probability. Boca Raton et al.:
Chapman & Hall/CRC.

Nelder, J. A., and R. W. M. Wedderburn. 1972. “Generalized Linear
Models.” *Journal of the Royal Statistical Society. Series A (General)*
135 (3): 370–84. <https://doi.org/10.2307/2344614>.

R Core Team. 2023. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.
