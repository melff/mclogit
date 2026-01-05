# The IWLS algorithm used to fit conditional logit models

The package “mclogit” fits conditional logit models using a maximum
likelihood estimator. It does this by maximizing the log-likelihood
function using an *iterative weighted least-squares* (IWLS) algorithm,
which follows the algorithm used by the
[`glm.fit()`](https://rdrr.io/r/stats/glm.html) function from the
“stats” package of *R* (Nelder and Wedderburn 1972; McCullagh and Nelder
1989; R Core Team 2023).

If \pi\_{ij} is the probability that individual i chooses alternative j
from his/her choice set \mathcal{S}\_i, where

\pi\_{ij}=\frac{\exp(\eta\_{ij})}{\sum\_{k\in\mathcal{S}\_i}\exp(\eta\_{ik})}

and if y\_{ij} is the dummy variable with equals 1 if individual i
chooses alternative j and equals 0 otherwise, the log-likelihood
function (given that the choices are identically independent distributed
given \pi\_{ij}) can be written as

\ell=\sum\_{i,j}y\_{ij}\ln\pi\_{ij}
=\sum\_{i,j}y\_{ij}\eta\_{ij}-\sum_i\ln\left(\sum_j\exp(\eta\_{ij})\right)

If the data are aggregated in the terms of counts such that n\_{ij} is
the number of individuals with the same choice set and the same choice
probabilities \pi\_{ij} that have chosen alternative j, the
log-likelihood is (given that the choices are identically independent
distributed given \pi\_{ij})

\ell=\sum\_{i,j}n\_{ij}\ln\pi\_{ij}
=\sum\_{i,j}n\_{ij}\eta\_{ij}-\sum_in\_{i+}\ln\left(\sum_j\exp(\eta\_{ij})\right)

where n\_{i+}=\sum\_{j\in\mathcal{S}\_i}n\_{ij}.

If

\eta\_{ij} =
\alpha_1x\_{1ij}+\cdots+\alpha_rx\_{rij}=\boldsymbol{x}\_{ij}'\boldsymbol{\alpha}

then the gradient of the log-likelihood with respect to the coefficient
vector \boldsymbol{\alpha} is

\frac{\partial\ell}{\partial\boldsymbol{\alpha}} = \sum\_{i,j}
\frac{\partial\eta\_{ij}}{\partial\boldsymbol{\alpha}}
\frac{\partial\ell}{\partial\eta\_{ij}} = \sum\_{i,j}
\boldsymbol{x}\_{ij} (n\_{ij}-n\_{i+}\pi\_{ij}) = \sum\_{i,j}
\boldsymbol{x}\_{ij} n\_{i+} (y\_{ij}-\pi\_{ij}) =
\boldsymbol{X}'\boldsymbol{N}(\boldsymbol{y}-\boldsymbol{\pi})

and the Hessian is

\frac{\partial^2\ell}{\partial\boldsymbol{\alpha}\partial\boldsymbol{\alpha}'}
= \sum\_{i,j} \frac{\partial\eta\_{ij}}{\partial\boldsymbol{\alpha}}
\frac{\partial\eta\_{ij}}{\partial\boldsymbol{\alpha}'}
\frac{\partial\ell^2}{\partial\eta\_{ij}^2} = - \sum\_{i,j,k}
\boldsymbol{x}\_{ij} n\_{i+} (\delta\_{jk}-\pi\_{ij}\pi\_{ik})
\boldsymbol{x}\_{ij}' = - \boldsymbol{X}'\boldsymbol{W}\boldsymbol{X}

Here y\_{ij}=n\_{ij}/n\_{i+}, while \boldsymbol{N} is a diagonal matrix
with diagonal elements n\_{i+}.

Newton-Raphson iterations then take the form

\boldsymbol{\alpha}^{(s+1)} = \boldsymbol{\alpha}^{(s)} - \left(
\frac{\partial^2\ell}{\partial\boldsymbol{\alpha}\partial\boldsymbol{\alpha}'}
\right)^{-1} \frac{\partial\ell}{\partial\boldsymbol{\alpha}} =
\boldsymbol{\alpha}^{(s)} + \left(
\boldsymbol{X}'\boldsymbol{W}\boldsymbol{X} \right)^{-1}
\boldsymbol{X}'\boldsymbol{N}(\boldsymbol{y}-\boldsymbol{\pi})

where \boldsymbol{\pi} and \boldsymbol{W} are evaluated at
\boldsymbol{\alpha}=\boldsymbol{\alpha}^{(s)}.

Multiplying by \boldsymbol{X}'\boldsymbol{W}\boldsymbol{X} gives

\boldsymbol{X}'\boldsymbol{W}\boldsymbol{X} \boldsymbol{\alpha}^{(s+1)}
= \boldsymbol{X}'\boldsymbol{W}\boldsymbol{X}
\boldsymbol{\alpha}^{(s)} +
\boldsymbol{X}'\boldsymbol{N}(\boldsymbol{y}-\boldsymbol{\pi}) =
\boldsymbol{X}'\boldsymbol{W}
\left(\boldsymbol{X}\boldsymbol{\alpha}^{(s)}+\boldsymbol{W}^-\boldsymbol{N}(\boldsymbol{y}-\boldsymbol{\pi})\right)
= \boldsymbol{X}'\boldsymbol{W}\boldsymbol{y}^\*

where \boldsymbol{W}^- is a generalized inverse of \boldsymbol{W} and
\boldsymbol{y}^\* is a “working response vector” with elements

y\_{ij}^\*=\boldsymbol{x}\_{ij}'\boldsymbol{\alpha}^{(s)}+\frac{y\_{ij}-\pi\_{ij}}{\pi\_{ij}}

The IWLS algorithm thus involves the following steps:

1.  Create some suitable starting values for \boldsymbol{\pi},
    \boldsymbol{W}, and \boldsymbol{y}^\*

2.  Construct the “working dependent variable” \boldsymbol{y}^\*

3.  Solve the equation

    \boldsymbol{X}'\boldsymbol{W}\boldsymbol{X} \boldsymbol{\alpha} =
    \boldsymbol{X}'\boldsymbol{W}\boldsymbol{y}^\*

    for \boldsymbol{\alpha}.

4.  Compute updated \boldsymbol{\eta}, \boldsymbol{\pi}, \boldsymbol{W},
    and \boldsymbol{y}^\*.

5.  Compute the updated value for the log-likelihood or the deviance

    d=2\sum\_{i,j}n\_{ij}\ln\frac{y\_{ij}}{\pi\_{ij}}

6.  If the decrease of the deviance (or the increase of the
    log-likelihood) is smaller than a given tolerance criterian
    (typically \Delta d \leq 10^{-7}) stop the algorighm and declare it
    as converged. Otherwise go back to step 2 with the updated value of
    \boldsymbol{\alpha}.

The starting values for the algorithm used by the *mclogit* package are
constructe as follows:

1.  Set

    \eta\_{ij}^{(0)} = \ln (n\_{ij}+\tfrac12) -
    \frac1{q_i}\sum\_{k\in\mathcal{S}\_i}\ln (n\_{ij}+\tfrac12)

    (where q_i is the size of the choice set \mathcal{S}\_i)

2.  Compute the starting values of the choice probabilities
    \pi\_{ij}^{(0)} according to the equation at the beginning of the
    page

3.  Compute intial values of the working dependent variable according to

    y\_{ij}^{\*(0)} =
    \eta\_{ij}^{(0)}+\frac{y\_{ij}-\pi\_{ij}^{(0)}}{\pi\_{ij}^{(0)}}

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
