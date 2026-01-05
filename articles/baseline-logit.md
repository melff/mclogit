# Baseline-category logit models

Multinomial baseline-category logit models are a generalisation of
logistic regression, that allow to model not only binary or dichotomous
responses, but also polychotomous responses. In addition, they allow to
model responses in the form of counts that have a pre-determined sum.
These models are described in Agresti (2002). Estimating these models is
also supported by the function
[`multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) in the *R*
package “nnet” (Venables and Ripley 2002). In the package “mclogit”, the
function to estimate these models is called
[`mblogit()`](https://melff.github.io/mclogit/reference/mblogit.md),
which uses the infrastructure for estimating conditional logit models,
exploiting the fact that baseline-category logit models can be
re-expressed as condigional logit models.

Baseline-category logit models are constructed as follows. Suppose a
categorical dependent variable or response with categories j=1,\ldots,q
is observed for individuals i=1,\ldots,n. Let \pi\_{ij} denote the
probability that the value of the dependent variable for individual i is
equal to j, then the baseline-category logit model takes the form:

\begin{aligned} \pi\_{ij} = \begin{cases}
\dfrac{\exp(\alpha\_{j0}+\alpha\_{j1}x\_{1i}+\cdots+\alpha\_{jr}x\_{ri})}
{1+\sum\_{k\>1}\exp(\alpha\_{k0}+\alpha\_{k1}x\_{1i}+\cdots+\alpha\_{kr}x\_{ri})}
& \text{for } j\>1\\\[20pt\] \dfrac{1}
{1+\sum\_{k\>1}\exp(\alpha\_{k0}+\alpha\_{k1}x\_{1i}+\cdots+\alpha\_{kr}x\_{ri})}
& \text{for } j=1 \end{cases} \end{aligned}

where the first category (j=1) is the baseline category.

Equivalently, the model can be expressed in terms of log-odds, relative
to the baseline-category:

\ln\frac{\pi\_{ij}}{\pi\_{i1}} =
\alpha\_{j0}+\alpha\_{j1}x\_{1i}+\cdots+\alpha\_{jr}x\_{ri}.

Here the relevant parameters of the model are the coefficients
\alpha\_{jk} which describe how the values of independent variables
(numbered k=1,\ldots,r) affect the relative chances of the response
taking a value j versus taking the value 1. Note that there is one
coefficient for each independent variable and *each response* other than
the baseline category.

## References

Agresti, Alan. 2002. *Categorical Data Analysis*. Second. New York:
Wiley.

Venables, W. N., and B. D. Ripley. 2002. *Modern Applied Statistics with
s*. Fourth. New York: Springer. <https://www.stats.ox.ac.uk/pub/MASS4/>.
