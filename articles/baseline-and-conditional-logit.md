# The relation between baseline logit and conditional logit models

Baseline-category logit models can be expressed as particular form of
conditional logit models. In a conditional logit model (without random
effects) the probability that individual i chooses alternative j from
choice set \mathcal{S}\_i is

\pi\_{ij} =
\frac{\exp(\eta\_{ij})}{\sum\_{k\in\mathcal{S}\_i}\exp(\eta\_{ik})}

where

\eta\_{ij} = \alpha_1x\_{1ij}+\cdots+\alpha_qx\_{qij}

In a baseline-category logit model, the set of alternatives is the same
for all individuals i that is \mathcal{S}\_i = {1,\ldots,q} and the
linear part of the model can be written like:

\eta\_{ij} = \beta\_{j0}+\beta\_{j1}x\_{i1}+\cdots+\beta\_{jr}x\_{ri}

where the coefficients in the equation for baseline category j are all
zero, i.e.

\beta\_{10} = \cdots = \beta\_{1r} = 0

After setting

\begin{aligned} x\_{(g\times(j-1))ij} = d\_{gj}, \quad
x\_{(g\times(j-1)+h)ij} = d\_{gj}x\_{hi}, \qquad \text{with }d\_{gj}=
\begin{cases} 0&\text{for } j\neq g\text{ or } j=g\text{ and } j=0\\
1&\text{for } j=g \text{ and } j\neq0\\ \end{cases} \end{aligned}

we have for the log-odds:

\begin{aligned} \begin{aligned} \ln\frac{\pi\_{ij}}{\pi\_{i1}}
&=\beta\_{j0}+\beta\_{ji}x\_{1i}+\cdots+\beta\_{jr}x\_{ri} \\
&=\sum\_{h}\beta\_{jh}x\_{hi}=\sum\_{g,h}\beta\_{jh}d\_{gj}x\_{hi}
=\sum\_{g,h}\alpha\_{g\times(j-1)+h}(d\_{gj}x\_{hi}-d\_{g1}x\_{hi})
=\sum\_{g,h}\alpha\_{g\times(j-1)+h}(x\_{(g\times(j-1)+h)ij}-x\_{(g\times(j-1)+h)i1})\\
&=\alpha\_{1}(x\_{1ij}-x\_{1i1})+\cdots+\alpha\_{s}(x\_{sij}-x\_{si1})
\end{aligned} \end{aligned}

where \alpha_1=\beta\_{21}, \alpha_2=\beta\_{22}, etc.

That is, the baseline-category logit model is translated into a
conditional logit model where the alternative-specific values of the
attribute variables are interaction terms composed of
alternativ-specific dummes and individual-specific values of
characteristics variables.

Analogously, the random-effects extension of the baseline-logit model
can be translated into a random-effects conditional logit model where
the random intercepts in the logit equations of the baseline-logit model
are translated into random slopes of category-specific dummy variables.
