# The relation between baseline logit and conditional logit models

Baseline-category logit models can be expressed as particular form of
conditional logit models. In a conditional logit model (without random
effects) the probability that individual $i$ chooses alternative $j$
from choice set $\mathcal{S}_{i}$ is

$$\pi_{ij} = \frac{\exp\left( \eta_{ij} \right)}{\sum\limits_{k \in \mathcal{S}_{i}}\exp\left( \eta_{ik} \right)}$$

where

$$\eta_{ij} = \alpha_{1}x_{1ij} + \cdots + \alpha_{q}x_{qij}$$

In a baseline-category logit model, the set of alternatives is the same
for all individuals $i$ that is $\mathcal{S}_{i} = {1,\ldots,q}$ and the
linear part of the model can be written like:

$$\eta_{ij} = \beta_{j0} + \beta_{j1}x_{i1} + \cdots + \beta_{jr}x_{ri}$$

where the coefficients in the equation for baseline category $j$ are all
zero, i.e.

$$\beta_{10} = \cdots = \beta_{1r} = 0$$

After setting

$$\begin{array}{r}
{x_{{(g \times {(j - 1)})}ij} = d_{gj},\quad x_{{(g \times {(j - 1)} + h)}ij} = d_{gj}x_{hi},\qquad{\text{with}\mspace{6mu}}d_{gj} = \begin{cases}
0 & {{\text{for}\mspace{6mu}}j \neq g{\mspace{6mu}\text{or}\mspace{6mu}}j = g{\mspace{6mu}\text{and}\mspace{6mu}}j = 0} \\
1 & {{\text{for}\mspace{6mu}}j = g{\mspace{6mu}\text{and}\mspace{6mu}}j \neq 0} \\
 & 
\end{cases}}
\end{array}$$

we have for the log-odds:

$$\begin{array}{r}
\begin{aligned}
{\ln\frac{\pi_{ij}}{\pi_{i1}}} & {= \beta_{j0} + \beta_{ji}x_{1i} + \cdots + \beta_{jr}x_{ri}} \\
 & {= \sum\limits_{h}\beta_{jh}x_{hi} = \sum\limits_{g,h}\beta_{jh}d_{gj}x_{hi} = \sum\limits_{g,h}\alpha_{g \times {(j - 1)} + h}\left( d_{gj}x_{hi} - d_{g1}x_{hi} \right) = \sum\limits_{g,h}\alpha_{g \times {(j - 1)} + h}\left( x_{{(g \times {(j - 1)} + h)}ij} - x_{{(g \times {(j - 1)} + h)}i1} \right)} \\
 & {= \alpha_{1}\left( x_{1ij} - x_{1i1} \right) + \cdots + \alpha_{s}\left( x_{sij} - x_{si1} \right)}
\end{aligned}
\end{array}$$

where $\alpha_{1} = \beta_{21}$, $\alpha_{2} = \beta_{22}$, etc.

That is, the baseline-category logit model is translated into a
conditional logit model where the alternative-specific values of the
attribute variables are interaction terms composed of
alternativ-specific dummes and individual-specific values of
characteristics variables.

Analogously, the random-effects extension of the baseline-logit model
can be translated into a random-effects conditional logit model where
the random intercepts in the logit equations of the baseline-logit model
are translated into random slopes of category-specific dummy variables.
