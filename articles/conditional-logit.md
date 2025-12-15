# Conditional logit models

Conditional logit models are motivated by a variety of considerations,
notably as a way to model binary panel data or responses in
case-control-studies. The variant supported by the package “mclogit” is
motivated by the analysis of discrete choices and goes back to McFadden
(1974). Here, a series of individuals $i = 1,\ldots,n$ is observed to
have made a choice (represented by a number $j$) from a choice set
$\mathcal{S}_{i}$, the set of alternatives at the individual’s disposal.
Each alternatives $j$ in the choice set can be described by the values
$x_{1ij},\ldots,x_{1ij}$ of $r$ attribute variables (where the variables
are enumerated as $i = 1,\ldots,r$). (Note that in contrast to the
baseline-category logit model, these values vary between choice
alternatives.) Conditional logit models then posit that individual $i$
chooses alternative $j$ from his or her choice set $\mathcal{S}_{i}$
with probability

$$\pi_{ij} = \frac{\exp\left( \alpha_{1}x_{1ij} + \cdots + \alpha_{r}x_{rij} \right)}{\sum\limits_{k \in \mathcal{S}_{i}}\exp\left( \alpha_{1}x_{1ik} + \cdots + \alpha_{r}x_{rik} \right)}.$$

It is worth noting that the conditional logit model does not require
that all individuals face the same choice sets. Only that the
alternatives in the choice sets can be distinguished from one another by
the attribute variables.

The similarities and differences of these models to baseline-category
logit model becomes obvious if one looks at the log-odds relative to the
first alternative in the choice set:

$$\ln\frac{\pi_{ij}}{\pi_{i1}} = \alpha_{1}\left( x_{1ij} - x_{1i1} \right) + \cdots + \alpha_{r}\left( x_{rij} - x_{ri1} \right).$$

Conditional logit models appear more parsimonious than baseline-category
logit models in so far as they have only one coefficient for each
independent variables.\[^1\] In the “mclogi" package, these models can
be estimated using the function
[`mclogit()`](https://melff.github.io/mclogit/reference/mclogit.md).

My interest in conditional logit models derives from my research into
the influence of parties' political positions on the patterns of voting.
Here, the political positions are the attributes of the alternatives and
the choice sets are the sets of parties that run candidates in a
countries at various points in time. For the application of the
conditional logit models, see Elff (2009).

## References

Elff, Martin. 2009. “Social Divisions, Party Positions, and Electoral
Behaviour.” *Electoral Studies* 28 (2): 297–308.
<https://doi.org/10.1016/j.electstud.2009.02.002>.

McFadden, Daniel. 1974. “Conditional Logit Analysis of Qualitative
Choice Behaviour.” In *Frontiers in Econometrics*, edited by Paul
Zarembka, 105–42. New York: Academic Press.
