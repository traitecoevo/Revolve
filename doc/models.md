
# Background


One-dimensional models

* @geritz_evolutionarily_1998 -- Evolution of seed size in multi-patch model with local adaptation
* @geritz_evolutionary_1999 --  Evolution of seed size with size-asymmetric competition for safe sites [results depend on starting conditions]
* @kisdi_evolutionary_1999 --  Size asymmetric competition
* @parvinen_evolutionary_2002 --  Evolutionary branching of dispersal strategies in structured metapopulations
* @parvinen_disturbance-generated_2009 --  Size structured meta population
* @bonsall_life_2004 --  Evolution of ...
* @dieckmann_origin_1999 --  Evolution under competition for a limiting resource
* @leimar_limiting_2013 --  Evolution of continuous traits with different competition kernels
* @dandrea_revising_2013 --  Evolution of seed size under fecundity-tolerance trade-off

### 2D

* @dieckmann_origin_1999 --  Evolution under competition for a limiting resource and assortative mating
- @ito_new_2007 -- Mulit-trait evolution with disruptive and directional selection

# Definitions

Consider a community of $N$ species potentially differing from one
another in $K$ quantitative traits. Let:

- $N$ be the number of species
- $K$ be the number of quantitative traits
- $x_{ij}$ be the value of the $j^{\textrm{th}}$ trait in species $i$
- $y_{i}$ be the abundance of individuals in species $i$ .

The phenotypic composition of the community is then described by the vector
of trait values
\begin{equation}x=\left(x_{11}, \ldots, x_{1K}, \ldots, x_{N1} , \ldots, x_{NK}\right)\end{equation}

and their corresponding abundances

\begin{equation}y=(y_{1},\ldots,y_{N}).\end{equation}

Define $\hat{f}(x^\prime,x, y)$ to be a scalar-valued function giving
the per capita rate of increase (fitness) for individuals with trait values
$x^\prime \in \Re^K$ growing in an environment shaped by resident
community $(x,y)$.

Let $\bar{y}$ be the value of $y$ satisfying $\hat{f}(x_{i},x,\bar{y})=0 \, \forall i$. Thus when $y = \bar{y}$, residents are at their demographic attractors (time-scale separation). We then write
$f(x^\prime,x) = \hat{f}(x^\prime,x, \bar{y})$ for the invasion
fitness of the mutant, i.e. per capita rate of increase in a
demographically stable resident population.

# Model catalogue

## Dieckmann and Doebeli 1999 -- Evolution under competition for a limiting resource

@dieckmann_origin_1999 present a general model of trait evolution where there is competition for a limiting resource and the strength of competition between individuals declines with phenotypic distance according to a Gaussian function:

Fitness is calculated as difference between births and deaths, given by

\begin{equation} \label{eq:DD99-fit1}  \hat{f}(x^\prime,x,y) = r\left(1- \frac{\sum_{i=1}^N y_i \, C(x^\prime,x_i)} {K(x^\prime)} \right),\end{equation}

where

- $x$ is a continuous phenotypic trait, such as beak size, which affects resource consumption.
- $r$ is the intrinsic birth rate
- $K(x^\prime)$ is the resource kernel which determines the equilibrium population density, or carrying capacity
- $C(x^\prime, x_i)$ is a competition kernel describing the relative change in death rate of individuals of type $x^\prime$ due to competition by individuals of type $x_i$.

Note that the birth rate $r$ is trait- and density- independent, so competition comes only via effects on death rate. The carrying capacity of type $x^\prime$ individuals is given by

\begin{equation} \label{eq:DD99-K}
K(x^\prime)=K_0 \exp \left(-\frac{x^2}{2\sigma_K^2}\right),\end{equation}

while competition kernel is given by

\begin{equation} \label{eq:DD99-C}
C(x^\prime,x)=\exp\left(-\frac{(x^\prime-x)^2}{2\sigma_C^2}\right).
\end{equation}

The parameters $\sigma_k$ and $\sigma_C$ are scaling factors for the resource
distribution and competition kernels respectively, while $K_0$ is a parameter for maximum population density.

### Parameters

The values used in the paper are $r=1$, $K_0=500$, $\sigma_K=1$, and $\sigma_C=0.4$, with $\sigma_K$ and $\sigma_C$ varied across a range from 0-2 in figure 4.

## Gertiz 1999 -- Evolution of seed size with size-asymmetric competition for safe sites

Relevant papers are @geritz_competition_1988, @geritz_evolutionary_1999, @geritz_evolutionarily_1995. These papers explore
dynamics leading to a stable polymorphism in seed size among plants
competing for safe sites (i.e. models falling with patch dynamics
frame-work). In brief:

-   @geritz_competition_1988 develops a model for competition and life-history
    evolution in safe-sites. Conditions for demographic equilibrium,
    invasion fitness and multi-specie coexistence are derived.

-   @geritz_evolutionarily_1995 found that small-scale spatial variation in seedling
    density favors the evolution of variation in seed size within
    individual plants if competition among seedlings is sufficiently
    asymmetric. The range of seed sizes is predicted to increase with
    adult resource status, reduction in juvenile mortality and more-even
    seed dispersal.

-   @geritz_evolutionary_1999 extends @geritz_competition_1988 and @geritz_evolutionarily_1995 models to incorporate differing     degrees of size-assymetric competition and juvenile mortality. Generalizes scenarios leading to progressively higher degrees of
    polymorphism.


### Safe-sites framework

Safe-sites are defined as an amount of vacant ground capable of
maintaining a single established adult plant. Occupation of a site gives
each adult access to an amount of resources $R$ which can be used for seed
production. Assume a large number of such sites, such that the fraction
of sites occupied by species i can be represented as a continuous
variable $P_i$. The models assume all occupied patches are vacated each
year (i.e. annual plants), although an extension to perennial case is
considered [see @geritz_evolutionarily_1995]. Except for seed size, different types are assumed to be identical in all other characteristics.

If all patches are linked by dispersal with no spatial patterning
(i.e., the island model) then the distribution of seeds among sites
will follow a Poisson process, with a rate parameter $y_i$ determined
by the mean number of seeds per patch produced by species i. The
probability of $n$ seeds landing in any given patch is then

\begin{equation}\label{eq:Poiss} p(n; y_i) = \frac{e^{-y_i}y_i^n}{n!} \sim \mathrm{Poisson}(y_i),\end{equation}

while the equilibrium number of patches occupied is given by

\begin{equation}P^* = 1 - e^{-y_i}.\end{equation}

(i.e., $1 - p(0;n_i)$).

### Fitness calculations

For seeds moving from establishment to maturity, Geritz assumes first a
period of frequency-independent growth, where offspring survival is
dependent on the seed size, followed by a period of competitive
interactions, where survival is determined by the relative seed-mass of
an individual compared to all other individuals in the same site (there
is no competition among sites). Let survival during the pre-competitive
and competitive phases be $s(x^\prime)$ and $g(x^\prime,x)$ respectively, where $x^\prime$ and $x$ are the offspring sizes of mutant and resident populations. Fitness is then given by

\begin{equation}\label{eq:fit1} \hat{f}(x^\prime,x, y) = \frac{R}{x^\prime}s(x^\prime)g(x^\prime,x, y).\end{equation}

Survival during the frequency-independent phase is modelled as

\begin{equation}s(x) = \max \{0, 1-2e^{-\beta x}\}\end{equation}

which gives as monotonic increasing curve with increasing $x$, approaching
towards 1 for large $x$.

Geritz assumes size-asymmetric competition among seedlings within each
patch is operating to determine survival during the competitive phase.

The probability of a seed with strategy $x^\prime$ establishing within a patch is given by

\begin{equation}\frac{c(x^\prime)}{c(x^\prime)+k_1c(x_1)+\cdots+k_Nc(x_N)},\end{equation}

where $k_1, \ldots, k_N$ is the number of seeds arriving and $c(m)$ is the competitive ability for each strategy. Competitive ability is given by the equation

\begin{equation}c(x) = \exp\left(\alpha \, x\right). \end{equation}


Thus to calculate $g(x^\prime,x)$, we need to sum the probability of
establishment over the entire range of possible resident densities and
their probability of occurrence. Using equation \ref{eq:Poiss} for the probability of drawing $x$ seeds of species $i$, given its
density $y_i$, we obtain

\begin{equation}\label{eq:g} g(x^\prime,x,y)= \sum_{k_1=0}^\infty\cdots \sum_{k_N=0}^\infty  \frac{c(x^\prime)}{c(x^\prime)+k_1c(x_1)+\cdots+k_Nc(x_N)} \times p(k_1; y_1) \times \cdots \times p(k_N; y_n).\end{equation}

### Parameter values

$\alpha \, R=6$, range - 0 -15

$\beta \, R=25$, range - 0-60

# References
