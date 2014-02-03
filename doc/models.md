% Revolve: An adaptive dynamics toolkit for R
% Daniel Falster; Rich FitzJohn

\newpage

# Background


The theory of adaptive dynamics offers a conceptual framework for understanding how phenotypic selection, driven by the selective
forces that arise through density- and frequency-dependent interactions
and population dynamics, shapes biological communities, in particular the distribution of quantitative traits [@geritz_evolutionarily_1998; @dieckmann_dynamical_1996; @dieckmann_can_1997; @dieckmann_adaptive_2007].  One direction for research in adaptive dynamics is to identify the ecological scenarios that allow for multiple species to coexist in a competitive context. The emphasis here is on the processes giving rise to fitness. A second direction for research is to describe the algorithms for modelling trait evolution. Several types of model can be identified , according to the different assumptions made about the evolutionary process, specifically the frequency and size of mutations. For the most part, these two directions are orthogonal: any given ecological model can be implemented under a range of assembly models.

The goal of this package is to provide functions for implementing a range of commonly used analysis techniques in adaptive dynamics research. These techniques are illustrated using a number of ecological models taken from published research. This includes some one-dimensional models (i.e. single trait evolving):

* @geritz_evolutionarily_1998 -- Evolution of seed size in multi-patch model with local adaptation
* @geritz_evolutionary_1999 --  Evolution of seed size with size-asymmetric competition for safe sites [results depend on starting conditions]
* @kisdi_evolutionary_1999 --  Size asymmetric competition
* @parvinen_evolutionary_2002 --  Evolutionary branching of dispersal strategies in structured metapopulations
* @parvinen_disturbance-generated_2009 --  Size structured meta population
* @bonsall_life_2004 --  Evolution of ...
* @dieckmann_origin_1999 --  Evolution under competition for a limiting resource
* @leimar_limiting_2013 --  Evolution of continuous traits with different competition kernels
* @dandrea_revising_2013 --  Evolution of seed size under fecundity-tolerance trade-off
* @scheffer_self-organized_2006 -- Evolution under competition for a limiting resource
* @leimar_limiting_2013 -- Limiting similarity, species packing, and the shape of competition kernels;

also some multi-dimensional models (two or more traits evolving):

* @dieckmann_origin_1999 --  Evolution under competition for a limiting resource with assortative mating
- @ito_new_2007 -- Multi-trait evolution with disruptive and directional selection.

Throughout all of what follows, we assume that a fitness function can be specified, giving the long-term rate of increase for a rare mutant in a resident community. Specifying this function implies two key assumptions.

1. Population sizes are sufficiently large and well mixed. The effect of this assumption is to remove stochastic elements from the population dynamics, which approach their deterministic limits in large populations.
2. The external (abitoic) environment is constant.

When stochastic population dynamics are of interest a finite size individual-based simulation model will be required.

# Definitions

Consider a community of $N$ species potentially differing from one another in $K$ quantitative traits. Let:

- $N$ be the number of species
- $K$ be the number of quantitative traits
- $x_{ij}$ be the value of the $j^{\textrm{th}}$ trait in species $i$
- $y_{i}$ be the abundance of individuals in species $i$ .

The phenotypic composition of the community is then described by the vector
of trait values
\begin{equation}x=\left(x_{11}, \ldots, x_{1K}, \ldots, x_{N1} , \ldots, x_{NK}\right)\end{equation}

and their corresponding abundances

\begin{equation}y=(y_{1},\ldots,y_{N}).\end{equation}

Define $f(x^\prime,x, y)$ to be a scalar-valued function giving
the per capita rate of increase (fitness) for individuals with trait values
$x^\prime \in \Re^K$ growing in an environment shaped by resident
community $(x,y)$.

Let $\bar{y}$ be the value of $y$ satisfying

\begin{equation}\label{ybar}  \hat{f}(x_{i},x,\bar{y})=0 \, \forall i. \end{equation}

Thus when $y = \bar{y}$, residents are at their demographic attractors (time scale separation). We then write
\begin{equation}\hat{f}(x^\prime,x) = f(x^\prime,x, y)\bigg|_{y=\bar{y}}.\end{equation}
for the invasion fitness of the mutant, i.e. per capita rate of increase in a
demographically stable resident population.

# Assembly algorithms

A number of assembly algorithms can be used with any given fitness function. The models are distinguished by the frequency and size of mutations, as illustrated in Fig. \ref{fig:D07-F1}. Here we adopt the terminology of @dieckmann_adaptive_2007.

If mutations occur frequently, evolution proceeds via **polymorphic and stochastic** model of trait evolution (Fig. [fig:D07-F1]a). In this model, polymorphic distributions of trait values stochastically drift and diffuse through selection and mutation.

![**Models of adaptive dynamics**
, from @dieckmann_adaptive_2007. Panel (a) illustrates the individual-based birth-death-mutation process (polymorphic and stochastic), panel (b) shows an evolutionary random walk (monomorphic and stochastic), panel (c) represents the gradient-ascent model (monomorphic and deterministic, described by the canonical equation of adaptive dynamics), and panel (d) depicts an evolutionary reaction- diffusion model (polymorphic and deterministic).\label{fig:D07-F1}](images/Dieckmann2007-Fig1.pdf)

When mutations are rare, ecological dynamics proceed much faster than evolutionary dynamics, allowing a formal time-scale separation between evolutionary and ecological processes (Fig. [fig:D07-F2]). Thus, whenever a successful mutation arises, it is allowed to spread before any additional mutations appear. This leads to the **monomorphic** model, so called because each species is represented by a single strategy. The term 'oligomorhpic' is also used to describe communities with several species (*oligo-*, containing a relatively small number of units).

Another assumption about the size of mutation then distinguishes between
two modes of adaptive evolution that can be handled using the
monomorphic model. If mutations are large, evolution proceeds as directed 'random walk' through trait space (panel b in Fig. \ref{fig:D07-F1}); the **monomorphic and stochastic** model.

On the other hand, if mutations are small, then evolution can be
modelled deterministically via gradient ascent, as a standard non-linear, dynamical system (panel c in Fig. \ref{fig:D07-F1}). This is the **monomorphic and deterministic** model.

![**Formal relations between the models of adaptive dynamics**
, from @dieckmann_adaptive_2007. The four classes of model are shown as rounded boxes, and the three derivations as arrows. Arrow labels highlight key assumptions.\label{fig:D07-F2}](images/Dieckmann2007-Fig2.pdf)

## Polymorphic and stochastic model

The polymorphic and stochastic model involves the fewest assumptions and code. Clearly, a stochastic approach is required whenever the key assumption of the deterministic model - that mutations are small - is violated. A stochastic approach may also be preferred in high dimensional systems, with many traits or
species, as solving of demographic attractors becomes more difficult. At a more fundamental level, stochastic approaches are inevitably more realistic, e.g. mutations are not rare in real communities.

### Numerical approach

The following algorithm can be used to step a stochastic system:

1.  Initialise the phenotypes $x=\left(x_{1}, \ldots, x_{i}, \ldots, x_{N}\right)$ for $N$ resident lineages at time $t=0$ (set $N=1$ for an initially monomorphic community). Define the extinction threshold $\epsilon$.

2.  \label{ps:Restart}
3.  ....
4.  .....ncrease $N$ by 1 and set $x_{N+1} = x_i^\prime$ .  Return to step \ref{ps:Restart}.

## Monomorphic systems

In many cases want to solve for the demographic attractor of the resident community, i.e. find values of $y$ satisfying equation \ref{ybar}. Biologically, we are assuming that mutations are rare and that a residents approach their deomgraphic attractor before a new mutation arises.

### Numerical approach
The following techniques may be applied to find  $\bar{y}$:

1.  **Analytical**: ideal, but only possible in some models and then only with single resident.
2.  **Multi-dimensional root solving**: This approach works well only if you have good initial guess for the solution. In systems with many species, obtaining a good initial guess is problematic.
3.  **Iteration**: Using a difference equation, $y_{i,t+1} = y_{i,t} \times 10^{\hat{f}(x_i,x, y_{t})}$, iterate the population until stable. This approach is fail proof, but a large number of iterations may be required to reach stability, especially when fitness is close to zero.
4.  **Solve system of ODEs**: Express problem as $\frac{\partial}{dt} y_{i} = y_{i} \hat{f}(x_i,x, y)$, then advance the system using an adaptive ODE solver. Depending on the system a stiff solver may perform better, this in turn requires calculation of jacobian.

For the root solving and ODE approach, it may be preferable to make $\log y$ the state variable, as this prevents negative population sizes being generated by the solver.

## Monomorphic and stochastic model

This model assumes infrequent but large mutations; evolution then follows a directed random walk through trait space [@dieckmann_dynamical_1996; @dieckmann_adaptive_2007].  Mathematically, the model is described by a master equation for the probability density $P(x,t)$ of realising a given phenotypic
distribution $x$ at time $t$ [see @dieckmann_dynamical_1996]. The shape of the fitness landscape influences the probability per unit time that a species transitions from its current trait value to another trait value in that direction. Evolution within each species then follows a Markovian *trait substitution sequence*. Mutational processes interact with the selective forces generated by ecological interaction to determine evolutionary trajectories.

### Numerical approach

The most efficient way to implement the monomorphic stochastic model is
using Gillespie’s minimal process method
[@gillespie_general_1976; @dieckmann_dynamical_1996]. Rather than
stepping the system over a given time interval and considering any
transitions that may occur in that time period, the system is stepped
between successive events by drawing the next invading mutant trait
value and the time at which the invasion occurs from appropriate
distributions. The method is founded on the idea that with a given
mutation rate, the distribution of waiting times follows an exponential
distribution. Appropriate algorithms are given by
@dieckmann_dynamical_1996, @ito_new_2007, @brannstrom_emergence_2010.

The following algorithm is adapted from @ito_new_2007:

1.  Initialise the phenotypes
    $x=\left(x_{1}, \ldots, x_{i}, \ldots, x_{N}\right)$ for $N$
    resident lineages at time $t=0$ (set $N=1$ for an initially
    monomorphic community). Define the extinction threshold $\epsilon$.

2.  \label{ms:Restart} Calculate equilibrium population sizes
    $y=\bar{y}$ satisfying
    $\hat{f}(x_{i},x,y)=0 \, \forall i$.

3.  \label{ms:extinct} Check whether $\bar{y}_{i} < \epsilon$ $\forall i$; if
    so, delete phenotype $x_i$ and decrease $N$ accordingly. Return to
    step \ref{ms:Restart}.

4.  \label{ms:Weight} Calculate the rate $w_i =\mu_i \bar{y}_i(x) $ for the
    emergence of a mutant from phenotype $x_i$, the immigration rate
    $w_{N+1} =I$, and the total mutation rate $w=\sum_{i=1}^{N+1} w_i$.

5.  \label{ms:newMutant} Choose lineage $i$ with probability $w_i/w$.

6.  Choose a new phenotype $x_i^\prime$ according to the mutation
    probability density $M(x_i^\prime, x_{i}, \sigma(x_{i}))$ if
    $i \in (1,N)$, or in the case of immigration, $M_I(x_i^\prime)$.
    Update time $t$ by adding $\Delta t =-(1/w)\ln p$, where
    $0\leq p \leq 1$ is a uniformly distributed random number.

7.  Calculate the invasion fitness of the new
    phenotype and check if it can invade. The probability of successful invasion depends on the magnitude of $f$ -- however, even if $f$ is positive the mutant may fail because of stochastic demographic effects. For models with continuous-time birth-death process, the probability the mutant succeeds is given $p = max(f,0)/b$ [@dieckmann_dynamical_1996]. Estimating a similar probability in structured populations is no trivial task. As a rough approximation, the probability the mutant succeeds might be estimated as $p = \max(\log(B),0)/log(R)$, where $R$ and $B$ are, respectively, the expected number of offspring and the expected number of offspring reaching maturity produced by an individual reaching maturity. Choose a uniformly distributed random number $0\leq p_r \leq 1$. If
    $p_r > p$ the mutant fails to invade, return to step
    \label{ms:newMutant}. Otherwise, increase $N$ by 1 and set $x_{N+1} = x_i^\prime$ .  Return to step \ref{ms:Restart}.

The above recipe estimates the waiting time between successive mutants,
and then asks whether such mutants can invade. It is also possible to
estimate the waiting time between successful mutants, using an extension
of this model [see @brannstrom_emergence_2010].
Although superior in many ways, this revised technique requires an
integration over the fitness landscape, which is costly in systems with
more than a single trait and is to be avoided.

## Monomorphic and deterministic model

Canonical equation....

## Other variants

- can extend to include genomes of diallelic loci [e.g. @dieckmann_origin_1999]see [issue here](https://github.com/dfalster/Revolve/issues/1)
- also allow for assortative mating of neutral marker trait [e.g. @dieckmann_origin_1999]

# A catalogue of ecological models

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
K(x^\prime)=K_0 \exp \left(-\frac{{x^\prime}^2}{2\sigma_K^2}\right),\end{equation}

while competition kernel is given by

\begin{equation} \label{eq:DD99-C}
C(x^\prime,x)=\exp\left(-\frac{(x^\prime-x)^2}{2\sigma_C^2}\right).
\end{equation}

The parameters $\sigma_k$ and $\sigma_C$ are scaling factors for the resource
distribution and competition kernels respectively, while $K_0$ is a parameter for maximum population density.

### Parameters

The values used in the paper are $r=1$, $K_0=500$, $\sigma_K=1$, and $\sigma_C=0.4$, with $\sigma_K$ and $\sigma_C$ varied across a range from 0-2 in figure 4. Mutations are generated according at rate $\mu = 0.001$, so total mutations per step is given by $r \, y_i \, \mu$. [NB: do we want to set mutation rate according to birth rate or total fitness ? Paper uses birth rate, but not all models have birth and death functions so perhaps take according to fitness?]

### History

According to @ito_new_2007, the model encapsulated in eqs. \ref{eq:DD99-fit1}-\ref{eq:DD99-C} has *"been investigated by many authors, including Christiansen and Fenchel (1977), Christiansen and Loeschcke (1980), Slatkin (1980), Case (1981), Seger (1985), Taper and Case (1985), Vincent et al. (1993), Metz et al. (1996), Doebeli (1996a),Drossel and McKane (2000), Day (2000), Ackermann and Doebeli (2004), and Doebeli et al. (2007)."*

## @ito_new_2007 -- Multi-trait evolution with competition for several limiting resources

@ito_new_2007 suggest an extension of the model from @dieckmann_origin_1999 to include a second resource and trait dimension. This model corresponds to the situation where there is [what shaped resource distribution - bivariate Gaussian?]. @vukics_speciation_2003 also investigate a similar extension. In principle these model could be extended to an arbitrary number of dimensions and resources.

Fitness in the multi-resource model is calculated as per eq. \ref{eq:DD99-fit1}, but note that $x^\prime$ and $x$ are now multi trait vectors. Multivariate extensions of eqs. \ref{eq:DD99-C}-\ref{eq:DD99-K} are then:

\begin{equation} \label{eq:DD99-K2}
K(x^\prime)=K_0 \exp \left(-\sum_{i}\frac{{{x_i}^\prime}^2}{2{\sigma_{Ki}}^2} \right),\end{equation}

and

\begin{equation} \label{eq:DD99-C2}
C(x^\prime,x)=\exp\left(-\sum_{i=1} \frac{({x_i}^\prime-x_i)^2}{2{\sigma_{Ci}}^2}\right), \end{equation}

where $i$ corresponds to the number of resources and trait dimensions.

### Parameters

Unknown because did not @ito_new_2007 did not implement.

## @ito_new_2007 -- Multi-trait evolution with directional selection

@ito_new_2007 extend the model from @dieckmann_origin_1999 to include a second trait dimension, but rather than treating the second trait in same fashion as first -- i.e. based on some stable resource distribution,  @ito_new_2007 set the second trait up to allow continuous directional selection.

In this model, fitness a new term corresponding to directional selection on either birth or death rate in the $j^{th}$ trait is added to the claulcation of fitness from eq. \ref{eq:DD99-fit1}:

\begin{equation} \label{eq:ID07-fit}  \hat{f}(x^\prime,x,y) = r\left(1- \frac{\sum_{i=1}^N y_i \, C(x^\prime,x_i)} {K(x^\prime)} \right) + d_0 \, D_j(x^\prime,x), \end{equation}

where $d_0$ is a constant influencing the strength of directional selection,

\begin{equation} \label{eq:ID07-D}   D_j(x^\prime,x) = x_j^\prime - \bar{x_j} \end{equation}

gives the difference between the $j^{th}$ trait and the community for mutant average:

\begin{equation} \label{eq:ID07-xj}  \bar{x_j} = \frac{\sum_{i=1}^N x_ij \, y_i}{\sum_{i=1}^N y_i} \end{equation}

In this formulation, directional selection in the   $j^{th}$ trait continues indefinitely.

### Parameters

General parameters used $K_0=100,000$, $r=1$,  $\mu=10^{-5}$.

Resource distribution:

- medium: $\sigma_k = 0.2, \sigma_C = 0.15$
- wide: $\sigma_k = 0.2, \sigma_C = 0.07$

Directional selection:

- absent: $d_0=0$
- weak: $d_0=0.5$
- strong: $d_0=10$ (did not show results)


### Results

The model produces three qualitatively different outcomes (see Fig 1), although the second type is only described verbally:

1. No directional selection: When $d_0=0$, the model produces a stable distribution of types, where number of types in each dimension depends on the width of resource distribution compared to the competition kernel
2. Strong directional selection: When directional selection in trait y is very strong, evolution in trait y is so swift that it prevents diversification in trait x. After the population has converged to x=0, it thus merely keeps evolving along this line in response to the directional selection pressure in y.
3. Weak directional selection: "When directional selection in trait y is finite and sufficiently weak, the initial branching in trait x is followed by a pattern of recurrent adaptive radiations and extinctions. This occurs because stochastic effects —- resulting from mutations as well as from the demography of finite populations -- cause the spontaneous breaking of the initial symmetry between the two diverging lineages. In particular, the population sizes of the two lineages will never be exactly equal. Since the more abundant lineage can evolve faster, it will move ahead in the race of responding to the weak directional selection pressure in trait y, thus increasing its relative fitness in terms of trait y......The evolutionary interplay between these two lineages causes the asymmetry in their population sizes to grow and their trait values y to diverge. This positive feedback continues until the posterior lineage becomes extinct. The positive feedback may be intensified by an additional effect: once the anterior lineage becomes sufficiently dominant, it experiences so little competition from the posterior lineage that it reverses its direction of gradual evolution in trait x, thus pushing the posterior lineage toward lower carrying capacity and accelerating its demise."

## Gertiz 1999 -- Evolution of seed size with size-asymmetric competition for safe sites

Relevant papers are @geritz_competition_1988, @geritz_evolutionary_1999, @geritz_evolutionarily_1995. These papers explore
dynamics leading to a stable polymorphism in seed size among plants
competing for safe sites (i.e. models falling with patch dynamics
frame-work). In brief:

-   @geritz_competition_1988 develops a model for competition and life-history evolution in safe-sites. Conditions for demographic equilibrium, invasion fitness and multi-specie coexistence are derived.

-   @geritz_evolutionarily_1995 found that small-scale spatial variation in seedling density favors the evolution of variation in seed size within individual plants if competition among seedlings is sufficiently asymmetric. The range of seed sizes is predicted to increase with adult resource status, reduction in juvenile mortality and more-even seed dispersal.

-   @geritz_evolutionary_1999 extends @geritz_competition_1988 and @geritz_evolutionarily_1995 models to incorporate differing     degrees of size-assymetric competition and juvenile mortality. Generalizes scenarios leading to progressively higher degrees of polymorphism.


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

# Extensions

Things that can be added to any model

- make quantitative traits depend on number of alleles, e.g. @dieckmann_origin_1999
- add neutral marker trait that influences mating probability, e.g. @dieckmann_origin_1999. This allows for assortative mating and linkage disequilibrium to develop; encourages disruptive selection.

# References
